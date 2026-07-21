# ============================================================
# Nechaev et al. 2010 (GSE18643) — remapping to dm6
#
# Paper processing: reads trimmed to 26nt, mapped with MOM
# (max oligonucleotide mapping), unique mappers only.
# We approximate with bowtie2 --very-sensitive -k 1, MAPQ>=20,
# and hard-trim reads to 26nt to match paper's read length.
#
# Samples:
#   1: SRR032455 — Pol II ChIP rep1 (41bp)
#   2: SRR032456 — Pol II ChIP rep2 (35bp)
#   3: SRR032457 — 5' short RNA rep1 (26bp)
#   4: SRR032458 — 5' short RNA rep2 (38bp)
#   5: SRR032459 — 3' short RNA rep1 (26bp)
#   6: SRR032460 — 3' short RNA rep2 (38bp)
# ============================================================

BOWTIE2="/shared/ucl/apps/bowtie2/bowtie2-2.2.5/bowtie2"
SAMTOOLS="/shared/ucl/apps/samtools/1.11/gnu-4.9.2/bin/samtools"
TRIM_GALORE="/shared/ucl/apps/trim_galore/0.6.10/TrimGalore-0.6.10/trim_galore"
JAVA="/shared/ucl/apps/temurin/17.0.2_8/jdk-17.0.2+8/bin/java"
SOFTWARE="/home/ucbtrm5/Scratch/software"
PICARD="$JAVA -jar $SOFTWARE/picard.jar"
SRA_SIF="$SOFTWARE/sra-tools_latest.sif"
BEDTOOLS="/shared/ucl/apps/bedtools/2.25.0/gnu-4.9.2/bin/bedtools"
CUTADAPT="/shared/ucl/apps/python/bundles/python39-6.0.0/venv/bin/cutadapt"

# ---- Samples ----
SRRS=(  _dummy      SRR032455     SRR032456     SRR032457     SRR032458     SRR032459     SRR032460  )
NAMES=( _dummy      polII_rep1    polII_rep2    scrna5pr_rep1 scrna5pr_rep2 scrna3pr_rep1 scrna3pr_rep2 )
TYPES=( _dummy      chipseq       chipseq       shortrna      shortrna      shortrna      shortrna )

SRR=${SRRS[$SGE_TASK_ID]}
NAME=${NAMES[$SGE_TASK_ID]}
TYPE=${TYPES[$SGE_TASK_ID]}

OUTDIR="nechaev_dm6"
INDEX="hsp90_chipseq_dm6/genome/dm6"
THREADS=16

mkdir -p $OUTDIR/{fastq,bam,qc,bedgraph}

set -euo pipefail
log() { echo "[$(date '+%H:%M:%S')] [task${SGE_TASK_ID}:${NAME}] $*"; }

log "Starting — SRR=${SRR} TYPE=${TYPE} THREADS=${THREADS}"

# ============================================================
# STEP 1: Download
# ============================================================
log "=== [1/4] Downloading ${SRR} ==="

if [ ! -f "$OUTDIR/fastq/${SRR}.fastq.gz" ]; then
    singularity exec $SRA_SIF fasterq-dump $SRR \
        --outdir $OUTDIR/fastq \
        --threads $THREADS \
        --progress
    gzip -f $OUTDIR/fastq/${SRR}.fastq
else
    log "Already downloaded, skipping"
fi

# ============================================================
# STEP 2: Trim
# Paper: hard-trimmed all reads to 26nt from 5' end
# We replicate this with --hardtrim5 26 for short RNA
# For Pol II ChIP: standard adapter trimming only
# ============================================================
log "=== [2/4] Trimming ==="

TRIMMED="$OUTDIR/qc/${SRR}_trimmed.fq.gz"

if [ ! -f "$TRIMMED" ]; then
    if [ "$TYPE" == "shortrna" ]; then
        # Hard trim to 26nt to match paper, then adapter trim
        $TRIM_GALORE \
        --path_to_cutadapt $CUTADAPT \
            --cores 4 \
            --fastqc \
            --fastqc_args "--outdir $OUTDIR/qc" \
            --output_dir $OUTDIR/qc \
            --hardtrim5 26 \
            $OUTDIR/fastq/${SRR}.fastq.gz
    else
        # Pol II ChIP: standard adapter trimming
        $TRIM_GALORE \
        --path_to_cutadapt $CUTADAPT \
            --cores 4 \
            --fastqc \
            --fastqc_args "--outdir $OUTDIR/qc" \
            --output_dir $OUTDIR/qc \
            $OUTDIR/fastq/${SRR}.fastq.gz
    fi
else
    log "Already trimmed, skipping"
fi

# ============================================================
# STEP 3: Align + sort + dedup (ChIP only)
# Paper: unique mappers only via MOM
# We use bowtie2 -k 1 + MAPQ>=20 to approximate this
# No dedup for short RNA — genuine signal, not PCR artefact
# ============================================================
log "=== [3/4] Aligning ==="

BAM="$OUTDIR/bam/${NAME}.sorted.bam"
DEDUP="$OUTDIR/bam/${NAME}.dedup.bam"
FINAL_BAM="$OUTDIR/bam/${NAME}.final.bam"

if [ ! -f "$FINAL_BAM" ]; then

    $BOWTIE2 \
        -x $INDEX \
        -U $TRIMMED \
        -p $THREADS \
        --no-unal \
        -k 1 \
        --very-sensitive \
        2> $OUTDIR/bam/${NAME}_bowtie2.log \
    | $SAMTOOLS view -bS -q 20 \
    | $SAMTOOLS sort -@ $THREADS -o $BAM

    $SAMTOOLS index $BAM
    log "Aligned: $($SAMTOOLS view -c $BAM) reads"

    if [ "$TYPE" == "chipseq" ]; then
        $PICARD MarkDuplicates \
            I=$BAM \
            O=$DEDUP \
            M=$OUTDIR/bam/${NAME}_dup_metrics.txt \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT \
            TMP_DIR=/tmp 2>/dev/null
        $SAMTOOLS index $DEDUP
        mv $DEDUP $FINAL_BAM
        mv ${DEDUP}.bai ${FINAL_BAM}.bai
        log "Deduped: $($SAMTOOLS view -c $BAM) -> $($SAMTOOLS view -c $FINAL_BAM) reads"
    else
        # No dedup for short RNA
        cp $BAM $FINAL_BAM
        cp ${BAM}.bai ${FINAL_BAM}.bai
        log "Short RNA (no dedup): $($SAMTOOLS view -c $FINAL_BAM) reads"
    fi

    $SAMTOOLS flagstat $FINAL_BAM > $OUTDIR/bam/${NAME}_flagstat.txt

else
    log "Already complete, skipping"
fi

# ============================================================
# STEP 4: Generate bedgraph (CPM normalised)
# Strand-specific for short RNA — we generate both fwd and
# rev and will determine strandedness from the data
# ============================================================
log "=== [4/4] Generating bedgraph ==="

MAPPED=$($SAMTOOLS view -c $FINAL_BAM)
SCALE=$(echo "1000000 / $MAPPED" | bc -l)
log "Library size: $MAPPED reads, scale factor: $SCALE"

if [ "$TYPE" == "shortrna" ]; then
    # Forward strand reads (flag -F 16: not reverse complement)
    $SAMTOOLS view -b -F 16 $FINAL_BAM \
        | $BEDTOOLS genomecov -ibam stdin -bg -scale $SCALE \
        | sort -k1,1 -k2,2n \
        > $OUTDIR/bedgraph/${NAME}_fwd.bedgraph

    # Reverse strand reads (flag -f 16: reverse complement)
    $SAMTOOLS view -b -f 16 $FINAL_BAM \
        | $BEDTOOLS genomecov -ibam stdin -bg -scale $SCALE \
        | sort -k1,1 -k2,2n \
        > $OUTDIR/bedgraph/${NAME}_rev.bedgraph

    FWD=$(wc -l < $OUTDIR/bedgraph/${NAME}_fwd.bedgraph)
    REV=$(wc -l < $OUTDIR/bedgraph/${NAME}_rev.bedgraph)
    log "Bedgraphs: ${FWD} fwd entries, ${REV} rev entries"
    log "Strand ratio fwd/rev: $(echo "$FWD / $REV" | bc -l)"
else
    # Pol II ChIP: unstranded
    $BEDTOOLS genomecov -ibam $FINAL_BAM -bg -scale $SCALE \
        | sort -k1,1 -k2,2n \
        > $OUTDIR/bedgraph/${NAME}.bedgraph

    log "Bedgraph: $(wc -l < $OUTDIR/bedgraph/${NAME}.bedgraph) entries"
fi

log "=== ${NAME} (${SRR}) done ==="
