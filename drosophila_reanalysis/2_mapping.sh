# ---- Full paths to all binaries ----
BOWTIE2="/shared/ucl/apps/bowtie2/bowtie2-2.2.5/bowtie2"
SAMTOOLS="/shared/ucl/apps/samtools/1.11/gnu-4.9.2/bin/samtools"
TRIM_GALORE="/shared/ucl/apps/trim_galore/0.6.10/TrimGalore-0.6.10/trim_galore"
FASTQC="/shared/ucl/apps/fastqc/0.11.8/FastQC/fastqc"
JAVA="/shared/ucl/apps/temurin/17.0.2_8/jdk-17.0.2+8/bin/java"
 
SOFTWARE="/home/ucbtrm5/Scratch/software"
PICARD="$JAVA -jar $SOFTWARE/picard.jar"
SRA_SIF="$SOFTWARE/sra-tools_latest.sif"
 
# ---- Samples ----
SRRS=( _dummy     SRR327604   SRR327605   SRR504969   SRR504970  )
NAMES=(_dummy     hsp90_rep1  hsp90_rep2  input_rep1  input_rep2 )
 
SRR=${SRRS[$SGE_TASK_ID]}
NAME=${NAMES[$SGE_TASK_ID]}
 
OUTDIR="hsp90_chipseq_dm6"
INDEX="$OUTDIR/genome/dm6"
THREADS=16
 
mkdir -p $OUTDIR/{fastq,bam,qc}
 
set -euo pipefail
log() { echo "[$(date '+%H:%M:%S')] [task${SGE_TASK_ID}:${NAME}] $*"; }
 
log "Starting — SRR=${SRR} THREADS=${THREADS}"
 
# ============================================================
# STEP 1: Download
# ============================================================
log "=== [1/3] Downloading ${SRR} ==="
 
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
# STEP 2: Trim + QC
# ============================================================
log "=== [2/3] Trimming ==="
 
TRIMMED="$OUTDIR/qc/${SRR}_trimmed.fq.gz"
 
if [ ! -f "$TRIMMED" ]; then
    $TRIM_GALORE \
        --cores 4 \
        --fastqc \
        --fastqc_args "--outdir $OUTDIR/qc" \
        --output_dir $OUTDIR/qc \
        $OUTDIR/fastq/${SRR}.fastq.gz
else
    log "Already trimmed, skipping"
fi
 
# ============================================================
# STEP 3: Align + coordinate-sort + dedup + name-sort
# ============================================================
log "=== [3/3] Aligning ==="
 
BAM="$OUTDIR/bam/${NAME}.sorted.bam"
DEDUP="$OUTDIR/bam/${NAME}.dedup.bam"
NAMESORTED="$OUTDIR/bam/${NAME}.dedup.nsorted.bam"
 
if [ ! -f "$NAMESORTED" ]; then
 
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
 
    $PICARD MarkDuplicates \
        I=$BAM \
        O=$DEDUP \
        M=$OUTDIR/bam/${NAME}_dup_metrics.txt \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=/tmp 2>/dev/null
 
    $SAMTOOLS index $DEDUP
    log "Deduped: $($SAMTOOLS view -c $BAM) → $($SAMTOOLS view -c $DEDUP) reads"
 
    $SAMTOOLS flagstat $DEDUP > $OUTDIR/bam/${NAME}_flagstat.txt
 
    $SAMTOOLS sort -n \
        -@ $THREADS \
        -o $NAMESORTED \
        $DEDUP
 
    log "Name-sorted BAM ready: $NAMESORTED"
 
else
    log "Already complete, skipping"
fi
 
log "=== ${NAME} (${SRR}) done ==="
fi

log "=== ${NAME} (${SRR}) done ==="
