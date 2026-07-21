# ============================================================
# Remap short RNA samples — files already trimmed to 26bp
# Trimmed files named: SRR03245X.26bp_5prime.fq.gz
# ============================================================

BOWTIE2="/shared/ucl/apps/bowtie2/bowtie2-2.2.5/bowtie2"
SAMTOOLS="/shared/ucl/apps/samtools/1.11/gnu-4.9.2/bin/samtools"
BEDTOOLS="/shared/ucl/apps/bedtools/2.25.0/gnu-4.9.2/bin/bedtools"

SRRS=(  _dummy      SRR032457     SRR032458     SRR032459     SRR032460  )
NAMES=( _dummy      scrna5pr_rep1 scrna5pr_rep2 scrna3pr_rep1 scrna3pr_rep2 )

SRR=${SRRS[$SGE_TASK_ID]}
NAME=${NAMES[$SGE_TASK_ID]}

OUTDIR="nechaev_dm6"
INDEX="hsp90_chipseq_dm6/genome/dm6"
THREADS=16

set -euo pipefail
log() { echo "[$(date '+%H:%M:%S')] [task${SGE_TASK_ID}:${NAME}] $*"; }

TRIMMED="$OUTDIR/qc/${SRR}.26bp_5prime.fq.gz"
BAM="$OUTDIR/bam/${NAME}.sorted.bam"
FINAL_BAM="$OUTDIR/bam/${NAME}.final.bam"

log "Trimmed file: $TRIMMED"
log "Output BAM: $FINAL_BAM"

# ============================================================
# Align
# ============================================================
log "=== Aligning ==="

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

# No dedup for short RNA
cp $BAM $FINAL_BAM
cp ${BAM}.bai ${FINAL_BAM}.bai

$SAMTOOLS flagstat $FINAL_BAM > $OUTDIR/bam/${NAME}_flagstat.txt

# ============================================================
# Generate strand-specific bedgraphs
# ============================================================
log "=== Generating bedgraphs ==="

MAPPED=$($SAMTOOLS view -c $FINAL_BAM)
SCALE=$(echo "1000000 / $MAPPED" | bc -l)

$SAMTOOLS view -b -F 16 $FINAL_BAM \
    | $BEDTOOLS genomecov -ibam stdin -bg -scale $SCALE \
    | sort -k1,1 -k2,2n \
    > $OUTDIR/bedgraph/${NAME}_fwd.bedgraph

$SAMTOOLS view -b -f 16 $FINAL_BAM \
    | $BEDTOOLS genomecov -ibam stdin -bg -scale $SCALE \
    | sort -k1,1 -k2,2n \
    > $OUTDIR/bedgraph/${NAME}_rev.bedgraph

log "Fwd: $(wc -l < $OUTDIR/bedgraph/${NAME}_fwd.bedgraph) entries"
log "Rev: $(wc -l < $OUTDIR/bedgraph/${NAME}_rev.bedgraph) entries"

log "=== ${NAME} (${SRR}) done ==="


#merge bams
samtools merge scrna3pr_rep1.sorted.bam scrna3pr_rep2.sorted.bam > 3pr_sorted.bam
samtools merge scrna5pr_rep1.sorted.bam scrna5pr_rep2.sorted.bam > 5pr_sorted.bam

#used bedtools to count in 150bp around tss
bedtools multicov -bams 3pr_sorted.bam 5pr_sorted.bam -bed dm6_tss_150bp_final.bed > scrna_tss_counts.bed

# Total genes
wc -l scrna_tss_counts.bed

# Paused genes (any signal in col 7 or 8)
awk '$7 > 0 || $8 > 0' scrna_tss_counts.bed | wc -l

# Hsp90 targets that are paused
awk '$7 > 0 || $8 > 0' scrna_tss_counts.bed | \
    cut -f4 | \
    grep -Fxf ../hsp90_promoter_genes.txt | wc -l
