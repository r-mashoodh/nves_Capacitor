# ============================================================
# Hsp90 ChIP-seq peak calling with MACS2
# Sawarkar et al. 2012 Cell — GSE31226 — dm6
#
# Parameters:
#   - Model building on (no --nomodel) — let MACS2 estimate
#     fragment length from cross-correlation
#   - Each ChIP rep matched to its own input
#   - --keep-dup all since Picard already removed duplicates
#   - --mfold 6 50 to guide model building as per paper
#   - -q 0.05 FDR threshold
# ============================================================

MACS2=~/.python3local/bin/macs2
BEDTOOLS=/shared/ucl/apps/bedtools/2.25.0/gnu-4.9.2/bin/bedtools
OUTDIR="hsp90_chipseq_dm6"
PEAKDIR="$OUTDIR/peaks_macs2_clean"

mkdir -p $PEAKDIR

set -euo pipefail
log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "=== [1/3] Calling peaks — rep1 ==="

$MACS2 callpeak \
    -t $OUTDIR/bam/hsp90_rep1.dedup.bam \
    -c $OUTDIR/bam/input_rep1.dedup.bam \
    -f BAM \
    -g 1.37e8 \
    -n hsp90_rep1 \
    --outdir $PEAKDIR \
    -q 0.05 \
    --keep-dup auto \
    --mfold 6 50 \
    --bdg \
    2> $PEAKDIR/hsp90_rep1_macs2.log

log "Rep1 peaks: $(wc -l < $PEAKDIR/hsp90_rep1_peaks.narrowPeak)"
log "Rep1 fragment length: $(grep 'predicted fragment length' $PEAKDIR/hsp90_rep1_macs2.log)"

log "=== [2/3] Calling peaks — rep2 ==="

$MACS2 callpeak \
    -t $OUTDIR/bam/hsp90_rep2.dedup.bam \
    -c $OUTDIR/bam/input_rep2.dedup.bam \
    -f BAM \
    -g 1.37e8 \
    -n hsp90_rep2 \
    --outdir $PEAKDIR \
    -q 0.05 \
    --keep-dup auto \
    --mfold 6 50 \
    --bdg \
    2> $PEAKDIR/hsp90_rep2_macs2.log

log "Rep2 peaks: $(wc -l < $PEAKDIR/hsp90_rep2_peaks.narrowPeak)"
log "Rep2 fragment length: $(grep 'predicted fragment length' $PEAKDIR/hsp90_rep2_macs2.log)"

log "=== [3/3] Intersecting replicates ==="

$BEDTOOLS intersect \
    -a $PEAKDIR/hsp90_rep1_peaks.narrowPeak \
    -b $PEAKDIR/hsp90_rep2_peaks.narrowPeak \
    -wa \
    > $PEAKDIR/hsp90_reproducible_peaks.bed

log "Reproducible peaks: $(wc -l < $PEAKDIR/hsp90_reproducible_peaks.bed)"
log "Done — peaks at $PEAKDIR/hsp90_reproducible_peaks.bed"
