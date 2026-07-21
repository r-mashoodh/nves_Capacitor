module load bowtie2/2.2.5

SOFTWARE="/home/ucbtrm5/Scratch/software"
OUTDIR="hsp90_chipseq_dm6"
THREADS=8

mkdir -p $OUTDIR/genome logs

set -euo pipefail
log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "Downloading dm6 FASTA ..."
wget -q -O $OUTDIR/genome/dm6.fa.gz \
    https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
gunzip $OUTDIR/genome/dm6.fa.gz

log "Downloading chrom sizes ..."
wget -q -O $OUTDIR/genome/dm6.chrom.sizes \
    https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes

log "Building Bowtie2 index ..."
bowtie2-build \
    $OUTDIR/genome/dm6.fa \
    $OUTDIR/genome/dm6

log "Done — index at $OUTDIR/genome/dm6"
