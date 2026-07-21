SAMTOOLS=/home/rm786/rds/rds-rmash-sHE3M4BjkMM/software/samtools-1.22.1/samtools
INDIR=/home/rm786/rds/rds-rmash-sHE3M4BjkMM/CutTag/aligned
OUTDIR=/home/rm786/rds/rds-rmash-sHE3M4BjkMM/CutTag/sort_n

mkdir -p $OUTDIR

for BAM in $INDIR/*.sorted.bam; do
    filename=$(basename $BAM | awk -F '.bam' '{print substr($1,1,11)}')
    echo "Name sorting: $filename"
    $SAMTOOLS sort -@ 8 -n $BAM -o $OUTDIR/${filename}.sorted.n.bam
done

