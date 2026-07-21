
## aligning h3k4me3 cut&tag samples to index

## build index
#module load bowtie2-2.3.5-intel-17.0.4-ypnhvtd

#RDS=/home/rm786/rds/rds-rmash-sHE3M4BjkMM
#GENOME=$RDS/new_hic_genome/Ni2.purged.folded.filtered_scaffolds_final.sm.renamed.fa
#INDEX_DIR=$RDS/bowtie2_index

#mkdir -p $INDEX_DIR
#bowtie2-build --threads 8 $GENOME $INDEX_DIR/Ni2

## align
module load bowtie2-2.3.5-intel-17.0.4-ypnhvtd

RDS=/home/rm786/rds/rds-rmash-sHE3M4BjkMM
INDEX=$RDS/CutTag/bowtie2_index/Ni2
INDIR=$RDS/CutTag/CutTag_trimmed
OUTDIR=$RDS/CutTag/aligned

SAMTOOLS=/home/rm786/rds/rds-rmash-sHE3M4BjkMM/software/samtools-1.22.1/samtools
BOWTIE2=/home/rm786/rds/rds-rmash-sHE3M4BjkMM/software/bowtie2-2.5.5-linux-x86_64/bowtie2

mkdir -p $OUTDIR logs

R1_FILES=($INDIR/*_val_1.fq.gz)
R1=${R1_FILES[$SLURM_ARRAY_TASK_ID]}
R2="${R1/_1_val_1.fq.gz/_2_val_2.fq.gz}"

SAMPLE=$(basename $R1 _1_val_1.fq.gz)
SHORT=$(echo $SAMPLE | cut -d'_' -f1-3)

echo "Processing: $SHORT"
echo "R1: $R1"
echo "R2: $R2"

#if [[ -f $OUTDIR/${SHORT}.sorted.bam.bai ]]; then
#    echo "Skipping $SHORT - already complete"
#    exit 0
#fi

$BOWTIE2 \
    -x $INDEX \
    -1 $R1 \
    -2 $R2 \
    --local --very-sensitive-local \
    --no-mixed --no-discordant \
    --phred33 \
    -I 10 -X 700 \
    -p 16 \
    2> $OUTDIR/${SHORT}_bowtie2.log \
| $SAMTOOLS view -bS -q 20 - \
| $SAMTOOLS sort -@ 8 -o $OUTDIR/${SHORT}.sorted.bam

$SAMTOOLS index $OUTDIR/${SHORT}.sorted.bam
echo "Done: $SHORT"
