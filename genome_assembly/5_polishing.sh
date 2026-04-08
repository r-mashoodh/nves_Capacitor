#!/bin/bash

## Polishing round

## round 2
## map to purged fa

TRIMMED='/mnt/home5/zoology/rm786/software/blobtoolkit/coverage/trimmed'
OUTPUT_DIR='aligned'
REFERENCE='~/software/purge_dups/Ni2_pacbio.asm/seqs/Ni2_pacbio.asm.purged.fa'
INDEX='./purged_index/purged.fa'

## index
#mkdir purged_index
#cd purged_index
#bwa index -p purged.fa $REFERENCE
#cd ..


# Ensure the output directory exists
mkdir -p $OUTPUT_DIR

for i in $TRIMMED/*_1_val_1.fq.gz; 
do 
    # Extract the base filename without path and without '_1_val_1.fq.gz'
    filename=$(basename "$i" "_1_val_1.fq.gz")
    
    # Define the corresponding read pair file
    read1=$TRIMMED/${filename}_1_val_1.fq.gz
    read2=$TRIMMED/${filename}_2_val_2.fq.gz
    
    # Run BWA mem to map reads and pipe to samtools to convert SAM to BAM
    bwa mem -t 32 $INDEX $read1 $read2 | samtools view -@30 -bS - > $OUTPUT_DIR/${filename}_bwa.bam
    
    # Sort BAM
    samtools sort -@30 -o $OUTPUT_DIR/${filename}_bwa_sorted.bam $OUTPUT_DIR/${filename}_bwa.bam

    # Index BAM
    samtools index $OUTPUT_DIR/${filename}_bwa_sorted.bam

    # Optional: Add logging or echo statements for progress monitoring
    echo "Processed: $filename"
done

## run pilon
java -Xmx140G -jar ~/software/pilon-1.24.jar --genome ~/software/purge_dups/Ni2_pacbio.asm/seqs/Ni2_pacbio.asm.purged.fa --bam ../aligned/GA_HLNF3DSX2_L4_bwa_sorted.bam --vcf --output Ni2_purged.pilon2.fa --diploid --changes --threads 32

## fasta stats
~/software/gfastats/build/bin/gfastats Ni2_purged.pilon2.fa.fasta --nstar-report --stats
