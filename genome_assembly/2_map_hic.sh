#!/bin/bash

## Mapping using arima hic script: 
#https://github.com/ArimaGenomics/mapping_pipeline/blob/master/arima_mapping_pipeline.sh

SRA='HiC_paired'
LABEL='NicV'
FILTER='./mapping_pipeline/filter_five_end.pl'
COMBINER='./mapping_pipeline/two_read_bam_combiner.pl'
PICARD='~/software/picard.jar'

# echo "### Step 0: Index reference"
# mkdir index
# cd index
# bwa index -p Ni2_pacbio.asm Ni2_pacbio.asm.fa

echo "### Step 1.A: FASTQ to BAM (1st)"
bwa mem -t 24 ./index/Ni2_pacbio.asm ../HiC/R1.fq.gz | samtools view -@ 24 -Sb - > HiC_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
bwa mem -t 24 Ni2_pacbio.asm ../HiC/R2.fq.gz | samtools view -@ 24 -Sb - > HiC_2.bam

echo "### Step 2.A: Filter 5' end (1st)"
samtools view -h HiC_1.bam | perl $FILTER | samtools view -Sb - > HiC_filt_1.bam

echo "### Step 2.B: Filter 5' end (2nd)"
samtools view -h HiC_2.bam | perl $FILTER | samtools view -Sb - > HiC_filt_2.bam

echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER HiC_filt_1.bam HiC_filt_2.bam samtools 10 | samtools view -bS -t Ni2_pacbio.asm.fa.fai - | samtools sort -@ 24 -o HiC_paired_tmp.bam -

echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=HiC_paired_tmp.bam OUTPUT=$SRA.bam ID=$SRA LB=$SRA SM=$LABEL PL=ILLUMINA PU=none

java -Xmx4G -Djava.io.tmpdir=temp/ -jar ~/software/picard.jar AddOrReplaceReadGroups INPUT=HiC_paired_tmp.bam OUTPUT=HiC_paired.bam ID=HiC_paired LB=HiC_paired SM=NicV PL=ILLUMINA PU=none
