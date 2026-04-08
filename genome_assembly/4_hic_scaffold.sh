## now map hic to purged.pilon2.fa.fasta

mkdir map_hic
cd map_hic
bwa index -p purged.pilon.fa ../Ni2_purged.pilon2.fa.fasta


# then scaffold with yahs
SRA='R'
LABEL='nves_hic'
IN_DIR='/mnt/home5/zoology/rm786/beetle_HiC/HiC'
REF='/mnt/home5/zoology/rm786/beetle_HiC/pilon/round2/Ni2_purged.pilon2.fa.fasta'
INDEX='/mnt/home5/zoology/rm786/beetle_HiC/pilon/round2/map_hic/purged.pilon.fa'
FAIDX='/mnt/home5/zoology/rm786/beetle_HiC/pilon/round2/Ni2_purged.pilon2.fa.fasta.fai'
OUT='HiC.purged.polished'

RAW_DIR='raw'
FILT_DIR='filt'
PAIR_DIR='paired'
TMP_DIR='tmp'

FILTER='/mnt/home5/zoology/rm786/beetle_HiC/map/mapping_pipeline/filter_five_end.pl'
COMBINER='/mnt/home5/zoology/rm786/beetle_HiC/map/mapping_pipeline/two_read_bam_combiner.pl'
STATS='/mnt/home5/zoology/rm786/beetle_HiC/map/mapping_pipeline/get_stats.pl'
PICARD='/mnt/home5/zoology/rm786/software/picard.jar'

MAPQ_FILTER=10
CPU=26

echo "### Step 0: Check output directories' existence & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR

#echo "### Make Index"
#bwa index -p purged.pilon.fa $REF

echo "### Step 1.A: FASTQ to BAM (1st)"
bwa mem -t $CPU $INDEX $IN_DIR/${SRA}1.fq.gz | samtools view -@ $CPU -Sb - > $RAW_DIR/${SRA}_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
bwa mem -t $CPU $INDEX $IN_DIR/${SRA}2.fq.gz | samtools view -@ $CPU -Sb - > $RAW_DIR/${SRA}_2.bam

echo "### Step 2.A: Filter 5' end (1st)"
samtools view -h $RAW_DIR/${SRA}_1.bam | perl $FILTER | samtools view -Sb - > $FILT_DIR/${SRA}_1.bam

echo "### Step 2.B: Filter 5' end (2nd)"
samtools view -h $RAW_DIR/${SRA}_2.bam | perl $FILTER | samtools view -Sb - > $FILT_DIR/${SRA}_2.bam

echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER $FILT_DIR/${SRA}_1.bam $FILT_DIR/${SRA}_2.bam samtools $MAPQ_FILTER | samtools view -bS -t $FAIDX - | samtools sort -@ $CPU -o $TMP_DIR/$SRA.bam -

#********** The command line looks like this in the new syntax:
#**********    AddOrReplaceReadGroups -INPUT tmp/R.bam -OUTPUT paired/HiC.purged.polished.bam -ID R -LB R -SM nves_hic -PL ILLUMINA -PU none

echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=$TMP_DIR/$SRA.bam OUTPUT=$PAIR_DIR/$OUT.bam ID=$SRA LB=$SRA SM=$LABEL PL=ILLUMINA PU=none

## scaffold with yahs
~/software/yahs/yahs $REF HiC.purged.polished.bam -e GATC,GANTC,CTNAG,TTAA -o Ni2.purged.polished

# then https://github.com/vgl-hub/gfastats
~/software/gfastats/build/bin/gfastats input.fa --nstar-report --stats
