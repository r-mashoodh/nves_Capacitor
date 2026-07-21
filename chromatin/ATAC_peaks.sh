# path to genrich executable
Genrich="/home/rm786/rds/rds-rmash-sHE3M4BjkMM/software/Genrich/Genrich"

# path to folders
PEAKS="/home/rm786/rds/rds-rmash-sHE3M4BjkMM/GRN_prediction/Genrich_peaks"
BAMS="/home/rm786/rds/rds-rmash-sHE3M4BjkMM/GRN_prediction/ATAC_data"

module load bedtools/2.20.1 

# Step 1: Name-sort all BAMs and index
for bam in "$BAMS"/*.sorted.renamed.bam; do
    echo "Sorting $bam by read name..."
    samtools sort -n -@ 8 -o "${bam%.bam}.n.bam" "$bam"
    
    echo "Indexing ${bam%.bam}.n.bam..."
    samtools index "${bam%.bam}.n.bam"
done

# Step 2: Run Genrich
#.sorted.renamed.n.bam
$Genrich -t "$BAMS/A3_FC.sorted.renamed.n.bam $BAMS/A4_FC.sorted.renamed.n.bam $BAMS/A5_FC.sorted.renamed.n.bam $BAMS/A6_FC.sorted.renamed.n.bam $BAMS/A3_NC.sorted.renamed.n.bam $BAMS/A4_NC.sorted.renamed.n.bam $BAMS/A5_NC.sorted.renamed.n.bam $BAMS/A6_NC.sorted.renamed.n.bam" \
-o $PEAKS/larvae.mapq13.narrowPeak -m 13 -q 0.05 -v -j -r -y -f $PEAKS/larvae.mapq13.outfile.noControl.log


# Step 3: Make BED file
cut -f1-6 "$PEAKS/larvae.mapq13.narrowPeak" > "$PEAKS/larvae.mapq13.bed"

# Step 4: Count reads in peaks
bedtools multicov \
-bams \
"$BAMS/A3_FC.sorted.renamed.mapq13.bam" \
"$BAMS/A4_FC.sorted.renamed.mapq13.bam" \
"$BAMS/A5_FC.sorted.renamed.mapq13.bam" \
"$BAMS/A6_FC.sorted.renamed.mapq13.bam" \
"$BAMS/A3_NC.sorted.renamed.mapq13.bam" \
"$BAMS/A4_NC.sorted.renamed.mapq13.bam" \
"$BAMS/A5_NC.sorted.renamed.mapq13.bam" \
"$BAMS/A6_NC.sorted.renamed.mapq13.bam" \
-bed "$PEAKS/larvae.mapq13.bed" \
> "$PEAKS/larvae.mapq13.counts"
