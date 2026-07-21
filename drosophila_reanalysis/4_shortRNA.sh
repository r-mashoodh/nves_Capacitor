

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
