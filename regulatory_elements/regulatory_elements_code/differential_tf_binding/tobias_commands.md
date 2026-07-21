## 1. Merge replicates and down-sample to 60M reads for each condition
samtools merge --threads 16 -o all_merged_NC_mapq13.bam out_beetle_atac/filtered_bams/A5_NC_mapq13.bam out_beetle_atac/filtered_bams/A6_NC_mapq13.bam out_beetle_atac/filtered_bams/A4_NC_mapq13.bam out_beetle_atac/filtered_bams/A3_NC_mapq13.bam

samtools index all_merged_NC_mapq13.bam

samtools idxstats  all_merged_NC_mapq13.bam | cut -f 1 | grep -v "chrMT" | xargs samtools view -b all_merged_NC_mapq13.bam > all_merged_NC_mapq13_nochrMT.bam                                                                     

samtools view -c -F 260 all_merged_NC_mapq13_nochrMT.bam   

samtools view -b -s 0.966 all_merged_NC_mapq13_nochrMT.bam > all_merged_NC_mapq13_nochrMT_60M.bam 
samtools view -b -s 0.495 all_merged_FC_mapq13_nochrMT.bam > all_merged_FC_mapq13_nochrMT_60M.bam 


## 2. Run footprinting with TOBIAS

conda activate tobias
conda install -c conda-forge pandas=2.2

TOBIAS ATACorrect --bam all_merged_NC_mapq13_nochrMT_60M.bam --genome Ni2.purged.folded.filtered_scaffolds_final.sm.renamed.nochrMT.fa --peaks larvae.noControl.peaks.rahia.genrich.q.05.mapq13.bed  --outdir ATACorrect_FC --cores 8 

TOBIAS ATACorrect --bam all_merged_FC_mapq13_nochrMT_60M.bam --genome Ni2.purged.folded.filtered_scaffolds_final.sm.renamed.nochrMT.fa --peaks larvae.noControl.peaks.rahia.genrich.q.05.mapq13.bed  --outdir ATACorrect_FC --cores 8 

TOBIAS FootprintScores --signal ATACorrect_NC/all_merged_NC_mapq13_nochrMT_60M_corrected.bw --regions larvae.noCo
ntrol.peaks.rahia.genrich.q.05.mapq13.bed --output ATACorrect_NC/all_merged_NC_mapq13_nochrMT_60M_footprints.bw --cores 
8 

TOBIAS FootprintScores --signal ATACorrect_FC/all_merged_FC_mapq13_nochrMT_60M_corrected.bw --regions larvae.noCo
ntrol.peaks.rahia.genrich.q.05.mapq13.bed --output ATACorrect_FC/all_merged_FC_mapq13_nochrMT_60M_footprints.bw --cores 
8 

Here we curated the JASPAR core insect motifs to build `JASPAR2026_CORE_insects_non-redundant_pfms_jasparCURATED.txt`, retaining motifs with an ortholog in *N. vespilloides*. For details on motifs to gene association, see `Ni2.JASPAR2026_insects.motif2factors.1motifpergene.txt`.

TOBIAS BINDetect --motifs JASPAR2026_CORE_insects_non-redundant_pfms_jasparCURATED.txt --signals ATACorrect_FC/all_merged_FC_mapq13_nochrMT_60M_footprints.bw ATACorrect_NC/all_merged_NC_mapq13_nochrMT_60M_footprints.bw --genome Ni2.purged.folded.filtered_scaffolds_final.sm.renamed.nochrMT.fa --peaks larvae.noControl.peaks.rahia.genrich.q.05.mapq13.noMT.bed --outdir BINDetect_output_v3 --cond_names FC NC --cores 8