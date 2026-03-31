# remap and quantify rna-seq samples with new genome

# load required modules
module load java/21.0.4

# nextflow command
nextflow run nf-core/rnaseq -profile singularity \
-w larvae_work --input samplesheet.csv --outdir rnaseq_results \
--fasta new_hic_genome/Ni2.purged.folded.filtered_scaffolds_final.sm.renamed.fa \
--gtf new_hic_genome/Nvi2.liftoff.mitohifi.final.ok.genename.for.cellranger.gtf \
--remove_ribo_rna -r 3.21.0 -c myriad.config -bg 
