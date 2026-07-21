# README: Genome Files

This folder contains the genome assembly and gene annotation files used for the *N. vespilloides* paper.

## Directory Structure
* **`genome/`**: Contains the genome sequence FASTA file (large file - hosted on zenodo).
* **`genes/`**: Contains predicted gene files (CDS, proteins) in GFF/FASTA formats.

## Reproducibility & Gene Names
For full reproducibility, we provide the specific GTF file used for running the CellRanger pipeline:
`genes/cellranger_gtf/Nvi2.liftoff.mitohifi.final.ok.genename.cellranger.gtf`

**Important Note on Nomenclature:** 
Gene names stored in this GTF and the resulting Seurat object differ from the final names used in the manuscript and figures (and other files located in 'genes'). To map and convert these identifiers, please use the gene master table in:
`genes/gene_models_allinfos_for_snRNAseq.tsv` (map `OLD_gene_name` to `gene name`).