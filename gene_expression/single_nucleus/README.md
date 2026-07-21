# README: Single-nucleus gene expression in larval heads of the burying beetle *N. vespilloides*

This folder contains the complete analysis pipeline for the single-nucleus RNA-sequencing (snRNA-seq) data presented in the manuscript.

## Key Information

* **Notebooks**: The analysis is broken down into **three numbered Jupyter notebooks** running R. Steps executed outside of R (e.g., Python scripts or Bash commands) are provided and explicitly referenced within the notebooks.
* **Reproducibility**: Large raw/processed input data required to run these notebooks are hosted in our associated **Zenodo repository** and/or NCBI/GEO repository.
* **Environment**: The exact computational environments are detailed in the `envs/` directory (Conda .yaml files).

## Important Note on Gene Names

Gene nomenclature in the Seurat object differs slightly from the final names used in the manuscript and figures. 

To map between them, please refer to the master table: `../../genome_assembly/genome_files/genes/gene_models_allinfos_for_snRNAseq.tsv`.
* **`OLD_gene_name`**: Name stored inside the Seurat object.
* **`gene name`**: Final name used in the paper.

The GTF file used for CellRanger can be found in `../../genome_assembly/genome_files/genes/genes/cellranger_gtf/Nvi2.liftoff.mitohifi.final.ok.genename.cellranger.gtf`.
