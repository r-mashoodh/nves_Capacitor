# README: Bulk chromatin accessibility (ATAC-seq) and H3K4me3 enrichment (CUT&Tag) in larval heads of the burying beetle *N. vespilloides*

This folder contains the code for the epigenomic data analysis presented in the manuscript.

## Key Information

* **Scripts**: The analysis is broken down into **several Python/R scripts, Bash commands, and a Snakemake pipeline**. Please see details in the "Code Details" section below.
* **Reproducibility**: Large raw/processed input data required to run these scripts are hosted in our associated **Zenodo repository** and/or NCBI/GEO repository.
* **Environment**: The exact computational environments are detailed in the `envs/` directory (Conda .yaml files).

## Code Details

The code is supplied in the `regulatory_elements_code/` directory:

* **Differential accessibility and binding analyses**: Found in `differential_access_activity/` (R scripts and Genrich commands).
* **Differential TF binding**: Found in `differential_access_activity/` (Python script and TOBIAS commands).
* **Promoters and enhancers definition & differential elements heatmap**: Executed via the provided Snakemake pipeline (`promoters_enhancers_final_plots.smk`).
* **GREAT region-to-gene association tests**: Found in `great_association_tests/` (Python scripts).
