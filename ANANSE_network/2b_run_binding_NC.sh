## remade the motifs2factors file using New_motifs2factors.R = Ni2.JASPAR2026_insects.pfm
## Ni2.JASPAR2026_insects.1motifpergene.pfm = elises trimmed list


module load miniconda/3
eval "$(conda shell.bash hook)"
conda activate /home/rm786/rds/rds-rmash-sHE3M4BjkMM/conda_envs/ananse

# --- Run ANANSE binding ---

ananse binding \
  -A ATAC_data/*NC.sorted.renamed.mapq13.bam \
  -r Genrich_peaks/larvae.mapq13.narrowPeak \
  -p manual_db/Ni2.JASPAR2026_insects.1motifpergene.pfm \
  -g Ni2 \
  -o NC.binding2

