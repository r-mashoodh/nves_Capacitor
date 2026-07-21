module load miniconda/3
eval "$(conda shell.bash hook)"
conda activate /home/rm786/rds/rds-rmash-sHE3M4BjkMM/conda_envs/ananse

# --- Run ANANSE network ---

ananse network FC.binding2/binding.h5 \
  -e RNAseq_data/FC_FC/*.quant.genes.sf \
  -g Ni2 \
  -n 4 \
  -o results2/FC_network.tsv

ananse network NC.binding2/binding.h5 \
  -e RNAseq_data/FC_NC/*.quant.genes.sf \
  -g Ni2 \
  -n 4 \
  -o results2/NC_network.tsv
