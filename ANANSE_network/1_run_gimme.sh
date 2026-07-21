module load miniconda/3
eval "$(conda shell.bash hook)"

## followed installation instructions here: https://github.com/vanheeringen-lab/ANANSE

# --- Environment setup ---
source ~/.bashrc

# Activate the ananse environment
conda activate /home/rm786/rds/rds-rmash-sHE3M4BjkMM/conda_envs/ananse

# Ensure your custom bin directory exists and contains the wrapper
mkdir -p /home/rm786/rds/rds-rmash-sHE3M4BjkMM/software/bin

# Create orthofinder wrapper if missing
if [ ! -f /home/rm786/rds/rds-rmash-sHE3M4BjkMM/software/bin/orthofinder ]; then
    cat << 'EOF' > /home/rm786/rds/rds-rmash-sHE3M4BjkMM/software/bin/orthofinder
#!/bin/bash
python /home/rm786/rds/rds-rmash-sHE3M4BjkMM/software/OrthoFinder_source/orthofinder.py "$@"
EOF
    chmod +x /home/rm786/rds/rds-rmash-sHE3M4BjkMM/software/bin/orthofinder
fi

# Add to PATH
export PATH=/home/rm786/rds/rds-rmash-sHE3M4BjkMM/software/bin:$PATH

# --- Run gimme motif2factors ---
cd /home/rm786/rds/rds-rmash-sHE3M4BjkMM/GRN_prediction

gimme motif2factors \
  --new-reference Ni2 \
  --database JASPAR2020_insects \
  --ortholog-references WBcel235 BDGP6.54 Tcas5.2 \
  --outdir ANANSE_insect_db2 \
  --tmpdir tmp \
  --threads 24 \
  --keep-intermediate
