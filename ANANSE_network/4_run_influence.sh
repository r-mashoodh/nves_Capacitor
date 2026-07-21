module load miniconda/3
eval "$(conda shell.bash hook)"
conda activate /home/rm786/rds/rds-rmash-sHE3M4BjkMM/conda_envs/ananse

## source: FC; target NC ... means you are getting the NC expressed/dominant network
## for this reason if looking at FC network then deseq2 log2foldchange should flip signs
## current Deseq2: ++ values = higher in NC, so need to flip this to get what is higher in FC
## used awk in RNAseq_data folder
# awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} {$2 = -$2; print}' gene_name_deseq2.tsv > gene_name_deseq2_flipped.tsv

# --- Run ANANSE influence ---

ananse influence  -s results2/FC_network.tsv \
                  -t results2/NC_network.tsv \
                  -d RNAseq_data/gene_name_deseq2.tsv \
                  -o results2/nc_influence.tsv \
                  -i 500_000 \
                  -n 8

# --- Plot ANANSE influence ---

ananse plot results2/nc_influence.tsv \
            --diff-network results2/nc_influence_diffnetwork.tsv \
            -o results2/nc_plot


## we can get the FC network by flipping

# --- Run ANANSE influence ---

ananse influence  -s results2/NC_network.tsv \
                  -t results2/FC_network.tsv \
                  -d RNAseq_data/gene_name_deseq2_flipped.tsv \
                  -o results2/fc_influence.tsv \
                  -i 500_000 \
                  -n 8


# --- Plot ANANSE influence ---

ananse plot results2/fc_influence.tsv \
            --diff-network results2/fc_influence_diffnetwork.tsv \
            -o results2/fc_plot
