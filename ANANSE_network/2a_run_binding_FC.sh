## remade the motifs2factors file using New_motifs2factors.R = Ni2.JASPAR2026_insects.pfm
## Ni2.JASPAR2026_insects.1motifpergene.pfm = trimmed list

module load miniconda/3
eval "$(conda shell.bash hook)"
conda activate /home/rm786/rds/rds-rmash-sHE3M4BjkMM/conda_envs/ananse

## using new bed file with less stringent q
## should I use a different bed file for each condition?
## manual db is assigned based on gimmemotifs, then manually curated using orthoFinder and rbh hits

# --- Run ANANSE binding ---
## note this is with mapq13 peaks
## manually made pfm file
ananse binding \
  -A ATAC_data/*FC.sorted.renamed.mapq13.bam \
  -r Genrich_peaks/larvae.mapq13.narrowPeak \
  -p manual_db/Ni2.JASPAR2026_insects.1motifpergene.pfm \
  -g Ni2 \
  -o FC.binding2
