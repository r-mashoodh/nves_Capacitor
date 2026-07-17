#!/bin/bash

## use liftoff to annotate gene models from original Vespilloides genome 

# https://github.com/agshumate/LiftoffTools
# conda install -c bioconda liftofftools

conda create --prefix /home/rm786/rds/rds-rmash-sHE3M4BjkMM/conda_envs/liftofftools -c bioconda liftofftools=0.4.4

#!/bin/bash

conda activate /home/rm786/rds/rds-rmash-sHE3M4BjkMM/conda_envs/liftofftools

GTF='/home/rm786/rds/rds-rmash-sHE3M4BjkMM/nves_genome/ref_Nicve_v1.0_top_level.gff3'
NEWGTF='/home/rm786/rds/rds-rmash-sHE3M4BjkMM/liftoff_assembly/Ni2.v2.liftoff.gff_polished'
REF_FA='/home/rm786/rds/rds-rmash-sHE3M4BjkMM/nves_genome/GCF_001412225.1_Nicve_v1.0_genomic.fna'
NEW_FA='/home/rm786/rds/rds-rmash-sHE3M4BjkMM/liftoff_assembly/Cunningham_new_asm/Ni2.purged.folded.filtered_scaffolds_final.fa'

liftofftools all -r $REF_FA -t $NEW_FA -rg $GTF -tg $NEWGTF
