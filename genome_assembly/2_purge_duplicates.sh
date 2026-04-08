#!/bin/bash

## purging duplicates

cd ~/software/purge_dups/

python3 ./scripts/pd_config.py -l Ni2 -n config.Ni2.json /mnt/home5/zoology/rm786/beetle_HiC/noHiC/Ni2_pacbio.asm.fa pb.fofn

python3 ./scripts/run_purge_dups.py config.Ni2.json src NVes2 -p bash

## get coverage plot
python3 ../../scripts/hist_plot.py -c cutoffs PB.stat PB.cov.png

## run busco
singularity exec --bind $(pwd):/busco_wd ~/software/busco.sif busco -i ~/software/purge_dups/Ni2_pacbio.asm/seqs/Ni2_pacbio.asm.purged.fa -m genome --cpu 16 -l arthropoda_odb10 -f
