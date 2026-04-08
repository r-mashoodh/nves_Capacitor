## Scaffolding with YAHS

### YAHS
~/software/yahs/yahs ../Ni2_pacbio.asm.fa ../../map/HiC_paired.bam

## mapping stats
python ~/software/fast_stats/fast_stats.py -n 50 -i yahs.out_scaffolds_final.fa

## run blobtoolkit

# run BUSCO
singularity exec --bind $(pwd):/busco_wd ~/software/busco.sif busco -i ../yahs.out_scaffolds_final.fa -m genome -cpu 32 -l arthropoda_odb10 -f
