#!/bin/bash
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH -o map_to_scaffold.out

# generate a wgs bam to check coverage of assembly using blobtoolkit
#bwa index -p purged.polished.scaffold ../scaffold/Ni2.purged.polished_scaffolds_final.fa
READ1='/mnt/home5/zoology/rm786/software/blobtoolkit/coverage/WA_HLNF3DSX2_L4_1_val_1.fq.gz'
READ2='/mnt/home5/zoology/rm786/software/blobtoolkit/coverage/WA_HLNF3DSX2_L4_2_val_2.fq.gz'
INDEX='purged.polished.scaffold'

bwa mem -t 32 $INDEX $READ1 $READ2 | samtools view -@30 -bS - > WA_purged.polished.scaffold.bam

samtools sort -@30 WA_purged.polished.scaffold.bam > WA_purged.polished.scaffold.sorted.bam
samtools index WA_purged.polished.scaffold.sorted.bam

### nextflow (myriad) -- BLOBTOOLKIT
### samplesheet_purged.csv
### WA_purged.polished.scaffold.bam

#sample,datatype,datafile
#sample1,hic,hic.cram
#sample2,illumina,illumina.cram
#sample2,illumina,illumina.cram

nextflow run sanger-tol/blobtoolkit -r main -resume -c myriad.config -bg -profile singularity --input samplesheet_purged.csv \
--outdir Ni2_purged.polished.scaffolded \
--fasta /home/ucbtrm5/Scratch/blobtoolkit/Ni2.purged.polished_scaffolds_final.fa \
--taxon 110193 --taxdump /home/ucbtrm5/Scratch/blobtoolkit/databases/taxdump \
--blastn /lustre/projects/BLASTdb/ \
--blastx /home/ucbtrm5/Scratch/blobtoolkit/databases/uniprot/reference_proteomes.dmnd \
--blastp /home/ucbtrm5/Scratch/blobtoolkit/databases/uniprot/reference_proteomes.dmnd 

## remove contaminants

## this gets only contigs with arthropod/no-hit contigs
## ie removes fungal
singularity exec docker://genomehubs/blobtoolkit:latest blobtools filter --param buscogenes_phylum--Keys=Ascomycota --fasta Ni2.purged.polished_scaffolds_final.fa --output filtered_funghi /home/ucbtrm5/Scratch/blobtoolkit/Ni2_purged.polished.scaffolded/blobtoolkit/Ni2.purged.polished_scaffolds_final/

## filtered now contains a folder where you can remake the plots 
singularity exec ../.singularity/pull/blobtoolkit.sf blobtools view --plot --out new_plots --view snail /home/ucbtrm5/Scratch/blobtoolkit/filtered_funghi
