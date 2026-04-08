## Assembly HiFi reads with hifiasm

## followed installation instructions here
## https://github.com/chhylp123/hifiasm

# conda create --prefix ~/home/rm786/rds/rds-rmash-sHE3M4BjkMM/my_conda_env
# conda install -c bioconda pbtk

### convert bam to fastq
# PacBio's bam2fastq
module load miniconda/3
#conda env create -f picrust2-env.yaml
eval "$(conda shell.bash hook)"
conda activate ~/home/rm786/rds/rds-rmash-sHE3M4BjkMM/my_conda_env
# generates Ni2.hifi_reads.fastq.gz
bam2fastq -o Ni2.hifi_reads Ni2.hifi_reads.bam


## run hifiasm
hifiasm -o NicV2/NicV2.asm -t36 --h1 HiC/R1.fq.gz --h2 HiC/R2.fq.gz PacBio/Ni2.hifi_reads.fastq.gz

#The FASTA file can be produced from GFA as follows:
awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa

# https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output

## HiFi reads only
## `prefix`.bp.p_ctg.gfa: assembly graph of primary contigs.
awk '/^S/{print ">"$2;print $3}' Ni2_pacbio.asm.bp.p_ctg.gfa > Ni2_pacbio.asm.fa

# get some assembly stats
python ~/software/fast_stats/fast_stats.py -n 50 -i Ni2_pacbio.asm.fa
~/software/quast/quast.py Ni2_pacbio.asm.fa --conserved-genes-finding
