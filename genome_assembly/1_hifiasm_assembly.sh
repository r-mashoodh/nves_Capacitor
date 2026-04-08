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
hifiasm -o NicV2/NicV2.asm -t36 PacBio/Ni2.hifi_reads.fastq.gz


