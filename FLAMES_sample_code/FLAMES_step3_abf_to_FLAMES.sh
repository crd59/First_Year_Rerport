#!/bin/bash
#SBATCH -A WARRIER-SL2-CPU
#SBATCH -a 1-2332%100
#SBATCH -J YG_GWAS_abf_to_FLAMES
#SBATCH -D /rds/user/CAMID/hpc-work/pgs_img_2/SSHscript/GWAS_FLAME/FLAMES_pipeline_v2
#SBATCH -o /rds/user/CAMID/hpc-work/pgs_img_2/SSHscript/GWAS_FLAME/FLAMES_pipeline_v2/step3_abf_to_FLAMES/abf_to_FLAMES_v2_%A_%a.log
#SBATCH -p sapphire                             # on sapphire partition
#SBATCH -t 30:00                                # HH:MM:SS with maximum 12:00:00 for SL3 or 36:00:00 for SL2
#SBATCH --mem=6G
#SBATCH --mail-type=ALL

# Source module environment
source /etc/profile.d/modules.sh

# Load R and samtools modules FIRST to prevent samtool environment overriding FLAMES 
module load R
module load samtools/1.14/gcc/v46lwk2d

# Set up micromamba AFTER loading modules
export MAMBA_EXE='/rds/user/CAMID/hpc-work/software_environments/bin/micromamba'
export MAMBA_ROOT_PREFIX='/rds/user/CAMID/hpc-work/software_environments/micromamba'
eval "$($MAMBA_EXE shell hook --shell bash)"

# Activate FLAMES environment
micromamba activate FLAMES

# Run the R script
Rscript FLAMES_step3_abf_to_FLAMES.R