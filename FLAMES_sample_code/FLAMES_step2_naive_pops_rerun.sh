#!/bin/bash
#SBATCH -A WARRIER-SL2-CPU
#SBATCH -a CHANGE_TO_ARRAYID_AS_NEEDED (Ex. 3,1,22,68)
#SBATCH -J YG_GWAS_PoPS
#SBATCH -D /rds/user/CAMID/hpc-work/pgs_img_2/SSHscript/GWAS_FLAME/FLAMES_pipeline_v2
#SBATCH -o /rds/user/CAMID/hpc-work/pgs_img_2/SSHscript/GWAS_FLAME/FLAMES_pipeline_v2/step2_naive_pops/naive_pops_rerun_v2_%A_%a.log
#SBATCH -p sapphire                             # on sapphire partition
#SBATCH -t 1:00:00                                # HH:MM:SS with maximum 12:00:00 for SL3 or 36:00:00 for SL2
#SBATCH --mem=32G
#SBATCH --mail-type=ALL

# Source module environment
source /etc/profile.d/modules.sh

# Set up micromamba exactly as in your .bashrc
export MAMBA_EXE='/rds/user/CAMID/hpc-work/software_environments/bin/micromamba'
export MAMBA_ROOT_PREFIX='/rds/user/CAMID/hpc-work/software_environments/micromamba'
eval "$($MAMBA_EXE shell hook --shell bash)"

# Activate your environment
micromamba activate FLAMES

# Load R module
module load R

# Run the R script
Rscript FLAMES_step2_naive_pops.R