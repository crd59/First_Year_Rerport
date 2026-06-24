#!/bin/bash

source /etc/profile.d/modules.sh

# Load R and samtools modules FIRST to prevent samtool environment overriding FLAMES 
module load R
module load samtools/1.14/gcc/v46lwk2d

export MAMBA_EXE='/rds/user/crd59/hpc-work/software_environments/bin/micromamba'
export MAMBA_ROOT_PREFIX='/rds/user/crd59/hpc-work/software_environments/micromamba'
eval "$($MAMBA_EXE shell hook --shell bash)"
micromamba activate FLAMES

python /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/FLAMES/FLAMES/FLAMES.py annotate -a /rds/user/crd59/hpc-work/FLAMES_example_code/Annotation_data/ -l /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/Cortical/LD_blocks_ordered_try.txt -sc FLAME_NAME -pc SNP.PP -id /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/Cortical/master_abf_result_file_paths.txt -p /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/Cortical/naive_path_Cortical.preds -m /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/Cortical/Cortical.genes.out -mt /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/Cortical/Cortical.gsa.out -o /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/Cortical