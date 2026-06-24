#!/bin/bash

source /etc/profile.d/modules.sh

# Load R and samtools modules FIRST to prevent samtool environment overriding FLAMES 
module load R
module load samtools/1.14/gcc/v46lwk2d

export MAMBA_EXE='/rds/user/crd59/hpc-work/software_environments/bin/micromamba'
export MAMBA_ROOT_PREFIX='/rds/user/crd59/hpc-work/software_environments/micromamba'
eval "$($MAMBA_EXE shell hook --shell bash)"
micromamba activate FLAMES

python /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/FLAMES/FLAMES/FLAMES.py FLAMES -id /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/Cortical/master_abf_result_file_paths.txt -o /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/Cortical
python /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/FLAMES/FLAMES/FLAMES.py FLAMES -id /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/LimbHeight/master_abf_result_file_paths.txt -o /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/LimbHeight
python /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/FLAMES/FLAMES/FLAMES.py FLAMES -id /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/TrunkWidths/master_abf_result_file_paths.txt -o /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/TrunkWidths