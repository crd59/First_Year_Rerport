#!/bin/bash

# Source module environment
source /etc/profile.d/modules.sh

# Set up micromamba exactly as in your .bashrc
export MAMBA_EXE='/rds/user/crd59/hpc-work/software_environments/bin/micromamba'
export MAMBA_ROOT_PREFIX='/rds/user/crd59/hpc-work/software_environments/micromamba'
eval "$($MAMBA_EXE shell hook --shell bash)"

# Activate your environment
micromamba activate FLAMES

python /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/pops/pops.py --num_feature_chunks 99 --magma_prefix /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/Cortical/Cortical --gene_annot_path /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/FLAMES/FLAMES_annotation_data/pops_features_pathway_naive/gene_annot.txt --feature_mat_prefix /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/FLAMES/FLAMES_annotation_data/pops_features_pathway_naive/munged_features/pops_features --control_features /rds/user/crd59/hpc-work/control.features --out_prefix /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/Cortical/naive_path_Cortical --verbose

python /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/pops/pops.py --num_feature_chunks 99 --magma_prefix /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/Cortical/Cortical --gene_annot_path /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/FLAMES/FLAMES_annotation_data/pops_features_pathway_naive/gene_annot.txt --feature_mat_prefix /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/FLAMES/FLAMES_annotation_data/pops_features_pathway_naive/munged_features/pops_features --control_features /rds/user/crd59/hpc-work/control.features --out_prefix /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/Cortical/naive_path_LimbHeight --verbose

python /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/pops/pops.py --num_feature_chunks 99 --magma_prefix /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/TrunkWidths/TrunkWidths --gene_annot_path /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/FLAMES/FLAMES_annotation_data/pops_features_pathway_naive/gene_annot.txt --feature_mat_prefix /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/FLAMES/FLAMES_annotation_data/pops_features_pathway_naive/munged_features/pops_features --control_features /rds/user/crd59/hpc-work/control.features --out_prefix /rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/TrunkWidths/naive_path_TrunkWidths --verbose
