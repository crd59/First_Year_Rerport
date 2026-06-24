
setwd("/rds/user/crd59/hpc-work/")

# Create dictionary with names and N 

mat <- matrix(ncol = 2, nrow =3)
mat <- as.data.frame(mat)
colnames(mat) <- c("GWAS", "N")
mat[,1] <-  c("Cortical", "LimbHeight", "TrunkWidths")
mat[,2] <- c(18037, 39556, 54881)

for(GWAS in 1:3){
  
  current_folder = paste0("/rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/", mat[GWAS,1])

  current_magma_output = paste(current_folder, mat[GWAS,1], sep = "/")
  
# Construct the Bash command: Naive_features (!! --num_feature_chunks 99 !!)
naive_path_pops_feature_folder = "/rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/FLAMES/FLAMES_annotation_data/pops_features_pathway_naive/"

naive_path_command <- paste("python /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/pops/pops.py", 
                            "--num_feature_chunks 99",
                            "--magma_prefix", current_magma_output,
                            "--gene_annot_path", paste0(naive_path_pops_feature_folder, "gene_annot.txt"),
                            "--feature_mat_prefix", paste0(naive_path_pops_feature_folder, "munged_features/pops_features"),
                            "--control_features", "./control.features",
                            "--out_prefix", paste0(current_folder, "/", "naive_path_",  mat[GWAS,1]),
                            "--verbose")

# Execute the Bash command
system(naive_path_command)

}