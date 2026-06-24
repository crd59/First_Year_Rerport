# Step 1a. Magma

setwd("/rds/user/crd59/hpc-work/")

# Create dictionary with names and N 

mat <- matrix(ncol = 2, nrow =3)
mat <- as.data.frame(mat)
colnames(mat) <- c("GWAS", "N")
mat[,1] <-  c("Cortical", "LimbHeight", "TrunkWidths")
mat[,2] <- c(18037, 39556, 54881)

for(GWAS in 1:3){

current_folder = paste0("/rds/user/crd59/hpc-work/FLAMES_Multivariate_GWAS/", mat[GWAS,1])

# !! create output folder if folder did not exist yet !!
if (!dir.exists(current_folder)) {
  dir.create(current_folder)
}

current_magma_output = paste(current_folder, mat[GWAS,1], sep = "/")

# Construct the Bash command
command <- paste("/rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/magma/magma", 
                 "--bfile /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/magma/g1000_eur", 
                 "--gene-annot /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/FLAMES/16_Feb_2021_magma_0kb.genes.annot",
                 "--pval", paste0("./formatted_GWAS/", mat[GWAS,1], "_GWAS.tsv"), paste0("N =", mat[GWAS,2]),
                 "--gene-model snp-wise=mean",
                 "--out", current_magma_output)

# Execute the Bash command
system(command)

# Step 1b. Magma Tissue
current_magma_tissue_input = paste0(current_magma_output, ".genes.raw")

# Construct the Bash command
command <- paste("/rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/magma/magma", 
                 "--gene-results", current_magma_tissue_input,
                 "--gene-covar /rds/user/crd59/rds-genetics_hpc-Nl99R8pHODQ/toolbox/FLAMES/FLAMES_annotation_data/gtex_v8_ts_avg_log2TPM.txt", 
                 "--out", current_magma_output)

# Execute the Bash command
system(command)

}
