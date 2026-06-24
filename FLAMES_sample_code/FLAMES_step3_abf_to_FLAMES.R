# Step 3a. Run Finemapping abf 
library(tidyverse)
library(coloc)
library(data.table, include.only = "fread")

setwd("~/Desktop/PhD /First_project/FLAMES_2.0/")

# Create dictionary with names and N 

mat <- matrix(ncol = 2, nrow =3)
mat <- as.data.frame(mat)
colnames(mat) <- c("GWAS", "N")
mat[,1] <-  c("Cortical", "LimbHeight", "TrunkWidths")
mat[,2] <- c(18037, 39556, 54881)

LD_blocks_new <- read.table("./../SuSiE_coloc/all_loci.tsv", sep = "\t", header = T)

for(GWAS in 1:3){
  
  full_GWAS <- read.table(paste0("./../formatted_GWAS/", mat[GWAS,1], "_GWAS.tsv.gz"), sep = "\t", header = T)
  
  current_folder = mat[GWAS,1]
  dir.create(current_folder)
  
  # set abf finemapping parameters
  
  munged_loci <- LD_blocks_new[grepl("trunk",LD_blocks_new$trait),]
  
  # !! create output folder if folder did not exist yet !!
  current_abf_folder = paste0(current_folder, "/", "abf/")
  
  if (!dir.exists(current_abf_folder)) {
    dir.create(current_abf_folder)
  }
  
  for (locus_trait in 1:nrow(munged_loci)) {

    start_pos <- munged_loci$start[locus_trait]
    end_pos <- munged_loci$end[locus_trait]
    chr <- munged_loci$chr[locus_trait]
    
    snps_in_locus <- full_GWAS[full_GWAS$BP >= start_pos &
                                 full_GWAS$BP <= end_pos &
                                 full_GWAS$CHR == chr,]
    
    snps_in_locus <- snps_in_locus[!duplicated(snps_in_locus$SNP),]
    
    current_beta = snps_in_locus$BETA
    current_varbeta = snps_in_locus$SE^2
    current_snp = snps_in_locus$SNP
    current_maf = snps_in_locus$MAF
    current_N = mat[GWAS,2]
    
    current_abf.df = list(beta = current_beta, varbeta = current_varbeta, N = current_N, 
                          snp = current_snp, type = "quant", MAF = current_maf)
    current_abf.result = finemap.abf(dataset=current_abf.df)
    
    # generate 95% credible set as requested by FLAMES
    credible_set = current_abf.result %>%
      arrange(desc(SNP.PP)) %>%
      mutate(cumsum_pp = cumsum(SNP.PP))
    credible_set = credible_set[1:min(which(credible_set$cumsum_pp >= 0.95)), ]
    
    credible_set_FLAMES = snps_in_locus %>% filter(SNP %in% credible_set$snp) %>% mutate(FLAME_NAME = paste(CHR, BP, A1, A2, sep = ":"))
    credible_set_FLAMES = full_join(credible_set, credible_set_FLAMES, join_by(snp == SNP)) %>% select(FLAME_NAME, SNP.PP)
    
    # write results
    credible_set_path = paste0(current_abf_folder, "/", munged_loci$locus[locus_trait], "_LD_block_", mat[GWAS,1], ".txt")
    write.table(credible_set_FLAMES, credible_set_path, row.names = FALSE, sep="\t", quote = FALSE) # FLAMES requested tab-separated columns
  }
  
  # Step 3b. generate master finemapping result path file for FLAMES annotate -id
  abf_results = list.files(current_abf_folder, full.names = TRUE)
  LD_id = sub("LD_(\\d+)_.*", "\\1", basename(abf_results))
  
  current_FLAMES_annotate_folder = paste0(current_folder, "/", "FLAMES_annotate/")
  
  if (!dir.exists(current_FLAMES_annotate_folder)) {
    dir.create(current_FLAMES_annotate_folder)
  }
  
  master.cs95.path.df = data.frame(Filename = abf_results,
                                   Annotfiles = paste0(current_FLAMES_annotate_folder, "FLAMES_annotated_", basename(abf_results), ".txt"))
  master.cs95.path.df$GenomicLocus <- gsub(x =master.cs95.path.df$Filename, pattern = paste0(mat[GWAS,1], "/abf//"), replacement= "")
  master.cs95.path.df$GenomicLocus <- gsub(x =master.cs95.path.df$GenomicLocus, pattern = paste0("_LD_block_", mat[GWAS,1],".txt"), replacement=  "")
  master.cs95.path.df$GenomicLocus <- as.integer(master.cs95.path.df$GenomicLocus)
  master.cs95.path.df <- master.cs95.path.df[order(master.cs95.path.df$GenomicLocus),]
  master.cs95.path.df <- master.cs95.path.df[,-3]
  
  master.cs95.path.name = paste0(current_folder, "/", "master_abf_result_file_paths.txt")
  write.table(master.cs95.path.df, master.cs95.path.name, row.names = FALSE, quote = FALSE, sep = "\t")
  
}

# Step 3c. Construct FLAMES annotate command
GWAS = 3
current_folder = mat[GWAS,1]
master.cs95.path.name = paste0(current_folder, "/", "master_abf_result_file_paths.txt")

## Sort out LD blcoks

LD_blocks <- read.table("./../SuSiE_coloc/all_loci.tsv", sep = "\t", header = T)
LD_blocks <- LD_blocks[grepl("trunk", LD_blocks$trait),]
LD_blocks <- LD_blocks[,c(2:5)]
colnames(LD_blocks)[1] <- "GenomicLocus"

LD_blocks_file <- paste0(mat[GWAS,1],"/LD_blocks.txt")
write.table(LD_blocks, paste0(mat[GWAS,1],"/LD_blocks.txt"), sep = "\t", row.names = F, col.names = T, quote = F)

## Sort out master file blocks

old_folder <- paste0("~/Desktop/PhD/First_project/FLAMES/FLAMES_Multivariate_GWAS/", mat[GWAS,1], "/",  mat[GWAS,1])

preds_file <- paste0("~/Desktop/PhD/First_project/FLAMES/FLAMES_Multivariate_GWAS/", mat[GWAS,1], "/naive_path_",  mat[GWAS,1],".preds")

  # Construct the Bash command
  command <- paste("python ~/Desktop/PhD/First_project/FLAMES/FLAMES/FLAMES.py annotate", 
                   "-a ~/Desktop/PhD/First_project/FLAMES/Annotation_data/",
                   "-l", LD_blocks_file, 
                   "-b GRCh37",
                   "-sc FLAME_NAME",
                   "-pc SNP.PP",
                   "-id", master.cs95.path.name,
                   "-p", preds_file,
                   "-m", paste0(old_folder, ".genes.out"),
                   "-mt", paste0(old_folder, ".gsa.out"),
                   "-o", current_folder)

  
  # FLAMES Result
  
  command <- paste("python ~/Desktop/PhD/First_project/FLAMES/FLAMES/FLAMES.py FLAMES", 
                   "-id", master.cs95.path.name,
                   "-o", current_folder)
  
  # Execute the Bash command
  system(command)
  
## Make sure there are no NA values
  GWAS = 3
  current_folder <- mat[GWAS,1]
  master.cs95.path.name = paste0(current_folder, "/", "master_abf_result_file_paths.txt")
  
  master <- read.table(master.cs95.path.name, header = T)
  
for(i in master$Filename){
  thing <- read.table(i, header = T)
  thing <- thing[complete.cases(thing),]
  write.table(thing, i, sep = "\t", col.names = T, row.names = F, quote = F)
}
  