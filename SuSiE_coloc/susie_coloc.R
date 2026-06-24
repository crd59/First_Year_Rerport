# Define functions

collapse_group <- function(group) {
  group_df <- do.call(rbind, group)
  
  # Row with smallest p drives the key SNP columns
  best_row <- group_df[which.min(group_df$p), ]
  
  data.frame(
    # From best p-value row
    uniqID        = best_row$uniqID,
    rsID          = best_row$rsID,
    chr           = best_row$chr,
    pos           = best_row$pos,
    p             = best_row$p,
    # Full merged coordinates
    start         = min(group_df$start),
    end           = max(group_df$end),
    # Summed columns
    nSNPs         = sum(group_df$nSNPs),
    nGWASSNPs     = sum(group_df$nGWASSNPs),
    nIndSigSNPs   = sum(group_df$nIndSigSNPs),
    nLeadSNPs     = sum(group_df$nLeadSNPs),
    # Pasted columns
    IndSigSNPs    = paste(group_df$IndSigSNPs, collapse = ";"),
    LeadSNPs      = paste(group_df$LeadSNPs,   collapse = ";"),
    stringsAsFactors = FALSE
  )
}

merge_overlapping_loci <- function(loci) {
  # Sort by chr and start position
  loci <- loci[order(loci$chr, loci$start), ]
  
  merged <- list()
  
  # Initialise the first group as a list of rows
  current_group <- list(loci[1, ])
  
  for (i in seq(2, nrow(loci))) {
    row <- loci[i, ]
    last <- current_group[[length(current_group)]]
    
    # If same chromosome and overlapping, add to current group
    if (row$chr == last$chr && row$start <= max(sapply(current_group, `[[`, "end"))) {
      current_group[[length(current_group) + 1]] <- row
    } else {
      # Collapse the current group and start a new one
      merged[[length(merged) + 1]] <- collapse_group(current_group)
      current_group <- list(row)
    }
  }
  # Collapse the final group
  merged[[length(merged) + 1]] <- collapse_group(current_group)
  
  return(do.call(rbind, merged))
}

# packages

library(susieR)
library(tidyverse)
library(coloc)
library(genetics.binaRies)
library(data.table)
library(Rfast)

# Genomic Risk Loci

setwd("~/Desktop/PhD/First_project/FUMA")
# 
 cortical_grl <- read.table("Cortical_FUMA/GenomicRiskLoci.txt", header = T)
 trunk_grl <- read.table("TrunkWidths_FUMA/GenomicRiskLoci.txt", header = T)
 limb_grl <-  read.table("LimbHeight_FUMA/GenomicRiskLoci.txt", header = T)
 
 femur_grl <- read.table("Femur_FUMA/SNP2GENE/GenomicRiskLoci.txt", header = T)
 FI_grl <- read.table("FI_FUMA/SNP2GENE/GenomicRiskLoci.txt", header = T)
 SA_grl <- read.table("SA_FUMA/SNP2GENE/GenomicRiskLoci.txt", header = T)
 IC_grl <- read.table("IC_FUMA/SNP2GENE/GenomicRiskLoci.txt", header = T)
 forearm_grl <- read.table("Forearm_FUMA/SNP2GENE/GenomicRiskLoci.txt", header = T)
 humerus_grl <- read.table("Humerus_FUMA/SNP2GENE/GenomicRiskLoci.txt", header = T)
 tibia_grl <- read.table("Tibia_FUMA/SNP2GENE/GenomicRiskLoci.txt", header = T)
 torso_grl <- read.table("Torso_FUMA/SNP2GENE/GenomicRiskLoci.txt", header = T)
 shoulders_grl <- read.table("Shoulder_width_FUMA/SNP2GENE/GenomicRiskLoci.txt", header = T)
 height_grl <- read.table("Height_FUMA/SNP2GENE/GenomicRiskLoci.txt", header = T)
 hip_grl <- read.table("Hip_width_FUMA/SNP2GENE/GenomicRiskLoci.txt", header = T)
 


# If the genomic span of a locus is less than 1 megabase, make it one 
# megabase by adding an equal number of bases on each end. 

traits <- c('cortical_grl', 'femur_grl', 'FI_grl', 'forearm_grl',
            'height_grl', 'hip_grl', 'humerus_grl', 'IC_grl',
            'limb_grl', 'SA_grl', 'shoulders_grl', 'tibia_grl',
            'torso_grl', 'trunk_grl')

 for(trait in traits){
 
 get_trait <- get(trait)  
   
 for(grl in 1:nrow(get_trait)){
 
 start <- get_trait$start[grl]
 end <- get_trait$end[grl]
 diff <- end - start
 
 # Extend to 1 Mb minimum
 
 if(diff < 1000000){
   
   to_1mb <- 1000000 - diff
   
   get_trait$start[grl] <- start - ceiling(to_1mb/2)
   get_trait$end[grl] <- end + ceiling(to_1mb/2)
 }
 }
 
 get_trait <- get_trait[,-1]
 get_trait <- merge_overlapping_loci(get_trait)
 get_trait[[paste0(trait, "_locus")]] <- 1:nrow(get_trait)
 
 # Give each locus a name
 assign(trait, get_trait) }
 
# 6## Define overlapping regions and run susie coloc
   # # Start only on overlap between the 3 multivariate GWAS
 
 trunk_grl_match_df <- as.data.frame(matrix(nrow = nrow(trunk_grl), ncol = 14))
 colnames(trunk_grl_match_df) <- traits
 trunk_grl_match_df$trunk_grl <- 1:nrow(trunk_grl)
 
 traits_trunk <- c('cortical_grl', 'femur_grl', 'FI_grl', 'forearm_grl',
             'height_grl', 'hip_grl', 'humerus_grl', 'IC_grl',
             'limb_grl', 'SA_grl', 'shoulders_grl', 'tibia_grl',
             'torso_grl')
 
 
 for(trait in traits_trunk){
   get_trait <- get(trait)
   
   for(i in 1:nrow(trunk_grl)){
     grl <- trunk_grl[i,]
     bps <- seq(grl$start, grl$end)
     
     if(nrow(get_trait[get_trait$chr == grl$chr & (get_trait$start %in% bps | get_trait$end %in% bps),]) > 0){
       trunk_grl_match_df[[trait]][i] <- paste0(get_trait[get_trait$chr == grl$chr &
                                                            (get_trait$start %in% bps |
                                                               get_trait$end %in% bps),][[paste0(trait, "_locus")]], collapse = ";")
     }
   }}


# ── 1. Combine all loci into one long data frame ──────────────────────────────
# 
 traits <- c('cortical_grl', 'femur_grl', 'FI_grl', 'forearm_grl',
             'height_grl', 'hip_grl', 'humerus_grl', 'IC_grl',
             'limb_grl', 'SA_grl', 'shoulders_grl', 'tibia_grl', 'torso_grl',
             'trunk_grl')
 
 all_loci <- bind_rows(lapply(traits, function(t) {
   df <- get(t)
   df$trait     <- t
   df$locus_id  <- df[[paste0(t, "_locus")]]
   df[, c("trait", "locus_id", "chr", "start", "end")]
 }))

 write.table(all_loci, "./../SuSiE_coloc/all_loci.tsv", sep = "\t", col.names = T)

all_loci <- read.table("./../SuSiE_coloc/all_loci.tsv", sep= "\t", header =T)

# ── 2. Get all unique trait pairs ─────────────────────────────────────────────
# 
 focal_traits <- c('trunk_grl', 'cortical_grl', 'limb_grl')
 focal_sans_limb <- c('trunk_grl', 'cortical_grl')
 other_traits <- setdiff(traits, focal_traits)
 
 trait_pairs <- c(
   # focal vs all others
   lapply(focal_traits, function(f) lapply(other_traits, function(o) c(f, o))),
   # focal vs focal (e.g. cortical vs limb)
   lapply('limb_grl', function(f) lapply(focal_sans_limb, function(o) c(f, o))),
   lapply('trunk_grl', function(f) lapply("cortical_grl", function(o) c(f, o)))
 ) %>% unlist(recursive = FALSE)
 
# # ── 3. Find overlapping loci for each pair ────────────────────────────────────
# 
 find_overlapping_pairs <- function(t1, t2, loci_df) {
   l1 <- filter(loci_df, trait == t1)
   l2 <- filter(loci_df, trait == t2)
   
   # Cross join then filter to overlapping intervals on same chromosome
   cross_join(l1, l2, suffix = c("_t1", "_t2")) %>%
     filter(chr_t1 == chr_t2,
            start_t1 <= end_t2,
            start_t2 <= end_t1) %>%
     transmute(
       trait1    = trait_t1,
       locus_id1 = locus_id_t1,
       trait2    = trait_t2,
       locus_id2 = locus_id_t2
     )
 }
 
 coloc_pairs <- bind_rows(lapply(trait_pairs, function(pair) {
   find_overlapping_pairs(pair[1], pair[2], all_loci)
 }))
 
# # ── 4. Inspect ────────────────────────────────────────────────────────────────
# 
 # How many pairs per trait combination?
 coloc_pairs_combinations <- coloc_pairs %>% count(trait1, trait2)
 
 write.table(coloc_pairs_combinations, "./../SuSiE_coloc/coloc_pairs_combinations.tsv",
             sep = "\t", col.names = T)
 
 write.table(coloc_pairs, "./../SuSiE_coloc/coloc_pairs.tsv",
             sep = "\t", col.names = T)

coloc_pairs <- read.table("./../SuSiE_coloc/coloc_pairs.tsv", sep = "\t", header = T)

#################################### COLOC #####################################

cortical_gwas <- read.table("./../formatted_GWAS/Cortical_GWAS.tsv.gz", header = T)
trunk_gwas <- read.table("./../formatted_GWAS/TrunkWidths_GWAS.tsv.gz", header = T)
limb_gwas <-  read.table("./../formatted_GWAS/LimbHeight_GWAS.tsv.gz", header = T)

cortical_gwas <- cortical_gwas[,c(1:3,5,6,4,7,8,9)]
trunk_gwas <- trunk_gwas[,c(1:3,5,6,4,7,8,9)]
limb_gwas <- limb_gwas[,c(1:3,5,6,4,7,8,9)]

femur_gwas <- read.table("./../formatted_GWAS/femur.tsv.gz", header = T)
FI_gwas <- read.table("./../formatted_GWAS/FI.tsv.gz", header = T)
SA_gwas <- read.table("./../formatted_GWAS/SA.tsv.gz", header = T)
IC_gwas <- read.table("./../formatted_GWAS/IC.tsv.gz", header = T)
forearm_gwas <- read.table("./../formatted_GWAS/forearm.tsv.gz", header = T)
humerus_gwas <- read.table("./../formatted_GWAS/humerus.tsv.gz", header = T)
tibia_gwas <- read.table("./../formatted_GWAS/tibia.tsv.gz", header = T)
torso_gwas <- read.table("./../formatted_GWAS/torso.tsv.gz", header = T)
shoulders_gwas <- read.table("./../formatted_GWAS/shoulder_width.tsv.gz", header = T)
height_gwas <- read.table("./../formatted_GWAS/height.tsv.gz", header = T)
hip_gwas <- read.table("./../formatted_GWAS/hip_width.tsv.gz", header = T)

cortical_gwas <- cortical_gwas[!duplicated(cortical_gwas$SNP),]
trunk_gwas <- trunk_gwas[!duplicated(trunk_gwas$SNP),]
limb_gwas <- limb_gwas[!duplicated(limb_gwas$SNP),]

femur_gwas <- femur_gwas[!duplicated(femur_gwas$SNP),]
FI_gwas <- FI_gwas[!duplicated(FI_gwas$SNP),]
SA_gwas <- SA_gwas[!duplicated(SA_gwas$SNP),]
IC_gwas <- IC_gwas[!duplicated(IC_gwas$SNP),]
forearm_gwas <- forearm_gwas[!duplicated(forearm_gwas$SNP),]
humerus_gwas <- humerus_gwas[!duplicated(humerus_gwas$SNP),]
tibia_gwas <- tibia_gwas[!duplicated(tibia_gwas$SNP),]
torso_gwas <- torso_gwas[!duplicated(torso_gwas$SNP),]
shoulders_gwas <- shoulders_gwas[!duplicated(shoulders_gwas$SNP),]
height_gwas <- height_gwas[!duplicated(height_gwas$SNP),]
hip_gwas <- hip_gwas[!duplicated(hip_gwas$SNP),]


# ── 0. Parameters ─────────────────────────────────────────────────────────────

KG_PLINK_PREFIX <- "~/Desktop/PhD/tools/1000G/"  # e.g. chr1.bed/bim/fam

# Named vector of sample sizes for each trait
sample_sizes <- c(
  cortical_grl  = 18037,
  femur_grl     = 53000,
  FI_grl        = 53000,
  forearm_grl   = 53000,
  height_grl    = 458303,
  hip_grl       = 53000,
  humerus_grl   = 53000,
  IC_grl        = 53000,
  limb_grl      = 39556,
  SA_grl        = 53000,
  shoulders_grl = 53000,
  tibia_grl     = 53000,
  torso_grl     = 53000,
  trunk_grl = 54881
)

# ── 1. Helper: extract SNPs in a window from a full GWAS summary stats df ─────
# Assumes gwas df has columns: SNP, CHR, BP, BETA, SE, P, A1, A2

extract_locus <- function(gwas_df, chr, start, end) {
  gwas_df[gwas_df$CHR == chr & gwas_df$BP >= start & gwas_df$BP <= end, ]
}

# ── 2. Helper: compute LD from 1000G for a set of SNPs ────────────────────────

compute_ld <- function(snps, chr) {
  snp_file <- tempfile(fileext = ".txt")
  plink_out <- tempfile()
  
  writeLines(snps, snp_file)
  
  system(paste(
    "plink",
    "--bfile",    paste0(KG_PLINK_PREFIX, "EUR_phase3_chr", chr),
    "--extract",  snp_file,
    "--r square",
    "--out",      plink_out,
    "--keep-allele-order",
    "--make-just-bim",
    "--memory 8000"
  ))
  
  ld_snps <- fread(paste0(plink_out, ".bim"))$V2
  ld_mat  <- as.matrix(fread(paste0(plink_out, ".ld"), header = FALSE))
  
  list(ld_mat = ld_mat, ld_snps = ld_snps)
}


# ── 3. Helper: run SuSiE on one trait in a window ─────────────────────────────

run_susie <- function(gwas_df, ld_snps, ld_mat, n) {
  # Align GWAS to PLINK SNP order
  gwas_df <- gwas_df[match(ld_snps, gwas_df$SNP), ]

  gwas_df$varbeta <- gwas_df$SE^2
  gwas_df <- gwas_df[,c(1,7,10,6,9)]
  colnames(gwas_df) <- c("snp", "beta", "varbeta", "MAF", "pvalues")
  gwas_df$position <- 1:nrow(gwas_df)
  
  gwas_df <- as.list(gwas_df)
  gwas_df$N <- n
  gwas_df$type <- "quant"
  

  colnames(ld_mat) <- gwas_df$snp
  rownames(ld_mat) <- gwas_df$snp
  
  gwas_df$LD <- ld_mat
  
  runsusie(gwas_df)
}


# ── 4. Main coloc loop ────────────────────────────────────────────────────────

coloc_results <- vector("list", nrow(coloc_pairs))

for (i in seq_len(nrow(coloc_pairs))) {
  
  t1 <- coloc_pairs$trait1[i];    t2 <- coloc_pairs$trait2[i]
  l1 <- coloc_pairs$locus_id1[i]; l2 <- coloc_pairs$locus_id2[i]
  
  cat(sprintf("Running coloc %d/%d: %s locus %s vs %s locus %s\n",
              i, nrow(coloc_pairs), t1, l1, t2, l2))
  
  # Get locus coordinates from the merged loci data frames
  locus1 <- all_loci[all_loci$locus_id == l1 & all_loci$trait == t1, ]
  locus2 <- all_loci[all_loci$locus_id == l2 & all_loci$trait == t2, ]
  
  # Define union window
  region_chr   <- locus1$chr
  region_start <- min(locus1$start, locus2$start)
  region_end   <- max(locus1$end,   locus2$end)
  
  # Extract summary stats for both traits in union window
  gwas1 <- extract_locus(get(paste0(gsub(x= t1, pattern = "_grl", replacement = ""), "_gwas")), region_chr, region_start, region_end)
  gwas2 <- extract_locus(get(paste0(gsub(x= t2, pattern = "_grl", replacement = ""), "_gwas")), region_chr, region_start, region_end)
  
  
  # Keep overlapping SNPs only
  shared_snps <- intersect(gwas1$SNP, gwas2$SNP)
  gwas1 <- gwas1[gwas1$SNP %in% shared_snps, ]
  gwas2 <- gwas2[gwas2$SNP %in% shared_snps, ]
  
  gwas2 <- merge(gwas2, gwas1[,c(1,6)], by = "SNP")
  gwas2 <- gwas2[,c(1:5,10,7,8,9)]
  colnames(gwas2)[6] <- "MAF"
  
  # Skip if too few SNPs
  if (length(shared_snps) < 50) {
    cat(sprintf("  Skipping: only %d shared SNPs\n", length(shared_snps)))
    next
  }
  
  # Compute LD
  ld <- tryCatch(
    compute_ld(shared_snps, region_chr),
    error = function(e) { cat("  LD computation failed:", e$message, "\n"); NULL }
  )
  if (is.null(ld)) next
  
  # Run SuSiE
  susie1 <- tryCatch(
    run_susie(gwas_df = gwas1, ld_snps = ld$ld_snps, ld_mat = ld$ld_mat, n = sample_sizes[t1]),
    error = function(e) { cat("  SuSiE failed for", t1, ":", e$message, "\n"); NULL }
  )
  susie2 <- tryCatch(
    run_susie(gwas_df = gwas2,ld_snps =  ld$ld_snps,ld_mat =  ld$ld_mat, n = sample_sizes[t2]),
    error = function(e) { cat("  SuSiE failed for", t2, ":", e$message, "\n"); NULL }
  )
  if (is.null(susie1) || is.null(susie2)) next
  
  # Run coloc
  res <- tryCatch(
    coloc.susie(susie1, susie2),
    error = function(e) { cat("  coloc failed:", e$message, "\n"); NULL }
  )
  if (is.null(res)) next
  
  # Store results with pair metadata
  coloc_results[[i]] <- res$summary %>%
    mutate(
      trait1    = t1, locus_id1 = l1,
      trait2    = t2, locus_id2 = l2,
      region    = sprintf("chr%s:%d-%d", region_chr, region_start, region_end),
      nSNPs     = length(shared_snps)
    )
}

# ── 5. Collate results ────────────────────────────────────────────────────────

coloc_results_df <- bind_rows(coloc_results)

# Significant coloc hits (PP.H4 >= 0.8)
coloc_hits <- coloc_results_df %>% filter(PP.H4.abf >= 0.8)

# Save
write.csv(coloc_results_df, "susie_coloc_results_all.csv",  row.names = FALSE)
write.csv(coloc_hits,       "susie_coloc_results_hits.csv", row.names = FALSE)

