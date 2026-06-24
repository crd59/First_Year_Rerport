# =============================================================================
#  LD-coloured locus-zoom plots for CDK6 and CENPW
#
#  For each gene: THREE plots (one per trait). In each plot every SNP is
#  coloured by its r^2 with that trait's lead SNP (computed from the gene's LD
#  matrix .rda), the lead SNP is a purple diamond, and an arrow on the lead
#  shows its direction of effect (up = positive beta, down = negative beta).
#
#  Each .rda loads an object called `ld`, a list with:
#     ld$ld_mat   -- the LD matrix in r values (NOT r^2)
#     ld$ld_snps  -- SNP names in the order they appear in ld_mat
#
#  CDK6 and CENPW each get their own clearly separated section below.
#  (Remove the CDK6 / CENPW blocks from your first script so they aren't drawn
#   twice in two different styles.)
# =============================================================================

setwd("~/Desktop/PhD/First_project/FLAMES_2.0/plot_files/")

# ============================== SETUP ========================================

library(data.table)
library(ggplot2)
library(ggrepel)

## ---- EDIT: paths to your three full summary-statistics files -----------------
sumstat_paths <- c(
  cortical_grl = "./../../formatted_GWAS/Cortical_GWAS.tsv.gz",
  limb_grl     = "./../../formatted_GWAS/LimbHeight_GWAS.tsv.gz",
  trunk_grl    = "./../../formatted_GWAS/TrunkWidths_GWAS.tsv.gz"
)

## ---- EDIT: column names AS THEY APPEAR in your sumstats ----------------------
COL_SNP  <- "SNP"    # rsID (must match the names in ld$ld_snps)
COL_CHR  <- "CHR"
COL_BP   <- "BP"
COL_P    <- "P"
COL_BETA <- "BETA"   # effect size; only its SIGN is used for the arrow.
## If you have OR instead of BETA, build a BETA column first: BETA <- log(OR).

# Read one trait's sumstats with read.table().
# (For a .gz file use read.table(gzfile(path), ...) instead.)
read_sumstats <- function(path) {
  raw <- read.table(path, header = TRUE, sep = "",
                    stringsAsFactors = FALSE, check.names = FALSE)
  d <- data.frame(
    SNP  = as.character(raw[[COL_SNP]]),
    CHR  = gsub("chr", "", as.character(raw[[COL_CHR]]), ignore.case = TRUE),
    BP   = as.numeric(raw[[COL_BP]]),
    P    = as.numeric(raw[[COL_P]]),
    BETA = as.numeric(raw[[COL_BETA]]),
    stringsAsFactors = FALSE
  )
  d$logP <- -log10(pmax(d$P, .Machine$double.xmin))
  d
}
ss <- lapply(sumstat_paths, read_sumstats)

# Load an LD matrix (.rda holding the `ld` list) and return an r matrix whose
# row/column names are the SNP IDs from ld$ld_snps.
load_ld <- function(path) {
  e <- new.env()
  load(path, envir = e)          # creates `ld` in environment e
  ld <- e$ld
  m <- as.matrix(ld$ld_mat)      # r values
  rownames(m) <- ld$ld_snps
  colnames(m) <- ld$ld_snps
  m
}

# "chr7:91753630-93169239" -> list(chr, start, end)
parse_region <- function(region) {
  chr <- sub("chr", "", sub(":.*$", "", region), ignore.case = TRUE)
  rng <- sub("^.*:", "", region)
  list(chr   = chr,
       start = as.numeric(sub("-.*$", "", rng)),
       end   = as.numeric(sub("^.*-", "", rng)))
}

# SNPs of one trait within a region.
region_subset <- function(trait, region) {
  r <- parse_region(region)
  d <- ss[[trait]]
  d[d$CHR == r$chr & d$BP >= r$start & d$BP <= r$end, , drop = FALSE]
}

# Look up a single SNP (position, logP, beta) in one trait.
lead_lookup <- function(trait, rsid) {
  d <- ss[[trait]]
  hit <- d[d$SNP == rsid, , drop = FALSE]
  if (nrow(hit) == 0) {
    warning(sprintf("Lead SNP '%s' not found in trait '%s'", rsid, trait))
    return(data.frame(SNP = rsid, BP = NA_real_, logP = NA_real_, BETA = NA_real_))
  }
  data.frame(SNP = rsid, BP = hit$BP[1], logP = hit$logP[1], BETA = hit$BETA[1])
}

## ---- Shared builder: one trait, coloured by LD, with effect-direction arrow --
make_ld_locuszoom <- function(gene, region, trait, lead_snp, ld_mat) {
  
  lead_snp <- unname(lead_snp)
  df <- region_subset(trait, region)
  
  # r^2 of every SNP with the lead SNP (ld_mat holds r, so square it)
  if (lead_snp %in% rownames(ld_mat)) {
    ldvec <- ld_mat[lead_snp, ]
    ld_df <- data.frame(SNP = names(ldvec), r2 = as.numeric(ldvec)^2)
    df <- merge(df, ld_df, by = "SNP", all.x = TRUE)
  } else {
    warning(sprintf("Lead '%s' not in LD matrix for %s/%s — all SNPs grey.",
                    lead_snp, gene, trait))
    df$r2 <- NA_real_
  }
  df <- df[order(df$r2, na.last = FALSE), ]   # low r^2 / no-LD underneath
  
  # lead SNP position + effect direction
  lead <- lead_lookup(trait, lead_snp)
  up   <- !is.na(lead$BETA) && lead$BETA > 0
  
  span <- diff(range(df$logP, na.rm = TRUE))
  if (!is.finite(span) || span == 0) span <- 1
  gap  <- 0.07 * span        # lift the arrow clear of the diamond
  alen <- 0.10 * span        # arrow length
  base <- lead$logP + gap
  arrow_df <- data.frame(
    x    = lead$BP / 1e6,
    y    = if (up) base else base + alen,   # arrow always sits above the point
    yend = if (up) base + alen else base    # only the direction flips
  )
  
  # nudge the rsID label to the side (away from the nearest edge) so its
  # connector line does not run up through the vertical arrow.
  rr       <- parse_region(region)
  width_mb <- (rr$end - rr$start) / 1e6
  mid_mb   <- (rr$start + rr$end) / 2 / 1e6
  nx       <- if ((lead$BP / 1e6) > mid_mb) -0.12 * width_mb else 0.12 * width_mb
  
  ggplot() +
    geom_point(data = df, aes(BP / 1e6, logP, fill = r2),
               shape = 21, size = 2, stroke = 0, alpha = 0.8) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed",
               colour = "grey50", linewidth = 0.3) +
    # lead SNP (yellow diamond matching r^2 = 1, with black outline)
    geom_point(data = lead, aes(BP / 1e6, logP),
               shape = 23, size = 3, fill = "yellow",
               colour = "black", stroke = 0.7) +
    # direction-of-effect arrow, sitting above the lead SNP (black)
    geom_segment(data = arrow_df, aes(x = x, xend = x, y = y, yend = yend),
                 arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                 linewidth = 0.7, colour = "black") +
    geom_label_repel(data = lead, aes(BP / 1e6, logP, label = SNP),
                     size = 3, fill = "white", alpha = 0.85, label.size = 0.15,
                     box.padding = 0.7, min.segment.length = 0, seed = 1,
                     nudge_x = nx, nudge_y = gap + alen) +
    scale_fill_viridis_c(option = "inferno", limits = c(0, 1),
                         na.value = "grey85",
                         name = expression(r^2~"with lead")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(title    = bquote(bolditalic(.(gene))~bold("\u2014")~bold(.(trait))),
         subtitle = region,
         x = paste0("Chromosome ", parse_region(region)$chr, " position (Mb)"),
         y = expression(-log[10](italic(P)))) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title    = element_text(hjust = 0.5),
          plot.subtitle = element_text(colour = "grey40", size = 9, hjust = 0.5),
          legend.key    = element_blank())
}

save_plot <- function(p, file, w = 7, h = 5, dpi = 320) {
  ggsave(file, p, width = w, height = h, dpi = dpi)
}


# ============================ CDK6 ===========================================
# All three traits referenced to the same lead SNP.
ld_CDK6      <- load_ld("LD_matCDK6.rda")
region_CDK6  <- "chr7:91753630-93169239"
leadSNP_CDK6 <- "rs42043"

p_CDK6_cortical <- make_ld_locuszoom("CDK6", region_CDK6, "cortical_grl",
                                     leadSNP_CDK6, ld_CDK6)
p_CDK6_limb     <- make_ld_locuszoom("CDK6", region_CDK6, "limb_grl",
                                     leadSNP_CDK6, ld_CDK6)
p_CDK6_trunk    <- make_ld_locuszoom("CDK6", region_CDK6, "trunk_grl",
                                     leadSNP_CDK6, ld_CDK6)

save_plot(p_CDK6_cortical, "locuszoom_CDK6_cortical_grl.png")
save_plot(p_CDK6_limb,     "locuszoom_CDK6_limb_grl.png")
save_plot(p_CDK6_trunk,    "locuszoom_CDK6_trunk_grl.png")


# ============================ CENPW ==========================================
# All three traits referenced to the same lead SNP.
ld_CENPW      <- load_ld("LD_matCENPW.rda")
region_CENPW  <- "chr6:126040435-127849292"
leadSNP_CENPW <- "rs1490384"

p_CENPW_cortical <- make_ld_locuszoom("CENPW", region_CENPW, "cortical_grl",
                                      leadSNP_CENPW, ld_CENPW)
p_CENPW_limb     <- make_ld_locuszoom("CENPW", region_CENPW, "limb_grl",
                                      leadSNP_CENPW, ld_CENPW)
p_CENPW_trunk    <- make_ld_locuszoom("CENPW", region_CENPW, "trunk_grl",
                                      leadSNP_CENPW, ld_CENPW)

save_plot(p_CENPW_cortical, "locuszoom_CENPW_cortical_grl.png")
save_plot(p_CENPW_limb,     "locuszoom_CENPW_limb_grl.png")
save_plot(p_CENPW_trunk,    "locuszoom_CENPW_trunk_grl.png")

# =============================================================================
# Six plot objects now exist (p_CDK6_cortical ... p_CENPW_trunk) and are saved
# as locuszoom_<GENE>_<trait>.png. Print any one interactively, e.g. p_CDK6_limb
# =============================================================================