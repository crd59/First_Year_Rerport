# =============================================================================
#  Locus-zoom plots for colocalisation results  (ggplot2)
#  One self-contained plot object per gene.
#
#  Structure:
#    1. SETUP  -- libraries, file paths, column config, helpers, palette, theme
#    2. PER-GENE BLOCKS -- each gene defines its own traits + highlighted leads
#                          and builds its own plot object (p_<GENE>) + saves a file
#
#  The mechanical parts (reading sumstats, subsetting a region, looking up a
#  lead SNP, the base aesthetics) live in shared helpers so they aren't copy-
#  pasted 12 times. Everything that actually DIFFERS between genes -- which
#  traits to show, which SNPs to highlight, and what colour category each
#  highlight gets -- is written out explicitly in that gene's block, so you can
#  edit any single plot without touching the others.
# =============================================================================

setwd("~/Desktop/PhD/First_project/FLAMES_2.0/plot_files")

# ============================== 1. SETUP =====================================

library(data.table)   # fast reading / subsetting of genome-wide sumstats
library(ggplot2)
library(ggrepel)      # non-overlapping SNP labels

## ---- EDIT THESE: paths to your three full summary-statistics files ----------
sumstat_paths <- c(
  cortical_grl = "./../../formatted_GWAS/Cortical_GWAS.tsv.gz",
  limb_grl     = "./../../formatted_GWAS/LimbHeight_GWAS.tsv.gz",
  trunk_grl    = "./../../formatted_GWAS/TrunkWidths_GWAS.tsv.gz"
)

## ---- EDIT THESE: the column names AS THEY APPEAR in your sumstats ------------
COL_SNP <- "SNP"   # rsID column
COL_CHR <- "CHR"   # chromosome
COL_BP  <- "BP"    # base-pair position (build must match your `region` strings)
COL_P   <- "P"     # p-value
## If you only have Z (or BETA/SE) rather than P, derive P before plotting, e.g.
##   dt[, P := 2 * pnorm(-abs(Z))]

# Read each file once, keep only the columns we need, standardise their names.
read_sumstats <- function(path) {
  dt <- fread(path)
  dt <- dt[, .(SNP = get(COL_SNP),
               CHR = as.character(get(COL_CHR)),
               BP  = as.numeric(get(COL_BP)),
               P   = as.numeric(get(COL_P)))]
  dt[, CHR := gsub("chr", "", CHR, ignore.case = TRUE)]
  # guard against P == 0 underflowing to Inf in -log10
  dt[, logP := -log10(pmax(P, .Machine$double.xmin))]
  setkey(dt, CHR, BP)
  dt[]
}
ss <- lapply(sumstat_paths, read_sumstats)   # named list of data.tables

# "chr7:91753630-93169239" -> list(chr="7", start=91753630, end=93169239)
parse_region <- function(region) {
  chr <- sub("chr", "", sub(":.*$", "", region), ignore.case = TRUE)
  rng <- sub("^.*:", "", region)
  list(chr   = chr,
       start = as.numeric(sub("-.*$", "", rng)),
       end   = as.numeric(sub("^.*-", "", rng)))
}

# All SNPs of one trait within a region, tagged with the trait name.
region_subset <- function(trait, region) {
  r <- parse_region(region)
  d <- ss[[trait]][CHR == r$chr & BP %between% c(r$start, r$end)]
  d[, trait := trait]
  as.data.frame(d)
}

# Look up a single lead SNP in one trait -> 1-row data.frame (BP, logP).
lead_lookup <- function(trait, rsid) {
  hit <- ss[[trait]][SNP == rsid]
  if (nrow(hit) == 0) {
    warning(sprintf("Lead SNP '%s' not found in trait '%s'", rsid, trait))
    return(data.frame(SNP = rsid, BP = NA_real_, logP = NA_real_))
  }
  data.frame(SNP = rsid, BP = hit$BP[1], logP = hit$logP[1])
}

## ---- Colour palettes --------------------------------------------------------
# Background association points: one muted colour per trait.
trait_pal <- c(
  cortical_grl = "#A6CEE3",   # light blue
  limb_grl     = "#B2DF8A",   # light green
  trunk_grl    = "#FDBF6F"    # light orange
)

# Highlighted (fine-mapped / colocalised) lead SNPs: vivid colour per category.
# The three single-trait colours are the darker version of that trait's
# background colour, so a non-colocalised causal variant visually matches its
# trait. Coloc pairs and the triple-coloc each get their own distinct colour.
highlight_pal <- c(
  "cortical_grl"           = "#1F78B4",   # cortical-only causal variant
  "limb_grl"               = "#33A02C",   # limb-only causal variant
  "trunk_grl"              = "#FF7F00",   # trunk-only causal variant
  "Coloc cortical-limb"    = "#6A3D9A",   # purple
  "Coloc cortical-trunk"   = "#A6761D",   # gold/brown
  "Coloc limb-trunk"       = "#E7298A",   # magenta
  "Triple coloc"           = "#E31A1C"    # red
)

## ---- Shared plot builder ----------------------------------------------------
# traits : character vector of traits to draw as background association points
# leads  : data.frame(trait, rsid, category)  where `category` is a key of
#          highlight_pal. Each row becomes one highlighted diamond + label.
make_locuszoom <- function(gene, region, traits, leads) {
  
  bg <- do.call(rbind, lapply(traits, region_subset, region = region))
  
  lead_df <- do.call(rbind, lapply(seq_len(nrow(leads)), function(i) {
    pos <- lead_lookup(leads$trait[i], leads$rsid[i])
    data.frame(SNP      = pos$SNP,
               BP       = pos$BP,
               logP     = pos$logP,
               trait    = leads$trait[i],
               category = leads$category[i])
  }))
  # keep category factor order consistent with the palette
  lead_df$category <- factor(lead_df$category, levels = names(highlight_pal))
  
  r <- parse_region(region)
  
  ggplot() +
    # background association signal, coloured by trait
    geom_point(data = bg,
               aes(x = BP / 1e6, y = logP, colour = trait),
               size = 1, alpha = 0.5) +
    # genome-wide significance reference line
    geom_hline(yintercept = -log10(5e-8),
               linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    # highlighted lead SNPs as bordered diamonds
    geom_point(data = lead_df,
               aes(x = BP / 1e6, y = logP, fill = category),
               shape = 23, size = 4.5, colour = "black", stroke = 0.7) +
    # rsID labels
    geom_label_repel(data = lead_df,
                     aes(x = BP / 1e6, y = logP, label = SNP),
                     size = 3, fill = "white", alpha = 0.85,
                     label.size = 0.15, box.padding = 0.6,
                     min.segment.length = 0, seed = 1) +
    scale_colour_manual(values = trait_pal, name = "Trait association",
                        drop = TRUE) +
    scale_fill_manual(values = highlight_pal, name = "Highlighted variant",
                      drop = TRUE) +
    guides(
      colour = guide_legend(order = 1, override.aes = list(size = 3, alpha = 1)),
      fill   = guide_legend(order = 2, override.aes = list(shape = 23, size = 4))
    ) +
    labs(title    = gene,
         subtitle = region,
         x = paste0("Chromosome ", r$chr, " position (Mb)"),
         y = expression(-log[10](italic(P)))) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold.italic", hjust = 0.5),
          plot.subtitle    = element_text(colour = "grey40", size = 9, hjust = 0.5),
          legend.key       = element_blank())
}

# Helper to save each plot (one file per gene, ready for PowerPoint).
save_plot <- function(p, gene, w = 8, h = 5, dpi = 320) {
  ggsave(sprintf("locuszoom_%s.png", gene), p,
         width = w, height = h, dpi = dpi)
}


# ============================ 2. PER-GENE BLOCKS =============================
# Each block: define traits + leads for that gene, build the plot, save it.
# Edit a block freely without affecting any other gene.


# -----------------------------------------------------------------------------
# SKP2  -- single causal variant (trunk only; no fine-mapped cortical variant)
#          row: cortical_grl vs trunk_grl, hit1 = NA, hit2 = rs7731023
# -----------------------------------------------------------------------------
p_SKP2 <- make_locuszoom(
  gene   = "SKP2",
  region = "chr5:35671147-36729847",
  traits = c("cortical_grl", "trunk_grl"),   # show both backgrounds
  leads  = data.frame(
    trait    = "trunk_grl",
    rsid     = "rs7731023",
    category = "trunk_grl"
  )
)
save_plot(p_SKP2, "SKP2")


# -----------------------------------------------------------------------------
# GRB10 -- single causal variant (cortical only; no fine-mapped trunk variant)
#          row: cortical_grl vs trunk_grl, hit1 = rs2237470, hit2 = NA
# -----------------------------------------------------------------------------
p_GRB10 <- make_locuszoom(
  gene   = "GRB10",
  region = "chr7:50244938-51289423",
  traits = c("cortical_grl", "trunk_grl"),
  leads  = data.frame(
    trait    = "cortical_grl",
    rsid     = "rs2237470",
    category = "cortical_grl"
  )
)
save_plot(p_GRB10, "GRB10")


# -----------------------------------------------------------------------------
# CDK6  -- TRIPLE colocalisation across cortical / limb / trunk
#          cortical~trunk (rs42039), cortical~limb (rs42039/rs42043),
#          limb~trunk (rs42039/rs42043).  Lead SNPs rs42039 & rs42043.
#          NOTE: adjust which rsID maps to which trait if you prefer.
# -----------------------------------------------------------------------------
p_CDK6 <- make_locuszoom(
  gene   = "CDK6",
  region = "chr7:91753630-93169239",
  traits = c("cortical_grl", "limb_grl", "trunk_grl"),
  leads  = data.frame(
    trait    = c("cortical_grl", "limb_grl",    "trunk_grl"),
    rsid     = c("rs42039",      "rs42043",     "rs42043"),
    category = c("Triple coloc", "Triple coloc","Triple coloc")
  )
)
save_plot(p_CDK6, "CDK6")


# -----------------------------------------------------------------------------
# CENPW -- TRIPLE colocalisation across cortical / limb / trunk
#          cortical~limb (rs1490384), cortical~trunk (rs1490384/rs9491652),
#          limb~trunk (rs9491652/rs1490384). Lead SNPs rs1490384 & rs9491652.
#          NOTE: adjust which rsID maps to which trait if you prefer.
# -----------------------------------------------------------------------------
p_CENPW <- make_locuszoom(
  gene   = "CENPW",
  region = "chr6:126040435-127849292",
  traits = c("cortical_grl", "limb_grl", "trunk_grl"),
  leads  = data.frame(
    trait    = c("cortical_grl", "limb_grl",     "trunk_grl"),
    rsid     = c("rs1490384",    "rs1490384",    "rs9491652"),
    category = c("Triple coloc", "Triple coloc", "Triple coloc")
  )
)
save_plot(p_CENPW, "CENPW")


# -----------------------------------------------------------------------------
# CTNNB1 -- pairwise colocalisation cortical~trunk (PP.H4 = 0.884)
#           hit1 = rs7630377 (cortical), hit2 = rs9311265 (trunk)
# -----------------------------------------------------------------------------
p_CTNNB1 <- make_locuszoom(
  gene   = "CTNNB1",
  region = "chr7:19925319-55281213",
  traits = c("cortical_grl", "trunk_grl"),
  leads  = data.frame(
    trait    = c("cortical_grl",         "trunk_grl"),
    rsid     = c("rs7630377",            "rs9311265"),
    category = c("Coloc cortical-trunk", "Coloc cortical-trunk")
  )
)
save_plot(p_CTNNB1, "CTNNB1")


# -----------------------------------------------------------------------------
# ITGB8 -- pairwise colocalisation cortical~limb (PP.H4 = 0.887)
#          hit1 = rs6973059 (cortical), hit2 = rs3757727 (limb)
# -----------------------------------------------------------------------------
p_ITGB8 <- make_locuszoom(
  gene   = "ITGB8",
  region = "chr7:19101819-20925320",
  traits = c("cortical_grl", "limb_grl"),
  leads  = data.frame(
    trait    = c("cortical_grl",        "limb_grl"),
    rsid     = c("rs6973059",           "rs3757727"),
    category = c("Coloc cortical-limb", "Coloc cortical-limb")
  )
)
save_plot(p_ITGB8, "ITGB8")


# -----------------------------------------------------------------------------
# TET2  -- pairwise colocalisation limb~trunk (PP.H4 = 0.9995)
#          hit1 = hit2 = rs2903385 (same lead SNP in both traits)
#          Third trait (cortical) shown as background only -- no fine-mapped
#          variant, so it is not annotated.
# -----------------------------------------------------------------------------
p_TET2 <- make_locuszoom(
  gene   = "TET2",
  region = "chr4:105613841-106680233",
  traits = c("cortical_grl", "limb_grl", "trunk_grl"),
  leads  = data.frame(
    trait    = c("limb_grl",         "trunk_grl"),
    rsid     = c("rs2903385",        "rs2903385"),
    category = c("Coloc limb-trunk", "Coloc limb-trunk")
  )
)
save_plot(p_TET2, "TET2")


# -----------------------------------------------------------------------------
# LCORL -- NO colocalisation cortical~limb (PP.H4 = 2.7e-08)
#          two separate causal variants, coloured by their own trait
#          hit1 = rs6813340 (cortical), hit2 = rs4144829 (limb)
#          Third trait (trunk) shown as background only -- no fine-mapped
#          variant, so it is not annotated.
# -----------------------------------------------------------------------------
p_LCORL <- make_locuszoom(
  gene   = "LCORL",
  region = "chr4:17338304-18416110",
  traits = c("cortical_grl", "limb_grl", "trunk_grl"),
  leads  = data.frame(
    trait    = c("cortical_grl", "limb_grl"),
    rsid     = c("rs6813340",    "rs4144829"),
    category = c("cortical_grl", "limb_grl")
  )
)
save_plot(p_LCORL, "LCORL")


# -----------------------------------------------------------------------------
# PIK3CD -- NO colocalisation cortical~limb (PP.H4 = 2.0e-05)
#           hit1 = rs79583983 (cortical), hit2 = rs9442571 (limb)
# -----------------------------------------------------------------------------
p_PIK3CD <- make_locuszoom(
  gene   = "PIK3CD",
  region = "chr1:8846137-10238997",
  traits = c("cortical_grl", "limb_grl"),
  leads  = data.frame(
    trait    = c("cortical_grl", "limb_grl"),
    rsid     = c("rs79583983",   "rs9442571"),
    category = c("cortical_grl", "limb_grl")
  )
)
save_plot(p_PIK3CD, "PIK3CD")


# -----------------------------------------------------------------------------
# PTEN  -- NO colocalisation cortical~trunk (PP.H4 = 8.7e-04)
#          hit1 = rs1234214 (cortical), hit2 = rs1044322 (trunk)
# -----------------------------------------------------------------------------
p_PTEN <- make_locuszoom(
  gene   = "PTEN",
  region = "chr10:89195453-90208342",
  traits = c("cortical_grl", "trunk_grl"),
  leads  = data.frame(
    trait    = c("cortical_grl", "trunk_grl"),
    rsid     = c("rs1234214",    "rs1044322"),
    category = c("cortical_grl", "trunk_grl")
  )
)
save_plot(p_PTEN, "PTEN")


# -----------------------------------------------------------------------------
# HMGA2 -- cortical~trunk COLOCALISE at rs7968682 (PP.H4 = 0.9975),
#          plus a third trait (limb) present in a separate row with its own
#          lead rs7959830 (no coloc) -> highlight + annotate it too.
# -----------------------------------------------------------------------------
p_HMGA2 <- make_locuszoom(
  gene   = "HMGA2",
  region = "chr12:65154559-66855026",
  traits = c("cortical_grl", "limb_grl", "trunk_grl"),
  leads  = data.frame(
    trait    = c("cortical_grl",         "trunk_grl",            "limb_grl"),
    rsid     = c("rs7968682",            "rs7968682",            "rs7959830"),
    category = c("Coloc cortical-trunk", "Coloc cortical-trunk", "limb_grl")
  )
)
save_plot(p_HMGA2, "HMGA2")


# -----------------------------------------------------------------------------
# PAPPA -- cortical~limb COLOCALISE (PP.H4 = 0.970), hit1 = rs10817865
#          (cortical), hit2 = rs3789283 (limb); plus a third trait (trunk)
#          present in a separate row with its own lead rs71505577 (no coloc)
#          -> highlight + annotate it too.
# -----------------------------------------------------------------------------
p_PAPPA <- make_locuszoom(
  gene   = "PAPPA",
  region = "chr9:117983738-119762704",
  traits = c("cortical_grl", "limb_grl", "trunk_grl"),
  leads  = data.frame(
    trait    = c("cortical_grl",        "limb_grl",            "trunk_grl"),
    rsid     = c("rs10817865",          "rs3789283",           "rs71505577"),
    category = c("Coloc cortical-limb", "Coloc cortical-limb", "trunk_grl")
  )
)
save_plot(p_PAPPA, "PAPPA")


# =============================================================================
# All plot objects (p_SKP2, p_GRB10, p_CDK6, ...) are now in your environment
# and saved as locuszoom_<GENE>.png. Print any one interactively, e.g.:  p_CDK6
# =============================================================================