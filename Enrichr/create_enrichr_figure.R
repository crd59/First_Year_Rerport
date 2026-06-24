
# =============================================================================
# Enrichment summary: per-trait proportion of simplified GO terms, by group/theme
# -----------------------------------------------------------------------------
# Each simplified term is tagged with the trait combination it is enriched in
# (cortical / limb / trunk, "+"-separated). We split those combinations so a
# "cortical+limb+trunk" term contributes one count to each of the three traits.
# For every (group, trait) we then take that count as a PROPORTION OF THE TRAIT
# TOTAL, so the three traits are on a comparable footing despite very different
# numbers of enriched terms. Bars are grouped by `group` and faceted by `theme`,
# with the theme strip on the right-hand side.
# =============================================================================

setwd("~/Desktop/PhD/first_year_report/results_figures/")

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(stringr)

# ---- config -----------------------------------------------------------------
infile   <- "./enrichment_plot.tsv"
outfile  <- "enrichment_barplot.pdf"
trait_levels  <- c("cortical", "limb", "trunk")
trait_palette <- c(cortical = "#2C7FB8",   # blue
                   limb     = "#E6892B",   # orange
                   trunk    = "#7B3294")   # purple

# ---- read -------------------------------------------------------------------
# The file was written with row names, so the first (unnamed) column is the
# index; read.delim picks it up as rownames automatically -> 4 data columns.
dat <- read.delim(infile, stringsAsFactors = FALSE)

# ---- reshape: one row per (term, individual trait) --------------------------
long <- dat |>
  mutate(traits = strsplit(traits, "\\+")) |>
  tidyr::unnest(traits) |>
  mutate(traits = factor(str_trim(traits), levels = trait_levels))

# ---- counts and per-trait proportions ---------------------------------------
plot_df <- long |>
  count(theme, group, traits, name = "n") |>
  group_by(traits) |>
  mutate(proportion = n / sum(n)) |>      # normalise within each trait
  ungroup()

# ---- ordering ---------------------------------------------------------------
# Order groups (within their theme band) by total term count; order themes the
# same way so the busiest blocks sit together.
group_tot <- plot_df |> count(theme, group, wt = n, name = "tot")
theme_tot <- group_tot |> group_by(theme) |>
  summarise(tot = sum(tot), .groups = "drop")

plot_df <- plot_df |>
  mutate(
    group = factor(group, levels = group_tot$group[order(group_tot$tot)]),
    theme = factor(theme, levels = theme_tot$theme[order(-theme_tot$tot)])
  )

# wrap long theme strip labels so the right-hand strips stay narrow
levels(plot_df$theme) <- str_wrap(levels(plot_df$theme), width = 18)

# ---- plot -------------------------------------------------------------------
p <- ggplot(plot_df, aes(x = proportion, y = group, fill = traits)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.1),
           width = 0.8) +
  facet_grid(rows = vars(theme), scales = "free_y", space = "free_y") +
  scale_fill_manual(values = trait_palette, name = "Trait",
                    drop = FALSE) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Proportion of trait's enriched terms",
       y = NULL,
       title = "Simplified GO terms by group",
       subtitle = "Each trait normalised to its own total") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.y.right = element_text(angle = 0, hjust = 0, size = 8),
    strip.background   = element_rect(fill = "grey92", colour = NA),
    panel.spacing.y    = unit(2, "pt"),
    legend.position    = "top",
    plot.title.position = "plot"
  )

# ---- save -------------------------------------------------------------------
ggsave(outfile, p, width = 9, height = 12, bg = "white", device = "pdf")
cat("wrote", outfile, "\n")
