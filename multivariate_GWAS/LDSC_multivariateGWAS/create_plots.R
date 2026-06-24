## Important traits H2 plots and RG plots

setwd("/rds/project/rds-Nl99R8pHODQ/UKB/Imaging_genetics/crd59/multivariate_GWAS/GWAS_height_new/LDSC/")

library(ggplot2)
library(tidyverse)

#### Create super long list 

long_list <- read.table("./Pairwise_rg.txt",header =T, strip.white = T)

long_list$p1 <- gsub(x = long_list$p1, pattern = "./", replacement = "")
long_list$p2 <- gsub(x = long_list$p2, pattern = "./../../../ldsc_analyses/important_sumstats/", replacement = "")
long_list$p1 <- gsub(x = long_list$p1, pattern = ".sumstats.gz", replacement = "")
long_list$p2 <- gsub(x = long_list$p2, pattern = ".sumstats.gz", replacement = "")
long_list$p2 <- gsub(x = long_list$p2, pattern = "./", replacement = "")

long_list$p1[long_list$p1 == "cortical_munged"] <- "Cortical"
long_list$p1[long_list$p1 == "limb_munged"] <-"LimbHeights"
long_list$p1[long_list$p1 == "trunk_munged"] <- "TrunkWidths"

long_list$p2[long_list$p2 == "torso_length"] <- "Torso_Length"
long_list$p2[long_list$p2 == "avg_forearm"] <- "Average_Forearm"
long_list$p2[long_list$p2 == "avg_femur"] <- "Average_Femur"
long_list$p2[long_list$p2 == "avg_humerus"] <- "Average_Humerus"
long_list$p2[long_list$p2 == "avg_tibia"] <- "Average_Tibia"
long_list$p2[long_list$p2 == "global_FI"] <- "FI"
long_list$p2[long_list$p2 == "new_height"] <- "Height"
long_list$p2[long_list$p2 == "hip_width"] <- "Hip_width"
long_list$p2[long_list$p2 == "global_IC"] <- "IC"
long_list$p2[long_list$p2 == "global_SA"] <- "SA"
long_list$p2[long_list$p2 == "shoulder_width"] <- "Shoulder_Width"
long_list$p2[long_list$p2 == "cortical_munged"] <- "Cortical"
long_list$p2[long_list$p2 == "limb_munged"] <-"LimbHeights"
long_list$p2[long_list$p2 == "trunk_munged"] <- "TrunkWidths"

long_list <- long_list %>%
  mutate(stars = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  ))

multivariate_traits <- long_list[c(1,2,14),]

rg_plot_multivariate <- ggplot(multivariate_traits, aes(x = p1, y = p2, fill = rg)) +
  geom_tile(color = "white") +
  geom_text(aes(label = stars), size = 5, color = "black") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(-1, 1),
    name = "rg"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  ) +
  labs(title = "Genetic Correlation Matrix")

ggsave(plot = rg_plot_multivariate, device = "png", filename = "./rg_plot_multivariate.png")

original_traits <- long_list[-c(1,2,14),]

traits<-c("FI", "IC", "SA", 
          "Average_Femur",
          "Average_Forearm", "Average_Humerus", "Average_Tibia", "Height", "Hip_width",
          "Shoulder_Width", "Torso_Length")

trait_order_p2 <- traits[c(1:11)]


original_traits$p2 <- factor(original_traits$p2, levels = trait_order_p2)

rg_plot <- ggplot(original_traits, aes(x = p1, y = p2, fill = rg)) +
  geom_tile(color = "white") +
  geom_text(aes(label = stars), size = 5, color = "black") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(-1, 1),
    name = "rg"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  ) +
  labs(title = "Genetic Correlation Matrix")

ggsave(plot = rg_plot, device = "png", filename = "./rg_plot_final.png")

write.table(long_list, "./pairwise_rg_important_traits.tsv", sep = "\t", col.names = T, row.names = F)


