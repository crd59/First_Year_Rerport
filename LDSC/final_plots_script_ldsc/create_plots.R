## Important traits H2 plots and RG plots

setwd("~/Desktop/PhD/First_project/final_plots_script_ldsc/")

library(ggplot2)
library(tidyverse)

individual <- read.table("./individual_heritabilities.txt", sep = "\t", header = T)

individual$Trait[11] <- "Height"

individual$Trait <- factor(individual$Trait, levels = individual$Trait)


  # Skeletal plot
  
  skeletal_h2_plot <- ggplot(individual, aes(x = Trait, y = h2)) +
  geom_col(color = "black", width = 0.8, fill = "steelblue") +
  geom_errorbar(aes(ymin = h2 - SE_h2, ymax = h2 + SE_h2),
                width = 0.2, linewidth = 0.6) +
  scale_fill_discrete(name = "Trait") +   # different colour for each bar
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",  # remove legend since each bar already labelled
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  ) +
  labs(
    x = "",
    y = "Heritability (h²)",
    title = "Trait Heritability"
  ) +
  expand_limits(y = 0)

ggsave(plot = skeletal_h2_plot, filename = "./important_h2_plot.pdf", device = "pdf",
       width = 7.85,
       height= 6.89)

#### Create super long list 

long_list <- read.table("./Pairwise_rg.txt",header =T, strip.white = T)

long_list$p1 <- gsub(x = long_list$p1, pattern = "./../important_sumstats/", replacement = "")
long_list$p2 <- gsub(x = long_list$p2, pattern = "./../important_sumstats/", replacement = "")
long_list$p1 <- gsub(x = long_list$p1, pattern = ".sumstats.gz", replacement = "")
long_list$p2 <- gsub(x = long_list$p2, pattern = ".sumstats.gz", replacement = "")
long_list$p1 <- gsub(x = long_list$p1, pattern = "./../new_GWAS_ldsc/", replacement = "")

long_list$p1[long_list$p1 == "avg_femur"] <- "Average_Femur"
long_list$p1[long_list$p1 == "avg_forearm"] <-"Average_Forearm"
long_list$p1[long_list$p1 == "avg_humerus"] <- "Average_Humerus"
long_list$p1[long_list$p1 == "avg_tibia"] <- "Average_Tibia"
long_list$p1[long_list$p1 == "global_FI"] <- "FI"
long_list$p1[long_list$p1 == "new_height"] <- "Height"
long_list$p1[long_list$p1 == "hip_width"] <- "Hip_width"
long_list$p1[long_list$p1 == "global_IC"] <- "IC"
long_list$p1[long_list$p1 == "global_SA"] <- "SA"
long_list$p1[long_list$p1 == "shoulder_width"] <- "Shoulder_Width"

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

traits<-c("FI", "IC",  "SA", 
          "Height", "Average_Femur",
          "Average_Forearm", "Average_Humerus", "Average_Tibia", "Hip_width",
          "Shoulder_Width", "Torso_Length")

trait_order <- traits[c(4,5,6,7,8,1,2,3,9,10)]

trait_order_p2 <- traits[c(11,10,9,3,2,1,8,7,6,5)]


long_list$p1 <- factor(long_list$p1, levels = trait_order)
long_list$p2 <- factor(long_list$p2, levels = trait_order_p2)

rg_plot <- ggplot(long_list, aes(x = p1, y = p2, fill = rg)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(-1, 1),
    name = "rg"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  ) +   
  geom_text(aes(label = rg),
                  size  = 6
  ) +
  labs(title = "Genetic Correlation Matrix")

ggsave(plot = rg_plot, device = "pdf", filename = "./rg_plot_final.pdf")

write.table(long_list, "./pairwise_rg_important_traits.tsv", sep = "\t", col.names = T, row.names = F)


