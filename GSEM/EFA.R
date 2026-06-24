## GSEM creation

setwd("~/Desktop/PhD/First_project/gsem/")

library(GenomicSEM)
library(EGAnet)
library(tidyverse)
library(reshape2)
library(psych)
library(lavaan) 
library(pls)

load("LDSCoutput.RData")

S <- LDSCoutput$S
V <- LDSCoutput$V
print("Genetic covariance matrix (S):")
print(round(S,4))
rG <- cov2cor(S)
print("Genetic correlations (cov2cor):")
print(round(rG,3))

# ---- 3. Exploratory inspection ----
# Heatmap of rG

rG_df <- melt(rG)
ggplot(rG_df, aes(Var1, Var2, fill = value)) +
  geom_tile() + theme_minimal() + labs(title="Genetic correlation matrix")

# Clustering (dendrogram) on absolute correlations
hc <- hclust(as.dist(1 - abs(rG)))

pdf("hc_exploration.pdf")
plot(hc, main="Hierarchical clustering (1 - |rG|)", xlab = "", )
dev.off()


paLDSC(S = rG, V = V, r = 1000, save.pdf = T, fa = T)


# Wanna look at item stability
# Calculate harmonic mean first

trait_ns <- c( N_SA = 53275.7, N_FI = 53225.09, N_IC = 53094.14, N_height = 458303,
               N_Hip_width = 53000, N_Shoulder_width= 53000, N_Torso_length= 53000, 
               N_Average_Femur = 53000, N_Average_Forearm = 53000, N_Average_Humerus = 53000,
               N_Average_Tibi = 53000)

harmonic_mean_n <- length(trait_ns) / sum(1 / trait_ns)
harmonic_mean_n

# median and mean for comparison
median_n <- median(trait_ns)
mean_n   <- mean(trait_ns)

cat("harmonic:", round(harmonic_mean_n), "\n",
    "median: ", median_n, "\n",
    "mean:   ", round(mean_n), "\n")

set.seed(123)
boot_ega <- bootEGA(rG, model = "glasso", n = round(harmonic_mean_n), iter = 1000)

boot_ega$stability$item.stability
boot_ega$stability$dimension.stability

# Want to try with different Ns Because n.obs is an approximation for rG matrices,
# check if conclusions change median_n, round(mean_n), min(trait_ns))

set.seed(123)
boot_ega_median <- bootEGA(rG, model = "glasso", n = median_n, iter = 1000)

boot_ega_median$stability$item.stability
boot_ega_median$stability$dimension.stability

boot_ega_mean <- bootEGA(rG, model = "glasso", n = round(mean_n), iter = 1000)

boot_ega_mean$stability$item.stability
boot_ega_mean$stability$dimension.stability

boot_ega_min <- bootEGA(rG, model = "glasso", n = min(trait_ns), iter = 1000)

boot_ega_min$stability$item.stability
boot_ega_min$stability$dimension.stability

pdf("boot_ega.pdf", width = 11.71, height = 5.85)
plot(boot_ega, plot.type = "structure")
dev.off()

## EFA with only odd numbered chromosomes

load("LDSCoutput_odd.rda")

S <- LDSCoutput_odd$S
V <- LDSCoutput_odd$V
print("Genetic covariance matrix (S):")
print(round(S,4))
rG <- cov2cor(S)
print("Genetic correlations (cov2cor):")
print(round(rG,3))

# Try EFA for 1:3 factors to see loadings (informal)
efa2 <- fa(rG, nfactors = 2, fm = "ml", n.obs = round(harmonic_mean_n), n.iter = 1000)
efa3 <- fa(rG, nfactors = 3, fm = "ml", n.obs = round(harmonic_mean_n), n.iter = 1000)
efa4 <- fa(rG, nfactors = 4, fm = "ml", n.obs = round(harmonic_mean_n), n.iter = 1000)

print(efa2$loadings)
print(efa3$loadings)
print(efa4$loadings)

save(efa2, file = "./EFA_odd_2.rda")
save(efa3, file = "./EFA_odd_3.rda")
save(efa4, file = "./EFA_odd_4.rda")


