## Run LDSC for the new GWAS

setwd("~/rds/rds-genetics_hpc-Nl99R8pHODQ/UKB/Imaging_genetics/crd59/multivariate_GWAS/GWAS_height_new/LDSC")

library(GenomicSEM)
library(tidyverse)

# Read all datasets

#vector of munged summary statisitcs
traits<-c(list.files("./../../../ldsc_analyses/important_sumstats/", full.names = T))
traits <- traits[-c(6,11)]

traits <- append(traits, "./trunk_munged.sumstats.gz")
traits <- append(traits, "./limb_munged.sumstats.gz")
traits <- append(traits, "./cortical_munged.sumstats.gz")

#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev<-c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,NA)

#vector of population prevalences
population.prev<-c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

#the folder of LD scores
ld<-"/home/crd59/rds/rds-genetics_hpc-Nl99R8pHODQ/ldsc/eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"/home/crd59/rds/rds-genetics_hpc-Nl99R8pHODQ/ldsc/eur_w_ld_chr/"

trait.names<-c("Average_Femur", "Average_Forearm", "Average_Humerus", "Average_Tibia",
               "FI", "IC", "SA", 
               "Hip_width", "Height", "Shoulder_width", "Torso_length", "Trunk_Widths",
               "Limb_Lengths", "Cortical")

#run LDSC
LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names)

#optional command to save the output as a .RData file for later use
save(LDSCoutput,file="./LDSCoutput_multivariate.RData")





