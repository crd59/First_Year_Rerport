## Create Sumstats File

library(GenomicSEM)
library(tidyverse)
library(psych)
library(lavaan) 

setwd("~/rds/rds-genetics_hpc-Nl99R8pHODQ/UKB/Imaging_genetics/crd59/multivariate_GWAS/GWAS_height_new/")

######### LOAD LDSC OUTPUT #######

# at this point use all chromosomes

load("./LDSCoutput_height.RData")

WASsumstats_files <- c('./../femur.tsv','./../forearm.tsv', './../humerus.tsv',
                       './../tibia.tsv','./../FI.tsv',
                       './../IC.tsv','./../SA.tsv',
                       './../hip_width.tsv', "./../height.tsv",
                       './../shoulder_width.tsv','./../torso.tsv')

# Reorder to fit LDSC order

trait.names<-c("Average_Femur", "Average_Forearm", "Average_Humerus", "Average_Tibia",
               "FI", "IC", "SA", 
               "Hip_width", "Height", "Shoulder_width", "Torso_length")

ref_snplist <- "./ref_rest.txt"  

## All continuous, no logistic scale

se.logit <- c(F,F,F, F, F, F, F, F, F,F,F)

OLS <- c(T,T,T,T,T,T,T,T,T,T,T)

linprob <- c(F,F,F,F,F,F,F,F,F,F,F)

## calculate N for cortical and height, all skeletal are 53000

N <- c(53000,53000,53000,53000, 53000,53000,53000, 53000, 458303,53000,53000)


sumstats <- sumstats(files=WASsumstats_files,
                                 ref=ref_snplist,
                                 trait.names=trait.names,
                                 se.logit=se.logit,
                                 OLS = OLS,
                                 linprob= linprob,
                                 N=N,
                                 betas=NULL,
                                 maf.filter=0.01,
                                 keep.indel=FALSE,
                                 parallel=TRUE)
write.table(sumstats_1_3, "./create_sumstat_folder/sumstats.tsv", sep = "\t", col.names = T)
