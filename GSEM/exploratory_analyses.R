library(GenomicSEM)
library(tidyverse)

setwd("~/Desktop/PhD/First_project/gsem/")

#vector of munged summary statisitcs
traits<-c(list.files("./../important_sumstats/", full.names = T))

#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev<-c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,NA)

#vector of population prevalences
population.prev<-c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,NA)

#the folder of LD scores
ld<-"~/Desktop/PhD/First_project/eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"~/Desktop/PhD/First_project/eur_w_ld_chr/"

trait.names<-c("Average_Femur", "Average_Forearm", "Average_Humerus", "Average_Tibia",
               "FI", "IC", "SA", 
               "Hip_width", "Height", "Shoulder_width", "Torso_length")

#run LDSC
LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names)

#optional command to save the output as a .RData file for later use
save(LDSCoutput,file="./LDSCoutput.RData")

### Split by chromosomes

#Even 

LDSCoutput_even<-ldsc(traits=traits,
                 sample.prev=sample.prev,
                 population.prev=population.prev,
                 ld=ld,
                 wld=wld,
                 trait.names=trait.names, 
                 select = "EVEN")

save(file = "./LDSCoutput_even.rda", LDSCoutput_even)

#Odd

LDSCoutput_odd<-ldsc(traits=traits,
                      sample.prev=sample.prev,
                      population.prev=population.prev,
                      ld=ld,
                      wld=wld,
                      trait.names=trait.names, select = "ODD")

save(file = "./LDSCoutput_odd.rda", LDSCoutput_odd)

