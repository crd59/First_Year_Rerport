## Run multlivariate GWAS

library(GenomicSEM)
library(tidyverse)
library(lavaan)


setwd("~/rds/rds-genetics_hpc-Nl99R8pHODQ/UKB/Imaging_genetics/crd59/multivariate_GWAS/GWAS_height_new/")

######### LOAD LDSC OUTPUT #######

# at this point use all chromosomes

load("./LDSCoutput_height.RData")

############# Run sumstats function ##########

merged_sumstats <- read.table("./create_sumstat_folder/merged_sumstats_hapmap_snps.tsv", sep = "\t", header = T)

#specify the model
model<-'
# Cortical factor
Cortical =~ FI + IC + SA

# Skeletal limbs
Skeletal_LimbHeight =~ Average_Femur + Average_Forearm + Average_Humerus + Average_Tibia + Height

# Skeletal trunk/widths
Skeletal_TrunkWidths =~ Hip_width + Shoulder_width + Torso_length + Height

# Allow factors to correlate
Cortical ~~ Skeletal_LimbHeight
Cortical ~~ Skeletal_TrunkWidths
Skeletal_LimbHeight ~~ Skeletal_TrunkWidths

Cortical ~ SNP
Skeletal_LimbHeight ~ SNP
Skeletal_TrunkWidths ~ SNP
'
#run the multivariate GWAS using parallel processing

# Do 3 thousand SNPs at a time. otherwise it doesn't run

for(iteration in 1:nrow(merged_sumstats)){

it_SNPs <- unique(merged_sumstats[1+((iteration-1)*5000):(iteration*5000),])
it_SNPs <- it_SNPs[complete.cases(it_SNPs),]

CorrelatedFactors <- userGWAS(covstruc = LDSCoutput,
                             SNPs = it_SNPs,
                             estimation = "DWLS",
                             model = model,
                             std.lv = T,
                             printwarn = TRUE,
                             sub=c("Cortical~SNP", "Skeletal_LimbHeight~SNP", "Skeletal_TrunkWidths~SNP"),
                             toler = FALSE,
                             SNPSE = FALSE,
                             parallel = F,
                             GC="standard",
                             MPI=FALSE,
                             smooth_check=TRUE,
                             fix_measurement=TRUE,
                             Q_SNP=TRUE)

write.table(CorrelatedFactors[[1]], paste0("./Cortical_GWAS/Cortical_GWAS_hapmap_part", iteration,".tsv"), sep = "\t", col.names = T)
write.table(CorrelatedFactors[[2]], paste0("./LimbHeight_GWAS/LimbHeight_GWAS_hapmap_part", iteration,".tsv"), sep = "\t", col.names = T)
write.table(CorrelatedFactors[[3]], paste0("./TrunkWidths_GWAS/TrunkWidths_GWAS_hapmap_part", iteration,".tsv"), sep = "\t", col.names = T)
}
