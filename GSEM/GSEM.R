# Actually do GenomicSEM

library(GenomicSEM)
library(EGAnet)
library(tidyverse)
library(reshape2)
library(psych)
library(lavaan) 
library(semPlot)

### ASSUMPTIONS (replace these)
# - ldsc_out exists (result of ldsc(...)) containing $S and $V
# - sumstats_files is a character vector of your 7 munged sumstats files (gz)
# - trait_names vector matches order used in ldsc_out
# - reference snplist path set for commonfactor

setwd("~/Desktop/PhD/First_project/gsem/")


#vector of munged summary statisitcs
traits<-c(list.files("./../important_sumstats/", full.names = T))

trait.names<-c("Average_Femur", "Average_Forearm", "Average_Humerus", "Average_Tibia",
               "FI", "IC", "SA", 
               "Hip_width", "Height", "Shoulder_width", "Torso_length")

# only even chrs
load("LDSCoutput_even.rda")

################### commonmodel #################

common_model <- commonfactor(covstruc = LDSCoutput_even, estimation = "DWLS")

########################### TWO-FACTOR MODEL ###########################

model_3f2 <- '
Skeletal =~ Average_Femur + Average_Forearm + Average_Humerus + Average_Tibia + Hip_width + Shoulder_width + Torso_length + Height
Cortical =~ FI + SA + IC
Skeletal ~~ Cortical
'

fit_3f2 <- usermodel(covstruc = LDSCoutput_even, model = model_3f2, estimation = "DWLS", std.lv = T)

save(x = fit_3f2, file = "./2_factor_GSEM.rda")

write.csv(fit_3f2$results, "./2_factor_model.csv")

############################ THREE FACTOR MODEL ###########################


model_3f31 <- '
# Cortical factor
Cortical =~ FI + IC + SA

# Skeletal - limbs + height
Skeletal_LimbHeight =~ Average_Femur + Average_Forearm + Average_Humerus + Average_Tibia + Height

# Skeletal - trunk/widths
Skeletal_TrunkWidths =~ Hip_width + Shoulder_width + Torso_length + Height

# Allow factors to correlate
Cortical ~~ Skeletal_LimbHeight
Cortical ~~ Skeletal_TrunkWidths
Skeletal_LimbHeight ~~ Skeletal_TrunkWidths
'

fit_3f31 <- usermodel(covstruc = LDSCoutput_even, model = model_3f31, estimation = "DWLS", std.lv = T)

save(x = fit_3f31, file = "./3_factor_GSEM.rda")

write.csv(fit_3f31$results, "./3_factor_model.csv")

############################ FOUR FACTOR MODEL ###########################


model_4f31 <- '
# Cortical factor
Cortical =~ FI + IC + SA

# Skeletal - limbs + height
Skeletal_LimbHeight =~ Average_Femur + Average_Tibia + Average_Forearm + Average_Humerus + Height
Skeletal_TrunkWidths =~ Hip_width + Shoulder_width + Torso_length + Height
Skeletal_residual =~ Average_Forearm + Shoulder_width + Height

# Allow factors to correlate
Cortical ~~ Skeletal_LimbHeight
Cortical ~~ Skeletal_TrunkWidths
Cortical ~~ Skeletal_residual

Skeletal_TrunkWidths ~~ Skeletal_LimbHeight
Skeletal_TrunkWidths ~~ Skeletal_residual

Skeletal_LimbHeight ~~ Skeletal_residual
'

fit_4f31 <- usermodel(covstruc = LDSCoutput_even, model = model_4f31, estimation = "DWLS", std.lv = T)

save(x = fit_3f31, file = "./3_factor_GSEM.rda")

write.csv(fit_3f31$results, "./3_factor_model.csv")

