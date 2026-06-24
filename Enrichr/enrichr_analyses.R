# Enrichr Analyses

setwd("~/Desktop/PhD/First_project/FLAMES_2.0")

library(devtools)
library(enrichR)

websiteLive <- getOption("enrichR.live")

if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
}

if (websiteLive) {
  dbs <- listEnrichrDbs()
  head(dbs)
}

dbs <- listEnrichrDbs()
dbs <- dbs[dbs$libraryName %in% c("OMIM_Disease",
                                  "GO_Biological_Process_2026",
                                  "GO_Cellular_Component_2026",
                                  "GO_Molecular_Function_2026",
                                  "DGIdb_Drug_Targets_2024",
                                  "KEGG_2026",
                                  "ClinVar_2025"
                                  ),]

## Without background genes

cortical_genes <- read.table("./Cortical/FLAMES_v2_wo_QSNP.tsv", header = T)
trunk_genes <- read.table("./TrunkWidths/FLAMES_v2_wo_QSNP.tsv", header = T)
limb_genes <- read.table("./LimbHeight/FLAMES_v2_wo_QSNP.tsv", header = T)

setwd("./../Enrichr/")

data(background)
length(background)

enriched_cortical_background <- enrichr(unique(cortical_genes$symbol), dbs$libraryName, background = background)
enriched_trunk_background <- enrichr(unique(trunk_genes$symbol), dbs$libraryName, background = background)
enriched_limb_background <- enrichr(unique(limb_genes$symbol), dbs$libraryName, background = background)

## Explore

#GOBP

cortical_BP <- enriched_cortical_background$GO_Biological_Process_2026
cortical_BP <- cortical_BP[cortical_BP$Adjusted.P.value <= 0.05,]
write.table(cortical_BP, "Cortical_GOBP.tsv", sep = "\t", col.names = T)

trunk_BP <- enriched_trunk_background$GO_Biological_Process_2026
trunk_BP <- trunk_BP[trunk_BP$Adjusted.P.value <= 0.05,]
write.table(trunk_BP, "TrunkWidths_GOBP.tsv", sep = "\t", col.names = T)

limb_BP <- enriched_limb_background$GO_Biological_Process_2026
limb_BP <- limb_BP[limb_BP$Adjusted.P.value <= 0.05,]
write.table(limb_BP, "LimbHeight_GOBP.tsv", sep = "\t", col.names = T)


#GOCC

cortical_CC <- enriched_cortical_background$GO_Cellular_Component_2026
cortical_CC <- cortical_CC[cortical_CC$Adjusted.P.value <= 0.05,]
write.table(cortical_CC, "Cortical_GOCC.tsv", sep = "\t", col.names = T)

trunk_CC <- enriched_trunk_background$GO_Cellular_Component_2026
trunk_CC <- trunk_CC[trunk_CC$Adjusted.P.value <= 0.05,]
write.table(trunk_CC, "TrunkWidths_GOCC.tsv", sep = "\t", col.names = T)

limb_CC <- enriched_limb_background$GO_Cellular_Component_2026
limb_CC <- limb_CC[limb_CC$Adjusted.P.value <= 0.05,]
write.table(limb_CC, "LimbHeight_GOCC.tsv", sep = "\t", col.names = T)

#GOMF

cortical_MF <- enriched_cortical_background$GO_Molecular_Function_2026
cortical_MF <- cortical_MF[cortical_MF$Adjusted.P.value <= 0.05,]
write.table(cortical_MF, "Cortical_GOMF.tsv", sep = "\t", col.names = T)

trunk_MF <- enriched_trunk_background$GO_Molecular_Function_2026
trunk_MF <- trunk_MF[trunk_MF$Adjusted.P.value <= 0.05,]
write.table(trunk_MF, "TrunkWidths_GOMF.tsv", sep = "\t", col.names = T)

limb_MF <- enriched_limb_background$GO_Molecular_Function_2026
limb_MF <- limb_MF[limb_MF$Adjusted.P.value <= 0.05,]
write.table(limb_MF, "LimbHeight_GOMF.tsv", sep = "\t", col.names = T)

#OMIM

cortical_OMIM <- enriched_cortical_background$OMIM_Disease
cortical_OMIM <- cortical_OMIM[cortical_OMIM$Adjusted.P.value <= 0.05,]
write.table(cortical_OMIM, "Cortical_OMIM.tsv", sep = "\t", col.names = T)

trunk_OMIM <- enriched_trunk_background$OMIM_Disease
trunk_OMIM <- trunk_OMIM[trunk_OMIM$Adjusted.P.value <= 0.05,]
write.table(trunk_OMIM, "TrunkWidths_OMIM.tsv", sep = "\t", col.names = T)

limb_OMIM <- enriched_limb_background$OMIM_Disease
limb_OMIM <- limb_OMIM[limb_OMIM$Adjusted.P.value <= 0.05,]
write.table(limb_OMIM, "LimbHeight_OMIM.tsv", sep = "\t", col.names = T)

#KEGG

cortical_KEGG <- enriched_cortical_background$KEGG_2026
cortical_KEGG <- cortical_KEGG[cortical_KEGG$Adjusted.P.value <= 0.05,]
write.table(cortical_KEGG, "Cortical_KEGG.tsv", sep = "\t", col.names = T)

trunk_KEGG <- enriched_trunk_background$KEGG_2026
trunk_KEGG <- trunk_KEGG[trunk_KEGG$Adjusted.P.value <= 0.05,]
write.table(trunk_KEGG, "TrunkWidths_KEGG.tsv", sep = "\t", col.names = T)

limb_KEGG <- enriched_limb_background$KEGG_2026
limb_KEGG <- limb_KEGG[limb_KEGG$Adjusted.P.value <= 0.05,]
write.table(limb_KEGG, "LimbHeight_KEGG.tsv", sep = "\t", col.names = T)

#DGIdb

cortical_DGIdb <- enriched_cortical_background$DGIdb_Drug_Targets_2024
cortical_DGIdb <- cortical_DGIdb[cortical_DGIdb$Adjusted.P.value <= 0.05,]
write.table(cortical_DGIdb, "Cortical_DGIdb.tsv", sep = "\t", col.names = T)

trunk_DGIdb <- enriched_trunk_background$DGIdb_Drug_Targets_2024
trunk_DGIdb <- trunk_DGIdb[trunk_DGIdb$Adjusted.P.value <= 0.05,]
write.table(trunk_DGIdb, "TrunkWidths_DGIdb.tsv", sep = "\t", col.names = T)

limb_DGIdb <- enriched_limb_background$DGIdb_Drug_Targets_2024
limb_DGIdb <- limb_DGIdb[limb_DGIdb$Adjusted.P.value <= 0.05,]
write.table(limb_DGIdb, "LimbHeight_DGIdb.tsv", sep = "\t", col.names = T)

#Clinvar

cortical_clinvar <- enriched_cortical_background$ClinVar_2025
cortical_clinvar <- cortical_clinvar[cortical_clinvar$Adjusted.P.value <= 0.05,]
write.table(cortical_clinvar, "Cortical_clinvar.tsv", sep = "\t", col.names = T)

trunk_clinvar <- enriched_trunk_background$ClinVar_2025
trunk_clinvar <- trunk_clinvar[trunk_clinvar$Adjusted.P.value <= 0.05,]
write.table(trunk_clinvar, "TrunkWidths_clinvar.tsv", sep = "\t", col.names = T)

limb_clinvar <- enriched_limb_background$ClinVar_2025
limb_clinvar <- limb_clinvar[limb_clinvar$Adjusted.P.value <= 0.05,]
write.table(limb_clinvar, "LimbHeight_clinvar.tsv", sep = "\t", col.names = T)

## Read again

cortical_BP <- read.table("./Cortical_GOBP.tsv", sep = "\t", header =T)
trunk_BP<- read.table("./TrunkWidths_GOBP.tsv", sep = "\t", header =T)
limb_BP <- read.table("./LimbHeight_GOBP.tsv", sep = "\t", header =T)

all_unique <- unique(c(cortical_BP$Term, trunk_BP$Term, limb_BP$Term))
writeLines(all_unique, "./all_sig_BP_terms.tsv")

groups <- read.table("GO_BP_term_group_mapping.tsv", sep = "\t", header = T)



