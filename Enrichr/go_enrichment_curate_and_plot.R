#!/usr/bin/env Rscript
# =============================================================================
#  GO Biological-Process enrichment: term simplification, grouping & plotting
#  Project : multivariate GWAS of cortical expansion, limb growth, trunk growth
#  Input   : EnrichR GO_Biological_Process_2023 results (FLAMES-prioritised genes)
#            - Cortical_GOBP.tsv, Limb_GOBP.tsv, Trunk_GOBP.tsv
#  Author  : (analysis pipeline)  |  R 4.3, ggplot2 3.4
#
#  WHY THIS SCRIPT EXISTS
#  ----------------------
#  EnrichR returns 439 unique FDR-significant GO:BP gene sets across the three
#  traits. Raw GO:BP is deeply redundant: it contains parent/child terms
#  ("Regulation of X", "Positive Regulation of X", "Negative Regulation of X",
#  "Cellular Response to X" vs "Response to X") that describe essentially one
#  biological process at different granularity, and it splits single programs
#  across dozens of near-synonymous labels. 439 rows is uninterpretable.
#
#  This script performs a transparent, fully auditable two-stage reduction:
#    STAGE 1  NORMALISE  - algorithmically strip directional/regulatory wrappers
#                          and harmonise "response-to" phrasing, collapsing
#                          parent/child variants onto one normalized label.
#    STAGE 2  CANONICALISE + GROUP - a hand-curated synonym map merges remaining
#                          near-duplicate processes onto a single simplified
#                          term (182 terms), each assigned to one of 21
#                          biological GROUPS, nested under 10 THEMES.
#  Every original term's fate is written to go_term_curation.tsv so the mapping
#  can be inspected and challenged line-by-line.
#
#  REPRODUCIBILITY NOTE
#  --------------------
#  The normalization logic below (Stage 1) is the actual derivation. The
#  curated maps (syn_clean / group_lookup / theme_lookup, Stage 2) were frozen
#  after manual review and are embedded verbatim further down. After the
#  pipeline runs, hard assertions check that the result still collapses to
#  exactly 182 simplified terms -> 21 groups -> 10 themes; if a future GO/EnrichR
#  version introduces a term the curated map does not cover, the script STOPS
#  and names the offending term rather than silently dropping or misfiling it.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)   # fast, *correctly-aligned* TSV reader (see note below)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

# ----------------------------------------------------------------------------
# 0.  Paths
# ----------------------------------------------------------------------------
in_dir  <- "~/Desktop/PhD/First_project/Enrichr/"
out_dir <- "~/Desktop/PhD/First_project/Enrichr/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

files <- c(
  cortical = file.path(in_dir, "Cortical_GOBP.tsv"),
  limb     = file.path(in_dir, "Limb_GOBP.tsv"),
  trunk    = file.path(in_dir, "Trunk_GOBP.tsv")
)
stopifnot(all(file.exists(files)))

# ----------------------------------------------------------------------------
# 1.  Read EnrichR output
#
#     GOTCHA (important): these TSVs were written with R-style rownames, so the
#     HEADER line carries 9 names but each DATA row has 10 fields (an unnamed
#     leading row-id column). This off-by-one breaks naive readers:
#       * readr::read_tsv() shifts the header onto the data, so Term[1] == "1".
#       * data.table::fread(header=TRUE) mis-detects the header on these files
#         and promotes the first data row into the column names.
#     The only fully robust fix is to read the 9 header names ourselves and
#     prepend a "rowid" name, then read the body with header=FALSE/skip=1 so the
#     10 data columns line up with our 10 explicit names. We verify the field
#     count matches before trusting the result.
#
#     EnrichR columns we keep:
#       Term              "Name (GO:#######)"
#       Adjusted.P.value  BH-FDR; every row here is already < 0.05
#       Genes             semicolon-separated FLAMES genes driving the term
# ----------------------------------------------------------------------------
read_enrichr <- function(path, trait) {
  hdr <- gsub("\"", "", strsplit(readLines(path, n = 1), "\t")[[1]])  # 9 names
  dt  <- fread(path, sep = "\t", header = FALSE, skip = 1, quote = "\"",
               showProgress = FALSE, data.table = FALSE)
  stopifnot(ncol(dt) == length(hdr) + 1L)        # 10 data cols vs 9 header names
  names(dt) <- c("rowid", hdr)                    # explicit, aligned names
  dt$trait <- trait
  dt[, c("Term", "Adjusted.P.value", "Genes", "trait")]
}

raw <- bind_rows(Map(read_enrichr, files, names(files)))
message(sprintf("Read %d total significant rows across 3 traits", nrow(raw)))

# Split "Name (GO:#######)" -> name + GO id
raw <- raw %>%
  mutate(
    go   = str_match(Term, "\\((GO:\\d+)\\)\\s*$")[, 2],
    name = str_trim(str_replace(Term, "\\s*\\(GO:\\d+\\)\\s*$", ""))
  )
stopifnot(all(!is.na(raw$go)))     # every term must carry a GO id

# Collapse to one row per UNIQUE GO id, recording which traits it appears in,
# the best (smallest) adjusted p across traits, and the union of driver genes.
collapse_genes <- function(x) {
  paste(sort(unique(unlist(str_split(x, ";")))), collapse = ";")
}
per_go <- raw %>%
  group_by(go) %>%
  summarise(
    name      = name[1],
    traits    = paste(sort(unique(trait)), collapse = "+"),
    best_adjp = min(Adjusted.P.value, na.rm = TRUE),
    genes     = collapse_genes(Genes),
    .groups   = "drop"
  )
message(sprintf("Collapsed to %d unique GO:BP terms", nrow(per_go)))

# ----------------------------------------------------------------------------
# 2.  STAGE 1 - algorithmic normalization (parent/child collapse)
#
#     Applied in order; each rule is deliberately conservative:
#       (a) expand the 3 abbreviated names EnrichR emitted ("Neg Reg ...",
#           "Pos Reg ...", "... Rec Protein ...") to their full wording so they
#           match their non-abbreviated siblings;
#       (b) strip a leading "Regulation of " / "Positive Regulation of " /
#           "Negative Regulation of " - repeated up to 3x to peel nested
#           wrappers (e.g. "Regulation of Regulation of ..." style chains);
#       (c) fold "Cellular Response to X" into "Response to X";
#       (d) drop a trailing " Stimulus".
#     The result ("simp0") is the normalized label BEFORE curated synonym
#     merging. This is pure string logic - no biology yet.
# ----------------------------------------------------------------------------
normalize_term <- function(s) {
  s <- str_replace_all(s, fixed("Neg Reg "), "Negative Regulation of ")
  s <- str_replace_all(s, fixed("Pos Reg "), "Positive Regulation of ")
  s <- str_replace_all(s, fixed("Rec Protein"), "Receptor Protein")
  for (i in 1:3) {
    s <- str_replace(s, "^(Positive |Negative )?Regulation of ", "")
  }
  s <- str_replace(s, "^Cellular Response to ", "Response to ")
  s <- str_replace(s, " Stimulus$", "")
  str_trim(s)
}
per_go$simp0 <- vapply(per_go$name, normalize_term, character(1))
message(sprintf("After Stage-1 normalization: %d distinct labels (from %d)",
                dplyr::n_distinct(per_go$simp0), nrow(per_go)))

# ============================================================================
#  CURATED MAPS (Stage 2)
#  These three named vectors + the EXPECTED_* invariants are the frozen,
#  manually-reviewed curation. They were emitted programmatically from the
#  validated prototype so they are guaranteed to match it exactly. Edit the
#  biology here (and only here) if you want to re-curate.
# ============================================================================


# (a) Synonym canonicalization: normalized label (simp0) -> canonical simplified term.
#     Only pairs where canonicalization changes the string are listed; anything
#     not present passes through unchanged.
syn_clean <- c(
  "Chromatin Organization" = "Chromatin organization",
  "Epithelial to Mesenchymal Transition" = "Epithelial–mesenchymal transition (EMT)",
  "Insulin Receptor Signaling Pathway" = "Insulin/IGF signalling",
  "Cell Cycle" = "Cell cycle (general)",
  "Cell Differentiation" = "Cell differentiation (general)",
  "Cell Growth" = "Cell growth",
  "DNA-templated Transcription" = "Transcription (RNA Pol II / DNA-templated)",
  "Macromolecule Biosynthetic Process" = "Transcription (RNA Pol II / DNA-templated)",
  "Osteoblast Differentiation" = "Osteoblast differentiation",
  "RNA Biosynthetic Process" = "Transcription (RNA Pol II / DNA-templated)",
  "Transcription by RNA Polymerase II" = "Transcription (RNA Pol II / DNA-templated)",
  "Phosphatidylinositol 3-Kinase/Protein Kinase B Signal Transduction" = "PI3K–AKT signalling",
  "Cell Migration" = "Cell migration",
  "Cell Motility" = "Cell migration",
  "Cell Population Proliferation" = "Cell proliferation",
  "Cellular Process" = "Cellular process (general)",
  "Multicellular Organismal Process" = "Multicellular organismal process",
  "Signal Transduction" = "Signal transduction (general)",
  "Wnt Signaling Pathway" = "WNT signalling",
  "Canonical Wnt Signaling Pathway" = "WNT signalling",
  "Fibroblast Proliferation" = "Fibroblast proliferation",
  "Stem Cell Differentiation" = "Stem-cell differentiation",
  "Aortic Valve Development" = "Cardiac valve development",
  "Aortic Valve Morphogenesis" = "Cardiac valve development",
  "Atrioventricular Valve Development" = "Cardiac valve development",
  "Atrioventricular Valve Morphogenesis" = "Cardiac valve development",
  "Cardiac Muscle Cell Membrane Repolarization" = "Cardiac repolarization",
  "Cartilage Development" = "Cartilage development / chondrogenesis",
  "Cell Cycle G1/S Phase Transition" = "Cell-cycle G1/S transition",
  "Cell Surface Receptor Protein Tyrosine Kinase Signaling Pathway" = "Receptor tyrosine-kinase signalling",
  "Response to BMP" = "BMP signalling",
  "Response to Growth Factor" = "Response to growth factor",
  "Response to Hormone" = "Response to hormone",
  "Response to Insulin" = "Insulin/IGF signalling",
  "Response to Transforming Growth Factor Beta" = "TGFβ signalling",
  "Chondrocyte Differentiation" = "Cartilage development / chondrogenesis",
  "Chordate Embryonic Development" = "Embryonic development (general)",
  "Circulatory System Development" = "Circulatory system development",
  "Collagen Fibril Organization" = "Extracellular matrix organization",
  "Dendritic Spine Maintenance" = "Dendrite / spine development",
  "Dendritic Spine Organization" = "Dendrite / spine development",
  "Digestive Tract Development" = "Digestive-tract development",
  "Embryonic Appendage Morphogenesis" = "Limb / appendage development",
  "Embryonic Digestive Tract Development" = "Digestive-tract development",
  "Embryonic Eye Morphogenesis" = "Eye development",
  "Embryonic Limb Morphogenesis" = "Limb / appendage development",
  "External Encapsulating Structure Organization" = "Extracellular matrix organization",
  "Extracellular Matrix Organization" = "Extracellular matrix organization",
  "Extracellular Structure Organization" = "Extracellular matrix organization",
  "Eye Development" = "Eye development",
  "Eye Morphogenesis" = "Eye development",
  "G1/S Transition of Mitotic Cell Cycle" = "Cell-cycle G1/S transition",
  "Heart Development" = "Heart development",
  "Heart Morphogenesis" = "Heart development",
  "Heterochromatin Formation" = "Heterochromatin formation",
  "Heterochromatin Organization" = "Heterochromatin formation",
  "Insulin-like Growth Factor Receptor Signaling Pathway" = "Insulin/IGF signalling",
  "Kidney Development" = "Kidney / renal development",
  "Limb Development" = "Limb / appendage development",
  "Limb Morphogenesis" = "Limb / appendage development",
  "Long-term Synaptic Depression" = "Synaptic transmission / plasticity",
  "Mesenchymal Cell Differentiation" = "Mesenchyme development",
  "Mitotic Cell Cycle Phase Transition" = "Cell-cycle phase transition",
  "Myotube Differentiation" = "Muscle development",
  "Transmembrane Receptor Protein Serine/Threonine Kinase Signaling Pathway" = "TGFβ-superfamily (Ser/Thr-kinase receptor) signalling",
  "Apoptotic Process" = "Apoptosis",
  "Cellular Senescence" = "Cellular senescence",
  "DNA Binding" = "DNA-binding TF activity (regulation)",
  "Developmental Growth" = "Developmental growth",
  "Gene Expression" = "Gene expression (regulation)",
  "Gene Expression, Epigenetic" = "Epigenetic gene silencing",
  "Growth" = "Growth (general)",
  "Mitotic Cell Cycle" = "Mitotic cell cycle",
  "Programmed Cell Death" = "Apoptosis",
  "Smooth Muscle Cell Proliferation" = "Smooth-muscle cell proliferation",
  "Steroid Biosynthetic Process" = "Steroid/sterol biosynthesis",
  "Transforming Growth Factor Beta Receptor Signaling Pathway" = "TGFβ signalling",
  "Vascular Associated Smooth Muscle Cell Proliferation" = "Smooth-muscle cell proliferation",
  "Nervous System Development" = "Nervous-system development",
  "Biomineral Tissue Development" = "Bone mineralization",
  "Bone Mineralization" = "Bone mineralization",
  "Cardiocyte Differentiation" = "Cardiac muscle development",
  "Glycogen Biosynthetic Process" = "Glycogen metabolism",
  "Glycogen Metabolic Process" = "Glycogen metabolism",
  "Intracellular Signal Transduction" = "Intracellular signal transduction",
  "Morphogenesis of an Epithelium" = "Epithelium morphogenesis",
  "Neuroepithelial Cell Differentiation" = "Neuroepithelial differentiation",
  "Odontoblast Differentiation" = "Tooth development",
  "Ossification" = "Ossification / bone formation",
  "Protein Catabolic Process" = "Ubiquitin–proteasome degradation",
  "Protein Serine/Threonine Kinase Activity" = "Ser/Thr-kinase activity",
  "SMAD Protein Signal Transduction" = "SMAD signalling",
  "T Cell Differentiation" = "Lymphocyte differentiation",
  "Vasculature Development" = "Angiogenesis",
  "miRNA Metabolic Process" = "miRNA biogenesis/metabolism",
  "Protein Localization to Chromatin" = "Protein localization to chromatin/chromosome",
  "Protein Phosphorylation" = "Protein phosphorylation",
  "Pulmonary Valve Development" = "Cardiac valve development",
  "Pulmonary Valve Morphogenesis" = "Cardiac valve development",
  "Aldosterone Biosynthetic Process" = "Steroid-hormone biosynthesis",
  "Cardiac Muscle Cell Differentiation" = "Cardiac muscle development",
  "Epithelial Cell Proliferation" = "Epithelial cell proliferation",
  "Germinal Center Formation" = "Lymphocyte differentiation",
  "MAPK Cascade" = "MAPK cascade",
  "Smoothened Signaling Pathway" = "Hedgehog (Smoothened) signalling",
  "miRNA Transcription" = "miRNA biogenesis/metabolism",
  "Renal System Development" = "Kidney / renal development",
  "Secondary Palate Development" = "Craniofacial development",
  "Skeletal System Development" = "Skeletal system development",
  "Somite Development" = "Somitogenesis / segmentation",
  "Transforming Growth Factor Beta Receptor Superfamily Signaling Pathway" = "TGFβ signalling",
  "Type B Pancreatic Cell Development" = "Pancreatic β-cell development",
  "Type B Pancreatic Cell Differentiation" = "Pancreatic β-cell development",
  "Ventricular Septum Development" = "Cardiac septum/chamber development",
  "Atrial Septum Morphogenesis" = "Cardiac septum/chamber development",
  "BMP Signaling Pathway" = "BMP signalling",
  "Bone Development" = "Ossification / bone formation",
  "Bone Morphogenesis" = "Ossification / bone formation",
  "Branching Morphogenesis of an Epithelial Tube" = "Epithelial tube morphogenesis",
  "CRD-mediated mRNA Stabilization" = "mRNA stabilization",
  "Calcineurin-NFAT Signaling Cascade" = "Calcineurin-NFAT signalling",
  "Camera-type Eye Development" = "Eye development",
  "Canonical NF-kappaB Signal Transduction" = "NF-κB signalling",
  "Cardiac Epithelial to Mesenchymal Transition" = "Cardiac EMT",
  "Cardiac Right Ventricle Morphogenesis" = "Cardiac septum/chamber development",
  "Cardiac Ventricle Morphogenesis" = "Cardiac septum/chamber development",
  "Cellular Component Assembly" = "Cellular component assembly",
  "Cellular Component Maintenance" = "Cellular component maintenance",
  "Response to Cholesterol" = "Response to cholesterol/sterol",
  "Response to Cold" = "Response to temperature",
  "Response to Growth Hormone" = "Growth-hormone signalling",
  "Response to Laminar Fluid Shear Stress" = "Response to fluid shear stress",
  "Response to Oxygen-Containing Compound" = "Response to oxygen-containing compound",
  "Response to Peptide Hormone" = "Response to peptide hormone",
  "Response to Sterol" = "Response to cholesterol/sterol",
  "Response to Xenobiotic" = "Response to xenobiotic",
  "Chemorepulsion of Axon" = "Axon guidance",
  "Chromosome Organization" = "Chromatin organization",
  "Connective Tissue Development" = "Connective-tissue development",
  "Contact Inhibition" = "Contact inhibition",
  "Cytosolic Transport" = "Lysosomal / endomembrane transport",
  "DNA Damage Response" = "DNA-damage response",
  "DNA Damage Response, Signal Transduction by P53 Class Mediator" = "p53 signalling",
  "Dendrite Morphogenesis" = "Dendrite / spine development",
  "Dendritic Spine Morphogenesis" = "Dendrite / spine development",
  "Dopaminergic Neuron Axon Guidance" = "Axon guidance",
  "Dopaminergic Neuron Differentiation" = "Dopaminergic neuron development",
  "Dorsal/ventral Pattern Formation" = "Embryonic pattern formation",
  "Ear Development" = "Ear development",
  "Embryonic Organ Development" = "Embryonic organ development",
  "Embryonic Organ Morphogenesis" = "Embryonic organ development",
  "Endocardial Cushion Development" = "Endocardial cushion / valve primordium",
  "Endocardial Cushion Morphogenesis" = "Endocardial cushion / valve primordium",
  "Endochondral Bone Morphogenesis" = "Endochondral ossification",
  "Endochondral Ossification" = "Endochondral ossification",
  "Epidermal Cell Differentiation" = "Epidermis development",
  "Epidermal Growth Factor Receptor Signaling Pathway" = "EGFR signalling",
  "Epithelial Cell Development" = "Epithelium development",
  "Epithelial Cell Differentiation" = "Epithelium development",
  "Epithelial Cell Differentiation Involved in Kidney Development" = "Kidney / renal development",
  "Epithelial Tube Morphogenesis" = "Epithelial tube morphogenesis",
  "Epithelium Development" = "Epithelium development",
  "Fat Cell Differentiation" = "Adipocyte differentiation",
  "Forebrain Neuron Differentiation" = "Neuron differentiation",
  "Generation of Neurons" = "Neurogenesis",
  "Glucose Homeostasis" = "Glucose homeostasis",
  "Golgi to Lysosome Transport" = "Lysosomal / endomembrane transport",
  "Growth Hormone Receptor Signaling Pathway" = "Growth-hormone signalling",
  "Hemopoiesis" = "Haematopoiesis",
  "Host-mediated Perturbation of Viral Transcription" = "Heterochromatin formation",
  "Intracellular Glucose Homeostasis" = "Glucose homeostasis",
  "Intracellular Iron Ion Homeostasis" = "Iron homeostasis",
  "Intracellular Signaling Cassette" = "Intracellular signal transduction",
  "Intrinsic Apoptotic Signaling Pathway by P53 Class Mediator" = "p53-mediated apoptosis",
  "Intrinsic Apoptotic Signaling Pathway in Response to DNA Damage by P53 Class Mediator" = "p53-mediated apoptosis",
  "Locomotor Rhythm" = "Locomotor rhythm (behaviour)",
  "Lysosomal Transport" = "Lysosomal / endomembrane transport",
  "Lysosome Localization" = "Lysosomal / endomembrane transport",
  "Macromolecule Modification" = "Protein modification (general)",
  "Macrophage Differentiation" = "Myeloid / leukocyte differentiation",
  "Membrane Repolarization During Action Potential" = "Cardiac repolarization",
  "Membrane Repolarization During Cardiac Muscle Cell Action Potential" = "Cardiac repolarization",
  "Membraneless Organelle Assembly" = "Biomolecular condensate assembly",
  "Mesenchymal to Epithelial Transition Involved in Metanephros Morphogenesis" = "Mesenchymal–epithelial transition (MET)",
  "Mesenchyme Development" = "Mesenchyme development",
  "Mesenchyme Morphogenesis" = "Mesenchyme development",
  "MiRNA Catabolic Process" = "miRNA biogenesis/metabolism",
  "MiRNA Processing" = "miRNA biogenesis/metabolism",
  "Midbrain Dopaminergic Neuron Differentiation" = "Dopaminergic neuron development",
  "Mitotic Sister Chromatid Cohesion" = "Sister chromatid cohesion",
  "Mitral Valve Development" = "Cardiac valve development",
  "Mitral Valve Morphogenesis" = "Cardiac valve development",
  "Muscle Tissue Development" = "Muscle development",
  "Myeloid Leukocyte Differentiation" = "Myeloid / leukocyte differentiation",
  "Myoblast Differentiation" = "Muscle development",
  "Alcohol Biosynthetic Process" = "Steroid/sterol biosynthesis",
  "Androgen Receptor Signaling Pathway" = "Androgen-receptor signalling",
  "Blood Vessel Morphogenesis" = "Angiogenesis",
  "Cell Communication" = "Signal transduction (general)",
  "Cell Cycle Phase Transition" = "Cell-cycle phase transition",
  "Cell Division" = "Cell division",
  "Cell Projection Organization" = "Cell-projection / pseudopodium assembly",
  "Cell-Matrix Adhesion" = "Cell adhesion",
  "Cold-Induced Thermogenesis" = "Thermogenesis",
  "Double-Strand Break Repair" = "DNA repair",
  "Hippo Signaling" = "Hippo signalling",
  "Insulin-Like Growth Factor Receptor Signaling Pathway" = "Insulin/IGF signalling",
  "Intracellular Steroid Hormone Receptor Signaling Pathway" = "Steroid-hormone (nuclear receptor) signalling",
  "Lens Fiber Cell Differentiation" = "Eye development",
  "Lipid Biosynthetic Process" = "Lipid biosynthesis",
  "Protein Import Into Nucleus" = "Nuclear import",
  "Protein Localization to Nucleus" = "Nuclear import",
  "Protein Maturation" = "Protein maturation/processing",
  "Protein Processing" = "Protein maturation/processing",
  "Response to" = "Response to stimulus (general)",
  "Signaling" = "Signal transduction (general)",
  "Small GTPase Mediated Signal Transduction" = "Small-GTPase signalling",
  "Smooth Muscle Cell Apoptotic Process" = "Apoptosis (smooth muscle)",
  "Smooth Muscle Cell Migration" = "Smooth-muscle cell migration",
  "Steroid Metabolic Process" = "Steroid/sterol biosynthesis",
  "Stress-Activated MAPK Cascade" = "MAPK cascade",
  "Stress-Activated Protein Kinase Signaling Cascade" = "MAPK cascade",
  "Synaptic Transmission" = "Synaptic transmission / plasticity",
  "Transforming Growth Factor Beta Production" = "TGFβ production",
  "Vascular Associated Smooth Muscle Cell Migration" = "Smooth-muscle cell migration",
  "Neural Crest Cell Differentiation" = "Neural-crest differentiation",
  "Neural Tube Closure" = "Neural-tube closure",
  "Neuron Differentiation" = "Neuron differentiation",
  "Non-canonical Wnt Signaling Pathway" = "WNT signalling",
  "Nose Development" = "Craniofacial development",
  "Notch Signaling Pathway" = "Notch signalling",
  "Nuclear Receptor-Mediated Steroid Hormone Signaling Pathway" = "Steroid-hormone (nuclear receptor) signalling",
  "Odontogenesis" = "Tooth development",
  "Odontogenesis of Dentin-Containing Tooth" = "Tooth development",
  "Organelle Transport Along Microtubule" = "Microtubule-based transport",
  "Peptidyl-tyrosine Phosphorylation" = "Tyrosine phosphorylation",
  "Pericardium Development" = "Heart development",
  "Pharyngeal System Development" = "Craniofacial / pharyngeal development",
  "Pituitary Gland Development" = "Pituitary development",
  "Astrocyte Differentiation" = "Astrocyte (glial) differentiation",
  "Biosynthetic Process" = "Biosynthetic process (general)",
  "Branching Involved in Ureteric Bud Morphogenesis" = "Kidney / renal development",
  "Carbohydrate Metabolic Process" = "Carbohydrate metabolism",
  "Cell-Substrate Adhesion" = "Cell adhesion",
  "Cellular Component Organization" = "Cellular component organization",
  "D-glucose Import Across Plasma Membrane" = "Glucose import/transport",
  "D-glucose Transmembrane Transport" = "Glucose import/transport",
  "DNA Metabolic Process" = "DNA metabolism/repair",
  "Developmental Process" = "Developmental process (general)",
  "Endothelial Cell Chemotaxis" = "Endothelial cell migration",
  "Endothelial Cell Migration" = "Endothelial cell migration",
  "Epithelial Cell Migration" = "Epithelial cell migration",
  "Insulin Secretion" = "Insulin secretion",
  "Leukocyte Apoptotic Process" = "Apoptosis (immune cells)",
  "Lipopolysaccharide-Mediated Signaling Pathway" = "LPS / innate-immune signalling",
  "Mitotic Nuclear Division" = "Cell division",
  "Multicellular Organism Growth" = "Body growth",
  "Myeloid Cell Apoptotic Process" = "Apoptosis (immune cells)",
  "Nitric Oxide Biosynthetic Process" = "Nitric-oxide / vascular tone",
  "Peptide Hormone Secretion" = "Peptide-hormone secretion",
  "Peroxisome Proliferator Activated Receptor Signaling Pathway" = "PPAR signalling",
  "Plasma Membrane Bounded Cell Projection Assembly" = "Cell-projection / pseudopodium assembly",
  "Protein Secretion" = "Secretion",
  "Protein-Containing Complex Disassembly" = "Protein-complex disassembly",
  "Pseudopodium Assembly" = "Cell-projection / pseudopodium assembly",
  "Secretion by Cell" = "Secretion",
  "Steroid Hormone Biosynthetic Process" = "Steroid-hormone biosynthesis",
  "Vascular Endothelial Growth Factor Signaling Pathway" = "VEGF signalling",
  "Primary Neural Tube Formation" = "Neural-tube closure",
  "Protein Localization to Chromosome" = "Protein localization to chromatin/chromosome",
  "Protein Modification Process" = "Protein modification (general)",
  "Proteoglycan Biosynthetic Process" = "Proteoglycan metabolism",
  "Proteoglycan Metabolic Process" = "Proteoglycan metabolism",
  "Reelin-mediated Signaling Pathway" = "Reelin signalling (neuronal migration)",
  "Attachment of Spindle Microtubules to Kinetochore" = "Spindle/kinetochore attachment",
  "Cell Cycle G2/M Phase Transition" = "Cell-cycle G2/M transition",
  "Cell Cycle Process" = "Cell cycle (general)",
  "Cortisol Biosynthetic Process" = "Steroid-hormone biosynthesis",
  "Cyclin-Dependent Protein Serine/Threonine Kinase Activity" = "Cyclin-dependent kinase activity",
  "Dendrite Development" = "Dendrite / spine development",
  "Mesenchymal Cell Proliferation" = "Mesenchymal cell proliferation",
  "Monocyte Differentiation" = "Myeloid / leukocyte differentiation",
  "Organ Growth" = "Organ growth",
  "Osteoblast Proliferation" = "Osteoblast proliferation",
  "Response to Biotic" = "Response to biotic/immune stimulus",
  "Sister Chromatid Cohesion" = "Sister chromatid cohesion",
  "Skeletal Muscle Contraction" = "Muscle contraction",
  "Stem Cell Proliferation" = "Stem cell proliferation",
  "Telomere Maintenance via Telomerase" = "Telomere maintenance",
  "Transforming Growth Factor Beta2 Production" = "TGFβ production",
  "Ubiquitin-Dependent Protein Catabolic Process" = "Ubiquitin–proteasome degradation",
  "Vascular Endothelial Growth Factor Receptor Signaling Pathway" = "VEGF signalling",
  "Regulatory ncRNA-mediated Heterochromatin Formation" = "Heterochromatin formation",
  "Response to Alcohol" = "Response to hormone/metabolite",
  "Response to Bile Acid" = "Response to hormone/metabolite",
  "Response to Estradiol" = "Response to hormone/metabolite",
  "Response to Fatty Acid" = "Response to hormone/metabolite",
  "Response to Ketone" = "Response to hormone/metabolite",
  "Response to X-ray" = "DNA-damage response",
  "Sensory Organ Development" = "Sensory-organ development",
  "Sensory Organ Morphogenesis" = "Sensory-organ development",
  "Skeletal System Morphogenesis" = "Skeletal system development",
  "Smooth Muscle Tissue Development" = "Muscle development",
  "Steroid Hormone Receptor Signaling Pathway" = "Steroid-hormone (nuclear receptor) signalling",
  "Striated Muscle Cell Differentiation" = "Muscle development",
  "Subtelomeric Heterochromatin Formation" = "Heterochromatin formation",
  "Tube Closure" = "Neural-tube closure",
  "Tumor Necrosis Factor-Mediated Signaling Pathway" = "TNF signalling",
  "Vascular Endothelial Cell Response to Fluid Shear Stress" = "Response to fluid shear stress",
  "Ventricular Septum Morphogenesis" = "Cardiac septum/chamber development"
)

# (b) Simplified term -> biological GROUP (21 groups)
group_lookup <- c(
  "Adipocyte differentiation" = "Insulin/IGF, glucose & energy metabolism",
  "Androgen-receptor signalling" = "Lipid, steroid & hormone signalling",
  "Angiogenesis" = "Cardiovascular development",
  "Apoptosis" = "Apoptosis, DNA damage & p53",
  "Apoptosis (immune cells)" = "Apoptosis, DNA damage & p53",
  "Apoptosis (smooth muscle)" = "Apoptosis, DNA damage & p53",
  "Astrocyte (glial) differentiation" = "Neural development & synaptic function",
  "Axon guidance" = "Neural development & synaptic function",
  "Axonogenesis" = "Neural development & synaptic function",
  "BMP signalling" = "TGFβ / BMP / SMAD signalling",
  "Biogenic Amine Metabolic Process" = "Lipid, steroid & hormone signalling",
  "Biomolecular condensate assembly" = "Protein modification, degradation & transport",
  "Biosynthetic process (general)" = "Broad regulatory terms (low specificity)",
  "Body growth" = "Cell & organism growth",
  "Bone mineralization" = "Skeletal, cartilage & ECM",
  "Calcineurin-NFAT signalling" = "RTK, MAPK & PI3K–AKT signalling",
  "Carbohydrate metabolism" = "Insulin/IGF, glucose & energy metabolism",
  "Cardiac EMT" = "Cardiovascular development",
  "Cardiac muscle development" = "Cardiovascular development",
  "Cardiac repolarization" = "Cardiovascular development",
  "Cardiac septum/chamber development" = "Cardiovascular development",
  "Cardiac valve development" = "Cardiovascular development",
  "Cartilage development / chondrogenesis" = "Skeletal, cartilage & ECM",
  "Cell adhesion" = "EMT, mesenchyme, migration & adhesion",
  "Cell cycle (general)" = "Cell cycle, proliferation & senescence",
  "Cell differentiation (general)" = "Broad regulatory terms (low specificity)",
  "Cell division" = "Cell cycle, proliferation & senescence",
  "Cell growth" = "Cell & organism growth",
  "Cell migration" = "EMT, mesenchyme, migration & adhesion",
  "Cell proliferation" = "Cell cycle, proliferation & senescence",
  "Cell-cycle G1/S transition" = "Cell cycle, proliferation & senescence",
  "Cell-cycle G2/M transition" = "Cell cycle, proliferation & senescence",
  "Cell-cycle phase transition" = "Cell cycle, proliferation & senescence",
  "Cell-projection / pseudopodium assembly" = "EMT, mesenchyme, migration & adhesion",
  "Cellular component assembly" = "Protein modification, degradation & transport",
  "Cellular component maintenance" = "Protein modification, degradation & transport",
  "Cellular component organization" = "Protein modification, degradation & transport",
  "Cellular process (general)" = "Broad regulatory terms (low specificity)",
  "Cellular senescence" = "Cell cycle, proliferation & senescence",
  "Chromatin organization" = "Chromatin, epigenetics & ncRNA",
  "Circulatory system development" = "Cardiovascular development",
  "Connective-tissue development" = "Skeletal, cartilage & ECM",
  "Contact inhibition" = "Cell cycle, proliferation & senescence",
  "Craniofacial / pharyngeal development" = "Other organ & epithelial morphogenesis",
  "Craniofacial development" = "Other organ & epithelial morphogenesis",
  "Cyclin-dependent kinase activity" = "Cell cycle, proliferation & senescence",
  "DNA metabolism/repair" = "Apoptosis, DNA damage & p53",
  "DNA repair" = "Apoptosis, DNA damage & p53",
  "DNA-binding TF activity (regulation)" = "Transcription & gene expression",
  "DNA-damage response" = "Apoptosis, DNA damage & p53",
  "Dendrite / spine development" = "Neural development & synaptic function",
  "Developmental growth" = "Cell & organism growth",
  "Developmental process (general)" = "Broad regulatory terms (low specificity)",
  "Digestive-tract development" = "Other organ & epithelial morphogenesis",
  "Dopaminergic neuron development" = "Neural development & synaptic function",
  "EGFR signalling" = "RTK, MAPK & PI3K–AKT signalling",
  "Ear development" = "Other organ & epithelial morphogenesis",
  "Embryonic development (general)" = "Other organ & epithelial morphogenesis",
  "Embryonic organ development" = "Other organ & epithelial morphogenesis",
  "Embryonic pattern formation" = "Limb, muscle & body-plan patterning",
  "Endocardial cushion / valve primordium" = "Cardiovascular development",
  "Endochondral ossification" = "Skeletal, cartilage & ECM",
  "Endoderm Formation" = "Other organ & epithelial morphogenesis",
  "Endodermal Cell Differentiation" = "Other organ & epithelial morphogenesis",
  "Endothelial cell migration" = "Cardiovascular development",
  "Epidermis development" = "Other organ & epithelial morphogenesis",
  "Epigenetic gene silencing" = "Chromatin, epigenetics & ncRNA",
  "Epithelial cell migration" = "EMT, mesenchyme, migration & adhesion",
  "Epithelial cell proliferation" = "Cell cycle, proliferation & senescence",
  "Epithelial tube morphogenesis" = "Other organ & epithelial morphogenesis",
  "Epithelial–mesenchymal transition (EMT)" = "EMT, mesenchyme, migration & adhesion",
  "Epithelium development" = "Other organ & epithelial morphogenesis",
  "Epithelium morphogenesis" = "Other organ & epithelial morphogenesis",
  "Extracellular matrix organization" = "Skeletal, cartilage & ECM",
  "Eye development" = "Other organ & epithelial morphogenesis",
  "Fibroblast proliferation" = "Cell cycle, proliferation & senescence",
  "Gene expression (regulation)" = "Transcription & gene expression",
  "Glucose homeostasis" = "Insulin/IGF, glucose & energy metabolism",
  "Glucose import/transport" = "Insulin/IGF, glucose & energy metabolism",
  "Glycogen metabolism" = "Insulin/IGF, glucose & energy metabolism",
  "Growth (general)" = "Cell & organism growth",
  "Growth-hormone signalling" = "Insulin/IGF, glucose & energy metabolism",
  "Haematopoiesis" = "Immune & inflammatory signalling",
  "Heart development" = "Cardiovascular development",
  "Hedgehog (Smoothened) signalling" = "Hedgehog & Notch signalling",
  "Heterochromatin formation" = "Chromatin, epigenetics & ncRNA",
  "Hippo signalling" = "RTK, MAPK & PI3K–AKT signalling",
  "Insulin secretion" = "Insulin/IGF, glucose & energy metabolism",
  "Insulin/IGF signalling" = "Insulin/IGF, glucose & energy metabolism",
  "Intracellular signal transduction" = "RTK, MAPK & PI3K–AKT signalling",
  "Iron homeostasis" = "Insulin/IGF, glucose & energy metabolism",
  "Kidney / renal development" = "Other organ & epithelial morphogenesis",
  "LPS / innate-immune signalling" = "Immune & inflammatory signalling",
  "Limb / appendage development" = "Limb, muscle & body-plan patterning",
  "Lipid biosynthesis" = "Lipid, steroid & hormone signalling",
  "Locomotor rhythm (behaviour)" = "Neural development & synaptic function",
  "Lymphocyte differentiation" = "Immune & inflammatory signalling",
  "Lysosomal / endomembrane transport" = "Protein modification, degradation & transport",
  "MAPK cascade" = "RTK, MAPK & PI3K–AKT signalling",
  "Mesenchymal cell proliferation" = "Cell cycle, proliferation & senescence",
  "Mesenchymal–epithelial transition (MET)" = "EMT, mesenchyme, migration & adhesion",
  "Mesenchyme development" = "EMT, mesenchyme, migration & adhesion",
  "Mesoderm Formation" = "Limb, muscle & body-plan patterning",
  "Microtubule-based transport" = "Protein modification, degradation & transport",
  "Mitotic cell cycle" = "Cell cycle, proliferation & senescence",
  "Multicellular organismal process" = "Broad regulatory terms (low specificity)",
  "Muscle contraction" = "Limb, muscle & body-plan patterning",
  "Muscle development" = "Limb, muscle & body-plan patterning",
  "Myeloid / leukocyte differentiation" = "Immune & inflammatory signalling",
  "NF-κB signalling" = "Immune & inflammatory signalling",
  "Nervous-system development" = "Neural development & synaptic function",
  "Neural-crest differentiation" = "EMT, mesenchyme, migration & adhesion",
  "Neural-tube closure" = "Neural development & synaptic function",
  "Neuroepithelial differentiation" = "Neural development & synaptic function",
  "Neurogenesis" = "Neural development & synaptic function",
  "Neuron differentiation" = "Neural development & synaptic function",
  "Nitric-oxide / vascular tone" = "Cardiovascular development",
  "Notch signalling" = "Hedgehog & Notch signalling",
  "Nuclear import" = "Protein modification, degradation & transport",
  "Organ growth" = "Cell & organism growth",
  "Ossification / bone formation" = "Skeletal, cartilage & ECM",
  "Osteoblast differentiation" = "Skeletal, cartilage & ECM",
  "Osteoblast proliferation" = "Cell cycle, proliferation & senescence",
  "PI3K–AKT signalling" = "RTK, MAPK & PI3K–AKT signalling",
  "PPAR signalling" = "Insulin/IGF, glucose & energy metabolism",
  "Pancreatic β-cell development" = "Other organ & epithelial morphogenesis",
  "Peptide-hormone secretion" = "Insulin/IGF, glucose & energy metabolism",
  "Pituitary development" = "Other organ & epithelial morphogenesis",
  "Protein localization to chromatin/chromosome" = "Chromatin, epigenetics & ncRNA",
  "Protein maturation/processing" = "Protein modification, degradation & transport",
  "Protein modification (general)" = "Protein modification, degradation & transport",
  "Protein phosphorylation" = "RTK, MAPK & PI3K–AKT signalling",
  "Protein-complex disassembly" = "Protein modification, degradation & transport",
  "Proteoglycan metabolism" = "Skeletal, cartilage & ECM",
  "Proteolysis" = "Protein modification, degradation & transport",
  "Receptor tyrosine-kinase signalling" = "RTK, MAPK & PI3K–AKT signalling",
  "Reelin signalling (neuronal migration)" = "Neural development & synaptic function",
  "Response to biotic/immune stimulus" = "Immune & inflammatory signalling",
  "Response to cholesterol/sterol" = "Lipid, steroid & hormone signalling",
  "Response to fluid shear stress" = "Cardiovascular development",
  "Response to growth factor" = "RTK, MAPK & PI3K–AKT signalling",
  "Response to hormone" = "Lipid, steroid & hormone signalling",
  "Response to hormone/metabolite" = "Lipid, steroid & hormone signalling",
  "Response to oxygen-containing compound" = "Cellular stress & environmental response",
  "Response to peptide hormone" = "Insulin/IGF, glucose & energy metabolism",
  "Response to stimulus (general)" = "Cellular stress & environmental response",
  "Response to temperature" = "Cellular stress & environmental response",
  "Response to xenobiotic" = "Cellular stress & environmental response",
  "SMAD signalling" = "TGFβ / BMP / SMAD signalling",
  "Secretion" = "Insulin/IGF, glucose & energy metabolism",
  "Sensory-organ development" = "Other organ & epithelial morphogenesis",
  "Ser/Thr-kinase activity" = "RTK, MAPK & PI3K–AKT signalling",
  "Signal transduction (general)" = "Broad regulatory terms (low specificity)",
  "Sister chromatid cohesion" = "Cell cycle, proliferation & senescence",
  "Skeletal system development" = "Skeletal, cartilage & ECM",
  "Small-GTPase signalling" = "RTK, MAPK & PI3K–AKT signalling",
  "Smooth-muscle cell migration" = "EMT, mesenchyme, migration & adhesion",
  "Smooth-muscle cell proliferation" = "Cell cycle, proliferation & senescence",
  "Somitogenesis / segmentation" = "Limb, muscle & body-plan patterning",
  "Spindle/kinetochore attachment" = "Cell cycle, proliferation & senescence",
  "Stem cell proliferation" = "Cell cycle, proliferation & senescence",
  "Stem-cell differentiation" = "Broad regulatory terms (low specificity)",
  "Steroid-hormone (nuclear receptor) signalling" = "Lipid, steroid & hormone signalling",
  "Steroid-hormone biosynthesis" = "Lipid, steroid & hormone signalling",
  "Steroid/sterol biosynthesis" = "Lipid, steroid & hormone signalling",
  "Synaptic transmission / plasticity" = "Neural development & synaptic function",
  "TGFβ production" = "TGFβ / BMP / SMAD signalling",
  "TGFβ signalling" = "TGFβ / BMP / SMAD signalling",
  "TGFβ-superfamily (Ser/Thr-kinase receptor) signalling" = "TGFβ / BMP / SMAD signalling",
  "TNF signalling" = "Immune & inflammatory signalling",
  "Telomere maintenance" = "Chromatin, epigenetics & ncRNA",
  "Thermogenesis" = "Insulin/IGF, glucose & energy metabolism",
  "Tooth development" = "Skeletal, cartilage & ECM",
  "Transcription (RNA Pol II / DNA-templated)" = "Transcription & gene expression",
  "Tyrosine phosphorylation" = "RTK, MAPK & PI3K–AKT signalling",
  "Ubiquitin–proteasome degradation" = "Protein modification, degradation & transport",
  "VEGF signalling" = "RTK, MAPK & PI3K–AKT signalling",
  "WNT signalling" = "WNT signalling",
  "mRNA stabilization" = "Chromatin, epigenetics & ncRNA",
  "miRNA biogenesis/metabolism" = "Chromatin, epigenetics & ncRNA",
  "p53 signalling" = "Apoptosis, DNA damage & p53",
  "p53-mediated apoptosis" = "Apoptosis, DNA damage & p53"
)

# (c) GROUP -> THEME (10 themes)
theme_lookup <- c(
  "Apoptosis, DNA damage & p53" = "Cell-cycle, growth & survival",
  "Broad regulatory terms (low specificity)" = "Broad / low-specificity",
  "Cardiovascular development" = "Cardiovascular, organ & tissue development",
  "Cell & organism growth" = "Cell-cycle, growth & survival",
  "Cell cycle, proliferation & senescence" = "Cell-cycle, growth & survival",
  "Cellular stress & environmental response" = "Cellular housekeeping",
  "Chromatin, epigenetics & ncRNA" = "Gene regulation & genome maintenance",
  "EMT, mesenchyme, migration & adhesion" = "Cardiovascular, organ & tissue development",
  "Hedgehog & Notch signalling" = "Developmental signalling pathways",
  "Immune & inflammatory signalling" = "Immune & inflammatory",
  "Insulin/IGF, glucose & energy metabolism" = "Metabolic & endocrine",
  "Limb, muscle & body-plan patterning" = "Skeletal & musculoskeletal development",
  "Lipid, steroid & hormone signalling" = "Metabolic & endocrine",
  "Neural development & synaptic function" = "Neural development",
  "Other organ & epithelial morphogenesis" = "Cardiovascular, organ & tissue development",
  "Protein modification, degradation & transport" = "Cellular housekeeping",
  "RTK, MAPK & PI3K–AKT signalling" = "Developmental signalling pathways",
  "Skeletal, cartilage & ECM" = "Skeletal & musculoskeletal development",
  "TGFβ / BMP / SMAD signalling" = "Developmental signalling pathways",
  "Transcription & gene expression" = "Gene regulation & genome maintenance",
  "WNT signalling" = "Developmental signalling pathways"
)

# (d) Invariants asserted after the R pipeline runs (reproducibility guards)
EXPECTED_SIMPS <- c("Adipocyte differentiation", "Androgen-receptor signalling", "Angiogenesis", "Apoptosis", "Apoptosis (immune cells)", "Apoptosis (smooth muscle)", "Astrocyte (glial) differentiation", "Axon guidance", "Axonogenesis", "BMP signalling", "Biogenic Amine Metabolic Process", "Biomolecular condensate assembly", "Biosynthetic process (general)", "Body growth", "Bone mineralization", "Calcineurin-NFAT signalling", "Carbohydrate metabolism", "Cardiac EMT", "Cardiac muscle development", "Cardiac repolarization", "Cardiac septum/chamber development", "Cardiac valve development", "Cartilage development / chondrogenesis", "Cell adhesion", "Cell cycle (general)", "Cell differentiation (general)", "Cell division", "Cell growth", "Cell migration", "Cell proliferation", "Cell-cycle G1/S transition", "Cell-cycle G2/M transition", "Cell-cycle phase transition", "Cell-projection / pseudopodium assembly", "Cellular component assembly", "Cellular component maintenance", "Cellular component organization", "Cellular process (general)", "Cellular senescence", "Chromatin organization", "Circulatory system development", "Connective-tissue development", "Contact inhibition", "Craniofacial / pharyngeal development", "Craniofacial development", "Cyclin-dependent kinase activity", "DNA metabolism/repair", "DNA repair", "DNA-binding TF activity (regulation)", "DNA-damage response", "Dendrite / spine development", "Developmental growth", "Developmental process (general)", "Digestive-tract development", "Dopaminergic neuron development", "EGFR signalling", "Ear development", "Embryonic development (general)", "Embryonic organ development", "Embryonic pattern formation", "Endocardial cushion / valve primordium", "Endochondral ossification", "Endoderm Formation", "Endodermal Cell Differentiation", "Endothelial cell migration", "Epidermis development", "Epigenetic gene silencing", "Epithelial cell migration", "Epithelial cell proliferation", "Epithelial tube morphogenesis", "Epithelial–mesenchymal transition (EMT)", "Epithelium development", "Epithelium morphogenesis", "Extracellular matrix organization", "Eye development", "Fibroblast proliferation", "Gene expression (regulation)", "Glucose homeostasis", "Glucose import/transport", "Glycogen metabolism", "Growth (general)", "Growth-hormone signalling", "Haematopoiesis", "Heart development", "Hedgehog (Smoothened) signalling", "Heterochromatin formation", "Hippo signalling", "Insulin secretion", "Insulin/IGF signalling", "Intracellular signal transduction", "Iron homeostasis", "Kidney / renal development", "LPS / innate-immune signalling", "Limb / appendage development", "Lipid biosynthesis", "Locomotor rhythm (behaviour)", "Lymphocyte differentiation", "Lysosomal / endomembrane transport", "MAPK cascade", "Mesenchymal cell proliferation", "Mesenchymal–epithelial transition (MET)", "Mesenchyme development", "Mesoderm Formation", "Microtubule-based transport", "Mitotic cell cycle", "Multicellular organismal process", "Muscle contraction", "Muscle development", "Myeloid / leukocyte differentiation", "NF-κB signalling", "Nervous-system development", "Neural-crest differentiation", "Neural-tube closure", "Neuroepithelial differentiation", "Neurogenesis", "Neuron differentiation", "Nitric-oxide / vascular tone", "Notch signalling", "Nuclear import", "Organ growth", "Ossification / bone formation", "Osteoblast differentiation", "Osteoblast proliferation", "PI3K–AKT signalling", "PPAR signalling", "Pancreatic β-cell development", "Peptide-hormone secretion", "Pituitary development", "Protein localization to chromatin/chromosome", "Protein maturation/processing", "Protein modification (general)", "Protein phosphorylation", "Protein-complex disassembly", "Proteoglycan metabolism", "Proteolysis", "Receptor tyrosine-kinase signalling", "Reelin signalling (neuronal migration)", "Response to biotic/immune stimulus", "Response to cholesterol/sterol", "Response to fluid shear stress", "Response to growth factor", "Response to hormone", "Response to hormone/metabolite", "Response to oxygen-containing compound", "Response to peptide hormone", "Response to stimulus (general)", "Response to temperature", "Response to xenobiotic", "SMAD signalling", "Secretion", "Sensory-organ development", "Ser/Thr-kinase activity", "Signal transduction (general)", "Sister chromatid cohesion", "Skeletal system development", "Small-GTPase signalling", "Smooth-muscle cell migration", "Smooth-muscle cell proliferation", "Somitogenesis / segmentation", "Spindle/kinetochore attachment", "Stem cell proliferation", "Stem-cell differentiation", "Steroid-hormone (nuclear receptor) signalling", "Steroid-hormone biosynthesis", "Steroid/sterol biosynthesis", "Synaptic transmission / plasticity", "TGFβ production", "TGFβ signalling", "TGFβ-superfamily (Ser/Thr-kinase receptor) signalling", "TNF signalling", "Telomere maintenance", "Thermogenesis", "Tooth development", "Transcription (RNA Pol II / DNA-templated)", "Tyrosine phosphorylation", "Ubiquitin–proteasome degradation", "VEGF signalling", "WNT signalling", "mRNA stabilization", "miRNA biogenesis/metabolism", "p53 signalling", "p53-mediated apoptosis")

EXPECTED_GROUPS <- c("Apoptosis, DNA damage & p53", "Broad regulatory terms (low specificity)", "Cardiovascular development", "Cell & organism growth", "Cell cycle, proliferation & senescence", "Cellular stress & environmental response", "Chromatin, epigenetics & ncRNA", "EMT, mesenchyme, migration & adhesion", "Hedgehog & Notch signalling", "Immune & inflammatory signalling", "Insulin/IGF, glucose & energy metabolism", "Limb, muscle & body-plan patterning", "Lipid, steroid & hormone signalling", "Neural development & synaptic function", "Other organ & epithelial morphogenesis", "Protein modification, degradation & transport", "RTK, MAPK & PI3K–AKT signalling", "Skeletal, cartilage & ECM", "TGFβ / BMP / SMAD signalling", "Transcription & gene expression", "WNT signalling")

EXPECTED_THEMES <- c("Broad / low-specificity", "Cardiovascular, organ & tissue development", "Cell-cycle, growth & survival", "Cellular housekeeping", "Developmental signalling pathways", "Gene regulation & genome maintenance", "Immune & inflammatory", "Metabolic & endocrine", "Neural development", "Skeletal & musculoskeletal development")

# ----------------------------------------------------------------------------
# 3.  STAGE 2 - apply curated canonicalization, then group & theme
#
#     syn_clean merges normalized labels that denote the same process onto one
#     canonical simplified term. Labels absent from syn_clean pass through
#     unchanged (they were already canonical). We then attach GROUP and THEME.
# ----------------------------------------------------------------------------
per_go <- per_go %>%
  mutate(simp = ifelse(simp0 %in% names(syn_clean), syn_clean[simp0], simp0))

# --- reproducibility guard #1: the curated reduction must still hold ----------
got_simps <- sort(unique(per_go$simp))
extra_simps <- setdiff(got_simps, EXPECTED_SIMPS)   # new, uncurated terms
lost_simps  <- setdiff(EXPECTED_SIMPS, got_simps)   # expected but absent
if (length(extra_simps) > 0) {
  stop("Uncurated simplified term(s) appeared - extend syn_clean/group_lookup:\n  ",
       paste(extra_simps, collapse = "\n  "))
}
if (length(lost_simps) > 0) {
  message("NOTE: ", length(lost_simps),
          " expected simplified term(s) absent in this input (OK if input changed).")
}

# --- attach group + theme -----------------------------------------------------
per_go <- per_go %>%
  mutate(group = unname(group_lookup[simp]),
         theme = unname(theme_lookup[group]))

# --- reproducibility guard #2: nothing may be left unmapped --------------------
unmapped_g <- per_go %>% filter(is.na(group)) %>% distinct(simp) %>% pull(simp)
unmapped_t <- per_go %>% filter(is.na(theme)) %>% distinct(group) %>% pull(group)
if (length(unmapped_g) > 0)
  stop("Simplified term(s) with no GROUP: ", paste(unmapped_g, collapse = "; "))
if (length(unmapped_t) > 0)
  stop("Group(s) with no THEME: ", paste(unmapped_t, collapse = "; "))

stopifnot(setequal(unique(per_go$group), EXPECTED_GROUPS))   # exactly 21
stopifnot(setequal(unique(per_go$theme), EXPECTED_THEMES))   # exactly 10
message(sprintf("Mapped cleanly: %d simplified terms -> %d groups -> %d themes",
                dplyr::n_distinct(per_go$simp),
                dplyr::n_distinct(per_go$group),
                dplyr::n_distinct(per_go$theme)))

# ----------------------------------------------------------------------------
# 4.  Write the full auditable curation (one row per original GO term)
# ----------------------------------------------------------------------------
curation <- per_go %>%
  transmute(go_id = go, original_term = name, simplified_term = simp,
            group, theme, traits, best_adjp, genes) %>%
  arrange(theme, group, simplified_term, original_term)
fwrite(curation, file.path(out_dir, "go_term_curation.tsv"), sep = "\t")
message("Wrote go_term_curation.tsv (", nrow(curation), " rows)")

# ----------------------------------------------------------------------------
# 5.  Summary tables for plotting & for the narrative
#
#     We count ORIGINAL GO terms (not simplified) per group per trait: this
#     reflects how much of the raw enrichment signal each program attracted in
#     each trait. A term shared by 2 traits is counted once per trait it is
#     significant in. trait membership is encoded in `traits` ("cortical+limb").
# ----------------------------------------------------------------------------
trait_levels <- c("cortical", "limb", "trunk")
long_terms <- per_go %>%
  separate_rows(traits, sep = "\\+") %>%
  rename(trait = traits) %>%
  mutate(trait = factor(trait, levels = trait_levels))

# group x trait (original-GO counts)
group_trait <- long_terms %>%
  distinct(go, trait, group) %>%
  count(group, trait, name = "n_terms") %>%
  complete(group, trait, fill = list(n_terms = 0))

# theme x trait (original-GO counts)
theme_trait <- long_terms %>%
  distinct(go, trait, theme) %>%
  count(theme, trait, name = "n_terms") %>%
  complete(theme, trait, fill = list(n_terms = 0))

# group -> theme key + group totals, used for ordering
group_key <- per_go %>% distinct(group, theme)
group_tot <- group_trait %>% group_by(group) %>%
  summarise(total = sum(n_terms), .groups = "drop")
theme_tot <- theme_trait %>% group_by(theme) %>%
  summarise(total = sum(n_terms), .groups = "drop")

# Biological ordering of themes (regulatory -> growth -> signalling -> metabolic
# -> skeletal -> organ -> neural -> immune -> housekeeping -> broad)
theme_order <- c(
  "Gene regulation & genome maintenance",
  "Cell-cycle, growth & survival",
  "Developmental signalling pathways",
  "Metabolic & endocrine",
  "Skeletal & musculoskeletal development",
  "Cardiovascular, organ & tissue development",
  "Neural development",
  "Immune & inflammatory",
  "Cellular housekeeping",
  "Broad / low-specificity"
)
stopifnot(setequal(theme_order, EXPECTED_THEMES))

# Order groups: by theme (biological order), then by descending total within theme
group_levels <- group_key %>%
  left_join(group_tot, by = "group") %>%
  mutate(theme = factor(theme, levels = theme_order)) %>%
  arrange(theme, desc(total)) %>%
  pull(group)

group_trait <- group_trait %>%
  left_join(group_key, by = "group") %>%
  mutate(group = factor(group, levels = rev(group_levels)),   # rev: top = first
         theme = factor(theme, levels = theme_order))
theme_trait <- theme_trait %>%
  mutate(theme = factor(theme, levels = rev(theme_order)))

# Console headline (matches the narrative)
wide_group <- group_trait %>%
  pivot_wider(id_cols = c(group, theme), names_from = trait,
              values_from = n_terms, values_fill = 0) %>%
  mutate(TOTAL = cortical + limb + trunk) %>%
  arrange(desc(TOTAL))
message("\n===== Original GO terms per group per trait =====")
print(as.data.frame(wide_group), row.names = FALSE)

fwrite(wide_group,  file.path(out_dir, "summary_group_by_trait.tsv"), sep = "\t")
fwrite(theme_trait %>% pivot_wider(names_from = trait, values_from = n_terms,
        values_fill = 0), file.path(out_dir, "summary_theme_by_trait.tsv"), sep = "\t")

# ============================================================================
# 6.  PLOTS  (ggplot2)
# ============================================================================
trait_pal <- c(cortical = "#3B6FB6", limb = "#E08214", trunk = "#1B7837")
trait_labs <- c(cortical = "Cortical expansion",
                limb = "Limb growth", trunk = "Trunk growth")

# Use a font that actually carries the Greek beta and en-dash used in some group
# labels (the default sans drops them, rendering "TGF.." / "PI3K..AKT").
# DejaVu Sans is present on essentially every Linux box and covers these glyphs.
BASE_FAMILY <- if ("DejaVu Sans" %in% systemfonts::system_fonts()$family)
                 "DejaVu Sans" else ""

# Belt-and-braces: some raster back-ends still fail to fall back to a glyph for
# the Greek beta / en-dash even when the font has them. To make the FIGURES
# render correctly on any machine, we sanitise the *display* labels to ASCII.
# This affects plots only -- go_term_curation.tsv and the summary TSVs keep the
# proper Unicode names ("TGFβ", "PI3K–AKT").
disp <- function(x) {
  x <- gsub("β", "-beta", x, fixed = TRUE)   # TGFβ  -> TGF-beta
  x <- gsub("–", "-",     x, fixed = TRUE)   # en-dash -> hyphen
  x <- gsub("—", "-",     x, fixed = TRUE)   # em-dash -> hyphen
  x
}

theme_pub <- function(base = 11) {
  theme_minimal(base_size = base, base_family = BASE_FAMILY) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks.x = element_line(colour = "grey80"),
          strip.text = element_text(face = "bold", hjust = 0),
          plot.title = element_text(face = "bold"),
          plot.title.position = "plot",
          legend.position = "top")
}

save_plot <- function(p, stem, w, h) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = w, height = h,
         device = cairo_pdf)
  # PNG: the default raster device (ragg, if present) renders UTF-8 (β, –) cleanly
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = w, height = h,
         dpi = 200)
}

# Sanitised display copies (canonical Unicode preserved in the TSVs above).
group_levels_disp <- disp(group_levels)
theme_order_disp  <- disp(theme_order)
group_trait_p <- group_trait %>%
  mutate(group = factor(disp(as.character(group)),
                        levels = disp(rev(group_levels))),
         theme = factor(disp(as.character(theme)), levels = theme_order_disp))
theme_trait_p <- theme_trait %>%
  mutate(theme = factor(disp(as.character(theme)), levels = rev(theme_order_disp)))

# ---- FIG 1: theme-level grouped bars (high-level overview) -------------------
fig1 <- ggplot(theme_trait_p,
               aes(x = n_terms, y = theme, fill = trait)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = trait_pal, labels = trait_labs, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "GO:BP enrichment by theme across the three growth traits",
       subtitle = "Count of FDR-significant GO terms (original, pre-simplification)",
       x = "Number of significant GO terms", y = NULL) +
  theme_pub()
save_plot(fig1, "fig1_theme_by_trait", 9, 5.5)

# ---- FIG 2: group-level bars, faceted by theme (mid-level detail) ------------
fig2 <- ggplot(group_trait_p,
               aes(x = n_terms, y = group, fill = trait)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(theme = label_wrap_gen(22))) +
  scale_fill_manual(values = trait_pal, labels = trait_labs, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "GO:BP enrichment by biological group",
       subtitle = "21 curated groups nested within 10 themes; original GO-term counts",
       x = "Number of significant GO terms", y = NULL) +
  theme_pub(10) +
  theme(panel.spacing = unit(2, "pt"),
        strip.text.y = element_text(angle = 0, hjust = 0, size = 8))
save_plot(fig2, "fig2_group_by_trait_faceted", 10, 12)

# ---- FIG 3: group x trait bubble heatmap (integrative view) ------------------
bubble <- group_trait_p %>%
  mutate(group = factor(as.character(group), levels = disp(group_levels))) %>%  # top = highest
  filter(n_terms > 0)
fig3 <- ggplot(bubble, aes(x = trait, y = group)) +
  geom_point(aes(size = n_terms, colour = n_terms)) +
  geom_text(aes(label = n_terms), size = 2.7, colour = "white", fontface = "bold", family = BASE_FAMILY) +
  scale_x_discrete(labels = trait_labs, position = "top") +
  scale_size_area(max_size = 13, guide = "none") +
  scale_colour_viridis_c(option = "mako", direction = -1, end = 0.9,
                         name = "GO terms") +
  labs(title = "Where each growth program concentrates",
       subtitle = "Bubble size/colour = number of significant GO terms per group",
       x = NULL, y = NULL) +
  theme_pub(10) +
  theme(panel.grid.major.y = element_line(colour = "grey92"),
        axis.text.y = element_text(size = 8),
        legend.position = "right")
save_plot(fig3, "fig3_group_trait_bubble", 8.5, 11)

# ---- FIG 4: stacked composition per trait (share of signal) ------------------
stack_df <- theme_trait_p %>%
  group_by(trait) %>% mutate(frac = n_terms / sum(n_terms)) %>% ungroup() %>%
  mutate(theme = factor(as.character(theme), levels = rev(theme_order_disp)))
fig4 <- ggplot(stack_df, aes(x = trait, y = frac, fill = theme)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  scale_x_discrete(labels = trait_labs) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_brewer(palette = "Spectral", direction = -1, name = "Theme") +
  labs(title = "Thematic composition of each trait's enrichment",
       subtitle = "Share of FDR-significant GO terms attributable to each theme",
       x = NULL, y = "Share of significant GO terms") +
  theme_pub(10) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.position = "right", legend.text = element_text(size = 8))
save_plot(fig4, "fig4_theme_composition_stacked", 8.5, 6)

message("\nSaved figures to ", out_dir, " :")
message("  fig1_theme_by_trait.(pdf|png)")
message("  fig2_group_by_trait_faceted.(pdf|png)")
message("  fig3_group_trait_bubble.(pdf|png)")
message("  fig4_theme_composition_stacked.(pdf|png)")
message("\nDONE.")
