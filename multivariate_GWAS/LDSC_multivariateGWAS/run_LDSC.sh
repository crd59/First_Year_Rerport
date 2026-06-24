#!/bin/bash

module load miniconda3

source activate ldsc

./../../../ldsc/ldsc.py \
--rg ./trunk_munged.sumstats.gz,./limb_munged.sumstats.gz,./cortical_munged.sumstats.gz,./../../../ldsc_analyses/important_sumstats/avg_femur.sumstats.gz,./../../../ldsc_analyses/important_sumstats/avg_forearm.sumstats.gz,./../../../ldsc_analyses/important_sumstats/avg_humerus.sumstats.gz,./../../../ldsc_analyses/important_sumstats/avg_tibia.sumstats.gz,./../../../ldsc_analyses/important_sumstats/global_FI.sumstats.gz,./../../../ldsc_analyses/important_sumstats/global_IC.sumstats.gz,./../../../ldsc_analyses/important_sumstats/global_SA.sumstats.gz,./../../../ldsc_analyses/important_sumstats/hip_width.sumstats.gz,./../../../ldsc_analyses/important_sumstats/shoulder_width.sumstats.gz,./../../../ldsc_analyses/important_sumstats/torso_length.sumstats.gz,./../../../ldsc_analyses/important_sumstats/new_height.sumstats.gz \
--ref-ld-chr ./../../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../../ldsc/eur_w_ld_chr/ \
--out trunk_ldsc

./../../../ldsc/ldsc.py \
--rg ./limb_munged.sumstats.gz,./cortical_munged.sumstats.gz,./../../../ldsc_analyses/important_sumstats/avg_femur.sumstats.gz,./../../../ldsc_analyses/important_sumstats/avg_forearm.sumstats.gz,./../../../ldsc_analyses/important_sumstats/avg_humerus.sumstats.gz,./../../../ldsc_analyses/important_sumstats/avg_tibia.sumstats.gz,./../../../ldsc_analyses/important_sumstats/global_FI.sumstats.gz,./../../../ldsc_analyses/important_sumstats/global_IC.sumstats.gz,./../../../ldsc_analyses/important_sumstats/global_SA.sumstats.gz,./../../../ldsc_analyses/important_sumstats/hip_width.sumstats.gz,./../../../ldsc_analyses/important_sumstats/shoulder_width.sumstats.gz,./../../../ldsc_analyses/important_sumstats/torso_length.sumstats.gz,./../../../ldsc_analyses/important_sumstats/new_height.sumstats.gz \
--ref-ld-chr ./../../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../../ldsc/eur_w_ld_chr/ \
--out limb_ldsc

./../../../ldsc/ldsc.py \
--rg ./cortical_munged.sumstats.gz,./../../../ldsc_analyses/important_sumstats/avg_femur.sumstats.gz,./../../../ldsc_analyses/important_sumstats/avg_forearm.sumstats.gz,./../../../ldsc_analyses/important_sumstats/avg_humerus.sumstats.gz,./../../../ldsc_analyses/important_sumstats/avg_tibia.sumstats.gz,./../../../ldsc_analyses/important_sumstats/global_FI.sumstats.gz,./../../../ldsc_analyses/important_sumstats/global_IC.sumstats.gz,./../../../ldsc_analyses/important_sumstats/global_SA.sumstats.gz,./../../../ldsc_analyses/important_sumstats/hip_width.sumstats.gz,./../../../ldsc_analyses/important_sumstats/shoulder_width.sumstats.gz,./../../../ldsc_analyses/important_sumstats/torso_length.sumstats.gz,./../../../ldsc_analyses/important_sumstats/new_height.sumstats.gz \
--ref-ld-chr ./../../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../../ldsc/eur_w_ld_chr/ \
--out cortical_ldsc