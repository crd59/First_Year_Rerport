#!/bin/bash

module load miniconda3

source activate ldsc

./../../ldsc/ldsc.py \
--rg ./../new_GWAS_ldsc/new_height.sumstats.gz,./../important_sumstats/avg_femur.sumstats.gz,./../important_sumstats/avg_forearm.sumstats.gz,./../important_sumstats/avg_humerus.sumstats.gz,./../important_sumstats/avg_tibia.sumstats.gz,./../important_sumstats/global_FI.sumstats.gz,./../important_sumstats/global_GMV.sumstats.gz,./../important_sumstats/global_IC.sumstats.gz,./../important_sumstats/global_SA.sumstats.gz,./../important_sumstats/height.sumstats.gz,./../important_sumstats/hip_width.sumstats.gz,./../important_sumstats/shoulder_width.sumstats.gz,./../important_sumstats/torso_length.sumstats.gz \
--ref-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--out new_height_rg
