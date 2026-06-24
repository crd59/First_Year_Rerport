#!/bin/bash

module load miniconda3

source activate ldsc


./../../ldsc/ldsc.py \
--rg ./../important_sumstats/avg_femur.sumstats.gz,./../important_sumstats/avg_forearm.sumstats.gz,./../important_sumstats/avg_humerus.sumstats.gz,./../important_sumstats/avg_tibia.sumstats.gz,./../important_sumstats/birth_length.sumstats.gz,./../important_sumstats/global_FI.sumstats.gz,./../important_sumstats/global_GMV.sumstats.gz,./../important_sumstats/global_IC.sumstats.gz,./../important_sumstats/global_SA.sumstats.gz,./../important_sumstats/HC.sumstats.gz,./../important_sumstats/height.sumstats.gz,./../important_sumstats/hip_width.sumstats.gz,./../important_sumstats/shoulder_width.sumstats.gz,./../important_sumstats/torso_length.sumstats.gz \
--ref-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--out avg_femur_rg

./../../ldsc/ldsc.py \
--rg ./../important_sumstats/avg_forearm.sumstats.gz,./../important_sumstats/avg_humerus.sumstats.gz,./../important_sumstats/avg_tibia.sumstats.gz,./../important_sumstats/birth_length.sumstats.gz,./../important_sumstats/global_FI.sumstats.gz,./../important_sumstats/global_GMV.sumstats.gz,./../important_sumstats/global_IC.sumstats.gz,./../important_sumstats/global_SA.sumstats.gz,./../important_sumstats/HC.sumstats.gz,./../important_sumstats/height.sumstats.gz,./../important_sumstats/hip_width.sumstats.gz,./../important_sumstats/shoulder_width.sumstats.gz,./../important_sumstats/torso_length.sumstats.gz \
--ref-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--out avg_forearm_rg

./../../ldsc/ldsc.py \
--rg ./../important_sumstats/avg_humerus.sumstats.gz,./../important_sumstats/avg_tibia.sumstats.gz,./../important_sumstats/birth_length.sumstats.gz,./../important_sumstats/global_FI.sumstats.gz,./../important_sumstats/global_GMV.sumstats.gz,./../important_sumstats/global_IC.sumstats.gz,./../important_sumstats/global_SA.sumstats.gz,./../important_sumstats/HC.sumstats.gz,./../important_sumstats/height.sumstats.gz,./../important_sumstats/hip_width.sumstats.gz,./../important_sumstats/shoulder_width.sumstats.gz,./../important_sumstats/torso_length.sumstats.gz \
--ref-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--out avg_humerus_rg

./../../ldsc/ldsc.py \
--rg ./../important_sumstats/avg_tibia.sumstats.gz,./../important_sumstats/birth_length.sumstats.gz,./../important_sumstats/global_FI.sumstats.gz,./../important_sumstats/global_GMV.sumstats.gz,./../important_sumstats/global_IC.sumstats.gz,./../important_sumstats/global_SA.sumstats.gz,./../important_sumstats/HC.sumstats.gz,./../important_sumstats/height.sumstats.gz,./../important_sumstats/hip_width.sumstats.gz,./../important_sumstats/shoulder_width.sumstats.gz,./../important_sumstats/torso_length.sumstats.gz \
--ref-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--out avg_tibia_rg

./../../ldsc/ldsc.py \
--rg ./../important_sumstats/birth_length.sumstats.gz,./../important_sumstats/global_FI.sumstats.gz,./../important_sumstats/global_GMV.sumstats.gz,./../important_sumstats/global_IC.sumstats.gz,./../important_sumstats/global_SA.sumstats.gz,./../important_sumstats/HC.sumstats.gz,./../important_sumstats/height.sumstats.gz,./../important_sumstats/hip_width.sumstats.gz,./../important_sumstats/shoulder_width.sumstats.gz,./../important_sumstats/torso_length.sumstats.gz \
--ref-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--out birth_length_rg

./../../ldsc/ldsc.py \
--rg ./../important_sumstats/global_FI.sumstats.gz,./../important_sumstats/global_GMV.sumstats.gz,./../important_sumstats/global_IC.sumstats.gz,./../important_sumstats/global_SA.sumstats.gz,./../important_sumstats/HC.sumstats.gz,./../important_sumstats/height.sumstats.gz,./../important_sumstats/hip_width.sumstats.gz,./../important_sumstats/shoulder_width.sumstats.gz,./../important_sumstats/torso_length.sumstats.gz \
--ref-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--out FI_rg

./../../ldsc/ldsc.py \
--rg ./../important_sumstats/global_GMV.sumstats.gz,./../important_sumstats/global_IC.sumstats.gz,./../important_sumstats/global_SA.sumstats.gz,./../important_sumstats/HC.sumstats.gz,./../important_sumstats/height.sumstats.gz,./../important_sumstats/hip_width.sumstats.gz,./../important_sumstats/shoulder_width.sumstats.gz,./../important_sumstats/torso_length.sumstats.gz \
--ref-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--out GMV_rg

./../../ldsc/ldsc.py \
--rg ./../important_sumstats/global_IC.sumstats.gz,./../important_sumstats/global_SA.sumstats.gz,./../important_sumstats/HC.sumstats.gz,./../important_sumstats/height.sumstats.gz,./../important_sumstats/hip_width.sumstats.gz,./../important_sumstats/shoulder_width.sumstats.gz,./../important_sumstats/torso_length.sumstats.gz \
--ref-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--out IC_rg

./../../ldsc/ldsc.py \
--rg ./../important_sumstats/global_SA.sumstats.gz,./../important_sumstats/HC.sumstats.gz,./../important_sumstats/height.sumstats.gz,./../important_sumstats/hip_width.sumstats.gz,./../important_sumstats/shoulder_width.sumstats.gz,./../important_sumstats/torso_length.sumstats.gz \
--ref-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--out SA_rg

./../../ldsc/ldsc.py \
--rg ./../important_sumstats/HC.sumstats.gz,./../important_sumstats/height.sumstats.gz,./../important_sumstats/hip_width.sumstats.gz,./../important_sumstats/shoulder_width.sumstats.gz,./../important_sumstats/torso_length.sumstats.gz \
--ref-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--out HC_rg

./../../ldsc/ldsc.py \
--rg ./../important_sumstats/height.sumstats.gz,./../important_sumstats/hip_width.sumstats.gz,./../important_sumstats/shoulder_width.sumstats.gz,./../important_sumstats/torso_length.sumstats.gz \
--ref-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--out height_rg

./../../ldsc/ldsc.py \
--rg ./../important_sumstats/hip_width.sumstats.gz,./../important_sumstats/shoulder_width.sumstats.gz,./../important_sumstats/torso_length.sumstats.gz \
--ref-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--out hip_width_rg

./../../ldsc/ldsc.py \
--rg ./../important_sumstats/shoulder_width.sumstats.gz,./../important_sumstats/torso_length.sumstats.gz \
--ref-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--w-ld-chr ./../../../../../ldsc/eur_w_ld_chr/ \
--out shoulder_width_rg

