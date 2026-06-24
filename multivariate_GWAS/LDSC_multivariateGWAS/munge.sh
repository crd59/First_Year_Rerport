#!/bin/bash

module load miniconda3

source activate ldsc

./../../../ldsc/munge_sumstats.py --sumstats ./../TrunkWidths_GWAS.tsv --out ./trunk_munged --merge-alleles ./../../../../../../ldsc/w_hm3.snplist --chunksize 500000 --p P --a1 A1 --a2 A2 --snp SNP --frq MAF --N 54881
./../../../ldsc/munge_sumstats.py --sumstats ./../LimbHeight_GWAS.tsv --out ./limb_munged --merge-alleles ./../../../../../../ldsc/w_hm3.snplist --chunksize 500000 --p P --a1 A1 --a2 A2 --snp SNP --frq MAF --N 39556
./../../../ldsc/munge_sumstats.py --sumstats ./../Cortical_GWAS.tsv --out ./cortical_munged --merge-alleles ./../../../../../../ldsc/w_hm3.snplist --chunksize 500000 --p P --a1 A1 --a2 A2 --snp SNP --frq MAF --N 18037
