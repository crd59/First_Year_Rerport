#!/bin/bash
#SBATCH -A WARRIER-SL2-CPU
#SBATCH -J GWAS_FLAMES
#SBATCH -D /rds/user/CAMID/rds-genetics_hpc-Nl99R8pHODQ/toolbox/CADD
#SBATCH -o /rds/user/CAMID/hpc-work/pgs_img_2/SSHscript/GWAS_FLAME/step0_CADD_download_v1_%A.log
#SBATCH -p icelake                           # on icelake partition
#SBATCH -t 12:00:00                            # HH:MM:SS with maximum 12:00:00 for SL3 or 36:00:00 for SL2
#SBATCH --mem=1G
#SBATCH --mail-type=ALL

wget -c https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz

echo "Download completed at $(date)"