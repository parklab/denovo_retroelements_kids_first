#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH --mem=32G
#SBATCH -o FM_N0D7XF4H_%j.out
#SBATCH --mail-type=NONE

REF=/n/data/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
SITES=HGDP_SGDP_1KG_germline_SVA_merged.vcf.orientation_site_level.txt
XTEA=/n/data/projects/XTEA/
OUT=cleavage_site_seqs

python ${XTEA}"x_cleavage_site.py" ${REF} ${SITES} ${OUT}
