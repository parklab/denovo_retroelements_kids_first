python extract_denovo_insertions.py \
  --vcf 1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz \
  --trios /n/no_backup2/dbmi/park/simon_chu/1kg_698/1kGP.3202_samples.pedigree_info.txt \
  --out denovo_insertion_candidates_singleton.tsv \
  --min-child-gq 0 \
  --min-parent-gq 0 \
  --max-cohort-af 0.01
