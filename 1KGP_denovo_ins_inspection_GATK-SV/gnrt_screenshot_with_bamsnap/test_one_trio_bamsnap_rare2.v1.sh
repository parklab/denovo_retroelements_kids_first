#!/usr/bin/env bash
set -euo pipefail

# Quick one-task test for the new rare2 insertion list.
bash ./trio_bamsnap_worker.rare2.v1.sh 1 \
  ./trio_bamsnap_tasks_rare2.lf.tsv \
  ./1000G_highcov_sample_map.tsv \
  ./REFERENCE/hg38_no_alt.fa \
  ./OneKGP/trio_cram_screenshot/gatk_sv_ins \
  ./OneKGP/trio_cram_screenshot/gatk_sv_ins/tmp \
  ./cram_list_child.txt
