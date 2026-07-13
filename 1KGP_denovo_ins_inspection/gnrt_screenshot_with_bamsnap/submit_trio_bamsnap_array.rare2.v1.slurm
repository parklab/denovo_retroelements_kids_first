#!/usr/bin/env bash
#SBATCH --job-name=trio_bamsnap_r2
#SBATCH --output=logs/trio_bamsnap_r2_%A_%a.out
#SBATCH --error=logs/trio_bamsnap_r2_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-2833%30

set -euo pipefail

TASKS_TSV="${TASKS_TSV:-./trio_bamsnap_tasks_rare2.lf.tsv}"
SAMPLE_MAP="${SAMPLE_MAP:-$PWD/1000G_highcov_sample_map.tsv}"
REF="${REF:-./REFERENCE/hg38_no_alt.fa}"
OUTDIR="${OUTDIR:-./OneKGP/trio_cram_screenshot/gatk_sv_ins}"
TMP_ROOT="${TMP_ROOT:-./OneKGP/trio_cram_screenshot/gatk_sv_ins/tmp}"
LOCAL_OVERRIDE_LIST="${LOCAL_OVERRIDE_LIST:-./cram_list_child.txt}"
WORKER="${WORKER:-$PWD/trio_bamsnap_worker.rare2.v1.sh}"

mkdir -p logs "$OUTDIR" "$TMP_ROOT"

bash "$WORKER" \
  "$SLURM_ARRAY_TASK_ID" \
  "$TASKS_TSV" \
  "$SAMPLE_MAP" \
  "$REF" \
  "$OUTDIR" \
  "$TMP_ROOT" \
  "$LOCAL_OVERRIDE_LIST"
