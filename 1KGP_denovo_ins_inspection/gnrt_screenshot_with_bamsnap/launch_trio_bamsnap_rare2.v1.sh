#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TASKS_TSV="${TASKS_TSV:-./trio_bamsnap_tasks_rare2.lf.tsv}"
REF="${REF:-./REFERENCE/hg38_no_alt.fa}"
OUTDIR="${OUTDIR:-./OneKGP/trio_cram_screenshot/gatk_sv_ins}"
TMP_ROOT="${TMP_ROOT:-./OneKGP/trio_cram_screenshot/gatk_sv_ins/tmp}"
LOCAL_OVERRIDE_LIST="${LOCAL_OVERRIDE_LIST:-./cram_list_child.txt}"

if [[ ! -s 1000G_highcov_sample_map.tsv ]]; then
  bash "${SCRIPT_DIR}/build_1000g_aws_sample_map.v2.sh" "$PWD"
else
  sed -i 's/\r$//' 1000G_highcov_sample_map.tsv
fi

[[ -s "$TASKS_TSV" ]] && sed -i 's/\r$//' "$TASKS_TSV"
[[ -s "$LOCAL_OVERRIDE_LIST" ]] && sed -i 's/\r$//' "$LOCAL_OVERRIDE_LIST"

ntasks=$(awk 'END{print NR-1}' "$TASKS_TSV")
echo "[INFO] Number of tasks: ${ntasks}"

echo "[INFO] Submitting array job ..."
sbatch \
  --array="1-${ntasks}%30" \
  --export=ALL,TASKS_TSV="$TASKS_TSV",SAMPLE_MAP="$PWD/1000G_highcov_sample_map.tsv",REF="$REF",OUTDIR="$OUTDIR",TMP_ROOT="$TMP_ROOT",LOCAL_OVERRIDE_LIST="$LOCAL_OVERRIDE_LIST",WORKER="$SCRIPT_DIR/trio_bamsnap_worker.rare2.v1.sh" \
  "${SCRIPT_DIR}/submit_trio_bamsnap_array.rare2.v1.slurm"
