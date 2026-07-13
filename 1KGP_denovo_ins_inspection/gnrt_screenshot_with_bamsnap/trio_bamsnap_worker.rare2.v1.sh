#!/usr/bin/env bash
set -euo pipefail

# Process one trio screenshot task.
# Compatible with bamsnap v1.
# This script:
#   1. reads one task from TASKS_TSV
#   2. resolves child/father/mother CRAMs from a local override list or the AWS sample map
#   3. extracts only the focal region using samtools
#   4. creates one 3-track bamsnap image
#
# Usage:
#   bash trio_bamsnap_worker.rare2.v1.sh TASK_ID TASKS_TSV SAMPLE_MAP REF OUTDIR TMP_ROOT [LOCAL_OVERRIDE_LIST]

TASK_ID="$1"
TASKS_TSV="$2"
SAMPLE_MAP="$3"
REF="$4"
OUTDIR="$5"
TMP_ROOT="$6"
LOCAL_OVERRIDE_LIST="${7:-}"

mkdir -p "$OUTDIR" "$TMP_ROOT"
WORKDIR="${TMP_ROOT}/task_${TASK_ID}"
mkdir -p "$WORKDIR"
trap 'rm -rf "$WORKDIR"' EXIT

strip_cr() {
  printf '%s' "$1" | tr -d '\r'
}

line=$(awk -F'\t' -v id="$TASK_ID" 'NR>1 && $1==id {gsub(/\r/,""); print; exit}' "$TASKS_TSV")
line=$(strip_cr "$line")

if [[ -z "$line" ]]; then
  echo "[ERROR] Task ID ${TASK_ID} not found in ${TASKS_TSV}" >&2
  exit 1
fi

IFS=$'\t' read -r task_id chrom pos start end child father mother output_prefix <<< "$line"

task_id=$(strip_cr "$task_id")
chrom=$(strip_cr "$chrom")
pos=$(strip_cr "$pos")
start=$(strip_cr "$start")
end=$(strip_cr "$end")
child=$(strip_cr "$child")
father=$(strip_cr "$father")
mother=$(strip_cr "$mother")
output_prefix=$(strip_cr "$output_prefix")

region="${chrom}:${start}-${end}"
out_png="${OUTDIR}/${output_prefix}.png"

if [[ -s "$out_png" ]]; then
  echo "[INFO] Output already exists, skip: $out_png"
  exit 0
fi

LOCAL_MAP="${WORKDIR}/local_override.map"
if [[ -n "$LOCAL_OVERRIDE_LIST" && -s "$LOCAL_OVERRIDE_LIST" ]]; then
  awk 'BEGIN{OFS="\t"} {
      gsub(/\r/,"");
      path=$0;
      n=split(path,a,"/");
      file=a[n];
      sub(/\.(final\.)?(cram|bam)$/,"",file);
      print file, path;
  }' "$LOCAL_OVERRIDE_LIST" > "$LOCAL_MAP"
fi

lookup_src() {
  local sample="$1"
  local src=""

  if [[ -s "$LOCAL_MAP" ]]; then
    src=$(awk -F'\t' -v s="$sample" '{
      gsub(/\r/,"");
      if ($1==s) {print $2; exit}
    }' "$LOCAL_MAP") || true
  fi

  if [[ -z "$src" ]]; then
    src=$(awk -F'\t' -v s="$sample" 'NR>1 {
      gsub(/\r/,"");
      if ($1==s) {print $4; exit}
    }' "$SAMPLE_MAP") || true
  fi

  src=$(strip_cr "$src")

  if [[ -z "$src" ]]; then
    echo "[ERROR] No CRAM/BAM source found for sample: $sample" >&2
    return 1
  fi

  printf '%s\n' "$src"
}

extract_one() {
  local sample="$1"
  local src="$2"
  local outbam="$3"

  sample=$(strip_cr "$sample")
  src=$(strip_cr "$src")
  outbam=$(strip_cr "$outbam")

  echo "[INFO] Extracting ${sample} from ${src} at ${region}"
  samtools view \
    -T "$REF" \
    -b \
    -o "$outbam" \
    "$src" \
    "$region"

  samtools index "$outbam"
}

child_src=$(lookup_src "$child")
father_src=$(lookup_src "$father")
mother_src=$(lookup_src "$mother")

child_src=$(strip_cr "$child_src")
father_src=$(strip_cr "$father_src")
mother_src=$(strip_cr "$mother_src")

child_bam="${WORKDIR}/${child}.${chrom}_${pos}.bam"
father_bam="${WORKDIR}/${father}.${chrom}_${pos}.bam"
mother_bam="${WORKDIR}/${mother}.${chrom}_${pos}.bam"

extract_one "$child" "$child_src" "$child_bam"
extract_one "$father" "$father_src" "$father_bam"
extract_one "$mother" "$mother_src" "$mother_bam"

echo "[INFO] Running bamsnap v1 for trio ${child}, ${father}, ${mother}"

bamsnap \
  --bam "$child_bam" \
  --bam "$father_bam" \
  --bam "$mother_bam" \
  --title "${child}_child" \
  --title "${father}_father" \
  --title "${mother}_mother" \
  --pos "$region" \
  --draw bamplot \
  --bamplot read \
  --margin 400 \
  --plot-margin-top 0 \
  --plot-margin-bottom 0 \
  --separator-height 150 \
  --show-soft-clipped \
  --read-color-by interchrom \
  --save-image-only \
  --ref "$REF" \
  --out "$out_png"

echo "[INFO] Finished: $out_png"
