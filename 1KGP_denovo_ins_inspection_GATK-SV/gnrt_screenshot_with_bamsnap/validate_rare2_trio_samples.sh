#!/usr/bin/env bash
set -euo pipefail

# Validate that trio samples in the task file can be found in either:
#   1. local override list, or
#   2. 1000G_highcov_sample_map.tsv
#
# Usage:
#   bash validate_rare2_trio_samples.sh [tasks.tsv] [sample_map.tsv] [local_override_list]

TASKS_TSV="${1:-./trio_bamsnap_tasks_rare2.lf.tsv}"
SAMPLE_MAP="${2:-./1000G_highcov_sample_map.tsv}"
LOCAL_OVERRIDE_LIST="${3:-./cram_list_child.txt}"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

awk -F'\t' 'NR>1 {
  gsub(/\r/,"");
  print $6; print $7; print $8
}' "$TASKS_TSV" | sort -u > "$tmpdir/needed.txt"

: > "$tmpdir/local_samples.txt"
if [[ -s "$LOCAL_OVERRIDE_LIST" ]]; then
  awk '{
    gsub(/\r/,"");
    path=$0;
    n=split(path,a,"/");
    file=a[n];
    sub(/\.(final\.)?(cram|bam)$/,"",file);
    print file;
  }' "$LOCAL_OVERRIDE_LIST" | sort -u > "$tmpdir/local_samples.txt"
fi

awk -F'\t' 'NR>1 {gsub(/\r/,""); print $1}' "$SAMPLE_MAP" | sort -u > "$tmpdir/aws_samples.txt"

cat "$tmpdir/local_samples.txt" "$tmpdir/aws_samples.txt" | sort -u > "$tmpdir/available.txt"

comm -23 "$tmpdir/needed.txt" "$tmpdir/available.txt" > missing_rare2_trio_samples.txt

echo "[INFO] Needed samples: $(wc -l < "$tmpdir/needed.txt")"
echo "[INFO] Available samples by local override/AWS map: $(wc -l < "$tmpdir/available.txt")"
echo "[INFO] Missing samples: $(wc -l < missing_rare2_trio_samples.txt)"
echo "[INFO] Missing sample list: missing_rare2_trio_samples.txt"
