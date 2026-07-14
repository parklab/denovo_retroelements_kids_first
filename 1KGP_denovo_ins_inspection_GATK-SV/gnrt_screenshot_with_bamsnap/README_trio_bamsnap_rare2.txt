Updated scripts for denovo_insertion_candidates_rare_2.tsv.

Input rows: 2833
Valid tasks written: 2833
Bad rows skipped: 0
Unique trio samples: 1676

Main files:
  trio_bamsnap_tasks_rare2.lf.tsv
  trio_bamsnap_worker.rare2.v1.sh
  submit_trio_bamsnap_array.rare2.v1.slurm
  launch_trio_bamsnap_rare2.v1.sh
  test_one_trio_bamsnap_rare2.v1.sh
  validate_rare2_trio_samples.sh
  trio_bamsnap_required_samples_rare2.txt
  trio_bamsnap_tasks_rare2_summary.tsv

Defaults:
  REF=./REFERENCE/hg38_no_alt.fa
  OUTDIR=./OneKGP/trio_cram_screenshot/gatk_sv_ins
  TMP_ROOT=./OneKGP/trio_cram_screenshot/gatk_sv_ins/tmp
  LOCAL_OVERRIDE_LIST=./cram_list_child.txt
  TASKS_TSV=./trio_bamsnap_tasks_rare2.lf.tsv

Recommended run:

  bash fix_crlf_for_trio_bamsnap.sh

  # Build or refresh the AWS sample map if needed:
  rm -f 1000G_highcov_sample_map.tsv
  bash build_1000g_aws_sample_map.v2.sh $PWD

  # Optional but recommended:
  bash validate_rare2_trio_samples.sh

  # Single-task smoke test:
  bash test_one_trio_bamsnap_rare2.v1.sh

  # Submit all tasks:
  bash launch_trio_bamsnap_rare2.v1.sh

The launcher dynamically submits:
  --array=1-2833%30

The .slurm file itself is also set to:
  #SBATCH --array=1-2833%30
