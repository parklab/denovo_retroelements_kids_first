bash fix_crlf_for_trio_bamsnap.sh
rm -f 1000G_highcov_sample_map.tsv
bash build_1000g_aws_sample_map.v2.sh $PWD
bash validate_rare2_trio_samples.sh
bash test_one_trio_bamsnap_rare2.v1.sh
bash launch_trio_bamsnap_rare2.v1.sh
