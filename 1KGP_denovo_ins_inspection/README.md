# Extract Candidate De Novo Insertions from a Multi-sample VCF

This repository identifies **candidate rare de novo insertions** in parent–child trios from a multi-sample VCF, such as the high-coverage 1000 Genomes structural-variant callset.

The script is intended for candidate discovery. Calls should be validated with read-level evidence and additional quality control before being treated as confirmed de novo events.

## Files

- `extract_denovo_insertions.py` — extraction and filtering program.
- `run_extract_v1_rare_ones.sh` — example command for the 1000 Genomes VCF.

## Requirements

- Python 3
- [`pysam`](https://pysam.readthedocs.io/)

Install the Python dependency with:

```bash
python3 -m pip install pysam
```

## Input files

### Multi-sample VCF

The VCF must contain genotypes for the children and both parents. Cohort allele frequency is read from:

1. `INFO/AF`, or
2. `INFO/AC` divided by `INFO/AN` when `AF` is unavailable.

Insertion alleles are recognized from symbolic ALT alleles such as `<INS>`, `<INS:ME:LINE1>`, and `<INS:ME:SVA>`, from `INFO/SVTYPE=INS`, or from explicit ALT sequences longer than the REF allele.

### Trio file

A whitespace-delimited file with at least the first three columns in this order:

```text
sampleID fatherID motherID sex
HG_CHILD HG_FATHER HG_MOTHER 1
```

The fourth column is optional and ignored. Rows with father or mother recorded as `0` are skipped. All three samples must be present in the VCF.

## Usage

Edit the input paths in `run_extract_v1_rare_ones.sh`, then run:

```bash
bash run_extract_v1_rare_ones.sh
```

Equivalent direct command (the output name uses `rare` because AF 0.01 is not a strict singleton threshold):

```bash
python3 extract_denovo_insertions.py \
  --vcf 1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz \
  --trios 1kGP.3202_samples.pedigree_info.txt \
  --out denovo_insertion_candidates_rare.tsv \
  --pass-only \
  --min-child-gq 20 \
  --min-parent-gq 20 \
  --max-cohort-af 0.01
```

## Candidate definition

For each insertion record and complete trio, the script requires:

1. The record passes the VCF filter when `--pass-only` is used.
2. The child has a fully called genotype and meets `--min-child-gq`.
3. The child carries an insertion ALT allele.
4. The insertion allele has cohort AF less than or equal to `--max-cohort-af`.
5. Both parents have fully called genotypes, meet `--min-parent-gq`, and lack the relevant ALT allele(s).

By default, each parent must have a reference genotype at the entire multiallelic record. With `--allele-specific`, a parent may carry another ALT allele but must not carry the same insertion allele found in the child.

Missing parental genotypes fail by default. `--allow-missing-parent-as-absent` changes this behavior but is not recommended for high-confidence analysis. Missing cohort AF also fails by default unless `--keep-missing-cohort-af` is specified.

## Output

The output is a tab-delimited file with one row per candidate variant–trio combination. It includes:

- Variant coordinates, ID, REF, ALT, insertion ALT, SV type, SV length, and FILTER status
- Cohort AF and its source (`INFO_AF` or `INFO_AC_AN`)
- Child, father, and mother sample IDs
- Genotypes and GQ values for all three samples

A processing summary is written to standard error, including the numbers of valid trios, scanned VCF records, insertion records, and candidate rows.

## Important notes

- `--max-cohort-af 0.01` means **AF ≤ 1%** and is a rare-variant filter, not a strict singleton filter. In 3,202 diploid samples, a single observed allele has AF of approximately `1 / 6404 = 0.000156`, although the exact value depends on `AN`.
- The script uses genotype and GQ filters but does not filter on read depth, allele balance, split reads, discordant pairs, or Mendelian likelihood.
- It does not verify pedigree relationships or sample identity.
- Explicit sequence insertions include any ALT allele longer than REF, so small insertion indels may also be retained.
