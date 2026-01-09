#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
from typing import Dict, Optional

import pysam

_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def open_text(path: str):
    """Open plain text or gzipped text for reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_info(info_str: str) -> Dict[str, Optional[str]]:
    """
    Parse VCF INFO string like 'END=123;STRAND=-;FLAG;K=V'
    into dict: {"END":"123", "STRAND":"-", "FLAG":None, "K":"V"}
    """
    d: Dict[str, Optional[str]] = {}
    if not info_str or info_str == ".":
        return d
    for part in info_str.split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        else:
            d[part] = None
    return d


def normalize_chrom(chrom: str, ref: pysam.FastaFile) -> str:
    """Try common contig naming fixes ('1' <-> 'chr1')."""
    if chrom in ref.references:
        return chrom
    if chrom.startswith("chr"):
        alt = chrom[3:]
        if alt in ref.references:
            return alt
    else:
        alt = "chr" + chrom
        if alt in ref.references:
            return alt
    return chrom  # may still fail


def fetch_seq_or_na(
    ref: pysam.FastaFile,
    chrom: str,
    start_1based: int,
    end_1based: int,
    strand: str,
) -> str:
    """
    Fetch reference sequence for [start_1based, end_1based] inclusive.
    Return 'NA' on any failure.
    """
    try:
        chrom2 = normalize_chrom(chrom, ref)

        # Convert to 0-based half-open
        start0 = start_1based - 1
        end0_excl = end_1based  # END inclusive -> exclusive

        if start_1based <= 0 or end_1based <= 0 or end_1based < start_1based:
            return "NA"

        seq = ref.fetch(chrom2, start0, end0_excl)

        if strand == "-":
            seq = revcomp(seq)
        elif strand != "+":
            return "NA"

        # pysam may return "" if range is weird; treat as NA
        return seq if seq else "NA"
    except Exception:
        return "NA"


def write_fasta_record(fh, header: str, seq: str, wrap: int = 60) -> None:
    fh.write(f">{header}\n")
    if seq == "NA":
        fh.write("NA\n")
        return
    for i in range(0, len(seq), wrap):
        fh.write(seq[i : i + wrap] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="Headerless VCF-like file (tab-delimited), .vcf or .vcf.gz")
    ap.add_argument("--ref", required=True, help="Reference FASTA (needs .fai index: samtools faidx ref.fa)")
    ap.add_argument("--out_fasta", required=True, help="Output FASTA (same order/count as input)")
    ap.add_argument("--out_tsv", required=True, help="Output TSV with 2 columns: header, sequence (same order/count)")
    ap.add_argument("--default_strand", default="+", choices=["+", "-"], help="Used if STRAND missing (default '+')")
    args = ap.parse_args()

    ref = pysam.FastaFile(args.ref)

    with open_text(args.vcf) as f_in, open(args.out_fasta, "w") as f_fa, open(args.out_tsv, "w") as f_tsv:
        # If you do NOT want a TSV header line, delete the next line.
        f_tsv.write("header\tsequence\n")

        rec_idx = 0
        for lineno, line in enumerate(f_in, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            rec_idx += 1
            fields = line.split("\t")

            # Need at least 8 columns for standard VCF-like layout
            # CHROM POS ID REF ALT QUAL FILTER INFO
            header = f"record{rec_idx}"
            seq = "NA"

            try:
                if len(fields) < 8:
                    raise ValueError("expected >= 8 columns")

                chrom = fields[0]
                start_1based = int(fields[1])
                info = parse_info(fields[7])

                end_str = info.get("END")
                if end_str is None:
                    raise ValueError("missing END in INFO")
                end_1based = int(end_str)

                strand = info.get("STRAND") or args.default_strand

                header = f"{chrom}:{start_1based}-{end_1based}({strand})"
                seq = fetch_seq_or_na(ref, chrom, start_1based, end_1based, strand)

            except Exception:
                # Keep header deterministic even on failure (still in-order)
                # Try to include chrom:pos if possible, else fallback record index
                if len(fields) >= 2:
                    try:
                        header = f"{fields[0]}:{fields[1]}-NA(?)"
                    except Exception:
                        header = f"record{rec_idx}"
                else:
                    header = f"record{rec_idx}"
                seq = "NA"

            # Outputs always get one record per input record (in same order)
            write_fasta_record(f_fa, header, seq, wrap=60)
            f_tsv.write(f"{header}\t{seq}\n")

    ref.close()


if __name__ == "__main__":
    main()

#python extract_segments_headerless.py \
#  --vcf segments.noheader.vcf.gz \
#  --ref ref.fa \
#  --out_fasta segments.fa \
#  --out_tsv segments.tsv
