#!/usr/bin/env python3

import argparse
import csv
import sys
from typing import List, Set, Tuple

import pysam


Trio = Tuple[str, str, str]  # child, father, mother


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Extract potential de novo rare insertions from a multi-sample VCF "
            "using trio metadata."
        )
    )

    p.add_argument(
        "--vcf",
        required=True,
        help="Input integrated VCF, e.g. 1KGP_integrated.vcf.gz",
    )

    p.add_argument(
        "--trios",
        required=True,
        help="Whitespace-delimited trio file: sampleID fatherID motherID sex",
    )

    p.add_argument(
        "--out",
        required=True,
        help="Output TSV of candidate de novo rare insertion calls",
    )

    p.add_argument(
        "--pass-only",
        action="store_true",
        help="Keep only records with FILTER=PASS or FILTER='.'",
    )

    p.add_argument(
        "--min-child-gq",
        type=float,
        default=0,
        help="Minimum child GQ required to call the insertion present. Default: 0",
    )

    p.add_argument(
        "--min-parent-gq",
        type=float,
        default=0,
        help="Minimum parent GQ required to call the insertion absent. Default: 0",
    )

    p.add_argument(
        "--max-cohort-af",
        type=float,
        default=0.01,
        help=(
            "Maximum cohort allele frequency allowed for the child's insertion allele. "
            "Default: 0.01, meaning keep AF <= 1%% and filter out AF > 1%%."
        ),
    )

    p.add_argument(
        "--keep-missing-cohort-af",
        action="store_true",
        help=(
            "Keep candidates when cohort AF cannot be read from INFO/AF or computed "
            "from AC/AN. By default, missing AF fails the AF filter."
        ),
    )

    p.add_argument(
        "--allow-missing-parent-as-absent",
        action="store_true",
        help=(
            "Treat missing parental GT as absent. Not recommended for high-confidence "
            "de novo calls."
        ),
    )

    p.add_argument(
        "--allele-specific",
        action="store_true",
        help=(
            "For multiallelic records, require parents to lack only the same insertion "
            "ALT allele carried by the child. By default, parents must be homozygous "
            "reference at the whole record."
        ),
    )

    return p.parse_args()


def load_trios(path: str) -> List[Trio]:
    trios: List[Trio] = []

    with open(path) as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 3:
                raise ValueError(
                    f"Line {line_no} in trio file has fewer than 3 columns: {line}"
                )

            # Skip header like:
            # sampleID fatherID motherID sex
            if fields[0].lower() in {
                "sampleid",
                "sample_id",
                "iid",
                "child",
                "childid",
            }:
                continue

            child = fields[0]
            father = fields[1]
            mother = fields[2]

            # Keep complete trios only.
            if father != "0" and mother != "0":
                trios.append((child, father, mother))

    return trios


def as_tuple(x):
    if x is None:
        return ()
    if isinstance(x, (tuple, list)):
        return x
    return (x,)


def to_float(x):
    if x is None:
        return None
    try:
        return float(x)
    except Exception:
        return None


def info_to_str(value):
    if value is None:
        return "."
    if isinstance(value, (tuple, list)):
        return ",".join(str(x) for x in value)
    return str(value)


def gt_to_str(call) -> str:
    gt = call.get("GT")

    if gt is None:
        return "."

    sep = "|" if getattr(call, "phased", False) else "/"

    return sep.join("." if a is None else str(a) for a in gt)


def gq_to_str(call) -> str:
    gq = call.get("GQ")
    return "." if gq is None else str(gq)


def gq_ok(call, min_gq: float) -> bool:
    if min_gq <= 0:
        return True

    gq = call.get("GQ")

    if gq is None:
        return False

    try:
        return float(gq) >= min_gq
    except Exception:
        return False


def gt_is_fully_called(call) -> bool:
    gt = call.get("GT")
    return gt is not None and len(gt) > 0 and all(a is not None for a in gt)


def alt_alleles_in_gt(call) -> Set[int]:
    """
    Return non-reference allele numbers in GT.

    Examples:
      0/0 -> set()
      0/1 -> {1}
      1/1 -> {1}
      1/2 -> {1, 2}
      ./1 -> {1}, but the caller should separately decide whether missing GT is valid.
    """

    gt = call.get("GT")

    if gt is None:
        return set()

    alleles = set()

    for a in gt:
        if a is not None and a > 0:
            alleles.add(a)

    return alleles


def insertion_alt_indices(record) -> Set[int]:
    """
    Identify which ALT allele numbers are insertions.

    Handles:
      ALT=<INS>
      ALT=<INS:ME>
      ALT=<INS:ME:ALU>
      ALT=<INS:ME:LINE1>
      ALT=<INS:ME:SVA>
      ALT=<INS:UNK>
      INFO/SVTYPE=INS
      explicit sequence insertions where len(ALT) > len(REF)
    """

    alts = record.alts or ()
    ins_indices: Set[int] = set()

    for i, alt in enumerate(alts, start=1):
        if alt is None:
            continue

        alt_upper = alt.upper()

        # Symbolic insertion alleles.
        if alt_upper.startswith("<INS") or alt_upper == "INS":
            ins_indices.add(i)
            continue

        # Explicit sequence insertion, if present.
        if (
            not alt_upper.startswith("<")
            and record.ref is not None
            and len(alt) > len(record.ref)
        ):
            ins_indices.add(i)

    svtype = record.info.get("SVTYPE")
    svtypes = set()

    if svtype is not None:
        if isinstance(svtype, (tuple, list)):
            svtypes = {str(x).upper() for x in svtype}
        else:
            svtypes = {str(svtype).upper()}

    # If record-level SVTYPE says insertion, treat all ALT alleles as insertion
    # alleles when ALT parsing did not already find them.
    if not ins_indices and any(x.startswith("INS") for x in svtypes):
        ins_indices = set(range(1, len(alts) + 1))

    return ins_indices


def estimate_cohort_af(record, alt_alleles: Set[int]):
    """
    Estimate cohort allele frequency for the child's insertion ALT allele(s).

    Priority:
      1. INFO/AF, Number=A
      2. INFO/AC divided by INFO/AN

    For multiallelic records, returns the maximum AF among the child's carried
    insertion ALT alleles.

    Returns:
      cohort_af, cohort_af_source

    cohort_af is None if it cannot be calculated.
    """

    af = as_tuple(record.info.get("AF"))

    if af:
        vals = []

        for alt_idx in alt_alleles:
            i = alt_idx - 1

            if 0 <= i < len(af):
                v = to_float(af[i])

                if v is not None:
                    vals.append(v)

        if vals:
            return max(vals), "INFO_AF"

    ac = as_tuple(record.info.get("AC"))
    an = to_float(record.info.get("AN"))

    if ac and an is not None and an > 0:
        vals = []

        for alt_idx in alt_alleles:
            i = alt_idx - 1

            if 0 <= i < len(ac):
                v = to_float(ac[i])

                if v is not None:
                    vals.append(v / an)

        if vals:
            return max(vals), "INFO_AC_AN"

    return None, "."


def passes_filter(record) -> bool:
    filters = set(record.filter.keys())
    return len(filters) == 0 or filters == {"PASS"}


def absent_for_alleles(
    call,
    forbidden_alleles: Set[int],
    min_gq: float,
    allow_missing: bool,
) -> bool:
    """
    Return True when the sample has no forbidden ALT alleles.

    By default, missing GT is not considered absent.
    """

    if not gt_is_fully_called(call):
        return allow_missing

    if not gq_ok(call, min_gq):
        return False

    return len(alt_alleles_in_gt(call) & forbidden_alleles) == 0


def filter_to_str(record) -> str:
    filters = list(record.filter.keys())

    if len(filters) == 0:
        return "."

    return ";".join(filters)


def alt_for_indices(record, indices: Set[int]) -> str:
    alts = record.alts or ()
    vals = []

    for idx in sorted(indices):
        if 1 <= idx <= len(alts):
            vals.append(alts[idx - 1])

    return ",".join(vals) if vals else "."


def main():
    args = parse_args()

    trios = load_trios(args.trios)

    if not trios:
        sys.exit("No complete trios found in trio metadata.")

    # First open the VCF only to inspect sample names.
    with pysam.VariantFile(args.vcf) as vf:
        vcf_samples = set(vf.header.samples)

    valid_trios: List[Trio] = []
    missing_trios: List[Trio] = []

    for child, father, mother in trios:
        if child in vcf_samples and father in vcf_samples and mother in vcf_samples:
            valid_trios.append((child, father, mother))
        else:
            missing_trios.append((child, father, mother))

    if not valid_trios:
        sys.exit("No complete trios have child/father/mother all present in the VCF.")

    trio_samples = sorted({sample for trio in valid_trios for sample in trio})

    n_records = 0
    n_ins_records = 0
    n_candidates = 0

    out_fields = [
        "chrom",
        "pos",
        "end",
        "id",
        "ref",
        "alt",
        "child_insertion_alt",
        "svtype",
        "svlen",
        "filter",
        "cohort_af",
        "cohort_af_source",
        "child",
        "father",
        "mother",
        "child_gt",
        "father_gt",
        "mother_gt",
        "child_gq",
        "father_gq",
        "mother_gq",
    ]

    with pysam.VariantFile(args.vcf) as vf:
        # Decode only trio sample columns. This is much faster than decoding all
        # ~3000 samples and is sufficient because the rare-frequency filter uses
        # INFO/AF or INFO/AC/AN.
        vf.subset_samples(trio_samples)

        with open(args.out, "w", newline="") as out_f:
            writer = csv.DictWriter(out_f, fieldnames=out_fields, delimiter="\t")
            writer.writeheader()

            for rec in vf:
                n_records += 1

                if args.pass_only and not passes_filter(rec):
                    continue

                ins_alleles = insertion_alt_indices(rec)

                if not ins_alleles:
                    continue

                n_ins_records += 1

                all_alt_alleles = set(range(1, len(rec.alts or ()) + 1))

                for child, father, mother in valid_trios:
                    child_call = rec.samples[child]
                    father_call = rec.samples[father]
                    mother_call = rec.samples[mother]

                    # Child must have a fully called genotype.
                    if not gt_is_fully_called(child_call):
                        continue

                    if not gq_ok(child_call, args.min_child_gq):
                        continue

                    # Child must carry one of the insertion ALT alleles.
                    child_ins_alleles = alt_alleles_in_gt(child_call) & ins_alleles

                    if not child_ins_alleles:
                        continue

                    # Rare insertion filter:
                    # keep only AF <= max_cohort_af.
                    cohort_af, cohort_af_source = estimate_cohort_af(
                        rec,
                        child_ins_alleles,
                    )

                    if args.max_cohort_af is not None:
                        if cohort_af is None:
                            if not args.keep_missing_cohort_af:
                                continue
                        else:
                            if cohort_af > args.max_cohort_af:
                                continue

                    # Default behavior:
                    # parents must be 0/0 at the whole record.
                    #
                    # With --allele-specific:
                    # parents only need to lack the same insertion ALT allele
                    # carried by the child.
                    if args.allele_specific:
                        forbidden_alleles = child_ins_alleles
                    else:
                        forbidden_alleles = all_alt_alleles

                    father_absent = absent_for_alleles(
                        father_call,
                        forbidden_alleles,
                        args.min_parent_gq,
                        args.allow_missing_parent_as_absent,
                    )

                    mother_absent = absent_for_alleles(
                        mother_call,
                        forbidden_alleles,
                        args.min_parent_gq,
                        args.allow_missing_parent_as_absent,
                    )

                    if not father_absent or not mother_absent:
                        continue

                    writer.writerow(
                        {
                            "chrom": rec.chrom,
                            "pos": rec.pos,
                            "end": (
                                rec.stop
                                if rec.stop is not None
                                else info_to_str(rec.info.get("END"))
                            ),
                            "id": rec.id or ".",
                            "ref": rec.ref or ".",
                            "alt": ",".join(rec.alts or []),
                            "child_insertion_alt": alt_for_indices(
                                rec,
                                child_ins_alleles,
                            ),
                            "svtype": info_to_str(rec.info.get("SVTYPE")),
                            "svlen": info_to_str(rec.info.get("SVLEN")),
                            "filter": filter_to_str(rec),
                            "cohort_af": (
                                "."
                                if cohort_af is None
                                else f"{cohort_af:.6g}"
                            ),
                            "cohort_af_source": cohort_af_source,
                            "child": child,
                            "father": father,
                            "mother": mother,
                            "child_gt": gt_to_str(child_call),
                            "father_gt": gt_to_str(father_call),
                            "mother_gt": gt_to_str(mother_call),
                            "child_gq": gq_to_str(child_call),
                            "father_gq": gq_to_str(father_call),
                            "mother_gq": gq_to_str(mother_call),
                        }
                    )

                    n_candidates += 1

    sys.stderr.write(f"Complete trios in VCF: {len(valid_trios)}\n")
    sys.stderr.write(
        f"Trios skipped because sample(s) missing from VCF: {len(missing_trios)}\n"
    )
    sys.stderr.write(f"VCF records scanned: {n_records}\n")
    sys.stderr.write(f"Insertion records scanned: {n_ins_records}\n")
    sys.stderr.write(f"Candidate de novo rare insertion rows written: {n_candidates}\n")


if __name__ == "__main__":
    main()
