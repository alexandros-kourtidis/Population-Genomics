#!/usr/bin/env python3
"""
Convert a biallelic SNP VCF + population map into TreeMix input format.

TreeMix format (gzipped):
  Header row: pop1 pop2 pop3 ...
  Each subsequent row: ref_count,alt_count for each pop (space-delimited)

Usage:
    python3 03_vcf_to_treemix.py \
        --vcf filtered.vcf.gz \
        --popmap population_map.txt \
        --out treemix.frq.gz
"""

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Convert a filtered VCF to TreeMix format")
    p.add_argument("--vcf", required=True, help="Filtered VCF or VCF.gz file")
    p.add_argument("--popmap", required=True, help="Two-column sample to population map")
    p.add_argument("--out", required=True, help="Output gzipped TreeMix file")
    return p.parse_args()


def load_popmap(path):
    """Return a dictionary mapping sample -> population."""
    popmap = {}
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                continue
            popmap[parts[0]] = parts[1]
    return popmap


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def open_output(path):
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    return gzip.open(out_path, "wt")


def main():
    args = parse_args()
    popmap = load_popmap(args.popmap)

    if not popmap:
        print("[ERROR] Population map is empty", file=sys.stderr)
        sys.exit(1)

    populations = sorted(set(popmap.values()))
    print(f"[INFO] Populations found: {', '.join(populations)}", file=sys.stderr)

    n_written = 0
    n_skipped = 0

    with open_text(args.vcf) as vcf_fh, open_output(args.out) as out_fh:
        sample_idx = []

        # --- header ---
        for line in vcf_fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip().split("\t")
                samples = fields[9:]
                missing = [s for s in samples if s not in popmap]
                if missing:
                    print(f"[WARN] {len(missing)} samples in VCF are missing from the popmap; skipping them", file=sys.stderr)

                sample_idx = [(i, popmap[s]) for i, s in enumerate(samples) if s in popmap]
                out_fh.write(" ".join(populations) + "\n")
                break

        if not sample_idx:
            print("[ERROR] No matching samples found between VCF and popmap", file=sys.stderr)
            sys.exit(1)

        # --- variants ---
        for line in vcf_fh:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            ref = fields[3]
            alt = fields[4]

            # Skip multi-allelic or indel sites; the filtering step should already remove them.
            if "," in alt or len(ref) != 1 or len(alt) != 1:
                n_skipped += 1
                continue

            fmt_fields = fields[8].split(":")
            if "GT" not in fmt_fields:
                n_skipped += 1
                continue

            gt_idx = fmt_fields.index("GT")
            counts = defaultdict(lambda: [0, 0])

            for col_i, pop in sample_idx:
                geno = fields[9 + col_i].split(":")[gt_idx]
                alleles = geno.replace("|", "/").split("/")
                for allele in alleles:
                    if allele in (".", ""):
                        continue
                    allele = int(allele)
                    if allele == 0:
                        counts[pop][0] += 1
                    elif allele == 1:
                        counts[pop][1] += 1

            if any(sum(counts[p]) == 0 for p in populations):
                n_skipped += 1
                continue

            row = " ".join(f"{counts[p][0]},{counts[p][1]}" for p in populations)
            out_fh.write(row + "\n")
            n_written += 1

    print(f"[INFO] SNPs written : {n_written}", file=sys.stderr)
    print(f"[INFO] SNPs skipped : {n_skipped}", file=sys.stderr)

    if n_written == 0:
        print("[ERROR] No SNPs written — check VCF and popmap", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
