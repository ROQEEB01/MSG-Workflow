#!/usr/bin/env python3

import argparse
import csv


def main():
    parser = argparse.ArgumentParser(
        description="Create a BED file from the biallelic SNP table for read-level SNP auditing."
    )
    parser.add_argument("--snps", required=True, help="TSV containing chr and pos columns.")
    parser.add_argument("--output", required=True, help="Output BED file.")
    args = parser.parse_args()

    with open(args.snps, newline="") as inp, open(args.output, "w") as out:
        reader = csv.DictReader(inp, delimiter="\t")
        if reader.fieldnames is None or not {"chr", "pos"}.issubset(reader.fieldnames):
            raise SystemExit("ERROR: SNP table must contain chr and pos columns.")

        n = 0
        for row in reader:
            chrom = row["chr"].strip()
            pos = int(row["pos"])
            if pos < 1:
                raise SystemExit(f"ERROR: Invalid 1-based position: {chrom}:{pos}")
            out.write(f"{chrom}\t{pos - 1}\t{pos}\n")
            n += 1

    print(f"SNP positions written: {n}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
