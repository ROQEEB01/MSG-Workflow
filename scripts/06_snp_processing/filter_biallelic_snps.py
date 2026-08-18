#!/usr/bin/env python3

import argparse
import csv


def main():
    parser = argparse.ArgumentParser(
        description="Retain biallelic SNPs and calculate global frequency metrics."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.input, newline="") as inp:
        reader = csv.DictReader(inp, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("ERROR: Input table has no header.")

        required = {
            "chr", "pos", "allele_count",
            "total_major_count", "total_minor_count"
        }
        missing = required - set(reader.fieldnames)
        if missing:
            raise SystemExit("ERROR: Missing columns: " + ", ".join(sorted(missing)))

        fields = list(reader.fieldnames) + [
            "global_minor_allele_frequency",
            "He_global"
        ]

        with open(args.output, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
            writer.writeheader()

            examined = retained = 0

            for row in reader:
                examined += 1
                if int(row["allele_count"]) != 2:
                    continue

                major = int(row["total_major_count"])
                minor = int(row["total_minor_count"])
                total = major + minor

                maf = minor / total if total > 0 else 0.0
                he = 2 * maf * (1 - maf)

                row["global_minor_allele_frequency"] = f"{maf:.8f}"
                row["He_global"] = f"{he:.8f}"

                writer.writerow(row)
                retained += 1

    print(f"Sites examined: {examined}")
    print(f"Biallelic SNPs retained: {retained}")


if __name__ == "__main__":
    main()
