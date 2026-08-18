#!/usr/bin/env python3

import argparse
import csv
import pysam


def main():
    parser = argparse.ArgumentParser(
        description="Extract reference sequence context and homopolymer status around SNPs."
    )
    parser.add_argument("--reference", required=True, help="Indexed reference FASTA.")
    parser.add_argument("--bed", required=True, help="BED file of SNP positions.")
    parser.add_argument("--output", required=True, help="Output TSV.")
    args = parser.parse_args()

    ref = pysam.FastaFile(args.reference)

    with open(args.bed) as inp, open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "chr", "pos", "ref_base", "context_11bp",
            "upstream_5", "downstream_5",
            "homopolymer_len", "in_homopolymer"
        ])

        for line in inp:
            if not line.strip():
                continue

            chrom, start, end = line.rstrip().split("\t")[:3]
            pos0 = int(start)
            pos1 = int(end)

            ref_len = ref.get_reference_length(chrom)
            context_start = max(0, pos0 - 5)
            context_end = min(ref_len, pos0 + 6)

            ref_base = ref.fetch(chrom, pos0, pos0 + 1).upper()
            context = ref.fetch(chrom, context_start, context_end).upper()
            upstream = ref.fetch(chrom, max(0, pos0 - 5), pos0).upper()
            downstream = ref.fetch(chrom, pos0 + 1, min(ref_len, pos0 + 6)).upper()

            hp_len = 1

            i = pos0 - 1
            while i >= 0 and ref.fetch(chrom, i, i + 1).upper() == ref_base:
                hp_len += 1
                i -= 1

            i = pos0 + 1
            while i < ref_len and ref.fetch(chrom, i, i + 1).upper() == ref_base:
                hp_len += 1
                i += 1

            writer.writerow([
                chrom, pos1, ref_base, context,
                upstream, downstream, hp_len, int(hp_len >= 3)
            ])

    ref.close()


if __name__ == "__main__":
    main()
