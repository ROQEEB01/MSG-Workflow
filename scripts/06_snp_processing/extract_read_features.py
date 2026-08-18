#!/usr/bin/env python3

import argparse
import csv
import pysam


def load_snps(path):
    snps = []
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, start, end = line.rstrip().split("\t")[:3]
            snps.append((chrom, int(start), int(end)))
    return snps


def main():
    parser = argparse.ArgumentParser(
        description="Extract per-read features at SNP positions from one BAM file."
    )
    parser.add_argument("--bam", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    snps = load_snps(args.bed)
    bam = pysam.AlignmentFile(args.bam, "rb")

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "chr", "pos", "read_base", "base_qual", "map_qual",
            "read_pos", "read_length", "rel_read_pos",
            "is_reverse", "n_mismatches", "alignment_score"
        ])

        for chrom, start, end in snps:
            for pileup_col in bam.pileup(
                chrom,
                start,
                end,
                min_mapping_quality=0,
                min_base_quality=0,
                truncate=True,
                stepper="nofilter"
            ):
                if pileup_col.pos != start:
                    continue

                for pileup_read in pileup_col.pileups:
                    if pileup_read.is_del or pileup_read.is_refskip:
                        continue

                    aln = pileup_read.alignment
                    qpos = pileup_read.query_position
                    if qpos is None or aln.query_sequence is None:
                        continue

                    read_length = aln.query_length or 0
                    base_qual = (
                        aln.query_qualities[qpos]
                        if aln.query_qualities is not None
                        else 0
                    )

                    writer.writerow([
                        chrom,
                        end,
                        aln.query_sequence[qpos].upper(),
                        base_qual,
                        aln.mapping_quality,
                        qpos,
                        read_length,
                        round(qpos / read_length, 6) if read_length else 0,
                        int(aln.is_reverse),
                        aln.get_tag("NM") if aln.has_tag("NM") else -1,
                        aln.get_tag("AS") if aln.has_tag("AS") else -1
                    ])

    bam.close()


if __name__ == "__main__":
    main()
