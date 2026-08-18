#!/usr/bin/env python3

import argparse
import csv


def parse_count_depth(value):
    parts = value.split("/")
    if len(parts) != 2:
        raise ValueError(f"Unexpected count/depth field: {value}")
    return int(parts[0]), int(parts[1])


def read_sample_order(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise SystemExit("ERROR: sample_order.tsv contains no samples.")
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Parse PoPoolation2 snp-frequency-diff _rc output."
    )
    parser.add_argument("--rc", required=True)
    parser.add_argument("--sample-order", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    samples = read_sample_order(args.sample_order)
    n_pools = len(samples)

    base_fields = [
        "chr", "pos", "ref", "allele_count", "allele_states",
        "deletion_sum", "snp_type", "major_alleles", "minor_alleles"
    ]
    pool_fields = []
    for i in range(1, n_pools + 1):
        pool_fields.extend([f"AB_{i}", f"depth_{i}", f"major_count_{i}", f"minor_count_{i}"])

    summary_fields = [
        "total_depth", "mean_depth_per_pool", "min_depth", "max_depth",
        "pools_with_data", "total_major_count", "total_minor_count"
    ]

    written = 0

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=base_fields + pool_fields + summary_fields,
            delimiter="\t"
        )
        writer.writeheader()

        with open(args.rc) as handle:
            for line_number, line in enumerate(handle, start=1):
                line = line.strip()
                if not line or line.startswith("##"):
                    continue

                fields = line.split()
                expected = 9 + 2 * n_pools
                if len(fields) != expected:
                    raise SystemExit(
                        f"ERROR: Expected {expected} columns on line {line_number}, "
                        f"found {len(fields)}."
                    )

                row = dict(zip(base_fields, fields[:9]))
                major_fields = fields[9:9+n_pools]
                minor_fields = fields[9+n_pools:9+2*n_pools]

                total_depth = 0
                total_major = 0
                total_minor = 0
                depths = []

                for i in range(n_pools):
                    major_count, major_depth = parse_count_depth(major_fields[i])
                    minor_count, minor_depth = parse_count_depth(minor_fields[i])

                    if major_depth != minor_depth:
                        raise SystemExit(
                            f"ERROR: Major/minor depth mismatch at "
                            f"{fields[0]}:{fields[1]}, pool {i+1}."
                        )

                    depth = major_depth
                    af = minor_count / depth if depth > 0 else None

                    row[f"AB_{i+1}"] = f"{af:.8f}" if af is not None else "NA"
                    row[f"depth_{i+1}"] = depth
                    row[f"major_count_{i+1}"] = major_count
                    row[f"minor_count_{i+1}"] = minor_count

                    depths.append(depth)
                    total_depth += depth
                    total_major += major_count
                    total_minor += minor_count

                row.update({
                    "total_depth": total_depth,
                    "mean_depth_per_pool": f"{sum(depths)/n_pools:.6f}",
                    "min_depth": min(depths),
                    "max_depth": max(depths),
                    "pools_with_data": sum(d > 0 for d in depths),
                    "total_major_count": total_major,
                    "total_minor_count": total_minor,
                })

                writer.writerow(row)
                written += 1

    print(f"Pools: {n_pools}")
    print(f"SNP rows written: {written}")


if __name__ == "__main__":
    main()
