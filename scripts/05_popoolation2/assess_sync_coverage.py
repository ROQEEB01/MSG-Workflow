#!/usr/bin/env python3

import argparse
import math


TIERS = [
    ("Tier 1", 5, 1.00),
    ("Tier 2", 5, 0.75),
    ("Tier 3", 5, 0.50),
    ("Tier 4", 3, 0.50),
    ("Tier 5", 10, 1.00),
]


def coverage(field):
    parts = field.split(":")
    if len(parts) < 4:
        raise ValueError(f"Unexpected sync field: {field}")
    return sum(int(x) for x in parts[:4])


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate coverage-retention tiers for a PoPoolation2 sync file."
    )
    parser.add_argument("--sync", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    counts = None
    n_pools = None
    total_sites = 0

    with open(args.sync) as handle:
        for line in handle:
            if not line.strip():
                continue

            fields = line.rstrip().split("\t")
            if n_pools is None:
                n_pools = len(fields) - 3
                if n_pools < 1:
                    raise SystemExit("ERROR: No population columns found.")
                counts = [0] * len(TIERS)

            total_sites += 1
            cov = [coverage(x) for x in fields[3:]]

            for i, (_, min_cov, fraction) in enumerate(TIERS):
                required = math.ceil(fraction * n_pools)
                if sum(c >= min_cov for c in cov) >= required:
                    counts[i] += 1

    if total_sites == 0:
        raise SystemExit("ERROR: Sync file is empty.")

    with open(args.output, "w") as out:
        out.write(
            "tier\tmin_coverage\trequired_pools\ttotal_pools\t"
            "sites_retained\tfraction_retained\n"
        )
        for (name, min_cov, fraction), retained in zip(TIERS, counts):
            required = math.ceil(fraction * n_pools)
            out.write(
                f"{name}\t{min_cov}\t{required}\t{n_pools}\t"
                f"{retained}\t{retained/total_sites:.6f}\n"
            )


if __name__ == "__main__":
    main()
