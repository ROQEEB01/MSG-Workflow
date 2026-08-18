#!/usr/bin/env python3

import argparse
import math


def coverage(field):
    parts = field.split(":")
    if len(parts) < 4:
        raise ValueError(f"Unexpected sync field: {field}")
    return sum(int(x) for x in parts[:4])


def main():
    parser = argparse.ArgumentParser(
        description="Filter PoPoolation2 sync positions by per-pool coverage."
    )
    parser.add_argument("--sync", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-coverage", type=int, required=True)
    parser.add_argument("--required-fraction", type=float, required=True)
    args = parser.parse_args()

    if not 0 < args.required_fraction <= 1:
        raise SystemExit("ERROR: required-fraction must be in (0, 1].")

    total = retained = 0
    n_pools = None

    with open(args.sync) as inp, open(args.output, "w") as out:
        for line in inp:
            if not line.strip():
                continue

            fields = line.rstrip().split("\t")

            if n_pools is None:
                n_pools = len(fields) - 3
                if n_pools < 1:
                    raise SystemExit("ERROR: No population columns found.")

            required = math.ceil(args.required_fraction * n_pools)
            cov = [coverage(x) for x in fields[3:]]

            total += 1
            if sum(c >= args.min_coverage for c in cov) >= required:
                out.write(line)
                retained += 1

    print(f"Pools: {n_pools}")
    print(f"Sites examined: {total}")
    print(f"Sites retained: {retained}")


if __name__ == "__main__":
    main()
