#!/usr/bin/env python3

import argparse
import csv
import os
import re
from collections import OrderedDict


def safe_filename(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


def main():
    parser = argparse.ArgumentParser(
        description="Prepare BAM lists for merging technical sequencing replicates."
    )
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    bam_list_dir = os.path.join(args.output_dir, "bam_lists")
    os.makedirs(bam_list_dir, exist_ok=True)

    samples = OrderedDict()

    required = {
        "sample_id", "pool_size", "ploidy", "site",
        "batch", "site_batch", "bam_path"
    }

    with open(args.metadata, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise SystemExit("ERROR: Metadata file has no header.")

        missing = required - set(reader.fieldnames)
        if missing:
            raise SystemExit(
                "ERROR: Missing metadata columns: " + ", ".join(sorted(missing))
            )

        for row in reader:
            sample_id = row["sample_id"].strip()
            if not sample_id:
                raise SystemExit("ERROR: Empty sample_id.")

            try:
                pool_size = int(row["pool_size"])
                ploidy = int(row["ploidy"])
            except ValueError:
                raise SystemExit(f"ERROR: Invalid pool_size/ploidy for {sample_id}")

            if pool_size < 1 or ploidy < 1:
                raise SystemExit(f"ERROR: Invalid pool_size/ploidy for {sample_id}")

            bam_path = os.path.abspath(os.path.expanduser(row["bam_path"].strip()))

            meta = {
                "pool_size": pool_size,
                "ploidy": ploidy,
                "site": row["site"].strip(),
                "batch": row["batch"].strip(),
                "site_batch": row["site_batch"].strip(),
            }

            if sample_id not in samples:
                samples[sample_id] = {**meta, "bams": []}
            else:
                for field, value in meta.items():
                    if samples[sample_id][field] != value:
                        raise SystemExit(
                            f"ERROR: Inconsistent {field} across technical replicates "
                            f"for {sample_id}"
                        )

            samples[sample_id]["bams"].append(bam_path)

    manifest_path = os.path.join(args.output_dir, "merge_manifest.tsv")
    used_names = set()

    with open(manifest_path, "w", newline="") as manifest:
        writer = csv.writer(manifest, delimiter="\t")
        writer.writerow([
            "sample_id", "pool_size", "ploidy", "chromosome_count",
            "site", "batch", "site_batch", "n_bams", "bam_list"
        ])

        for sample_id, info in samples.items():
            filename = safe_filename(sample_id)
            if filename in used_names:
                raise SystemExit(f"ERROR: Duplicate safe filename: {filename}")
            used_names.add(filename)

            bam_list_path = os.path.join(bam_list_dir, f"{filename}.bamlist")
            with open(bam_list_path, "w") as out:
                for bam in info["bams"]:
                    out.write(f"{bam}\n")

            writer.writerow([
                sample_id,
                info["pool_size"],
                info["ploidy"],
                info["pool_size"] * info["ploidy"],
                info["site"],
                info["batch"],
                info["site_batch"],
                len(info["bams"]),
                os.path.abspath(bam_list_path),
            ])

    print(f"Biological pools: {len(samples)}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
