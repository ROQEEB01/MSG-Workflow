#!/usr/bin/env python3

import argparse
import csv
import glob
import os


def read_stats(stats_dir):
    records = {}
    for path in sorted(glob.glob(os.path.join(stats_dir, "*_mapping_summary.tsv"))):
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            row = next(reader, None)
            if row:
                records[row["sample_id"]] = row
    return records


def main():
    parser = argparse.ArgumentParser(
        description="Apply biological-pool quality filters and define final pool order."
    )
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--stats-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--min-mapped-reads", type=int, default=50000)
    parser.add_argument("--min-mapping-rate", type=float, default=50.0)
    parser.add_argument("--min-pool-size", type=int, default=5)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    stats = read_stats(args.stats_dir)

    qc_path = os.path.join(args.output_dir, "sample_qc.tsv")
    selected_path = os.path.join(args.output_dir, "selected_samples.tsv")
    bam_list_path = os.path.join(args.output_dir, "final_bams.txt")
    order_path = os.path.join(args.output_dir, "sample_order.tsv")

    fields = [
        "sample_id", "pool_size", "ploidy", "chromosome_count",
        "site", "batch", "site_batch",
        "total_reads", "mapped_reads", "unmapped_reads",
        "mapping_rate", "bam_path", "status", "exclusion_reason"
    ]

    with open(args.manifest, newline="") as manifest, \
         open(qc_path, "w", newline="") as qc_out, \
         open(selected_path, "w", newline="") as selected_out, \
         open(bam_list_path, "w") as bam_out, \
         open(order_path, "w", newline="") as order_out:

        reader = csv.DictReader(manifest, delimiter="\t")
        qc_writer = csv.DictWriter(qc_out, delimiter="\t", fieldnames=fields)
        selected_writer = csv.DictWriter(
            selected_out, delimiter="\t", fieldnames=fields
        )
        order_writer = csv.writer(order_out, delimiter="\t")

        qc_writer.writeheader()
        selected_writer.writeheader()
        order_writer.writerow([
            "pool_index", "sample_id", "pool_size", "ploidy",
            "chromosome_count", "site", "batch", "site_batch", "bam_path"
        ])

        pool_index = 0

        for manifest_row in reader:
            sample_id = manifest_row["sample_id"]
            if sample_id not in stats:
                raise SystemExit(f"ERROR: Missing mapping statistics for {sample_id}")

            stat = stats[sample_id]

            pool_size = int(stat["pool_size"])
            ploidy = int(stat["ploidy"])
            chromosome_count = int(stat["chromosome_count"])
            mapped = int(stat["mapped_reads"])
            unmapped = int(stat["unmapped_reads"])
            total = int(stat["total_reads"])
            mapping_rate = float(stat["mapping_rate"])

            if chromosome_count != pool_size * ploidy:
                raise SystemExit(
                    f"ERROR: chromosome_count mismatch for {sample_id}"
                )

            reasons = []
            if mapped < args.min_mapped_reads:
                reasons.append(f"mapped_reads<{args.min_mapped_reads}")
            if mapping_rate < args.min_mapping_rate:
                reasons.append(f"mapping_rate<{args.min_mapping_rate}")
            if pool_size < args.min_pool_size:
                reasons.append(f"pool_size<{args.min_pool_size}")

            status = "PASS" if not reasons else "FAIL"

            record = {
                "sample_id": sample_id,
                "pool_size": pool_size,
                "ploidy": ploidy,
                "chromosome_count": chromosome_count,
                "site": stat["site"],
                "batch": stat["batch"],
                "site_batch": stat["site_batch"],
                "total_reads": total,
                "mapped_reads": mapped,
                "unmapped_reads": unmapped,
                "mapping_rate": f"{mapping_rate:.4f}",
                "bam_path": stat["bam_path"],
                "status": status,
                "exclusion_reason": ";".join(reasons) if reasons else "-"
            }

            qc_writer.writerow(record)

            if status == "PASS":
                selected_writer.writerow(record)
                bam_out.write(f"{stat['bam_path']}\n")

                pool_index += 1
                order_writer.writerow([
                    pool_index,
                    sample_id,
                    pool_size,
                    ploidy,
                    chromosome_count,
                    stat["site"],
                    stat["batch"],
                    stat["site_batch"],
                    stat["bam_path"],
                ])

    print(f"QC report: {qc_path}")
    print(f"Selected samples: {selected_path}")
    print(f"Final BAM list: {bam_list_path}")
    print(f"Sample order: {order_path}")


if __name__ == "__main__":
    main()
