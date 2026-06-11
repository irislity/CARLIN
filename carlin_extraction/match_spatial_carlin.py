"""
match_spatial_carlin.py
-----------------------
Join spatial barcodes and CARLIN sequences by read ID.

Inputs
------
  --spatial   <sample>_edit2_spatial.tsv   (read_id → spatial_barcode + bc1+bc2 + umi)
  --carlin    <sample>_carlin.tsv          (read_id → carlin_seq, from extract_carlin.py)
  --out       <sample>_spatial_carlin.tsv  (joined output)

Output columns
--------------
  read_name       read name (before first space)
  umi             9 bp UMI
  spatial_barcode Visium HD spot ID (e.g. s_002um_01944_02002-1)
  barcode         corrected bc1+bc2 (28 bp)
  carlin_seq      extracted CARLIN amplicon sequence
  carlin_len      length of CARLIN amplicon
  carlin_status   ok / no_primer3 / no_primer5 / bad_length

Only reads that are valid (valid==1) in the spatial TSV are included.
Reads with no CARLIN match are kept with carlin_seq = NA.

Usage
-----
    python3 carlin_extraction/match_spatial_carlin.py \\
        --spatial spatial_barcode_edit2/results/676-1_edit2_spatial.tsv \\
        --carlin  carlin_extraction/results/676-1_carlin.tsv \\
        --out     carlin_extraction/results/676-1_spatial_carlin.tsv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_carlin(carlin_path: Path) -> dict[str, dict]:
    """Return {read_name: {carlin_seq, carlin_len, carlin_status}}."""
    carlin: dict[str, dict] = {}
    with open(carlin_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            carlin[row["read_name"]] = {
                "carlin_seq":    row["carlin_seq"],
                "carlin_len":    row["carlin_len"],
                "carlin_status": row["status"],
            }
    return carlin


def main():
    p = argparse.ArgumentParser(description="Join spatial barcodes with CARLIN sequences")
    p.add_argument("--spatial", required=True, type=Path, help="edit2_spatial TSV")
    p.add_argument("--carlin",  required=True, type=Path, help="carlin TSV (from extract_carlin.py)")
    p.add_argument("--out",     required=True, type=Path, help="Output TSV")
    args = p.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading CARLIN sequences from {args.carlin.name} ...")
    carlin = load_carlin(args.carlin)
    print(f"  {len(carlin):,} CARLIN records loaded")

    total = valid = matched = unmatched = 0

    with open(args.spatial, newline="") as fin, open(args.out, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, delimiter="\t", fieldnames=[
            "read_name", "umi", "spatial_barcode", "barcode",
            "carlin_seq", "carlin_len", "carlin_status"
        ])
        writer.writeheader()

        for row in reader:
            total += 1
            if row["valid"] != "1":
                continue
            valid += 1

            read_name = row["read_id"].split(" ")[0]
            c = carlin.get(read_name)

            if c and c["carlin_status"] == "ok":
                matched += 1
            else:
                unmatched += 1

            writer.writerow({
                "read_name":       read_name,
                "umi":             row["umi"],
                "spatial_barcode": row["spatial_barcode"],
                "barcode":         row["barcode"],
                "carlin_seq":      c["carlin_seq"] if c else "NA",
                "carlin_len":      c["carlin_len"] if c else "NA",
                "carlin_status":   c["carlin_status"] if c else "no_carlin",
            })

    print(f"  Total reads in spatial TSV : {total:,}")
    print(f"  Valid spatial barcodes     : {valid:,}")
    print(f"  Matched to CARLIN (ok)     : {matched:,} ({100*matched/valid:.1f}%)")
    print(f"  No/failed CARLIN           : {unmatched:,} ({100*unmatched/valid:.1f}%)")
    print(f"  Output: {args.out}")


if __name__ == "__main__":
    main()
