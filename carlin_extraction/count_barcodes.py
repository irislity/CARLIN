"""
count_barcodes.py
-----------------
Count unique spatial barcodes, unique CARLIN sequences, and unique
spatial+CARLIN combinations from a spatial_carlin TSV.

Only rows with carlin_status == "ok" are counted for CARLIN and combination
statistics. All valid rows (spatial_barcode != NA) are counted for spatial.

Usage
-----
    python3 carlin_extraction/count_barcodes.py \
        --input carlin_extraction/results/676-1_spatial_carlin.tsv \
        --out   carlin_extraction/results/676-1_counts.tsv
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description="Count unique barcodes in spatial_carlin TSV")
    p.add_argument("--input", required=True, type=Path, help="spatial_carlin TSV")
    p.add_argument("--out",   required=True, type=Path, help="Output summary TSV")
    args = p.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)

    spatial_bcs   = set()
    carlin_bcs    = set()
    combinations  = set()

    spatial_counts  = Counter()              # reads per spatial barcode
    carlin_counts   = Counter()              # reads per CARLIN sequence
    combo_counts    = Counter()              # reads per (spatial, carlin) pair
    spot_carlins: dict[str, set] = {}        # unique CARLINs per spot

    total = ok = 0

    with open(args.input, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            total += 1
            spot = row["spatial_barcode"]
            carlin = row["carlin_seq"]
            status = row["carlin_status"]

            if spot != "NA":
                spatial_bcs.add(spot)
                spatial_counts[spot] += 1

            if status == "ok":
                ok += 1
                carlin_bcs.add(carlin)
                carlin_counts[carlin] += 1
                combinations.add((spot, carlin))
                combo_counts[(spot, carlin)] += 1
                spot_carlins.setdefault(spot, set()).add(carlin)

    # CARLINs per spot stats
    carlins_per_spot = {spot: len(cs) for spot, cs in spot_carlins.items()}
    cps_values = sorted(carlins_per_spot.values())

    # Summary stats
    print(f"  Total reads              : {total:,}")
    print(f"  Reads with ok CARLIN     : {ok:,} ({100*ok/total:.1f}%)")
    print(f"  Unique spatial barcodes  : {len(spatial_bcs):,}")
    print(f"  Unique CARLIN sequences  : {len(carlin_bcs):,}")
    print(f"  Unique spot+CARLIN combos: {len(combinations):,}")
    print()

    # CARLINs per spot summary
    n_spots = len(cps_values)
    print(f"  CARLINs per spot:")
    print(f"    Min    : {cps_values[0]}")
    print(f"    Median : {cps_values[n_spots//2]}")
    print(f"    Mean   : {sum(cps_values)/n_spots:.1f}")
    print(f"    Max    : {cps_values[-1]}")
    print()

    # Top 10 most common CARLIN sequences
    print("  Top 10 CARLIN sequences by read count:")
    for seq, cnt in carlin_counts.most_common(10):
        print(f"    {cnt:>8,}  {seq}")

    # Write summary TSV
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["total_reads",              total])
        w.writerow(["reads_ok_carlin",          ok])
        w.writerow(["unique_spatial_barcodes",  len(spatial_bcs)])
        w.writerow(["unique_carlin_sequences",  len(carlin_bcs)])
        w.writerow(["unique_spot_carlin_combos",len(combinations)])
        w.writerow(["mean_carlins_per_spot",    f"{sum(cps_values)/n_spots:.2f}"])
        w.writerow(["median_carlins_per_spot",  cps_values[n_spots//2]])
        w.writerow(["max_carlins_per_spot",     cps_values[-1]])

    # Write per-combination count table
    combo_out = args.out.with_name(args.out.stem + "_combos.tsv")
    with open(combo_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["spatial_barcode", "carlin_seq", "read_count"])
        for (spot, carlin), cnt in combo_counts.most_common():
            w.writerow([spot, carlin, cnt])

    # Write per-spot CARLIN count table
    spot_out = args.out.with_name(args.out.stem + "_per_spot.tsv")
    with open(spot_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["spatial_barcode", "unique_carlins", "total_reads"])
        for spot, n_carl in sorted(carlins_per_spot.items(), key=lambda x: -x[1]):
            w.writerow([spot, n_carl, spatial_counts[spot]])

    print(f"\n  Summary written    : {args.out}")
    print(f"  Combo table        : {combo_out}")
    print(f"  CARLINs per spot   : {spot_out}")


if __name__ == "__main__":
    main()
