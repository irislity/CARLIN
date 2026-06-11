"""
top_carlin_sequences.py
-----------------------
Report the most abundant CARLIN sequences across samples by summing
read counts from the combos TSV produced by count_barcodes.py.

Usage
-----
    python3 carlin_extraction/top_carlin_sequences.py \
        --input carlin_extraction/results/676-1_counts_combos.tsv \
                carlin_extraction/results/676-4_counts_combos.tsv \
        --labels 676-1 676-4 \
        --top    10 \
        --out    carlin_extraction/results/top_carlin_sequences.tsv
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description="Report most abundant CARLIN sequences")
    p.add_argument("--input",  nargs="+", required=True, type=Path)
    p.add_argument("--labels", nargs="+", required=True)
    p.add_argument("--top",    type=int, default=10, help="Number of top sequences (default: 10)")
    p.add_argument("--out",    required=True, type=Path)
    args = p.parse_args()

    assert len(args.input) == len(args.labels)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    # Collect per-sample counts
    sample_counts: dict[str, Counter] = {}
    for path, label in zip(args.input, args.labels):
        counts: Counter = Counter()
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                counts[row["carlin_seq"]] += int(row["read_count"])
        sample_counts[label] = counts
        print(f"  {label}: {sum(counts.values()):,} total reads, {len(counts):,} unique CARLIN sequences")

    # Determine ranking by total reads across all samples
    combined: Counter = Counter()
    for counts in sample_counts.values():
        combined.update(counts)

    top_seqs = [seq for seq, _ in combined.most_common(args.top)]

    # Print to console
    for label, counts in sample_counts.items():
        total = sum(counts.values())
        print(f"\n=== {label} (top {args.top}) ===")
        print(f"{'Rank':<5} {'Reads':>12} {'% of total':>10}  CARLIN sequence")
        for i, seq in enumerate(top_seqs, 1):
            cnt = counts[seq]
            print(f"{i:<5} {cnt:>12,} {100*cnt/total:>9.1f}%  {seq}")

    # Write TSV
    fieldnames = ["rank", "carlin_seq", "carlin_len"] + \
                 [f"{l}_reads" for l in args.labels] + \
                 [f"{l}_pct" for l in args.labels]

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for i, seq in enumerate(top_seqs, 1):
            row = {"rank": i, "carlin_seq": seq, "carlin_len": len(seq)}
            for label, counts in sample_counts.items():
                total = sum(counts.values())
                cnt = counts[seq]
                row[f"{label}_reads"] = cnt
                row[f"{label}_pct"]   = f"{100*cnt/total:.2f}"
            w.writerow(row)

    print(f"\nSaved: {args.out}")


if __name__ == "__main__":
    main()
