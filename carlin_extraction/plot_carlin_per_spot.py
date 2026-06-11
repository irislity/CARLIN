"""
plot_carlin_per_spot.py
-----------------------
Plot the distribution of unique CARLIN sequences per spot.

Usage
-----
    python3 carlin_extraction/plot_carlin_per_spot.py \
        --input carlin_extraction/results/676-1_counts_per_spot.tsv \
                carlin_extraction/results/676-4_counts_per_spot.tsv \
        --labels 676-1 676-4 \
        --out    carlin_extraction/results/carlin_per_spot_distribution.png
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


def load_per_spot(path: Path) -> list[int]:
    values = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            values.append(int(row["unique_carlins"]))
    return values


def main():
    p = argparse.ArgumentParser(description="Plot CARLIN per spot distribution")
    p.add_argument("--input",  nargs="+", required=True, type=Path)
    p.add_argument("--labels", nargs="+", required=True)
    p.add_argument("--out",    required=True, type=Path)
    args = p.parse_args()

    assert len(args.input) == len(args.labels), "Must supply one label per input file"

    n = len(args.input)
    fig, axes = plt.subplots(n, 2, figsize=(13, 5 * n))
    if n == 1:
        axes = [axes]   # ensure 2D indexing
    colors = ["#2196F3", "#E91E63", "#4CAF50", "#FF9800"]

    for row_axes, path, label, color in zip(axes, args.input, args.labels, colors):
        values = load_per_spot(path)
        arr = np.array(values)

        median = int(np.median(arr))
        mean   = arr.mean()
        mx     = int(arr.max())

        # Left — linear x, clipped at 50
        row_axes[0].hist(arr[arr <= 50], bins=50, color=color, alpha=0.8)
        row_axes[0].set_xlabel("Unique CARLINs per spot")
        row_axes[0].set_ylabel("Number of spots")
        row_axes[0].set_title(f"{label} — zoomed (1–50)\nmedian={median}  mean={mean:.1f}")
        row_axes[0].xaxis.set_major_locator(ticker.MultipleLocator(5))
        row_axes[0].grid(axis="y", alpha=0.3)

        # Right — log x, full range
        bins_log = np.logspace(0, np.log10(mx + 1), 60)
        row_axes[1].hist(arr, bins=bins_log, color=color, alpha=0.8)
        row_axes[1].set_xscale("log")
        row_axes[1].set_xlabel("Unique CARLINs per spot (log scale)")
        row_axes[1].set_ylabel("Number of spots")
        row_axes[1].set_title(f"{label} — full range (log scale)\nmax={mx}")
        row_axes[1].grid(axis="y", alpha=0.3)

    fig.suptitle("Unique CARLIN sequences per spatial spot", fontsize=13, fontweight="bold")
    plt.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
