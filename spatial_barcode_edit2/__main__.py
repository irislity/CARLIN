"""
CLI for spatial_barcode_edit2.

Usage
-----
    python -m spatial_barcode_edit2 \
        --fastq      data/676-1_R1.fastq.gz \
        --whitelist  data/hd_6.5mm_bc_list.txt \
        --output     results/676-1_edit2.tsv \
        --map-cache  data/bc_correction_maps.pkl

Arguments
---------
--fastq       R1 FASTQ file (.fastq or .fastq.gz).
--whitelist   Combined Visium HD whitelist (hd_6.5mm_bc_list.txt, 29–31 bp entries).
--output      Output TSV path. Defaults to stdout if omitted.
--map-cache   Optional path to save/load precomputed correction maps (.pkl).
              Building the maps takes ~5–15 min and ~1 GB RAM each.
              Saving to cache avoids rebuilding on repeated runs.

Output columns
--------------
read_id, umi, bc1, bc2, barcode, bc1_edit_distance, bc2_edit_distance,
bc1_corrected, bc2_corrected, valid

  barcode = bc1 + bc2 (28 bp), or "NA" if correction failed.
  bc*_edit_distance: 0 = exact match, 1–2 = corrected, -1 = uncorrectable.
  valid: 1 if both bc1 and bc2 resolved, 0 otherwise.
"""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

from .whitelist import load_combined_whitelist, split_combined_whitelist
from .correction_map import build_correction_map
from .matcher import process_fastq, write_tsv


def build_or_load_maps(
    whitelist_path: Path,
    cache_path: Path | None,
) -> tuple[set, set, dict, dict]:
    """
    Return (wl1, wl2, map1, map2), loading from cache if available.

    If cache_path is given and the file exists, maps are loaded from it.
    Otherwise the maps are built from scratch and saved to cache_path.
    """
    if cache_path is not None and Path(cache_path).exists():
        print(f"[edit2] Loading correction maps from cache: {cache_path}", file=sys.stderr)
        with open(cache_path, "rb") as f:
            result = pickle.load(f)
        wl1, wl2, map1, map2 = result
        return wl1, wl2, map1, map2

    print(f"[edit2] Loading whitelist: {whitelist_path}", file=sys.stderr)
    combined = load_combined_whitelist(whitelist_path)
    wl1, wl2 = split_combined_whitelist(combined)

    print(f"[edit2] Building bc1 correction map ({len(wl1):,} sequences)…", file=sys.stderr)
    map1 = build_correction_map(wl1)
    print(f"[edit2] Building bc2 correction map ({len(wl2):,} sequences)…", file=sys.stderr)
    map2 = build_correction_map(wl2)

    if cache_path is not None:
        print(f"[edit2] Saving correction maps to: {cache_path}", file=sys.stderr)
        with open(cache_path, "wb") as f:
            pickle.dump((wl1, wl2, map1, map2), f)

    return wl1, wl2, map1, map2


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Edit-distance-2 Visium HD barcode matching.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--fastq", required=True, metavar="FILE",
                        help="R1 FASTQ file (.fastq or .fastq.gz)")
    parser.add_argument("--whitelist", required=True, metavar="FILE",
                        help="Combined Visium HD whitelist (hd_6.5mm_bc_list.txt)")
    parser.add_argument("--output", default=None, metavar="TSV",
                        help="Output TSV path (default: stdout)")
    parser.add_argument("--map-cache", default=None, metavar="PKL",
                        help="Path to save/load precomputed correction maps (.pkl)")

    args = parser.parse_args()

    wl1, wl2, map1, map2 = build_or_load_maps(
        Path(args.whitelist),
        Path(args.map_cache) if args.map_cache else None,
    )

    print(f"[edit2] Processing: {args.fastq}", file=sys.stderr)
    results = process_fastq(args.fastq, wl1, wl2, map1, map2)
    _tracked = write_tsv(results, args.output)

    total    = _tracked["total"]
    valid    = _tracked["valid"]
    corrected = _tracked["corrected"]
    pct_valid     = 100 * valid / total if total else 0.0
    pct_corrected = 100 * corrected / total if total else 0.0
    print(
        f"[edit2] Done. {total:,} reads | "
        f"{valid:,} valid ({pct_valid:.1f}%) | "
        f"{corrected:,} corrected ({pct_corrected:.1f}%)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
