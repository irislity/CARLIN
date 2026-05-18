"""
convert_whitelist.py
--------------------
Convert the entire Visium HD whitelist (hd_6.5mm_bc_list.txt) to a
barcode → spatial_barcode mapping TSV.

Each line N (0-based) in the whitelist maps to:
    row = N // HD_6_5MM_NCOLS
    col = N % HD_6_5MM_NCOLS
    spatial_barcode = s_002um_{row:05}_{col:05}-1

Output TSV columns:
    barcode    spatial_barcode

Usage
-----
    python convert_whitelist.py \\
        --whitelist  ../spatial_barcode_match/data/hd_6.5mm_bc_list.txt \\
        --output     whitelist_spatial.tsv

    # custom gem group
    python convert_whitelist.py \\
        --whitelist  ../spatial_barcode_match/data/hd_6.5mm_bc_list.txt \\
        --output     whitelist_spatial.tsv \\
        --gem-group  2
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

from spatial_barcode import SquareBinIndex

HD_6_5MM_NCOLS: int = 3_350
HD_6_5MM_PITCH_UM: int = 2


def convert_whitelist(
    whitelist_path: str | Path,
    output_path: str | Path,
    ncols: int = HD_6_5MM_NCOLS,
    size_um: int = HD_6_5MM_PITCH_UM,
    gem_group: int = 1,
) -> int:
    """Convert every entry in the whitelist to a spatial barcode.

    Parameters
    ----------
    whitelist_path : str | Path
        Path to hd_6.5mm_bc_list.txt (plain text or .gz).
    output_path : str | Path
        Destination TSV path.
    ncols : int
        Number of columns in the slide grid (default 3350 for 6.5 mm).
    size_um : int
        Spot pitch in µm (default 2).
    gem_group : int
        Gem group suffix (default 1).

    Returns
    -------
    int
        Total number of entries written.
    """
    whitelist_path = Path(whitelist_path)
    output_path = Path(output_path)

    opener = gzip.open if whitelist_path.suffix == ".gz" else open

    print(f"[convert_whitelist] Reading: {whitelist_path}", file=sys.stderr)
    print(f"[convert_whitelist] Writing: {output_path}", file=sys.stderr)

    count = 0
    with (
        opener(whitelist_path, "rb") as in_fh,
        open(output_path, "w", newline="") as out_fh,
    ):
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["barcode", "spatial_barcode"])

        for i, raw_line in enumerate(in_fh):
            seq = raw_line.rstrip(b"\n\r")
            if not seq or seq.startswith(b"#"):
                continue

            row = i // ncols
            col = i % ncols
            spatial_bc = SquareBinIndex(row=row, col=col, size_um=size_um).with_gem_group(gem_group)

            writer.writerow([seq.decode(), spatial_bc])
            count += 1

            if count % 1_000_000 == 0:
                print(f"[convert_whitelist]   {count:>12,} entries written...", file=sys.stderr)

    print(f"[convert_whitelist] Done: {count:,} entries written to {output_path}", file=sys.stderr)
    return count


def main() -> None:
    p = argparse.ArgumentParser(
        description="Convert entire Visium HD whitelist to barcode → spatial_barcode TSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--whitelist", required=True, metavar="FILE",
                   help="hd_6.5mm_bc_list.txt (plain text or .gz)")
    p.add_argument("--output", required=True, metavar="TSV",
                   help="Output TSV path")
    p.add_argument("--ncols", type=int, default=HD_6_5MM_NCOLS,
                   help=f"Slide grid columns (default: {HD_6_5MM_NCOLS})")
    p.add_argument("--size-um", type=int, default=HD_6_5MM_PITCH_UM,
                   help=f"Spot pitch in µm (default: {HD_6_5MM_PITCH_UM})")
    p.add_argument("--gem-group", type=int, default=1,
                   help="Gem group suffix (default: 1)")

    args = p.parse_args()
    convert_whitelist(
        whitelist_path=args.whitelist,
        output_path=args.output,
        ncols=args.ncols,
        size_um=args.size_um,
        gem_group=args.gem_group,
    )


if __name__ == "__main__":
    main()
