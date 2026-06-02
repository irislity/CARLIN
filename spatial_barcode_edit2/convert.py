"""
convert.py  (spatial_barcode_edit2 version)
--------------------------------------------
Convert edit-distance-corrected Visium HD barcodes to spatial barcodes.

Difference from spatial_barcode_conversion/convert.py
------------------------------------------------------
The original pipeline (spatial_barcode_match) outputs a single combined
barcode sequence of 29–31 bp, which can be looked up directly in the
whitelist.

The edit-distance pipeline (spatial_barcode_edit2) outputs bc1 and bc2
as SEPARATE 14 bp sequences (after correction).  The whitelist contains
bc1 + junction (1–3 bp) + bc2 = 29–31 bp entries, so a direct lookup
of the 28 bp bc1+bc2 string will always miss.

Solution
--------
Build a (bc1, bc2) → (row, col) index by splitting every whitelist entry:
    bc1 = entry[0:14]
    bc2 = entry[-14:]
    row = line_number // 3350
    col = line_number % 3350

Each (bc1, bc2) pair is unique across the 11,222,500-entry whitelist
(one spot per row/col combination), so this mapping is unambiguous.

Input TSV columns (from spatial_barcode_edit2/__main__.py)
----------------------------------------------------------
    read_id  umi  bc1  bc2  barcode
    bc1_edit_distance  bc2_edit_distance
    bc1_corrected  bc2_corrected  valid

Output
------
Same TSV + one extra column: spatial_barcode  (e.g. s_002um_00042_00137-1)

Usage
-----
    python -m spatial_barcode_edit2.convert \\
        --whitelist  data/hd_6.5mm_bc_list.txt \\
        --input      results/676-1_edit2.tsv \\
        --output     results/676-1_edit2_spatial.tsv

    # multiple files
    python -m spatial_barcode_edit2.convert \\
        --whitelist  data/hd_6.5mm_bc_list.txt \\
        --input      results/676-1_edit2.tsv results/676-4_edit2.tsv \\
        --outdir     results/
"""

from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).parent.parent))
from spatial_barcode_conversion.spatial_barcode import SquareBinIndex

HD_6_5MM_NCOLS: int = 3_350
HD_6_5MM_PITCH_UM: int = 2
BC1_LENGTH: int = 14
BC2_LENGTH: int = 14


def load_bc1bc2_index(
    path: str | Path,
    *,
    progress: bool = True,
) -> dict[tuple[bytes, bytes], tuple[int, int]]:
    """
    Build a (bc1, bc2) → (row, col) index from the combined whitelist.

    Each line N (0-based) maps to row = N // HD_6_5MM_NCOLS, col = N % HD_6_5MM_NCOLS.
    bc1 = first 14 bp, bc2 = last 14 bp of each whitelist entry.
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    index: dict[tuple[bytes, bytes], tuple[int, int]] = {}
    skipped = 0

    if progress:
        print(f"[convert] Loading whitelist index: {path}", file=sys.stderr)

    with opener(path, "rb") as fh:
        for i, line in enumerate(fh):
            entry = line.strip()
            if not entry or entry.startswith(b"#"):
                continue
            if len(entry) < BC1_LENGTH + BC2_LENGTH:
                skipped += 1
                continue
            bc1 = entry[:BC1_LENGTH]
            index[(bc1, entry[-BC2_LENGTH:])] = (i // HD_6_5MM_NCOLS, i % HD_6_5MM_NCOLS)
            if progress and (i + 1) % 1_000_000 == 0:
                print(f"[convert]   {i+1:>12,} entries...", file=sys.stderr)

    if progress:
        print(f"[convert] Index built: {len(index):,} entries", file=sys.stderr)
    if skipped:
        print(f"[convert] WARNING: skipped {skipped} short entries", file=sys.stderr)
    return index


def bc1bc2_to_spatial(
    bc1: str | bytes,
    bc2: str | bytes,
    index: dict[tuple[bytes, bytes], tuple[int, int]],
    size_um: int = HD_6_5MM_PITCH_UM,
    gem_group: int = 1,
) -> Optional[str]:
    """
    Convert a (bc1, bc2) pair to a spatial barcode string.

    Returns None if the pair is not in the whitelist index.
    """
    if isinstance(bc1, str):
        bc1 = bc1.encode()
    if isinstance(bc2, str):
        bc2 = bc2.encode()

    result = index.get((bc1, bc2))
    if result is None:
        return None

    row, col = result
    return SquareBinIndex(row=row, col=col, size_um=size_um).with_gem_group(gem_group)


def convert_tsv(
    input_path: str | Path,
    output_path: str | Path,
    index: dict[tuple[bytes, bytes], tuple[int, int]],
    size_um: int = HD_6_5MM_PITCH_UM,
    gem_group: int = 1,
) -> dict[str, int]:
    """
    Add a spatial_barcode column to one edit2 TSV file.

    Input columns: read_id umi bc1 bc2 barcode bc1_edit_distance bc2_edit_distance
                   bc1_corrected bc2_corrected valid
    Output:        same + spatial_barcode
    """
    input_path  = Path(input_path)
    output_path = Path(output_path)
    counts = {"total": 0, "valid": 0, "converted": 0, "not_in_whitelist": 0}

    with (
        open(input_path, newline="") as in_fh,
        open(output_path, "w", newline="") as out_fh,
    ):
        reader = csv.DictReader(in_fh, delimiter="\t")
        fieldnames = (reader.fieldnames or []) + ["spatial_barcode"]
        writer = csv.DictWriter(
            out_fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()

        for row in reader:
            counts["total"] += 1
            spatial_bc = "NA"

            if row.get("valid") == "1":
                counts["valid"] += 1
                spatial_bc_result = bc1bc2_to_spatial(
                    row.get("bc1", ""),
                    row.get("bc2", ""),
                    index,
                    size_um=size_um,
                    gem_group=gem_group,
                )
                if spatial_bc_result is not None:
                    spatial_bc = spatial_bc_result
                    counts["converted"] += 1
                else:
                    counts["not_in_whitelist"] += 1

            row["spatial_barcode"] = spatial_bc
            writer.writerow(row)

            if counts["total"] % 1_000_000 == 0:
                print(
                    f"[convert]   {counts['total']:>12,} rows processed "
                    f"({counts['converted']:,} converted)...",
                    file=sys.stderr,
                )

    pct = 100 * counts["converted"] / counts["valid"] if counts["valid"] else 0.0
    print(
        f"[convert] {input_path.name}: "
        f"{counts['converted']:,}/{counts['valid']:,} valid barcodes converted "
        f"({pct:.1f}%)  →  {output_path}",
        file=sys.stderr,
    )
    return counts


def _parse_args():
    import argparse

    p = argparse.ArgumentParser(
        description=(
            "Convert spatial_barcode_edit2 TSV output to spatial barcode format.\n"
            "Input TSV must have columns: read_id umi bc1 bc2 barcode ... valid"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--whitelist", required=True, metavar="FILE",
                   help="hd_6.5mm_bc_list.txt (plain text or .gz)")
    p.add_argument("--input", required=True, nargs="+", metavar="TSV",
                   help="One or more edit2 TSV files to convert")

    out_group = p.add_mutually_exclusive_group()
    out_group.add_argument("--output", default=None, metavar="TSV",
                           help="Output file (single input only)")
    out_group.add_argument("--outdir", default=None, metavar="DIR",
                           help="Output directory; files written as <stem>_spatial.tsv")

    p.add_argument("--size-um", type=int, default=HD_6_5MM_PITCH_UM,
                   help=f"Spot pitch in µm (default: {HD_6_5MM_PITCH_UM})")
    p.add_argument("--gem-group", type=int, default=1,
                   help="Gem group suffix (default: 1)")
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    inputs = [Path(p) for p in args.input]

    if args.output and len(inputs) > 1:
        print("Error: --output can only be used with a single input file.", file=sys.stderr)
        sys.exit(1)

    if args.outdir:
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        outputs = [outdir / (p.stem + "_spatial.tsv") for p in inputs]
    elif args.output:
        outputs = [Path(args.output)]
    else:
        outputs = [p.parent / (p.stem + "_spatial.tsv") for p in inputs]

    index = load_bc1bc2_index(args.whitelist)

    total_stats: dict[str, int] = {"total": 0, "valid": 0, "converted": 0, "not_in_whitelist": 0}
    for inp, out in zip(inputs, outputs):
        stats = convert_tsv(inp, out, index, size_um=args.size_um, gem_group=args.gem_group)
        for k in total_stats:
            total_stats[k] += stats[k]

    if len(inputs) > 1:
        pct = 100 * total_stats["converted"] / total_stats["valid"] if total_stats["valid"] else 0
        print(
            f"[convert] Total: {total_stats['converted']:,}/{total_stats['valid']:,} "
            f"converted ({pct:.1f}%) across {len(inputs)} files",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
