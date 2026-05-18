"""
convert.py
----------
Convert already-matched Visium HD sequence barcodes to spatial barcodes.

Input
-----
TSV files produced by the spatial_barcode_match pipeline:
    read_id  umi  barcode  corrected  valid

  - ``barcode``   : bc1+bc2 combined sequence (29-31 bp) when valid=1, else "NA"
  - ``corrected`` : 0/1 flag (was a 1-mismatch correction applied?)
  - ``valid``     : 0/1

Whitelist
---------
``hd_6.5mm_bc_list.txt`` from 10X Genomics — one bc1+bc2 sequence per line.
The file has 11,222,500 entries arranged in row-major order on a 3350 × 3350
grid at 2 µm pitch.  Line N (0-based) → row = N // 3350, col = N % 3350.

Output
------
Same TSV with one extra column appended:
    spatial_barcode  e.g. s_002um_00042_00137-1

Usage
-----
As a library:
    from convert import load_whitelist_index, seq_to_spatial, convert_tsv

As a CLI:
    # single file
    python convert.py \\
        --whitelist  hd_6.5mm_bc_list.txt \\
        --input      676-1_barcodes.tsv \\
        --output     676-1_spatial.tsv

    # multiple files at once
    python convert.py \\
        --whitelist  hd_6.5mm_bc_list.txt \\
        --input      676-1_barcodes.tsv 676-4_barcodes.tsv \\
        --outdir     ./results
"""

from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path

from spatial_barcode import SquareBinIndex


# ---------------------------------------------------------------------------
# Grid constants for the Visium HD 6.5 mm slide
# ---------------------------------------------------------------------------

#: Total number of spots on a 6.5 mm HD slide (3350 × 3350).
HD_6_5MM_TOTAL_SPOTS: int = 11_222_500

#: Number of columns (spots per row) on the 6.5 mm HD slide.
HD_6_5MM_NCOLS: int = 3_350

#: Base spot pitch in micrometers.
HD_6_5MM_PITCH_UM: int = 2


# ---------------------------------------------------------------------------
# Whitelist index
# ---------------------------------------------------------------------------

def load_whitelist_index(
    path: str | Path,
    *,
    progress: bool = True,
) -> dict[bytes, int]:
    """Load the combined whitelist into a sequence → 0-based-index dict.

    The 0-based line number IS the spatial index:
        row = index // HD_6_5MM_NCOLS
        col = index % HD_6_5MM_NCOLS

    Parameters
    ----------
    path : str | Path
        Path to ``hd_6.5mm_bc_list.txt`` (plain text or .gz).
    progress : bool
        Print a progress message every 1 M entries (default True).

    Returns
    -------
    dict[bytes, int]
        Maps each barcode sequence (bytes) to its 0-based line number.
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    index: dict[bytes, int] = {}

    if progress:
        print(f"[convert] Loading whitelist: {path}", file=sys.stderr)

    with opener(path, "rb") as fh:
        for i, line in enumerate(fh):
            seq = line.rstrip(b"\n\r")
            if seq and not seq.startswith(b"#"):
                index[seq] = i
            if progress and (i + 1) % 1_000_000 == 0:
                print(f"[convert]   {i+1:>12,} entries loaded...", file=sys.stderr)

    if progress:
        print(f"[convert] Whitelist loaded: {len(index):,} entries", file=sys.stderr)
    return index


# ---------------------------------------------------------------------------
# Single-sequence conversion
# ---------------------------------------------------------------------------

def seq_to_spatial(
    seq: str | bytes,
    whitelist_index: dict[bytes, int],
    ncols: int = HD_6_5MM_NCOLS,
    size_um: int = HD_6_5MM_PITCH_UM,
    gem_group: int = 1,
) -> str | None:
    """Convert one bc1+bc2 sequence to a spatial barcode string.

    Parameters
    ----------
    seq : str | bytes
        Matched barcode sequence from the ``barcode`` column.
    whitelist_index : dict[bytes, int]
        Output of :func:`load_whitelist_index`.
    ncols : int
        Number of columns in the slide grid (default 3350 for 6.5 mm).
    size_um : int
        Spot pitch in µm (default 2).
    gem_group : int
        Gem group suffix (default 1).

    Returns
    -------
    str | None
        e.g. ``"s_002um_00042_00137-1"``, or ``None`` if seq not in whitelist.
    """
    if isinstance(seq, str):
        seq = seq.encode()

    idx = whitelist_index.get(seq)
    if idx is None:
        return None

    row = idx // ncols
    col = idx % ncols
    return SquareBinIndex(row=row, col=col, size_um=size_um).with_gem_group(gem_group)


# ---------------------------------------------------------------------------
# TSV conversion
# ---------------------------------------------------------------------------

def convert_tsv(
    input_path: str | Path,
    output_path: str | Path,
    whitelist_index: dict[bytes, int],
    ncols: int = HD_6_5MM_NCOLS,
    size_um: int = HD_6_5MM_PITCH_UM,
    gem_group: int = 1,
) -> dict[str, int]:
    """Add a ``spatial_barcode`` column to one barcodes TSV.

    Reads columns: read_id, umi, barcode, corrected, valid
    Writes same columns + spatial_barcode (or "NA" when valid=0 or not found).

    Parameters
    ----------
    input_path : str | Path
        Path to a ``*_barcodes.tsv`` file from spatial_barcode_match.
    output_path : str | Path
        Destination TSV path.
    whitelist_index : dict[bytes, int]
        Output of :func:`load_whitelist_index`.
    ncols, size_um, gem_group :
        Grid / format parameters (see :func:`seq_to_spatial`).

    Returns
    -------
    dict[str, int]
        Summary counts: total, valid, converted, not_in_whitelist.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    counts = {"total": 0, "valid": 0, "converted": 0, "not_in_whitelist": 0}

    with (
        open(input_path, newline="") as in_fh,
        open(output_path, "w", newline="") as out_fh,
    ):
        reader = csv.DictReader(in_fh, delimiter="\t")
        fieldnames = (reader.fieldnames or []) + ["spatial_barcode"]
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()

        for row in reader:
            counts["total"] += 1
            spatial_bc = "NA"

            if row.get("valid") == "1":
                counts["valid"] += 1
                barcode_seq = row.get("barcode", "NA")
                if barcode_seq != "NA":
                    result = seq_to_spatial(
                        barcode_seq, whitelist_index,
                        ncols=ncols, size_um=size_um, gem_group=gem_group,
                    )
                    if result is not None:
                        spatial_bc = result
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


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
    import argparse

    p = argparse.ArgumentParser(
        description=(
            "Convert matched Visium HD barcodes (TSV) to spatial barcode format.\n"
            "Input TSV must have columns: read_id  umi  barcode  corrected  valid"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--whitelist", required=True, metavar="FILE",
        help="hd_6.5mm_bc_list.txt (plain text or .gz)",
    )
    p.add_argument(
        "--input", required=True, nargs="+", metavar="TSV",
        help="One or more *_barcodes.tsv files to convert",
    )

    out_group = p.add_mutually_exclusive_group()
    out_group.add_argument(
        "--output", default=None, metavar="TSV",
        help="Output file (only valid when --input has a single file)",
    )
    out_group.add_argument(
        "--outdir", default=None, metavar="DIR",
        help=(
            "Output directory; each input file is written as "
            "<stem>_spatial.tsv inside this directory"
        ),
    )

    p.add_argument("--ncols", type=int, default=HD_6_5MM_NCOLS,
                   help=f"Slide grid columns (default: {HD_6_5MM_NCOLS})")
    p.add_argument("--size-um", type=int, default=HD_6_5MM_PITCH_UM,
                   help=f"Spot pitch in µm (default: {HD_6_5MM_PITCH_UM})")
    p.add_argument("--gem-group", type=int, default=1,
                   help="Gem group suffix (default: 1)")
    return p.parse_args()


def main() -> None:
    args = _parse_args()

    # Resolve output paths
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
        # Default: write next to input with _spatial suffix
        outputs = [p.parent / (p.stem + "_spatial.tsv") for p in inputs]

    # Load whitelist once
    wl_index = load_whitelist_index(args.whitelist)

    # Convert each file
    total_stats: dict[str, int] = {"total": 0, "valid": 0, "converted": 0, "not_in_whitelist": 0}
    for inp, out in zip(inputs, outputs):
        stats = convert_tsv(
            inp, out, wl_index,
            ncols=args.ncols,
            size_um=args.size_um,
            gem_group=args.gem_group,
        )
        for k in total_stats:
            total_stats[k] += stats[k]

    if len(inputs) > 1:
        pct = 100 * total_stats["converted"] / total_stats["valid"] if total_stats["valid"] else 0
        print(
            f"[convert] Total: {total_stats['converted']:,}/{total_stats['valid']:,} "
            f"valid barcodes converted ({pct:.1f}%) across {len(inputs)} files",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
