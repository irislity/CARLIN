"""
Command-line interface for spatial_barcode_match.

Usage — combined whitelist (hd_6.5mm_bc_list.txt from 10X website)
--------------------------------------------------------------------
python -m spatial_barcode_match \\
    --r1         reads_R1.fastq.gz \\
    --whitelist  hd_6.5mm_bc_list.txt \\
    --output     results.tsv \\
    --workers    8

Output TSV columns:
    read_id  umi  barcode  corrected  reverse_complemented  valid

"barcode" is the full matched sequence (29-31 bp) or "NA" when the read
could not be matched even after 1-mismatch correction.
"""

from __future__ import annotations

import argparse
import os
import sys

from .whitelist import Whitelist
from .matcher import process_fastq_combined, write_tsv_combined


def _parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="python -m spatial_barcode_match",
        description="Match SPATIAL-HD-v1-3P barcodes in R1 FASTQ reads.",
    )
    p.add_argument("--r1", required=True, metavar="FASTQ",
                   help="R1 FASTQ file (.fastq or .fastq.gz)")
    p.add_argument("--whitelist", required=True, metavar="FILE",
                   help="Combined barcode whitelist (hd_6.5mm_bc_list.txt)")
    p.add_argument("--output", default=None, metavar="TSV",
                   help="output TSV file (default: stdout)")
    p.add_argument("--workers", type=int, default=1, metavar="N",
                   help=("number of parallel worker processes "
                         "(default: 1; use -1 for all CPU cores)"))
    p.add_argument("--chunk-size", type=int, default=50_000, metavar="N",
                   help="reads per worker batch (default: 50000)")
    p.add_argument("--min-offset", type=int, default=8,
                   help="smallest barcode start offset to search (default: 8)")
    p.add_argument("--max-offset", type=int, default=12,
                   help="largest barcode start offset to search (default: 12)")
    p.add_argument("--confidence", type=float, default=0.975,
                   help="posterior threshold for 1-mismatch correction (default: 0.975)")
    p.add_argument("--reverse-complement", action="store_true", default=False,
                   help="also try the reverse complement of each barcode window")
    return p.parse_args(argv)


def main(argv=None) -> None:
    args = _parse_args(argv)

    workers = os.cpu_count() if args.workers == -1 else args.workers

    print(f"[spatial_barcode_match] Loading whitelist: {args.whitelist}",
          file=sys.stderr)
    wl = Whitelist.from_file(args.whitelist, progress=True)
    print(f"[spatial_barcode_match] {len(wl):,} barcodes  "
          f"(length {wl.min_seq_len}–{wl.max_seq_len} bp)", file=sys.stderr)
    print(f"[spatial_barcode_match] workers={workers}  "
          f"chunk={args.chunk_size:,}  confidence={args.confidence}  "
          f"revcomp={args.reverse_complement}",
          file=sys.stderr)

    results = process_fastq_combined(
        fastq_path=args.r1,
        whitelist=wl,
        min_offset=args.min_offset,
        max_offset=args.max_offset,
        confidence_threshold=args.confidence,
        reverse_complement=args.reverse_complement,
        num_workers=workers,
        chunk_size=args.chunk_size,
    )

    write_tsv_combined(results, out_path=args.output)

    if args.output:
        print(f"[spatial_barcode_match] Done → {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
