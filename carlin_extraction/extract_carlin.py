"""
extract_carlin.py
-----------------
Extract CARLIN amplicon sequences from R2_valid FASTQ files.

For each read, finds Primer5 and Primer3 and extracts the sequence between
them. Outputs a TSV with the read name, extraction status, and CARLIN sequence.

OriginalCARLIN amplicon structure (between primers):
  [prefix: CGCCG][segment1][PAM][segment2][PAM]...[segment10][postfix: TGGGAGCT]
  Total unedited length: ~230 bp

Usage
-----
    python3 pipeline/extract_carlin.py \
        --r2  spatial_barcode_edit2/results/676-1_R2_valid.fastq.gz \
        --out spatial_barcode_edit2/results/676-1_carlin.tsv

    python3 pipeline/extract_carlin.py \
        --r2  spatial_barcode_edit2/results/676-4_R2_valid.fastq.gz \
        --out spatial_barcode_edit2/results/676-4_carlin.tsv
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

# OriginalCARLIN primer sequences
PRIMER5 = "GAGCTGTACAAGTAAGCGGC"
PRIMER3 = "CGACTGTGCCTTCTAGTTGC"

# Expected amplicon length range:
# - Unedited: ~275 bp
# - Min: prefix(5) + postfix(8) = 13 bp (extreme full deletion)
# - Large deletions commonly seen at ~30-70 bp
AMPLICON_MIN =  10
AMPLICON_MAX =  300


def find_primer(seq: str, primer: str, max_mismatch: int = 2) -> int:
    """
    Return the start position of primer in seq, allowing up to max_mismatch
    mismatches. Returns -1 if not found.
    Searches only the first 60 bp for Primer5, full read for Primer3.
    """
    plen = len(primer)
    for i in range(len(seq) - plen + 1):
        mismatches = sum(a != b for a, b in zip(seq[i:i+plen], primer))
        if mismatches <= max_mismatch:
            return i
    return -1


def extract_carlin(r2_path: Path, out_path: Path) -> None:
    total = found = 0
    no_p5 = no_p3 = bad_len = 0

    with gzip.open(r2_path, "rt") as fin, open(out_path, "w") as fout:
        fout.write("read_name\tstatus\tcarlin_seq\tcarlin_len\n")

        while True:
            header = fin.readline()
            if not header:
                break
            seq  = fin.readline().rstrip()
            fin.readline()   # +
            fin.readline()   # qual
            total += 1

            read_name = header[1:].split(" ")[0].rstrip()

            # Find Primer5 (search first 60 bp only)
            p5_pos = find_primer(seq[:60], PRIMER5)
            if p5_pos == -1:
                no_p5 += 1
                fout.write(f"{read_name}\tno_primer5\t\t\n")
                continue

            p5_end = p5_pos + len(PRIMER5)

            # Find Primer3 (search after Primer5)
            p3_pos = find_primer(seq[p5_end:], PRIMER3)
            if p3_pos == -1:
                no_p3 += 1
                fout.write(f"{read_name}\tno_primer3\t{seq[p5_end:]}\t{len(seq[p5_end:])}\n")
                continue

            carlin = seq[p5_end : p5_end + p3_pos]

            if not (AMPLICON_MIN <= len(carlin) <= AMPLICON_MAX):
                bad_len += 1
                fout.write(f"{read_name}\tbad_length\t{carlin}\t{len(carlin)}\n")
                continue

            found += 1
            fout.write(f"{read_name}\tok\t{carlin}\t{len(carlin)}\n")

    print(f"  Total reads     : {total:,}")
    print(f"  Extracted (ok)  : {found:,} ({100*found/total:.1f}%)")
    print(f"  No Primer5      : {no_p5:,}")
    print(f"  No Primer3      : {no_p3:,}")
    print(f"  Bad length      : {bad_len:,}")
    print(f"  Output          : {out_path}")


def main():
    p = argparse.ArgumentParser(description="Extract CARLIN amplicon from R2_valid FASTQ")
    p.add_argument("--r2",  required=True, type=Path, help="R2_valid FASTQ (.fastq.gz)")
    p.add_argument("--out", required=True, type=Path, help="Output TSV")
    p.add_argument("--max-mismatch", type=int, default=2,
                   help="Max mismatches allowed in primer matching (default: 2)")
    args = p.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    extract_carlin(args.r2, args.out)


if __name__ == "__main__":
    main()
