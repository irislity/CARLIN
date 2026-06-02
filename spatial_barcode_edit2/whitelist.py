"""
Load the combined Visium HD whitelist and extract separate bc1 / bc2 sets.

The combined whitelist (e.g. hd_6.5mm_bc_list.txt) contains full
bc1 + junction + bc2 sequences of variable length (29–31 bp):

    [ bc1: 14 bp ][ junction: 1–3 bp ][ bc2: 14 bp ]
    └──────────────────────────────────────────────────┘
               total 29, 30, or 31 bp

bc1 is always the FIRST 14 bp of each entry.
bc2 is always the LAST  14 bp of each entry.

Splitting this way yields:
  ~3,350 unique bc1 sequences  (one per grid row)
  ~3,350 unique bc2 sequences  (one per grid column)
  3,350 × 3,350 = 11,222,500 total spots on the 6.5 mm slide

Source:
  chemistry_defs.json  SPATIAL-HD-v1-3P:
    bc1 length=14, offset=11
    bc2 length=14, offset=25
  preflight.py:
    VISIUM_HD_SLIDE_STATS = SlideStats(nrows=3350, ncols=3350, spot_pitch=2)
"""

from __future__ import annotations

import gzip
from pathlib import Path

BC1_LENGTH: int = 14
BC2_LENGTH: int = 14


def load_combined_whitelist(path: str | Path) -> list[bytes]:
    """
    Load the combined whitelist file and return all sequences as a list.

    Skips blank lines and comment lines (starting with '#').
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    seqs: list[bytes] = []
    with opener(path, "rb") as fh:
        for line in fh:
            seq = line.strip()
            if seq and not seq.startswith(b"#"):
                seqs.append(seq)
    return seqs


def split_combined_whitelist(
    combined: list[bytes],
) -> tuple[set[bytes], set[bytes]]:
    """
    Split combined bc1+junction+bc2 sequences into separate bc1 and bc2 sets.

    bc1 = first 14 bp, bc2 = last 14 bp.
    Entries shorter than 28 bp are skipped (counted in skipped).
    """
    import sys

    bc1_set: set[bytes] = set()
    bc2_set: set[bytes] = set()
    skipped = 0

    for entry in combined:
        if len(entry) < BC1_LENGTH + BC2_LENGTH:
            skipped += 1
            continue
        bc1_set.add(entry[:BC1_LENGTH])
        bc2_set.add(entry[-BC2_LENGTH:])

    if skipped:
        print(
            f"[whitelist] WARNING: skipped {skipped} entries shorter than "
            f"{BC1_LENGTH + BC2_LENGTH} bp",
            file=sys.stderr,
        )

    return bc1_set, bc2_set
