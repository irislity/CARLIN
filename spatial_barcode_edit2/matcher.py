"""
JointBc1Bc2 barcode matching for SPATIAL-HD-v1-3P with edit-distance-2 correction.

Mirrors Space Ranger's barcode_correction.rs JointBc1Bc2 logic, replacing
the Posterior (1-mismatch Hamming) corrector with the precomputed
edit-distance-2 map corrector from corrector.py.

R1 read layout (SPATIAL-HD-v1-3P)
----------------------------------
  pos  0–8   : UMI            (9 bp)
  pos  9–10  : spacer         (2 bp, ignored)
  pos  11–24 : bc1            (14 bp, canonical offset = 11)
  pos  25–38 : bc2            (14 bp)
  pos  39+   : polyT tail     (ignored)

The junction between bc1 and bc2 can shift by ±1–2 bp due to ligation
variability, making the observed bc1+bc2 region 29–31 bp total.
This is handled by the offset scan (offsets 8–12) and the SEARCH_PADDING
on bc2 length during correction.

Four correction cases (mirrors barcode_correction.rs lines 163–270)
--------------------------------------------------------------------
  (bc1 valid, bc2 valid)   → pass through, no correction.
  (bc1 valid, bc2 invalid) → anchor bc1 end, search lengths for bc2.
  (bc1 invalid, bc2 valid) → anchor bc2 start, scan offsets for bc1.
  (bc1 invalid, bc2 invalid) → joint search over all (offset, len1, len2).

Source:
  lib/rust/cr_lib/src/stages/barcode_correction.rs  lines 144–270
  lib/rust/cr_types/src/rna_read.rs                 JointBc1Bc2 extraction
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

from .corrector import correct_segment, CorrectionResult
from .whitelist import BC1_LENGTH, BC2_LENGTH

# ---------------------------------------------------------------------------
# Constants — mirror barcode_correction.rs / rna_read.rs
# ---------------------------------------------------------------------------

MIN_OFFSET: int = 8
MAX_OFFSET: int = 12
BC1_DEFAULT_OFFSET: int = 11
UMI_OFFSET: int = 0
UMI_LENGTH: int = 9
SEARCH_PADDING: int = 1
LARGE_DISTANCE: int = 1000


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------

@dataclass
class MatchResult:
    """Full per-read result from match_barcode."""
    read_id: str
    umi: bytes
    bc1: bytes
    bc2: bytes
    barcode: Optional[bytes]          # bc1 + bc2 (28 bp) or None
    bc1_edit_distance: int            # 0=exact, 1–2=corrected, -1=uncorrectable
    bc2_edit_distance: int
    bc1_corrected: bool
    bc2_corrected: bool
    valid: bool


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _no_correction(bc: bytes) -> CorrectionResult:
    """Exact whitelist hit — edit distance 0."""
    return CorrectionResult(corrected=bc, edit_distance=0)


def _failed() -> Optional[CorrectionResult]:
    """Uncorrectable segment."""
    return None


# ---------------------------------------------------------------------------
# Core matching
# ---------------------------------------------------------------------------

def match_barcode(
    r1_seq: bytes,
    wl1: set[bytes],
    wl2: set[bytes],
    map1: dict,
    map2: dict,
    min_offset: int = MIN_OFFSET,
    max_offset: int = MAX_OFFSET,
    bc1_len: int = BC1_LENGTH,
    bc2_len: int = BC2_LENGTH,
) -> tuple[Optional[CorrectionResult], Optional[CorrectionResult], int]:
    """
    Find the best (bc1_result, bc2_result, offset) for one R1 read.

    Implements the four-case JointBc1Bc2 correction strategy using the
    precomputed edit-distance-2 maps instead of Bayesian 1-mismatch.

    Returns (bc1_result, bc2_result, best_offset).
    Either result may be None if uncorrectable.
    """
    read_len = len(r1_seq)

    best_bc1: Optional[CorrectionResult] = None
    best_bc2: Optional[CorrectionResult] = None
    best_total_dist = LARGE_DISTANCE * 2 + 1
    best_offset = BC1_DEFAULT_OFFSET

    for offset in range(min_offset, max_offset + 1):
        for l1 in range(bc1_len - SEARCH_PADDING, bc1_len + SEARCH_PADDING + 1):
            for l2 in range(bc2_len - SEARCH_PADDING, bc2_len + SEARCH_PADDING + 1):
                if l1 <= 0 or l2 <= 0:
                    continue
                if offset + l1 + l2 > read_len:
                    continue

                s1 = r1_seq[offset: offset + l1]
                s2 = r1_seq[offset + l1: offset + l1 + l2]

                r1 = correct_segment(s1, wl1, map1)
                r2 = correct_segment(s2, wl2, map2)

                d1 = r1.edit_distance if r1 is not None else LARGE_DISTANCE
                d2 = r2.edit_distance if r2 is not None else LARGE_DISTANCE
                total = d1 + d2

                if total < best_total_dist:
                    best_total_dist = total
                    best_bc1 = r1
                    best_bc2 = r2
                    best_offset = offset

    return best_bc1, best_bc2, best_offset


# ---------------------------------------------------------------------------
# FASTQ processing
# ---------------------------------------------------------------------------

def process_fastq(
    fastq_path: str | Path,
    wl1: set[bytes],
    wl2: set[bytes],
    map1: dict,
    map2: dict,
    min_offset: int = MIN_OFFSET,
    max_offset: int = MAX_OFFSET,
) -> Iterator[MatchResult]:
    """
    Parse a FASTQ file and yield a MatchResult for every R1 read.
    """
    from .fastq import read_fastq

    min_read_len = min_offset + BC1_LENGTH + BC2_LENGTH

    for rec in read_fastq(fastq_path):
        if len(rec.seq) < min_read_len:
            print(
                f"[spatial_barcode_edit2] WARNING: {rec.header!r} is only "
                f"{len(rec.seq)} bp (need ≥ {min_read_len}), skipping.",
                file=sys.stderr,
            )
            continue

        umi = rec.seq[UMI_OFFSET: UMI_OFFSET + UMI_LENGTH]

        bc1_res, bc2_res, _ = match_barcode(
            rec.seq, wl1, wl2, map1, map2,
            min_offset=min_offset,
            max_offset=max_offset,
        )

        bc1 = bc1_res.corrected if bc1_res is not None else b"N" * BC1_LENGTH
        bc2 = bc2_res.corrected if bc2_res is not None else b"N" * BC2_LENGTH
        bc1_ed = bc1_res.edit_distance if bc1_res is not None else -1
        bc2_ed = bc2_res.edit_distance if bc2_res is not None else -1
        bc1_ok = bc1_res is not None
        bc2_ok = bc2_res is not None
        valid = bc1_ok and bc2_ok

        yield MatchResult(
            read_id=rec.header,
            umi=umi,
            bc1=bc1,
            bc2=bc2,
            barcode=(bc1 + bc2) if valid else None,
            bc1_edit_distance=bc1_ed,
            bc2_edit_distance=bc2_ed,
            bc1_corrected=bc1_ed > 0 if bc1_ok else False,
            bc2_corrected=bc2_ed > 0 if bc2_ok else False,
            valid=valid,
        )


def write_tsv(
    results: Iterator[MatchResult],
    out_path: Optional[str | Path] = None,
) -> dict[str, int]:
    """
    Write MatchResult objects to a TSV file (or stdout).

    Output columns:
        read_id, umi, bc1, bc2, barcode,
        bc1_edit_distance, bc2_edit_distance,
        bc1_corrected, bc2_corrected, valid

    Returns summary counts dict.
    """
    fh = open(out_path, "w", newline="") if out_path else sys.stdout
    try:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "read_id", "umi", "bc1", "bc2", "barcode",
            "bc1_edit_distance", "bc2_edit_distance",
            "bc1_corrected", "bc2_corrected", "valid",
        ])
        counts = {"total": 0, "valid": 0, "corrected": 0}
        for read_id, m in ((m.read_id, m) for m in results):
            counts["total"] += 1
            if m.valid:
                counts["valid"] += 1
            if m.bc1_corrected or m.bc2_corrected:
                counts["corrected"] += 1
            writer.writerow([
                read_id,
                m.umi.decode(),
                m.bc1.decode(),
                m.bc2.decode(),
                m.barcode.decode() if m.barcode else "NA",
                m.bc1_edit_distance,
                m.bc2_edit_distance,
                int(m.bc1_corrected),
                int(m.bc2_corrected),
                int(m.valid),
            ])
    finally:
        if out_path:
            fh.close()
    return counts
