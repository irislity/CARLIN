"""
barcode_correction.py
---------------------
Python port of lib/rust/cr_lib/src/stages/barcode_correction.rs
(JointBc1Bc2 branch, lines 144-287)

Implements joint extraction and correction of the two-part Visium HD barcode
(bc1 + bc2) from an R1 read sequence.

Chemistry: SPATIAL-HD-v1-3P
  - bc1: 14 bp, nominal start offset 11 (flexible: 8–12)
  - bc2: 14 bp, immediately following bc1
  - Both corrected jointly: up to 1 mismatch allowed per segment
  - Correction target: the provided whitelist sets for bc1 and bc2
"""

from __future__ import annotations

from typing import Iterable

# ---------------------------------------------------------------------------
# Hamming distance helpers
# ---------------------------------------------------------------------------

def hamming_distance(a: bytes | str, b: bytes | str) -> int:
    """Compute Hamming distance between two equal-length sequences."""
    if len(a) != len(b):
        return len(a)  # treat length mismatch as maximum distance
    if isinstance(a, str):
        a, b = a.encode(), b.encode()
    return sum(x != y for x, y in zip(a, b))


def closest_in_whitelist(
    seq: bytes,
    whitelist: set[bytes],
    max_dist: int = 1,
) -> tuple[bytes | None, int]:
    """Find the closest whitelist entry within max_dist Hamming distance.

    Returns (best_match, distance), or (None, max_dist+1) if no match found.
    """
    best: bytes | None = None
    best_dist = max_dist + 1
    for entry in whitelist:
        if len(entry) != len(seq):
            continue
        d = hamming_distance(seq, entry)
        if d < best_dist:
            best_dist = d
            best = entry
            if d == 0:
                break
    return best, best_dist


# ---------------------------------------------------------------------------
# Barcode extraction parameters (from chemistry_defs.json / mod.rs)
# ---------------------------------------------------------------------------

class JointBc1Bc2Params:
    """Extraction parameters for SPATIAL-HD-v1-3P.

    Attributes
    ----------
    bc1_length : int
        Expected length of bc1 segment (14 bp).
    bc2_length : int
        Expected length of bc2 segment (14 bp).
    min_offset : int
        Minimum start position of bc1 in the R1 read (inclusive). Default 8.
    max_offset : int
        Maximum start position of bc1 in the R1 read (inclusive). Default 12.
    search_padding : int
        Extra ±1 bp tried around nominal lengths during correction. Default 1.
    """

    def __init__(
        self,
        bc1_length: int = 14,
        bc2_length: int = 14,
        min_offset: int = 8,
        max_offset: int = 12,
        search_padding: int = 1,
    ):
        self.bc1_length = bc1_length
        self.bc2_length = bc2_length
        self.min_offset = min_offset
        self.max_offset = max_offset
        self.search_padding = search_padding


# Default parameters for SPATIAL-HD-v1-3P
HD_V1_3P_PARAMS = JointBc1Bc2Params(
    bc1_length=14,
    bc2_length=14,
    min_offset=8,
    max_offset=12,
)


# ---------------------------------------------------------------------------
# Core correction logic
# ---------------------------------------------------------------------------

def _try_correct(
    read: bytes,
    offset: int,
    length: int,
    whitelist: set[bytes],
    max_dist: int = 1,
) -> tuple[bytes | None, int]:
    """Extract a subsequence and attempt whitelist correction.

    Returns (corrected_seq, distance) or (None, large_distance).
    """
    end = offset + length
    if end > len(read) or offset < 0:
        return None, 9999
    seq = read[offset:end]
    return closest_in_whitelist(seq, whitelist, max_dist=max_dist)


def correct_joint_bc1_bc2(
    r1_seq: bytes | str,
    bc1_whitelist: set[bytes],
    bc2_whitelist: set[bytes],
    params: JointBc1Bc2Params = HD_V1_3P_PARAMS,
    max_dist: int = 1,
) -> tuple[bytes | None, bytes | None, int, int]:
    """Extract and correct bc1 + bc2 from an R1 read.

    Mirrors the JointBc1Bc2 branch of correct_barcode_in_read() in
    lib/rust/cr_lib/src/stages/barcode_correction.rs (lines 144–270).

    Parameters
    ----------
    r1_seq : bytes | str
        Full R1 read sequence.
    bc1_whitelist : set[bytes]
        Set of valid bc1 sequences (14-mers).
    bc2_whitelist : set[bytes]
        Set of valid bc2 sequences (14-mers).
    params : JointBc1Bc2Params
        Extraction geometry for this chemistry. Defaults to HD_V1_3P_PARAMS.
    max_dist : int
        Maximum Hamming distance allowed for correction (default 1).

    Returns
    -------
    bc1 : bytes | None
        Corrected bc1 sequence, or None if uncorrectable.
    bc2 : bytes | None
        Corrected bc2 sequence, or None if uncorrectable.
    bc1_offset : int
        Start position of bc1 in the read (–1 if not found).
    bc2_offset : int
        Start position of bc2 in the read (–1 if not found).
    """
    if isinstance(r1_seq, str):
        r1_seq = r1_seq.encode()

    p = params
    pad = p.search_padding
    read_len = len(r1_seq)

    # --- Try nominal positions first (offset 11 as in chemistry_defs.json) ---
    nominal_offset = 11  # default mid-point of [min_offset, max_offset] window
    bc1_nom, d1_nom = _try_correct(r1_seq, nominal_offset, p.bc1_length, bc1_whitelist, max_dist)
    bc2_nom_offset = nominal_offset + p.bc1_length
    bc2_nom, d2_nom = _try_correct(r1_seq, bc2_nom_offset, p.bc2_length, bc2_whitelist, max_dist)

    # Case 1: both valid at nominal positions
    if bc1_nom is not None and bc2_nom is not None and d1_nom == 0 and d2_nom == 0:
        return bc1_nom, bc2_nom, nominal_offset, bc2_nom_offset

    # General search across allowed offset range
    length_range_1 = range(
        max(1, p.bc1_length - pad),
        min(read_len, p.bc1_length + pad) + 1,
    )
    length_range_2 = range(
        max(1, p.bc2_length - pad),
        min(read_len, p.bc2_length + pad) + 1,
    )

    # -- Determine initial validity at default positions --
    bc1_default, d1_default = _try_correct(
        r1_seq, nominal_offset, p.bc1_length, bc1_whitelist, max_dist
    )
    bc2_default_offset = nominal_offset + p.bc1_length
    bc2_default, d2_default = _try_correct(
        r1_seq, bc2_default_offset, p.bc2_length, bc2_whitelist, max_dist
    )

    bc1_valid = bc1_default is not None
    bc2_valid = bc2_default is not None

    LARGE_DIST = 9999

    # Case 2: bc1 valid, bc2 invalid → search bc2 immediately after bc1
    if bc1_valid and not bc2_valid:
        bc1_end = nominal_offset + p.bc1_length
        best = (None, LARGE_DIST, -1, -1)
        for len2 in length_range_2:
            bc2_cand, d2 = _try_correct(r1_seq, bc1_end, len2, bc2_whitelist, max_dist)
            if bc2_cand is not None and d2 < best[1]:
                best = (bc2_cand, d2, nominal_offset, bc1_end)
        bc2_result, _, bc1_off, bc2_off = best
        return bc1_default, bc2_result, bc1_off, bc2_off

    # Case 3: bc2 valid, bc1 invalid → search bc1 ending right before bc2 start
    if bc2_valid and not bc1_valid:
        bc2_start = bc2_default_offset
        max_bc1_len = p.bc1_length + pad
        best = (None, LARGE_DIST, -1)
        for offset in range(p.min_offset, p.max_offset + 1):
            bc1_len = min(bc2_start - offset, max_bc1_len)
            if bc1_len <= 0:
                continue
            bc1_cand, d1 = _try_correct(r1_seq, offset, bc1_len, bc1_whitelist, max_dist)
            if bc1_cand is not None and d1 < best[1]:
                best = (bc1_cand, d1, offset)
        bc1_result, _, bc1_off = best
        return bc1_result, bc2_default, bc1_off, bc2_default_offset

    # Case 4: neither valid → exhaustive search over all offsets × lengths
    best_total = LARGE_DIST
    best_bc1 = best_bc2 = None
    best_bc1_off = best_bc2_off = -1

    for offset in range(p.min_offset, p.max_offset + 1):
        for len1 in length_range_1:
            bc2_start = offset + len1
            for len2 in length_range_2:
                if bc2_start + len2 > read_len:
                    continue
                bc1_cand, d1 = _try_correct(r1_seq, offset, len1, bc1_whitelist, max_dist)
                bc2_cand, d2 = _try_correct(r1_seq, bc2_start, len2, bc2_whitelist, max_dist)
                d1 = d1 if bc1_cand is not None else LARGE_DIST
                d2 = d2 if bc2_cand is not None else LARGE_DIST
                total = d1 + d2
                if total < best_total:
                    best_total = total
                    best_bc1 = bc1_cand
                    best_bc2 = bc2_cand
                    best_bc1_off = offset
                    best_bc2_off = bc2_start

    return best_bc1, best_bc2, best_bc1_off, best_bc2_off


# ---------------------------------------------------------------------------
# Whitelist loading helper
# ---------------------------------------------------------------------------

def load_whitelist(path: str) -> set[bytes]:
    """Load a barcode whitelist text file (one barcode per line, # = comment).

    Parameters
    ----------
    path : str
        Path to a plain-text or gzip-compressed whitelist (.txt or .txt.gz).
    """
    import gzip

    opener = gzip.open if path.endswith(".gz") else open
    barcodes: set[bytes] = set()
    with opener(path, "rb") as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith(b"#"):
                barcodes.add(line)
    return barcodes


# ---------------------------------------------------------------------------
# Quick test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # Synthetic test: exact match at nominal positions
    bc1_wl = {b"ACGTACGTACGTAC"}  # 14-mer
    bc2_wl = {b"TTGGTTGGTTGGTT"}  # 14-mer

    # Build R1 with bc1 at offset 11, bc2 immediately after
    prefix = b"N" * 11
    r1 = prefix + b"ACGTACGTACGTAC" + b"TTGGTTGGTTGGTT" + b"NNNNNNNN"

    bc1, bc2, off1, off2 = correct_joint_bc1_bc2(r1, bc1_wl, bc2_wl)
    assert bc1 == b"ACGTACGTACGTAC", bc1
    assert bc2 == b"TTGGTTGGTTGGTT", bc2
    print(f"bc1={bc1!r} at offset {off1}, bc2={bc2!r} at offset {off2}")

    # Test with 1 mismatch in bc2
    r1_mut = prefix + b"ACGTACGTACGTAC" + b"TTGGTTGGTTGGTA" + b"NNNNNNNN"  # last base A→T
    bc1, bc2, off1, off2 = correct_joint_bc1_bc2(r1_mut, bc1_wl, bc2_wl)
    assert bc1 == b"ACGTACGTACGTAC"
    assert bc2 == b"TTGGTTGGTTGGTT"  # corrected back
    print("1-mismatch correction passed.")
    print("All barcode_correction.py tests passed.")
