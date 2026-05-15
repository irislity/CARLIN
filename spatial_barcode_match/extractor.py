"""
JointBc1Bc2 barcode extraction for SPATIAL-HD-v1-3P.

The Visium HD slide encodes each spot's position as two concatenated 14-bp
oligo sequences (bc1, bc2) that are ligated onto the R1 read.  Because the
ligation junction is not at a perfectly fixed position, Spaceranger scans a
small window of starting offsets and picks whichever position maximises exact
whitelist hits before falling back to the canonical positions.

Source references
-----------------
* Chemistry definition (offsets, lengths, extraction method, offset window):
    lib/python/cellranger/chemistry_defs.json  lines 1705-1754
* extract_barcode(), JointBc1Bc2 branch:
    lib/rust/cr_types/src/rna_read.rs  lines 316-368
* BarcodeExtraction::JointBc1Bc2 struct:
    lib/rust/cr_types/src/chemistry/mod.rs  lines 1238-1264
"""

from __future__ import annotations

from dataclasses import dataclass

from .whitelist import Whitelist

# ---------------------------------------------------------------------------
# SPATIAL-HD-v1-3P chemistry constants
# Source: lib/python/cellranger/chemistry_defs.json  lines 1705-1754
# ---------------------------------------------------------------------------

#: Canonical start offset of bc1 inside R1.
BC1_DEFAULT_OFFSET: int = 11
#: Fixed length of the bc1 segment.
BC1_LENGTH: int = 14
#: Canonical start offset of bc2 inside R1 (= BC1_DEFAULT_OFFSET + BC1_LENGTH).
BC2_DEFAULT_OFFSET: int = 25
#: Fixed length of the bc2 segment.
BC2_LENGTH: int = 14

#: UMI is the first 9 bases of R1.
UMI_OFFSET: int = 0
UMI_LENGTH: int = 9

# ---------------------------------------------------------------------------
# JointBc1Bc2 extraction parameters
# Source: chemistry_defs.json  "barcode_extraction": {"min_offset": 8, "max_offset": 12}
# ---------------------------------------------------------------------------

#: Smallest bc1 start offset tried during extraction.
MIN_OFFSET: int = 8
#: Largest bc1 start offset tried during extraction.
MAX_OFFSET: int = 12


@dataclass
class ExtractionResult:
    """Raw output of :func:`extract_barcode`."""

    bc1_seq: bytes
    bc2_seq: bytes
    umi_seq: bytes

    #: Actual offset chosen for bc1 (may differ from BC1_DEFAULT_OFFSET).
    bc1_offset: int
    bc1_len: int

    #: Actual start of bc2 (= bc1_offset + bc1_len).
    bc2_offset: int
    bc2_len: int

    #: Whether bc1 / bc2 were found verbatim in their whitelists.
    bc1_in_whitelist: bool
    bc2_in_whitelist: bool


def extract_barcode(
    r1: bytes,
    wl1: Whitelist,
    wl2: Whitelist,
    min_offset: int = MIN_OFFSET,
    max_offset: int = MAX_OFFSET,
    bc1_len: int = BC1_LENGTH,
    bc2_len: int = BC2_LENGTH,
) -> ExtractionResult:
    """
    Extract bc1, bc2, and UMI from a raw R1 read.

    Algorithm (JointBc1Bc2)
    -----------------------
    For each start offset in ``[min_offset, max_offset]``:

    * Extract ``bc1 = R1[offset : offset+bc1_len]``
    * Extract ``bc2 = R1[offset+bc1_len : offset+bc1_len+bc2_len]``
    * Count how many of the two segments are an exact whitelist hit.

    The canonical positions (``offset=11, len=14``) are appended **last** and
    win all ties.  This faithfully mirrors the Rust ``max_by_key`` which
    returns the **last** maximum in the iterator — the default (chained last)
    therefore wins whenever no better-scoring candidate exists.

    Source: lib/rust/cr_types/src/rna_read.rs  lines 336-367

    Parameters
    ----------
    r1:
        R1 nucleotide sequence as bytes (e.g. ``b"ACGT..."``).
    wl1, wl2:
        Whitelists for bc1 and bc2.
    min_offset, max_offset:
        Search window for bc1 start position.
    bc1_len, bc2_len:
        Expected segment lengths (both 14 for ``SPATIAL-HD-v1-3P``).
    """
    read_len = len(r1)
    default = (BC1_DEFAULT_OFFSET, bc1_len, bc2_len)

    # Build candidate list: all offsets in window, then default appended last.
    # Rust: (*min_offset..=*max_offset).flat_map(...).chain(once(default))
    candidates: list[tuple[int, int, int]] = []
    for offset in range(min_offset, max_offset + 1):
        if offset + bc1_len + bc2_len <= read_len:
            candidates.append((offset, bc1_len, bc2_len))
    candidates.append(default)

    # Pick the candidate with the highest whitelist-hit count.
    # Ties go to the **later** candidate (>= not >) so that default wins
    # zero-hit ties, matching the Rust last-wins max_by_key behaviour.
    # Source: rna_read.rs  lines 358-366
    best_score = -1
    best = default
    for off, l1, l2 in candidates:
        end1 = off + l1
        end2 = end1 + l2
        if end2 > read_len:
            continue
        score = wl1.contains(r1[off:end1]) + wl2.contains(r1[end1:end2])
        if score >= best_score:
            best_score = score
            best = (off, l1, l2)

    off, l1, l2 = best
    end1 = off + l1
    end2 = end1 + l2
    b1 = r1[off:end1]
    b2 = r1[end1:end2]
    umi = r1[UMI_OFFSET: UMI_OFFSET + UMI_LENGTH]

    return ExtractionResult(
        bc1_seq=b1,
        bc2_seq=b2,
        umi_seq=umi,
        bc1_offset=off,
        bc1_len=l1,
        bc2_offset=end1,
        bc2_len=l2,
        bc1_in_whitelist=wl1.contains(b1),
        bc2_in_whitelist=wl2.contains(b2),
    )
