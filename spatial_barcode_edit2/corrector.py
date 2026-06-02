"""
Runtime barcode correction using a precomputed edit-distance-2 map.

Replaces Space Ranger's internal SymSpellMatch corrector for Visium HD.

Lookup priority
---------------
1. Exact whitelist match (edit distance 0) → return as-is, no correction flag.
2. Precomputed map lookup (edit distance 1 or 2) → corrected if unambiguous.
3. No match → return None (read will not get a valid barcode).

Ambiguity handling
------------------
If the observed sequence is equidistant from two or more valid barcodes,
the correction_map stores None for that entry and the read is discarded.
This mirrors the Bayesian posterior approach in the 1-mismatch corrector
where low-confidence corrections are also discarded.

Source:
  barcode/src/corrector.rs — Posterior::correct_barcode (logic adapted)
  stubs/mod.rs line 43     — SymSpellMatch::with_whitelist(..., EditMetric::Levenshtein)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass
class CorrectionResult:
    corrected: bytes   # corrected barcode sequence
    edit_distance: int  # 0 = exact, 1–2 = corrected


def correct_segment(
    observed: bytes,
    whitelist: set[bytes],
    correction_map: dict,
) -> Optional[CorrectionResult]:
    """
    Attempt to correct *observed* onto the whitelist.

    Returns a CorrectionResult on success, or None if uncorrectable.
    """
    # Exact whitelist match — no correction needed
    if observed in whitelist:
        return CorrectionResult(corrected=observed, edit_distance=0)

    # Precomputed map lookup
    entry = correction_map.get(observed)
    if entry is None:
        return None  # not in map or ambiguous marker missing

    correct_bc, dist = entry
    if correct_bc is None:
        return None  # ambiguous

    return CorrectionResult(corrected=correct_bc, edit_distance=dist)
