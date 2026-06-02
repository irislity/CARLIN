"""
Build a precomputed edit-distance-2 correction map for a barcode whitelist.

Mirrors Space Ranger's MAKE_CORRECTION_MAP stage (make_correction_map.rs).

Algorithm
---------
For every valid barcode in the whitelist:
  1. Enumerate all sequences within edit distance 1 and 2  (via neighbors.py).
  2. Skip neighbors that are themselves valid barcodes (no correction needed).
  3. Record: observed_sequence → (correct_barcode, edit_distance).
  4. On collision (two valid barcodes equidistant): mark as ambiguous (None).

The resulting map is queried at runtime in corrector.py.

Memory estimate (Python dict)
------------------------------
  bc1 or bc2 alone: ~3,350 whitelist seqs × ~6,700 neighbors each
  ≈ 22 M map entries × ~50 bytes/entry ≈ ~1 GB peak during build.
  Final map (after dedup) is significantly smaller.

Time estimate
-------------
  ~5–15 minutes in Python for one segment (bc1 or bc2).
  Use --map-cache in the CLI to avoid rebuilding on repeated runs.

Source:
  make_correction_map.rs — CorrectionMapBuilder::new / build
  stubs/mod.rs line 43 — SymSpellMatch::with_whitelist(&probe.0, EditMetric::Levenshtein)
"""

from __future__ import annotations

from typing import Optional

from .neighbors import edit2_neighbors


def build_correction_map(
    whitelist: list[bytes] | set[bytes],
) -> dict[bytes, Optional[tuple[bytes, int]]]:
    """
    Build and return a correction map for the given whitelist.

    Parameters
    ----------
    whitelist:
        Iterable of valid barcode sequences (bc1 or bc2 set).

    Returns
    -------
    dict mapping observed_seq → (correct_bc, edit_distance) or None.
        None means the observed sequence is equidistant from two or more
        valid barcodes (ambiguous — discard the read).
    """
    whitelist_set = set(whitelist)
    correction_map: dict[bytes, Optional[tuple[bytes, int]]] = {}

    for i, valid_bc in enumerate(whitelist):
        for neighbor, dist in edit2_neighbors(valid_bc):
            if neighbor in whitelist_set:
                continue  # already valid — no correction needed

            if neighbor in correction_map:
                existing_bc, existing_dist = correction_map[neighbor] or (None, dist)
                if correction_map[neighbor] is None:
                    pass  # already ambiguous
                elif existing_dist == dist:
                    correction_map[neighbor] = None   # ambiguous
                elif dist < existing_dist:
                    correction_map[neighbor] = (valid_bc, dist)
                # else existing is closer — keep it
            else:
                correction_map[neighbor] = (valid_bc, dist)

    return correction_map
