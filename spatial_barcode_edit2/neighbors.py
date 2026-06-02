"""
Enumerate all DNA sequences within edit distance 1 or 2 of a given sequence.

Edit operations allowed:
  - Substitution : change one base to another (A→C, A→G, A→T, etc.)
  - Deletion     : remove one base          (length n → n-1)
  - Insertion    : insert one base anywhere  (length n → n+1)

For a 14 bp sequence:
  Edit distance 1 : 14×3 (sub) + 14 (del) + 15×4 (ins) = 116 neighbors
  Edit distance 2 : ~6,700 unique neighbors (after dedup)

Source reference:
  Space Ranger make_correction_map.rs comment:
  "~6GB memory is needed to build all the sequences two edits away
   from a whitelist of size 5k with 19 bases each in memory."
"""

from __future__ import annotations

from typing import Iterator

BASES: bytes = b"ACGT"


def edit1_neighbors(seq: bytes) -> Iterator[tuple[bytes, int]]:
    """Yield (neighbor, 1) for every sequence within edit distance 1 of seq."""
    n = len(seq)
    a = bytearray(seq)

    # Substitutions
    for i in range(n):
        orig = a[i]
        for b in BASES:
            if b != orig:
                a[i] = b
                yield bytes(a), 1
        a[i] = orig

    # Deletions
    for i in range(n):
        yield seq[:i] + seq[i + 1:], 1

    # Insertions
    for i in range(n + 1):
        for b in BASES:
            yield seq[:i] + bytes([b]) + seq[i:], 1


def edit2_neighbors(seq: bytes) -> Iterator[tuple[bytes, int]]:
    """
    Yield (neighbor, distance) for every unique sequence within edit distance 2.

    Yields distance-1 neighbors first, then distance-2 neighbors.
    The original sequence is never yielded.
    """
    seen: set[bytes] = {seq}

    d1_neighbors = list(edit1_neighbors(seq))

    for n1, _ in d1_neighbors:
        if n1 not in seen:
            seen.add(n1)
            yield n1, 1

        for n2, _ in edit1_neighbors(n1):
            if n2 not in seen:
                seen.add(n2)
                yield n2, 2
