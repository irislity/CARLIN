"""
Whitelist loading and membership checking.

Supports two usage modes:

1. **Split whitelist** (separate bc1 / bc2 files, 14 bp each)
   — used when the slide design is available as two plain-text files.

2. **Combined whitelist** (``hd_6.5mm_bc_list.txt`` from the 10X website)
   — all valid bc1+bc2 combinations pre-concatenated, 29-31 bp each.
   This is what you get from https://www.10xgenomics.com when downloading
   the Visium HD barcode list.  The variable length arises because the
   bc1 and bc2 oligos on the slide are not all exactly 14 bp.

Source references
-----------------
* Whitelist enum (Plain / Trans / SpatialHd):
    lib/rust/barcode/src/whitelist.rs  lines 508-565
* WhitelistSource::SlideFile iteration (oligos → BcSegSeq):
    lib/rust/barcode/src/whitelist.rs  lines 312-316
* load_oligos() stub (proprietary binary reader, not available publicly):
    lib/rust/slide_design/src/stubs/mod.rs  line 22
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path


class Whitelist:
    """
    An in-memory set of valid barcode sequences.

    Works for both the **split** case (one file per segment, fixed 14 bp)
    and the **combined** case (``hd_6.5mm_bc_list.txt``, variable 29-31 bp).

    Corresponds to ``Whitelist::Plain`` in the Rust source.

    Source: lib/rust/barcode/src/whitelist.rs  lines 512-565
    """

    def __init__(self, sequences: set[bytes]) -> None:
        self._seqs: set[bytes] = sequences
        if sequences:
            lengths = {len(s) for s in sequences}
            self._min_len: int = min(lengths)
            self._max_len: int = max(lengths)
        else:
            self._min_len = 0
            self._max_len = 0

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_file(cls, path: str | Path, *, progress: bool = False) -> "Whitelist":
        """
        Load from a plain-text or gzip-compressed file (one barcode per line).

        For large files like ``hd_6.5mm_bc_list.txt`` (~11 M lines) loading
        takes ~30-60 s and uses ~1-2 GB of RAM.  Pass ``progress=True`` to
        print a line-count update every million entries.

        Mirrors the text-file iteration path in WhitelistSource::TxtFile:
            lib/rust/barcode/src/whitelist.rs  lines 282-295
        """
        path = Path(path)
        opener = gzip.open if path.suffix == ".gz" else open
        seqs: set[bytes] = set()
        n = 0
        with opener(path, "rb") as fh:
            for raw in fh:
                bc = raw.strip()
                if bc and not bc.startswith(b"#"):
                    seqs.add(bc)
                    n += 1
                    if progress and n % 1_000_000 == 0:
                        print(f"  loaded {n:,} barcodes…", file=sys.stderr)
        if progress:
            print(f"  done — {n:,} barcodes loaded.", file=sys.stderr)
        return cls(seqs)

    @classmethod
    def from_sequences(cls, sequences) -> "Whitelist":
        """Build directly from an iterable of byte strings."""
        return cls({bytes(s) for s in sequences})

    # ------------------------------------------------------------------
    # Membership
    # ------------------------------------------------------------------

    def contains(self, seq: bytes) -> bool:
        """
        Exact-match lookup.

        Mirrors: Whitelist::contains  lib/rust/barcode/src/whitelist.rs  line 571
        """
        return seq in self._seqs

    # ------------------------------------------------------------------
    # Sequence-length helpers (used by the combined-whitelist matcher)
    # ------------------------------------------------------------------

    @property
    def min_seq_len(self) -> int:
        """Shortest sequence in this whitelist."""
        return self._min_len

    @property
    def max_seq_len(self) -> int:
        """Longest sequence in this whitelist."""
        return self._max_len

    @property
    def is_variable_length(self) -> bool:
        """True when the whitelist contains sequences of different lengths."""
        return self._min_len != self._max_len

    def __len__(self) -> int:
        return len(self._seqs)

    def __repr__(self) -> str:
        if self.is_variable_length:
            return (f"Whitelist({len(self._seqs):,} sequences, "
                    f"len {self._min_len}–{self._max_len} bp)")
        return f"Whitelist({len(self._seqs):,} sequences, {self._min_len} bp)"
