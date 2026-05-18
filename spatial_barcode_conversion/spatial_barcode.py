"""
spatial_barcode.py
------------------
Python port of lib/rust/barcode/src/binned.rs

Implements SquareBinIndex: parses and formats Visium HD spatial barcodes.

Spatial barcode format:
    s_{size_um:03}um_{row:05}_{col:05}-{gem_group}

Example:
    s_008um_00042_00137-1  →  row=42, col=137, size_um=8
"""

from __future__ import annotations
import re

SQUARE_BIN_PREFIX = "s"
MICROMETER = "um"

# Regex for parsing: s_008um_00042_00137 (gem group suffix optional)
_PATTERN = re.compile(
    r"^s_(\d{3})um_(\d{5})_(\d{5})(?:-\d+)?$"
)


class SquareBinIndex:
    """
    Represents a single Visium HD spatial barcode position.

    Attributes
    ----------
    row : int
        Zero-based row index on the slide grid.
    col : int
        Zero-based column index on the slide grid.
    size_um : int
        Bin size in micrometers (e.g., 2, 8, 16, 48).
    """

    def __init__(self, row: int, col: int, size_um: int):
        self.row = row
        self.col = col
        self.size_um = size_um

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_barcode_string(cls, barcode: str) -> "SquareBinIndex":
        """Parse a spatial barcode string (gem group suffix is ignored).

        Parameters
        ----------
        barcode : str
            e.g. "s_008um_00042_00137-1" or "s_008um_00042_00137"

        Raises
        ------
        ValueError
            If the string does not match the expected format.
        """
        m = _PATTERN.match(barcode.strip())
        if m is None:
            raise ValueError(f"Unable to parse '{barcode}' as a spatial barcode.")
        size_um = int(m.group(1))
        row = int(m.group(2))
        col = int(m.group(3))
        return cls(row=row, col=col, size_um=size_um)

    @classmethod
    def from_bytes(cls, b: bytes) -> "SquareBinIndex":
        """Parse from bytes (e.g. values read directly from a matrix .bcs array)."""
        return cls.from_barcode_string(b.decode())

    # ------------------------------------------------------------------
    # Formatting
    # ------------------------------------------------------------------

    def __str__(self) -> str:
        """Format without gem group suffix."""
        return f"{SQUARE_BIN_PREFIX}_{self.size_um:03}{MICROMETER}_{self.row:05}_{self.col:05}"

    def __repr__(self) -> str:
        return f"SquareBinIndex(row={self.row}, col={self.col}, size_um={self.size_um})"

    def with_gem_group(self, gem_group: int = 1) -> str:
        """Format with gem group suffix (standard output form)."""
        return f"{self}-{gem_group}"

    # ------------------------------------------------------------------
    # Binning
    # ------------------------------------------------------------------

    def binned(self, bin_scale: int) -> "SquareBinIndex":
        """Derive a coarser-resolution barcode by aggregating bin_scale x bin_scale tiles.

        Parameters
        ----------
        bin_scale : int
            Factor to multiply size_um and divide row/col by.
            E.g. bin_scale=4 converts 2µm → 8µm bins.
        """
        return SquareBinIndex(
            row=self.row // bin_scale,
            col=self.col // bin_scale,
            size_um=self.size_um * bin_scale,
        )

    def scale(self, pitch_um: int) -> int:
        """Number of pitch-sized spots this bin spans along one axis."""
        return self.size_um // pitch_um

    # ------------------------------------------------------------------
    # Comparison / hashing
    # ------------------------------------------------------------------

    def __eq__(self, other) -> bool:
        if not isinstance(other, SquareBinIndex):
            return NotImplemented
        return (self.row, self.col, self.size_um) == (other.row, other.col, other.size_um)

    def __hash__(self) -> int:
        return hash((self.row, self.col, self.size_um))

    def __lt__(self, other) -> bool:
        return (self.size_um, self.row, self.col) < (other.size_um, other.row, other.col)


# ---------------------------------------------------------------------------
# Convenience helpers
# ---------------------------------------------------------------------------

def parse_spatial_barcode(barcode: str | bytes) -> tuple[int, int]:
    """Return (row, col) from a spatial barcode string or bytes object.

    Parameters
    ----------
    barcode : str | bytes
        e.g. b"s_008um_00042_00137-1" or "s_008um_00042_00137-1"
    """
    if isinstance(barcode, bytes):
        barcode = barcode.decode()
    idx = SquareBinIndex.from_barcode_string(barcode)
    return idx.row, idx.col


def make_spatial_barcode(row: int, col: int, size_um: int = 8, gem_group: int = 1) -> str:
    """Build a spatial barcode string from coordinates.

    Parameters
    ----------
    row : int
        Grid row (zero-based).
    col : int
        Grid column (zero-based).
    size_um : int
        Bin size in micrometers (default 8).
    gem_group : int
        Gem group (default 1).
    """
    return SquareBinIndex(row=row, col=col, size_um=size_um).with_gem_group(gem_group)


# ---------------------------------------------------------------------------
# Quick test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # Round-trip test (mirrors test_square_bin in binned.rs)
    b = SquareBinIndex(row=1, col=2, size_um=16)
    assert str(b) == "s_016um_00001_00002", str(b)
    assert b.with_gem_group() == "s_016um_00001_00002-1"
    assert SquareBinIndex.from_barcode_string("s_016um_00001_00002") == b
    assert SquareBinIndex.from_barcode_string("s_016um_00001_00002-1") == b

    # Binning test: 2µm base → 8µm (scale=4)
    base = SquareBinIndex(row=42, col=137, size_um=2)
    binned = base.binned(4)
    assert binned.size_um == 8
    assert binned.row == 10
    assert binned.col == 34

    print("All spatial_barcode.py tests passed.")
    print(f"Example: {make_spatial_barcode(42, 137, size_um=8)}")
