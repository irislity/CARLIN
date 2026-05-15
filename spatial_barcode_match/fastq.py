"""
Minimal FASTQ parser (4-line records, optionally gzip-compressed).

Only the sequence (line 2) and quality (line 4) are used downstream;
the header (line 1) is preserved for traceability.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


@dataclass
class FastqRecord:
    header: str      # full @… line, without the leading '@'
    seq: bytes       # nucleotide sequence
    qual: bytes      # Phred+33 quality bytes, same length as seq


def read_fastq(path: str | Path) -> Iterator[FastqRecord]:
    """
    Yield :class:`FastqRecord` objects from a FASTQ file.

    Accepts plain ``.fastq`` / ``.fq`` or gzip-compressed
    ``.fastq.gz`` / ``.fq.gz`` files.

    Raises
    ------
    ValueError
        If a record is malformed (missing lines, seq/qual length mismatch).
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open

    with opener(path, "rt") as fh:
        while True:
            header_line = fh.readline()
            if not header_line:
                break                       # clean EOF
            if not header_line.startswith("@"):
                raise ValueError(f"Expected '@' header line, got: {header_line!r}")

            seq_line  = fh.readline().rstrip("\n")
            plus_line = fh.readline()       # '+' separator — ignored
            qual_line = fh.readline().rstrip("\n")

            if not qual_line:
                raise ValueError("Truncated FASTQ record")

            seq  = seq_line.encode()
            qual = qual_line.encode()

            if len(seq) != len(qual):
                raise ValueError(
                    f"Seq/qual length mismatch for {header_line.strip()!r}: "
                    f"{len(seq)} vs {len(qual)}"
                )

            yield FastqRecord(
                header=header_line[1:].rstrip("\n"),  # strip leading '@' and newline
                seq=seq,
                qual=qual,
            )
