"""
Minimal FASTQ parser. Mirrors spatial_barcode_match/fastq.py.
Handles plain .fastq and gzip-compressed .fastq.gz files.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


@dataclass
class FastqRecord:
    header: str    # full @… line, without the leading '@'
    seq: bytes     # nucleotide sequence
    qual: bytes    # Phred+33 quality bytes, same length as seq


def read_fastq(path: str | Path) -> Iterator[FastqRecord]:
    """
    Yield FastqRecord objects from a FASTQ file.

    Accepts plain .fastq / .fq or gzip-compressed .fastq.gz / .fq.gz files.
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open

    with opener(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq  = fh.readline().rstrip("\n").encode()
            fh.readline()  # '+' line — ignored
            qual = fh.readline().rstrip("\n").encode()
            yield FastqRecord(
                header=header[1:].rstrip("\n"),
                seq=seq,
                qual=qual,
            )
