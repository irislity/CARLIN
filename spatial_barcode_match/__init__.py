"""
spatial_barcode_match — standalone SPATIAL-HD-v1-3P barcode matching.

Reimplements the Spaceranger JointBc1Bc2 extraction + Posterior correction
pipeline in pure Python, referencing the original Rust source throughout.

Two usage modes
---------------
**Combined whitelist** (``hd_6.5mm_bc_list.txt`` from the 10X website)::

    from spatial_barcode_match import Whitelist, process_fastq_combined, write_tsv_combined
    import os

    wl = Whitelist.from_file("hd_6.5mm_bc_list.txt", progress=True)
    results = process_fastq_combined(
        "reads_R1.fastq.gz", wl,
        num_workers=os.cpu_count(),   # parallel processing
    )
    write_tsv_combined(results, "output.tsv")

**Split whitelist** (separate bc1 / bc2 files, if available)::

    from spatial_barcode_match import Whitelist, process_fastq, write_tsv

    wl1 = Whitelist.from_file("bc1.txt")
    wl2 = Whitelist.from_file("bc2.txt")
    results = process_fastq("reads_R1.fastq.gz", wl1, wl2)
    write_tsv(results, "output.tsv")
"""

from .whitelist import Whitelist
from .extractor import extract_barcode, ExtractionResult
from .corrector import correct_segment
from .matcher import (
    match_barcode,
    MatchResult,
    process_fastq,
    write_tsv,
    ReadMatchResult,
    match_barcode_combined,
    CombinedReadMatchResult,
    process_fastq_combined,
    write_tsv_combined,
)
from .fastq import read_fastq, FastqRecord

__all__ = [
    # Whitelist
    "Whitelist",
    # Low-level
    "extract_barcode",
    "ExtractionResult",
    "correct_segment",
    # Split-whitelist API
    "match_barcode",
    "MatchResult",
    "process_fastq",
    "write_tsv",
    "ReadMatchResult",
    # Combined-whitelist API  ← use this with hd_6.5mm_bc_list.txt
    "match_barcode_combined",
    "CombinedReadMatchResult",
    "process_fastq_combined",
    "write_tsv_combined",
    # FASTQ I/O
    "read_fastq",
    "FastqRecord",
]
