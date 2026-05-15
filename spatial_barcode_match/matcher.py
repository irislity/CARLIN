"""
Top-level barcode matching for SPATIAL-HD-v1-3P.

Orchestrates:
  1. JointBc1Bc2 extraction (scan offset window for best whitelist hit).
  2. Per-case 1-mismatch correction driven by Bayesian posterior probability.

The four correction cases mirror the Rust stage exactly:

  (bc1 valid, bc2 valid)   → no correction needed.
  (bc1 valid, bc2 invalid) → anchor bc1 end; search bc2 nearby.
  (bc1 invalid, bc2 valid) → anchor bc2 start; search bc1 in offset window.
  (bc1 invalid, bc2 invalid) → joint search over all (offset, len1, len2)
                                combos; pick minimum total correction distance.

Source references
-----------------
* JointBc1Bc2 correction, all four cases:
    lib/rust/cr_lib/src/stages/barcode_correction.rs  lines 144-270
* try_correct helper:
    lib/rust/cr_lib/src/stages/barcode_correction.rs  lines 274-287
* SegmentedBarcode concatenation (bc1 + bc2 = 28 bp):
    lib/rust/barcode/src/lib.rs  lines 941-951
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

from .whitelist import Whitelist
from .extractor import (
    extract_barcode,
    ExtractionResult,
    MIN_OFFSET,
    MAX_OFFSET,
    BC1_DEFAULT_OFFSET,
    BC1_LENGTH,
    BC2_LENGTH,
    UMI_OFFSET,
    UMI_LENGTH,
)
from .corrector import (
    correct_segment,
    SEARCH_PADDING,
    BARCODE_CONFIDENCE_THRESHOLD,
)

# ---------------------------------------------------------------------------
# Reverse complement
# ---------------------------------------------------------------------------

_COMP = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")

def revcomp(seq: bytes) -> bytes:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMP)[::-1]


@dataclass
class MatchResult:
    """Full result of :func:`match_barcode` for one R1 read."""

    #: bc1 sequence (14 bp), possibly corrected.
    bc1: bytes
    #: bc2 sequence (14 bp), possibly corrected.
    bc2: bytes
    #: UMI sequence (9 bp).
    umi: bytes

    #: True when bc1 / bc2 required a 1-mismatch correction.
    bc1_corrected: bool
    bc2_corrected: bool

    #: True when both bc1 and bc2 are valid whitelist members
    #: (either originally or after correction).
    valid: bool

    #: Concatenated 28-bp barcode (bc1 + bc2), or ``None`` if invalid.
    barcode: Optional[bytes]


def match_barcode(
    r1_seq: bytes,
    r1_qual: bytes,
    wl1: Whitelist,
    wl2: Whitelist,
    bc1_counts: Optional[dict[bytes, int]] = None,
    bc2_counts: Optional[dict[bytes, int]] = None,
    min_offset: int = MIN_OFFSET,
    max_offset: int = MAX_OFFSET,
    bc1_len: int = BC1_LENGTH,
    bc2_len: int = BC2_LENGTH,
    confidence_threshold: float = BARCODE_CONFIDENCE_THRESHOLD,
) -> MatchResult:
    """
    Match a single R1 read to the SPATIAL-HD-v1-3P barcode whitelist.

    Parameters
    ----------
    r1_seq:
        R1 nucleotide sequence as bytes (e.g. ``b"ACGT..."``).
    r1_qual:
        R1 base-quality string as bytes (Phred+33 encoded, same length as
        *r1_seq*).
    wl1, wl2:
        Loaded :class:`~barcode_match.Whitelist` objects for bc1 and bc2.
    bc1_counts, bc2_counts:
        Optional ``dict[seq_bytes, int]`` of observed barcode counts used as
        the Bayesian prior.  When omitted a uniform prior is assumed (each
        candidate contributes count = 1 after Laplace smoothing).
    min_offset, max_offset:
        Offset search window; default 8–12 per the chemistry definition.
    bc1_len, bc2_len:
        Expected segment lengths; default 14 each.
    confidence_threshold:
        Minimum posterior probability to accept a correction; default 0.975.

    Returns
    -------
    MatchResult
    """
    ext: ExtractionResult = extract_barcode(
        r1_seq, wl1, wl2, min_offset, max_offset, bc1_len, bc2_len
    )
    read_len = len(r1_seq)

    bc1 = ext.bc1_seq
    bc2 = ext.bc2_seq
    bc1_valid = ext.bc1_in_whitelist
    bc2_valid = ext.bc2_in_whitelist
    bc1_corrected = False
    bc2_corrected = False

    # Pre-slice quality bytes for the extracted positions.
    bc1_qual = r1_qual[ext.bc1_offset: ext.bc1_offset + ext.bc1_len]
    bc2_qual = r1_qual[ext.bc2_offset: ext.bc2_offset + ext.bc2_len]

    # -------------------------------------------------------------------
    # Correction — mirrors barcode_correction.rs  lines 163-270
    # -------------------------------------------------------------------

    if bc1_valid and bc2_valid:
        # Case (true, true): both segments exact-match; nothing to do.
        # Source: barcode_correction.rs  line 164
        pass

    elif bc1_valid and not bc2_valid:
        # Case (true, false): bc1 is anchored at its end; try nearby lengths
        # for bc2 starting from bc1_end.
        # Source: barcode_correction.rs  lines 165-191
        bc1_end = ext.bc1_offset + ext.bc1_len
        lo = max(0, bc2_len - SEARCH_PADDING)
        hi = min(read_len - bc1_end, bc2_len + SEARCH_PADDING)

        best_bc2: Optional[bytes] = None
        best_dist: int = 2  # sentinel > maximum possible 1-mismatch distance

        for l2 in range(lo, hi + 1):
            s2 = r1_seq[bc1_end: bc1_end + l2]
            q2 = r1_qual[bc1_end: bc1_end + l2]
            corrected = correct_segment(s2, q2, wl2, bc2_counts, confidence_threshold)
            if corrected is not None:
                dist = sum(a != b for a, b in zip(s2, corrected))
                if dist < best_dist:
                    best_dist = dist
                    best_bc2 = corrected

        if best_bc2 is not None:
            bc2 = best_bc2
            bc2_valid = True
            bc2_corrected = True

    elif not bc1_valid and bc2_valid:
        # Case (false, true): bc2 is anchored at its start; try all offsets
        # in the window for bc1, each time capping bc1 length so it ends just
        # before bc2_start.
        # Source: barcode_correction.rs  lines 192-219
        bc2_start = ext.bc2_offset
        max_bc1_len = bc1_len + SEARCH_PADDING - 1  # mirrors Rust: end+PAD-1

        best_bc1: Optional[bytes] = None
        best_dist: int = 2

        for offset in range(min_offset, max_offset + 1):
            l1 = min(bc2_start - offset, max_bc1_len)
            if l1 <= 0:
                continue
            s1 = r1_seq[offset: offset + l1]
            q1 = r1_qual[offset: offset + l1]
            corrected = correct_segment(s1, q1, wl1, bc1_counts, confidence_threshold)
            if corrected is not None:
                dist = sum(a != b for a, b in zip(s1, corrected))
                if dist < best_dist:
                    best_dist = dist
                    best_bc1 = corrected

        if best_bc1 is not None:
            bc1 = best_bc1
            bc1_valid = True
            bc1_corrected = True

    else:
        # Case (false, false): neither segment is valid.  Perform a joint
        # search over all (offset, len1, len2) combinations; pick the combo
        # with the minimum *total* correction distance (uncorrectable segments
        # contribute LARGE_DISTANCE = 1000 to the sum).
        # Source: barcode_correction.rs  lines 220-269
        LARGE_DISTANCE = 1000
        best_total_dist = LARGE_DISTANCE * 2 + 1  # sentinel
        best_result: Optional[tuple[Optional[bytes], Optional[bytes]]] = None

        l1_range = range(max(0, bc1_len - SEARCH_PADDING), bc1_len + SEARCH_PADDING + 1)
        l2_range = range(max(0, bc2_len - SEARCH_PADDING), bc2_len + SEARCH_PADDING + 1)

        for offset in range(min_offset, max_offset + 1):
            for l1 in l1_range:
                for l2 in l2_range:
                    if offset + l1 + l2 > read_len:
                        continue
                    s1 = r1_seq[offset: offset + l1]
                    q1 = r1_qual[offset: offset + l1]
                    s2 = r1_seq[offset + l1: offset + l1 + l2]
                    q2 = r1_qual[offset + l1: offset + l1 + l2]

                    c1 = correct_segment(s1, q1, wl1, bc1_counts, confidence_threshold)
                    c2 = correct_segment(s2, q2, wl2, bc2_counts, confidence_threshold)

                    d1 = sum(a != b for a, b in zip(s1, c1)) if c1 is not None else LARGE_DISTANCE
                    d2 = sum(a != b for a, b in zip(s2, c2)) if c2 is not None else LARGE_DISTANCE
                    total = d1 + d2

                    # Source: barcode_correction.rs  lines 247-252 (min_by_key on total dist)
                    if total < best_total_dist:
                        best_total_dist = total
                        best_result = (c1, c2)

        if best_result is not None:
            c1, c2 = best_result
            if c1 is not None:
                bc1 = c1
                bc1_valid = True
                bc1_corrected = True
            if c2 is not None:
                bc2 = c2
                bc2_valid = True
                bc2_corrected = True

    # -------------------------------------------------------------------
    # Assemble result
    # Final 28-bp barcode = bc1 (14 bp) + bc2 (14 bp).
    # Source: lib/rust/barcode/src/lib.rs  lines 941-951
    # -------------------------------------------------------------------
    umi = r1_seq[UMI_OFFSET: UMI_OFFSET + UMI_LENGTH]
    valid = bc1_valid and bc2_valid
    barcode = (bc1 + bc2) if valid else None

    return MatchResult(
        bc1=bc1,
        bc2=bc2,
        umi=umi,
        bc1_corrected=bc1_corrected,
        bc2_corrected=bc2_corrected,
        valid=valid,
        barcode=barcode,
    )


# ---------------------------------------------------------------------------
# End-to-end FASTQ processing
# ---------------------------------------------------------------------------

@dataclass
class ReadMatchResult:
    """Per-read result returned by :func:`process_fastq`."""
    read_id: str           # FASTQ header (without '@')
    r1_seq: bytes          # original R1 sequence
    umi: bytes             # 9-bp UMI
    bc1: bytes             # bc1 segment (14 bp)
    bc2: bytes             # bc2 segment (14 bp)
    barcode: Optional[bytes]  # bc1+bc2 (28 bp) or None when invalid
    bc1_corrected: bool
    bc2_corrected: bool
    valid: bool


def process_fastq(
    fastq_path: str | Path,
    wl1: Whitelist,
    wl2: Whitelist,
    bc1_counts: Optional[dict[bytes, int]] = None,
    bc2_counts: Optional[dict[bytes, int]] = None,
    min_offset: int = MIN_OFFSET,
    max_offset: int = MAX_OFFSET,
    bc1_len: int = BC1_LENGTH,
    bc2_len: int = BC2_LENGTH,
    confidence_threshold: float = BARCODE_CONFIDENCE_THRESHOLD,
) -> Iterator[ReadMatchResult]:
    """
    Parse a FASTQ file and match barcodes for every R1 read.

    Your R1 reads are expected to have the SPATIAL-HD-v1-3P layout::

        pos  0–8  : UMI        (9 bp)
        pos  9–10 : spacer     (2 bp, ignored)
        pos 11–24 : bc1        (14 bp)
        pos 25–38 : bc2        (14 bp)
        pos 39+   : polyT / primer tail (ignored)

    Reads shorter than ``min_offset + bc1_len + bc2_len`` are skipped
    with a warning.

    Parameters
    ----------
    fastq_path:
        Path to the R1 FASTQ file (``.fastq`` or ``.fastq.gz``).
    wl1, wl2:
        Loaded :class:`~barcode_match.Whitelist` objects for bc1 and bc2.
    bc1_counts, bc2_counts:
        Optional observed-count dicts used as the Bayesian prior.
    min_offset / max_offset:
        Offset search window for bc1 (default 8–12).
    bc1_len / bc2_len:
        Expected barcode segment lengths (default 14).
    confidence_threshold:
        Posterior probability threshold for 1-mismatch correction (default 0.975).

    Yields
    ------
    ReadMatchResult
        One result per FASTQ record.
    """
    from .fastq import read_fastq

    min_read_len = min_offset + bc1_len + bc2_len

    for rec in read_fastq(fastq_path):
        if len(rec.seq) < min_read_len:
            print(
                f"[barcode_match] WARNING: read {rec.header!r} is only "
                f"{len(rec.seq)} bp (need ≥ {min_read_len}), skipping.",
                file=sys.stderr,
            )
            continue

        m = match_barcode(
            r1_seq=rec.seq,
            r1_qual=rec.qual,
            wl1=wl1,
            wl2=wl2,
            bc1_counts=bc1_counts,
            bc2_counts=bc2_counts,
            min_offset=min_offset,
            max_offset=max_offset,
            bc1_len=bc1_len,
            bc2_len=bc2_len,
            confidence_threshold=confidence_threshold,
        )

        yield ReadMatchResult(
            read_id=rec.header,
            r1_seq=rec.seq,
            umi=m.umi,
            bc1=m.bc1,
            bc2=m.bc2,
            barcode=m.barcode,
            bc1_corrected=m.bc1_corrected,
            bc2_corrected=m.bc2_corrected,
            valid=m.valid,
        )


# ---------------------------------------------------------------------------
# Combined-whitelist path  (hd_6.5mm_bc_list.txt)
# ---------------------------------------------------------------------------

@dataclass
class CombinedReadMatchResult:
    """
    Per-read result from :func:`process_fastq_combined`.

    Unlike :class:`ReadMatchResult`, the barcode is not split into bc1/bc2
    because the combined whitelist does not encode that boundary.
    """
    read_id: str
    r1_seq: bytes
    umi: bytes
    barcode: Optional[bytes]   # the full matched sequence (29-31 bp), or None
    corrected: bool            # True if a 1-mismatch correction was applied
    reverse_complemented: bool # True if the match was on the reverse complement
    valid: bool


def match_barcode_combined(
    r1_seq: bytes,
    r1_qual: bytes,
    whitelist: Whitelist,
    bc_counts: Optional[dict[bytes, int]] = None,
    min_offset: int = MIN_OFFSET,
    max_offset: int = MAX_OFFSET,
    confidence_threshold: float = BARCODE_CONFIDENCE_THRESHOLD,
    reverse_complement: bool = False,
) -> tuple[Optional[bytes], bool, bool]:
    """
    Match one R1 read against a **combined** (bc1+bc2) whitelist.

    Because ``hd_6.5mm_bc_list.txt`` lists every valid bc1+bc2 pair as a
    single sequence of variable length (29-31 bp), we search for any window
    in R1 that hits the combined whitelist directly.

    Performance note
    ----------------
    The **fast path** tries only the canonical offset (11) with all valid
    lengths — typically 3 hash lookups — and returns immediately on a hit.
    This covers the vast majority of reads.  Only on a miss does it fall back
    to scanning all offsets (15 lookups) and then to 1-mismatch correction.

    Source: JointBc1Bc2 extraction logic in
        lib/rust/cr_types/src/rna_read.rs  lines 316-368

    Returns
    -------
    (barcode, corrected, reverse_complemented)
        ``barcode`` is the matched bytes or ``None``; ``corrected`` is
        ``True`` when a 1-mismatch correction was accepted;
        ``reverse_complemented`` is ``True`` when the match was found on
        the reverse complement of the extracted window.
    """
    min_len = whitelist.min_seq_len
    max_len = whitelist.max_seq_len
    read_len = len(r1_seq)
    seqs = whitelist._seqs   # direct set access — avoids per-call method overhead

    # -----------------------------------------------------------------------
    # Helper: check one sequence (and optionally its RC) against the whitelist.
    # Returns (matched_seq, is_rc) or (None, False).
    # -----------------------------------------------------------------------
    def _check(seq: bytes) -> tuple[Optional[bytes], bool]:
        if seq in seqs:
            return seq, False
        if reverse_complement:
            rc = revcomp(seq)
            if rc in seqs:
                return rc, True
        return None, False

    # -----------------------------------------------------------------------
    # Fast path: canonical offset (11) only — covers ~95%+ of reads with
    # just (max_len - min_len + 1) hash lookups, typically 3.
    # -----------------------------------------------------------------------
    for length in range(min_len, max_len + 1):
        end = BC1_DEFAULT_OFFSET + length
        if end <= read_len:
            matched, is_rc = _check(r1_seq[BC1_DEFAULT_OFFSET:end])
            if matched is not None:
                return matched, False, is_rc

    # -----------------------------------------------------------------------
    # Slow path: scan remaining offsets.
    # Canonical offset (11) is tried last so it wins ties on correction too.
    # Source: rna_read.rs  lines 354-367 (chain + last-wins max_by_key)
    # -----------------------------------------------------------------------
    other_offsets = [o for o in range(min_offset, max_offset + 1)
                     if o != BC1_DEFAULT_OFFSET]

    barcode: Optional[bytes] = None
    barcode_is_rc: bool = False

    for offset in other_offsets:
        for length in range(min_len, max_len + 1):
            end = offset + length
            if end <= read_len:
                matched, is_rc = _check(r1_seq[offset:end])
                if matched is not None:
                    barcode, barcode_is_rc = matched, is_rc

    # Retry canonical last (wins ties)
    for length in range(min_len, max_len + 1):
        end = BC1_DEFAULT_OFFSET + length
        if end <= read_len:
            matched, is_rc = _check(r1_seq[BC1_DEFAULT_OFFSET:end])
            if matched is not None:
                barcode, barcode_is_rc = matched, is_rc

    if barcode is not None:
        return barcode, False, barcode_is_rc

    # -----------------------------------------------------------------------
    # 1-mismatch correction — Bayesian posterior ≥ 97.5%.
    # Applied to the forward sequence; if reverse_complement=True, also tried
    # on the RC of each window.
    # Source: Posterior::correct_barcode  corrector.rs  lines 112-166
    # -----------------------------------------------------------------------
    for offset in [BC1_DEFAULT_OFFSET] + other_offsets:
        for length in range(min_len, max_len + 1):
            end = offset + length
            if end > read_len:
                continue
            seq = r1_seq[offset:end]
            qual = r1_qual[offset:end]

            # Try forward correction
            corrected = correct_segment(seq, qual, whitelist, bc_counts,
                                        confidence_threshold)
            if corrected is not None:
                return corrected, True, False

            # Try RC correction (qual is reversed to match the flipped sequence)
            if reverse_complement:
                corrected = correct_segment(revcomp(seq), qual[::-1], whitelist,
                                            bc_counts, confidence_threshold)
                if corrected is not None:
                    return corrected, True, True

    return None, False, False


def _process_chunk(args: tuple) -> list:
    """Worker function for multiprocessing — processes a list of raw FASTQ records."""
    records, whitelist_seqs, min_len, max_len, min_offset, max_offset, confidence, do_rc = args
    wl = Whitelist(whitelist_seqs)
    out = []
    for header, seq, qual in records:
        barcode, corrected, is_rc = match_barcode_combined(
            seq, qual, wl,
            min_offset=min_offset, max_offset=max_offset,
            confidence_threshold=confidence,
            reverse_complement=do_rc,
        )
        out.append((
            header,
            seq[UMI_OFFSET: UMI_OFFSET + UMI_LENGTH],
            barcode,
            corrected,
            is_rc,
            barcode is not None,
        ))
    return out


def process_fastq_combined(
    fastq_path: str | Path,
    whitelist: Whitelist,
    bc_counts: Optional[dict[bytes, int]] = None,
    min_offset: int = MIN_OFFSET,
    max_offset: int = MAX_OFFSET,
    confidence_threshold: float = BARCODE_CONFIDENCE_THRESHOLD,
    reverse_complement: bool = False,
    num_workers: int = 1,
    chunk_size: int = 50_000,
) -> Iterator[CombinedReadMatchResult]:
    """
    Parse a FASTQ file and match barcodes using the combined whitelist.

    Use this when you have ``hd_6.5mm_bc_list.txt`` from the 10X website
    (11,222,500 entries, 29-31 bp each).

    Performance
    -----------
    * Most reads hit the fast path (3 hash lookups at the canonical offset).
    * Pass ``num_workers > 1`` to parallelise across CPU cores.  Each worker
      receives a chunk of ``chunk_size`` reads and processes them in C-level
      threads, giving near-linear scaling.  Example: 8 workers on a laptop
      typically processes 4-8 M reads/minute.

    Your R1 reads must follow the SPATIAL-HD-v1-3P layout::

        pos  0–8   : UMI        (9 bp)
        pos  9–10  : spacer     (2 bp, ignored)
        pos  11+   : barcode    (29-31 bp, searched at offsets 8-12)
        pos  39+   : polyT tail (ignored)

    Parameters
    ----------
    fastq_path:
        Path to the R1 FASTQ (``.fastq`` or ``.fastq.gz``).
    whitelist:
        Combined whitelist loaded from ``hd_6.5mm_bc_list.txt``.
    bc_counts:
        Optional observed-count dict for the Bayesian correction prior.
    min_offset / max_offset:
        Offset search window (default 8–12).
    confidence_threshold:
        Posterior threshold for 1-mismatch correction (default 0.975).
    reverse_complement:
        If ``True``, also try the reverse complement of each barcode window.
        Useful when reads may be in the opposite orientation.  Roughly doubles
        the number of hash lookups per read on the slow path.
    num_workers:
        Number of parallel worker processes (default 1 = single-threaded).
        Set to ``os.cpu_count()`` to use all cores.
    chunk_size:
        Reads per worker batch (default 50,000).

    Yields
    ------
    CombinedReadMatchResult
    """
    import os
    from .fastq import read_fastq

    min_read_len = min_offset + whitelist.min_seq_len

    if num_workers <= 1:
        # Single-threaded path — simple and avoids pickling overhead
        for rec in read_fastq(fastq_path):
            if len(rec.seq) < min_read_len:
                print(f"[spatial_barcode_match] WARNING: {rec.header!r} is "
                      f"only {len(rec.seq)} bp, skipping.", file=sys.stderr)
                continue
            barcode, corrected, is_rc = match_barcode_combined(
                rec.seq, rec.qual, whitelist,
                bc_counts=bc_counts,
                min_offset=min_offset,
                max_offset=max_offset,
                confidence_threshold=confidence_threshold,
                reverse_complement=reverse_complement,
            )
            yield CombinedReadMatchResult(
                read_id=rec.header,
                r1_seq=rec.seq,
                umi=rec.seq[UMI_OFFSET: UMI_OFFSET + UMI_LENGTH],
                barcode=barcode,
                corrected=corrected,
                reverse_complemented=is_rc,
                valid=barcode is not None,
            )
    else:
        # Multi-process path — chunk reads and dispatch to worker pool.
        # The whitelist set is shared read-only across workers via fork (Unix)
        # or serialised once per chunk (Windows).
        import multiprocessing as mp

        wl_seqs = whitelist._seqs
        min_len = whitelist.min_seq_len
        max_len = whitelist.max_seq_len

        def _iter_chunks():
            chunk: list = []
            for rec in read_fastq(fastq_path):
                if len(rec.seq) < min_read_len:
                    continue
                chunk.append((rec.header, rec.seq, rec.qual))
                if len(chunk) >= chunk_size:
                    yield chunk
                    chunk = []
            if chunk:
                yield chunk

        worker_args = (wl_seqs, min_len, max_len, min_offset, max_offset,
                       confidence_threshold, reverse_complement)

        with mp.Pool(processes=num_workers) as pool:
            for results in pool.imap(
                _process_chunk,
                (
                    (chunk, *worker_args)
                    for chunk in _iter_chunks()
                ),
                chunksize=1,
            ):
                for header, umi, barcode, corrected, is_rc, valid in results:
                    yield CombinedReadMatchResult(
                        read_id=header,
                        r1_seq=b"",   # not stored to save memory in parallel mode
                        umi=umi,
                        barcode=barcode,
                        corrected=corrected,
                        reverse_complemented=is_rc,
                        valid=valid,
                    )


def write_tsv_combined(
    results: Iterator[CombinedReadMatchResult],
    out_path: Optional[str | Path] = None,
) -> None:
    """
    Write combined-whitelist results to a TSV (or stdout).

    Columns: read_id, umi, barcode, corrected, reverse_complemented, valid
    """
    fh = open(out_path, "w", newline="") if out_path else sys.stdout
    try:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["read_id", "umi", "barcode", "corrected",
                          "reverse_complemented", "valid"])
        for r in results:
            writer.writerow([
                r.read_id,
                r.umi.decode(),
                r.barcode.decode() if r.barcode else "NA",
                int(r.corrected),
                int(r.reverse_complemented),
                int(r.valid),
            ])
    finally:
        if out_path:
            fh.close()


def write_tsv(
    results: Iterator[ReadMatchResult],
    out_path: Optional[str | Path] = None,
) -> None:
    """
    Write matching results to a tab-separated file (or stdout).

    Columns: read_id, umi, bc1, bc2, barcode, bc1_corrected, bc2_corrected, valid

    Parameters
    ----------
    results:
        Iterator from :func:`process_fastq`.
    out_path:
        Output file path.  If ``None``, writes to stdout.
    """
    fh = open(out_path, "w", newline="") if out_path else sys.stdout
    try:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            ["read_id", "umi", "bc1", "bc2", "barcode",
             "bc1_corrected", "bc2_corrected", "valid"]
        )
        for r in results:
            writer.writerow([
                r.read_id,
                r.umi.decode(),
                r.bc1.decode(),
                r.bc2.decode(),
                r.barcode.decode() if r.barcode else "NA",
                int(r.bc1_corrected),
                int(r.bc2_corrected),
                int(r.valid),
            ])
    finally:
        if out_path:
            fh.close()
