"""
Unit tests for barcode_match.

Tests mirror the scenarios described in the Spaceranger source and cover:
  - exact match at default offsets
  - exact match found at a shifted offset (offset-scan logic)
  - 1-mismatch correction of bc2 when bc1 is valid
  - 1-mismatch correction of bc1 when bc2 is valid
  - both segments invalid and uncorrectable → valid=False
  - quality-weighted rejection (low-confidence correction not applied)
"""

from __future__ import annotations

import sys
import os

# Allow running directly: python test_matcher.py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spatial_barcode_match import Whitelist, match_barcode
from spatial_barcode_match.extractor import (
    BC1_DEFAULT_OFFSET,
    BC1_LENGTH,
    BC2_DEFAULT_OFFSET,
    BC2_LENGTH,
    UMI_LENGTH,
    extract_barcode,
)
from spatial_barcode_match.corrector import correct_segment, phred_error_prob

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

BC1 = b"AAAACCCCGGGGTTTT"[:BC1_LENGTH]   # 14-bp bc1
BC2 = b"TTTTGGGGCCCCAAAA"[:BC2_LENGTH]   # 14-bp bc2
UMI = b"ACGTACGTA"                        # 9-bp UMI (bases 0-8 of R1)

HIGH_QUAL = bytes([40 + 33] * 60)         # Q40 everywhere → P(err) ≈ 1e-4
LOW_QUAL  = bytes([2  + 33] * 60)         # Q2  everywhere → P(err) ≈ 0.63


def _make_r1(umi=UMI, pad_before_bc1=b"NN", bc1=BC1, bc2=BC2, pad_after=b"") -> bytes:
    """
    Assemble an R1 with the canonical layout:
      [0:9]  UMI
      [9:11] 2-byte pad (to reach offset 11)
      [11:25] bc1
      [25:39] bc2
    """
    before = umi + pad_before_bc1
    assert len(before) == BC1_DEFAULT_OFFSET, f"pre-bc1 region must be {BC1_DEFAULT_OFFSET} bytes"
    return before + bc1 + bc2 + pad_after


def _make_shifted_r1(shift: int) -> tuple[bytes, bytes]:
    """R1 where bc1 starts at (BC1_DEFAULT_OFFSET + shift) instead of 11."""
    extra = b"N" * shift
    before = UMI + b"NN" + extra          # UMI(9) + pad(2) + extra_shift
    r1 = before + BC1 + BC2 + b"NNNN"
    qual = HIGH_QUAL[:len(r1)]
    return r1, qual


# ---------------------------------------------------------------------------
# Whitelist / whitelist.py
# ---------------------------------------------------------------------------

class TestWhitelist:
    def test_contains_exact(self):
        wl = Whitelist.from_sequences([BC1, BC2])
        assert wl.contains(BC1)
        assert wl.contains(BC2)
        assert not wl.contains(b"GGGGGGGGGGGGGG")

    def test_len(self):
        wl = Whitelist.from_sequences([BC1, BC2, b"CCCCGGGGTTTTAAAA"[:BC1_LENGTH]])
        assert len(wl) == 3


# ---------------------------------------------------------------------------
# Extractor / extractor.py
# ---------------------------------------------------------------------------

class TestExtractBarcode:
    def _wl(self):
        return Whitelist.from_sequences([BC1]), Whitelist.from_sequences([BC2])

    def test_default_offset_exact_hit(self):
        r1 = _make_r1()
        qual = HIGH_QUAL[:len(r1)]
        wl1, wl2 = self._wl()
        ext = extract_barcode(r1, wl1, wl2)

        assert ext.bc1_seq == BC1
        assert ext.bc2_seq == BC2
        assert ext.umi_seq == UMI
        assert ext.bc1_in_whitelist
        assert ext.bc2_in_whitelist
        assert ext.bc1_offset == BC1_DEFAULT_OFFSET

    def test_shifted_offset_preferred(self):
        """
        If bc1 starts at offset 12 (shift=+1), the scanner should find it
        even though the canonical offset is 11.

        Mirrors: rna_read.rs  lines 336-367 (max_by_key on whitelist hits).
        """
        # Build R1 where bc1 starts at offset 12
        r1 = UMI + b"NNN" + BC1 + BC2 + b"NNNN"   # offset 12 = 9+3
        qual = HIGH_QUAL[:len(r1)]
        wl1, wl2 = self._wl()
        ext = extract_barcode(r1, wl1, wl2)

        # Should detect shift and choose offset=12 for score=(True,True)
        assert ext.bc1_seq == BC1
        assert ext.bc2_seq == BC2
        assert ext.bc1_in_whitelist
        assert ext.bc2_in_whitelist
        assert ext.bc1_offset == 12

    def test_no_hit_falls_back_to_default(self):
        """When no candidate hits the whitelist, default positions are used."""
        r1 = b"A" * 50
        qual = HIGH_QUAL[:50]
        wl1, wl2 = self._wl()
        ext = extract_barcode(r1, wl1, wl2)

        assert ext.bc1_offset == BC1_DEFAULT_OFFSET
        assert not ext.bc1_in_whitelist
        assert not ext.bc2_in_whitelist


# ---------------------------------------------------------------------------
# Corrector / corrector.py
# ---------------------------------------------------------------------------

class TestCorrector:
    def test_exact_not_corrected(self):
        wl = Whitelist.from_sequences([BC1])
        # BC1 is already in whitelist; corrector only runs on invalid barcodes
        # (the caller checks validity first), but we can still verify it won't
        # erroneously produce a different answer for a sequence 1 base off.
        mutant = bytearray(BC1)
        mutant[0] = ord("C") if BC1[0] != ord("C") else ord("A")
        result = correct_segment(bytes(mutant), HIGH_QUAL[:BC1_LENGTH], wl)
        assert result == BC1

    def test_low_quality_rejected(self):
        """
        A 1-mismatch barcode with very low quality (Q2) produces many
        competing candidates; posterior for any single one stays below 97.5%.
        With only one whitelist entry the posterior IS 100%, so we make the
        whitelist large to dilute it.
        """
        wl_seqs = [BC1]
        # Add 200 decoy barcodes differing from BC1 at position 0
        bases = b"ACGT"
        for i in range(200):
            decoy = bytearray(BC1)
            decoy[0] = bases[i % 4]
            wl_seqs.append(bytes(decoy))
        wl = Whitelist.from_sequences(wl_seqs)

        mutant = bytearray(BC1)
        mutant[0] = bases[(bases.index(BC1[0:1]) + 1) % 4]  # one substitution
        # Low quality: many candidates tie → posterior < 0.975 → no correction
        result = correct_segment(bytes(mutant), LOW_QUAL[:BC1_LENGTH], wl)
        # With 200 decoys, the posterior for the correct BC1 entry is ~1/201
        assert result is None

    def test_phred_conversion(self):
        assert abs(phred_error_prob(43) - 10 ** (-(43 - 33) / 10)) < 1e-12


# ---------------------------------------------------------------------------
# Matcher / matcher.py
# ---------------------------------------------------------------------------

class TestMatchBarcode:
    def _wls(self):
        return Whitelist.from_sequences([BC1]), Whitelist.from_sequences([BC2])

    def _run(self, r1, qual=None):
        wl1, wl2 = self._wls()
        if qual is None:
            qual = HIGH_QUAL[:len(r1)]
        return match_barcode(r1, qual, wl1, wl2)

    # --- case (true, true) -------------------------------------------------
    def test_both_exact(self):
        r1 = _make_r1()
        res = self._run(r1)
        assert res.valid
        assert res.barcode == BC1 + BC2
        assert not res.bc1_corrected
        assert not res.bc2_corrected

    # --- case (true, false) — correct bc2 ----------------------------------
    def test_bc1_valid_bc2_one_mismatch(self):
        """bc2 has one substitution at position 3; should be corrected."""
        bad_bc2 = bytearray(BC2)
        bad_bc2[3] ^= 0x01 | 0x10      # flip some bits → guaranteed different base
        # Make sure it's a valid base
        bad_bc2[3] = [b for b in b"ACGT" if b != BC2[3]][0]
        bad_bc2 = bytes(bad_bc2)

        r1 = _make_r1(bc2=bad_bc2)
        res = self._run(r1)

        assert res.valid, "bc2 should be correctable"
        assert res.bc2 == BC2
        assert res.bc2_corrected
        assert not res.bc1_corrected

    # --- case (false, true) — correct bc1 ----------------------------------
    def test_bc2_valid_bc1_one_mismatch(self):
        bad_bc1 = bytearray(BC1)
        bad_bc1[5] = [b for b in b"ACGT" if b != BC1[5]][0]
        bad_bc1 = bytes(bad_bc1)

        r1 = _make_r1(bc1=bad_bc1)
        res = self._run(r1)

        assert res.valid, "bc1 should be correctable"
        assert res.bc1 == BC1
        assert res.bc1_corrected
        assert not res.bc2_corrected

    # --- both invalid, uncorrectable ---------------------------------------
    def test_both_invalid_uncorrectable(self):
        """Random bytes → neither segment is correctable → valid=False."""
        r1 = b"AAAAAAAAANNGGGGGGGGGGGGGGTTTTTTTTTTTTTT"
        res = self._run(r1)
        assert not res.valid
        assert res.barcode is None

    # --- UMI ---------------------------------------------------------------
    def test_umi_extracted(self):
        r1 = _make_r1()
        res = self._run(r1)
        assert res.umi == UMI


# ---------------------------------------------------------------------------
# FASTQ parsing / process_fastq
# ---------------------------------------------------------------------------

class TestFastq:
    """Tests using the real-looking 44-bp reads from the user's data."""

    # Five reads taken verbatim from the user's sample FASTQ.
    # Layout: [0:9] UMI  [9:11] spacer  [11:25] bc1  [25:39] bc2  [39:44] tail
    SAMPLE_READS = [
        ("SH00291:64:BXA23111-2810:1:1101:1015:1015 1:N:0:TAACGCGTGA+CCCTAACTTC",
         b"ATGGTTGGCGTTCATCCACACATATATTCAGCATCATGCGTTTT",
         b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG9"),
        ("SH00291:64:BXA23111-2810:1:1101:1025:1015 1:N:0:TAACGCGTGA+CCCTAACTTC",
         b"TTTGACTACGTCAGAATCATCTGTCCATCTGCATGAATGGGTTT",
         b"GGGGGGGGG9GGGG9GGGGGGGGGGGGGGGGGGGGG9GGGGGGG"),
        ("SH00291:64:BXA23111-2810:1:1101:1065:1015 1:N:0:TAACGCGTGA+CCCTAACTTC",
         b"GGGTCATCGCCATGCTTCGCTCAGCAGGAGCAAGTCGGATTTTT",
         b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"),
    ]

    def _write_fastq(self, path, records):
        with open(path, "w") as fh:
            for header, seq, qual in records:
                fh.write(f"@{header}\n{seq.decode()}\n+\n{qual.decode()}\n")

    def test_fastq_reader(self):
        import tempfile, os
        from spatial_barcode_match.fastq import read_fastq

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fastq", delete=False) as fh:
            path = fh.name
            for header, seq, qual in self.SAMPLE_READS:
                fh.write(f"@{header}\n{seq.decode()}\n+\n{qual.decode()}\n")
        try:
            recs = list(read_fastq(path))
            assert len(recs) == 3
            assert recs[0].seq == self.SAMPLE_READS[0][1]
            assert recs[0].qual == self.SAMPLE_READS[0][2]
            assert recs[1].seq == self.SAMPLE_READS[1][1]
        finally:
            os.unlink(path)

    def test_process_fastq_extracts_umi_and_positions(self):
        """
        Verify that process_fastq correctly slices UMI, bc1, bc2 from
        the 44-bp reads at the expected offsets, regardless of whether
        the sequences are in the (empty) whitelists.
        """
        import tempfile, os
        from spatial_barcode_match.matcher import process_fastq
        from spatial_barcode_match.extractor import (
            UMI_OFFSET, UMI_LENGTH,
            BC1_DEFAULT_OFFSET, BC1_LENGTH,
            BC2_DEFAULT_OFFSET, BC2_LENGTH,
        )

        # Empty whitelists — all reads will be invalid, but positions must be right
        wl1 = Whitelist.from_sequences([])
        wl2 = Whitelist.from_sequences([])

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fastq", delete=False) as fh:
            path = fh.name
            for header, seq, qual in self.SAMPLE_READS:
                fh.write(f"@{header}\n{seq.decode()}\n+\n{qual.decode()}\n")
        try:
            results = list(process_fastq(path, wl1, wl2))
            assert len(results) == 3
            for res, (_, seq, _) in zip(results, self.SAMPLE_READS):
                # UMI must come from [0:9]
                assert res.umi == seq[UMI_OFFSET: UMI_OFFSET + UMI_LENGTH]
                # bc1 and bc2 must come from their canonical positions
                # (no whitelist hit → default offsets are used)
                assert res.bc1 == seq[BC1_DEFAULT_OFFSET: BC1_DEFAULT_OFFSET + BC1_LENGTH]
                assert res.bc2 == seq[BC2_DEFAULT_OFFSET: BC2_DEFAULT_OFFSET + BC2_LENGTH]
                # All invalid since whitelists are empty
                assert not res.valid
                assert res.barcode is None
        finally:
            os.unlink(path)

    def test_process_fastq_valid_when_whitelisted(self):
        """
        Plant the exact bc1/bc2 sequences from a read into the whitelists;
        that read must come back valid=True with the correct 28-bp barcode.
        """
        import tempfile, os
        from spatial_barcode_match.matcher import process_fastq
        from spatial_barcode_match.extractor import BC1_DEFAULT_OFFSET, BC1_LENGTH, BC2_DEFAULT_OFFSET, BC2_LENGTH

        header, seq, qual = self.SAMPLE_READS[0]
        bc1_seq = seq[BC1_DEFAULT_OFFSET: BC1_DEFAULT_OFFSET + BC1_LENGTH]
        bc2_seq = seq[BC2_DEFAULT_OFFSET: BC2_DEFAULT_OFFSET + BC2_LENGTH]

        wl1 = Whitelist.from_sequences([bc1_seq])
        wl2 = Whitelist.from_sequences([bc2_seq])

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fastq", delete=False) as fh:
            path = fh.name
            fh.write(f"@{header}\n{seq.decode()}\n+\n{qual.decode()}\n")
        try:
            results = list(process_fastq(path, wl1, wl2))
            assert len(results) == 1
            r = results[0]
            assert r.valid
            assert r.barcode == bc1_seq + bc2_seq
            assert not r.bc1_corrected
            assert not r.bc2_corrected
        finally:
            os.unlink(path)


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_all():
    suites = [TestWhitelist, TestExtractBarcode, TestCorrector, TestMatchBarcode, TestFastq]
    passed = 0
    failed = 0
    for cls in suites:
        obj = cls()
        for name in dir(obj):
            if not name.startswith("test_"):
                continue
            try:
                getattr(obj, name)()
                print(f"  PASS  {cls.__name__}.{name}")
                passed += 1
            except Exception as exc:
                print(f"  FAIL  {cls.__name__}.{name}  →  {exc}")
                failed += 1
    print(f"\n{passed} passed, {failed} failed")
    return failed == 0


if __name__ == "__main__":
    ok = run_all()
    sys.exit(0 if ok else 1)
