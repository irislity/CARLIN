"""
Bayesian posterior barcode correction (up to 1 mismatch).

For each invalid barcode segment the corrector enumerates every single-base
substitution, checks whitelist membership, and computes a posterior
probability weighted by the sequencing quality score and a barcode-count
prior.  A correction is accepted only when the best candidate's posterior
exceeds 97.5 %.

Source references
-----------------
* Constants (BC_MAX_QV, BASE_OPTS, BARCODE_CONFIDENCE_THRESHOLD):
    lib/rust/barcode/src/corrector.rs  lines 9-11, 84
* Posterior::correct_barcode (core algorithm):
    lib/rust/barcode/src/corrector.rs  lines 112-166
* fn probability (Phred → error probability):
    lib/rust/barcode/src/corrector.rs  lines 168-172
* SEARCH_PADDING (length padding for joint correction):
    lib/rust/cr_lib/src/stages/barcode_correction.rs  line 156
"""

from __future__ import annotations

from typing import Optional

# ---------------------------------------------------------------------------
# Constants — mirrors corrector.rs
# ---------------------------------------------------------------------------

#: Cap on Phred quality value before conversion to probability.
#: Source: corrector.rs  line 9
BC_MAX_QV: int = 66

#: The four possible bases for substitution trials.
#: Source: corrector.rs  line 10
BASE_OPTS: bytes = b"ACGT"

#: Minimum posterior probability required to accept a 1-mismatch correction.
#: Source: corrector.rs  line 84
BARCODE_CONFIDENCE_THRESHOLD: float = 0.975

#: Extra length positions searched on either side of the nominal length when
#: both bc1 and bc2 are invalid.
#: Source: barcode_correction.rs  line 156
SEARCH_PADDING: int = 1


# ---------------------------------------------------------------------------
# Core helpers
# ---------------------------------------------------------------------------

def phred_error_prob(qv: int) -> float:
    """
    Convert a raw quality byte (Phred+33 encoded) to an error probability.

    ``P(error) = 10 ^ (-(qv - 33) / 10)``

    Source: fn probability  lib/rust/barcode/src/corrector.rs  lines 168-172
    """
    q = min(qv, BC_MAX_QV) - 33
    return 10.0 ** (-q / 10.0)


def correct_segment(
    seq: bytes,
    qual: bytes,
    whitelist: "Whitelist",
    bc_counts: Optional[dict[bytes, int]] = None,
    confidence_threshold: float = BARCODE_CONFIDENCE_THRESHOLD,
) -> Optional[bytes]:
    """
    Try to correct *seq* onto the whitelist using a single-base substitution.

    For every position × alternative base:

    1. Substitute the base and check whitelist membership.
    2. Compute ``likelihood = P(error at qv) × (1 + observed_count)``
       (Laplace smoothing: adding 1 ensures unobserved barcodes still
       contribute a small prior).
    3. Track the candidate with the highest likelihood and the total
       likelihood across all candidates.
    4. Accept the best candidate if
       ``best_likelihood / total_likelihood >= confidence_threshold``.

    Returns the corrected sequence on success, or ``None`` if no single
    substitution passes the threshold.

    Source: Posterior::correct_barcode
        lib/rust/barcode/src/corrector.rs  lines 112-166
    """
    from .whitelist import Whitelist  # avoid circular import at module level

    a = bytearray(seq)
    best_bc: Optional[bytes] = None
    best_likelihood: float = 0.0
    total_likelihood: float = 0.0

    for pos in range(len(a)):
        # Quality byte: use BC_MAX_QV when no quality provided.
        # Source: corrector.rs  line 127
        qv = qual[pos] if qual else BC_MAX_QV
        p_err = phred_error_prob(qv)

        existing = a[pos]
        for alt in BASE_OPTS:
            if alt == existing:
                continue
            a[pos] = alt
            candidate = bytes(a)

            if whitelist.contains(candidate):
                # Laplace (additive) smoothing.
                # Source: corrector.rs  lines 137-141
                raw_count = bc_counts.get(candidate, 0) if bc_counts else 0
                bc_count = 1 + raw_count
                likelihood = p_err * bc_count
                total_likelihood += likelihood
                if likelihood > best_likelihood:
                    best_likelihood = likelihood
                    best_bc = candidate

        a[pos] = existing  # restore before next position

    if best_bc is None or total_likelihood == 0.0:
        return None

    # Accept only if posterior >= threshold.
    # Source: corrector.rs  lines 153-163
    if best_likelihood / total_likelihood >= confidence_threshold:
        return best_bc
    return None
