# P8: finite outer MS conversion for the heavy-vertex family

Original P8 remains OPEN. This follows the
[complete quartic chord primitives](assessment-2026-09-10-p8-fermion-complete-quartic-chords.md).
No user intervention is needed for the current research.

## The missing outer reference is now bounded

[S6.143](../problems/P8/s6/continuation/s6_143/FORMULATION.md)
completes the finite outer MS-reference conversion for the full
S6.138 heavy-vertex fermion-insertion family. It retains the complete
dimensionally regulated inner physical on-shell subtraction before
performing the outer subtraction.

The regulated spectral normalization is derived with the same
four-component Dirac trace, d=4-2epsilon and loop measure as the
inherited scheme. Its leading outer reference has Laurent coefficients

    1/(2epsilon^2) -7/(6epsilon) +47/18 +pi^2/12.

A separate direct radial integral at epsilon=1/2 checks the
normalization independently. The exact nonzero mass-ratio correction
has a nonzero simple pole: that pole and the finite prefactor products
are retained. The remaining compact finite integral is bounded by

    |delta F|/(C/Q) <=r(20-2log r), r=1/(4mF^2), C=2NY/Q.

The endpoint estimates are uniform in the declared domain. An exact
dyadic logarithm bound avoids floating-point fits. A literal local
parent parameter map fixes the finite interaction conversion and its
forward coefficient, including the heavy-pole dependence. Locality
in the parent is not confused with a constant four-light amplitude.

At the shared reference, |F_MS| is bounded by approximately
1.49134255240484e-208. The conversion changes b2 by less than
3.41341221406136e-413 relative to the tree coefficient.
Adding it once to the earlier full-family error leaves the family
bound below 1e-622 in the common MS interaction scheme.

## Exact scope

Only scalar_Phi4_W2_F0 advances to the common-MS paired-family status.
Five other primitive rows remain wholly unevaluated at this checkpoint.
The two complete quartic chord primitives and older quadratic paired
sector retain their existing scopes. Other finite field and parameter
conversions, first-order-shift squares and canonical contributions
remain separately owned. This is not full two-loop matching.

## Independent verification

The report pins 18 source/proof/test files and 20 fields: 54 named
identities, 54 scalar entries, 28 proof gates, 9 controls and
53 rejected inputs. Final private science passed 177 tests in
24.03 seconds; fresh repository science passed 177 in 24.06 seconds.
Ordinary replay passed 202 tests in 2610.04 seconds.
Independent native command-line replay passed.
The complete P8 snapshot passed 21660 tests in 3277.18 seconds.

All 539 captured test files were present and unchanged. Path-list SHA-256:

    b0c12f9957c62f6a976bdb7d2e2a4bbba19731b5d41a19a0cf6cb0a847ef923b

The adapter passed 128 original tuple comparisons. Full-run counters:
34024 domain fallbacks, 7340 exact descents and 88 mixed fallbacks.
The snapshot predates S6.144. Only full regression uses the audited
exact GCD adapter. Native, ordinary, CLI and direct science retain
unmodified SymPy with interpreter-only allowances.

The 33009-character native report was transferred losslessly, and
all eighteen source hashes were independently verified. Report SHA-256:

    4ec6b42374a9d0bf3832d32389cd73dcb3a9e0de4ce1cdd1fbdd2a52f7215026

Independent tests include complex-regulator Laurent extraction, the
nonzero mass-ratio pole and finite part, the half-dimensional direct
radial integral and frontier mutations. The analytic arguments are
written proofs, not formalization or independent peer review.
Exact 22-file staging excludes later continuations and unrelated P4/P9 work.

## Continuing work and nonclosure

The remaining mixed quartic primitive has passed final private and
repository science and native ancestry replay. Its ordinary, CLI and
full independent replays are running; it is not yet counted as a
published independent checkpoint here. A private quadratic calculation
now accounts for both overlapping Yukawa vertex subtractions and the
local whole-fermion-cycle term. Finite MS quadratic references remain
distinct from that proposed nonlocal on-shell bound.

Remaining primitive/local-reference work, matching and canonical
contributions, the full two-loop pole/error and truncation estimates
remain open. So do V contours, finite-gravity G and common-parent B.
No all-orders UV construction is added to the adopted finite-EFT and
necessary-positivity contract. Scoped P8(a) and A.20-A.23 are unchanged.
Original P8(b) and original P8 remain OPEN.
