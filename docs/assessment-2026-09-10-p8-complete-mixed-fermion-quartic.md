# P8: the remaining mixed quartic primitive is bounded

Original P8 remains OPEN. This follows the
[finite outer MS conversion](assessment-2026-09-10-p8-fermion-outer-ms-conversion.md).
No user intervention is needed for the current research.

## All four primitive quartic rows now have bounds

[S6.144](../problems/P8/s6/continuation/s6_144/FORMULATION.md)
bounds the complete scalar_Phi4_W1_F2 family. The actual nonlocal
tree Hessian, all six cyclic fermion boxes, both assignments in
each of three channels and their one-half weights are retained.

The proper box is split exactly into a zero-soft kernel and its
full momentum-dependent difference. A joint massive/massive/massless
integral bounds the latter. Only the fermion lines receive the
large soft-scaling Cauchy circle: the physical light and heavy
propagators remain outside that operation.

A direct complex-disc bound on the product of the two light
denominators is used. It does not assume that an estimate for an
individual real-momentum propagator remains valid after continuation.
All chord-diagonal and unequal-momentum regions are integrable.

The zero-soft kernel is a fixed-scale second fermion-mass derivative.
Its signed spectral density is not treated as a positive exact
spectral measure. Explicit bubble differences and heavy-triangle
integrals bound its physical momentum dependence.

## Entire MS reference, not only momentum dependence

The proper quartic MS counterterm is paired before the remaining
outer reference is subtracted. Both regulated soft and hard pieces
are retained. In units CY/Q, their combined leading finite
reference is

    -46-32 log(mF^2), C=2NY/Q.

The finite mass-ratio correction is retained and bounded by
r[300+100 log(1/r)], r=1/(4mF^2). Its difference integral is
convergent before the regulator limit. A literal parent parameter
map supplies the finite outer MS conversion, whose relative b2
coefficient is minus that reference.

The complete mixed primitive error is bounded by approximately
1.84389805055489e-610, or 4.60974512638723e-11 relative to the
tree forward coefficient. Only this row advances. All four
quartic primitive rows now have bounds, while four vacuum/quadratic
rows remain wholly unevaluated at this checkpoint. Other matching
and canonical contributions are not assigned zero.

## Independent verification

The report pins 18 source/proof/test files and 20 fields:
78 named identities, 78 scalar entries, 25 proof gates, 9 controls
and 90 rejected inputs. Final private science passed 285 tests
in 24.12 seconds; fresh repository science passed 285 in
24.13 seconds. Ordinary replay passed 310 tests in 2514.79 seconds.
Independent native command-line replay passed.
The complete P8 snapshot passed 21970 tests in 3246.76 seconds.

All 541 captured test files were present and unchanged. Path-list SHA-256:

    1c874594b12763206b03faeea542dc788360fe05f8abbcb86ff9f66953d35c06

The adapter passed 128 original tuple comparisons. Full-run counters:
34022 domain fallbacks, 7340 exact descents and 88 mixed fallbacks.
The snapshot predates S6.145. Only full regression uses the
audited exact GCD adapter. Native, ordinary, CLI and direct science
retain unmodified SymPy with interpreter-only allowances.

The 44061-character native report was transferred losslessly.
All eighteen source hashes were independently checked. Report SHA-256:

    32208fdf95cbaa1bfdef6e3d43f84450b0b2d41e9ec043fc544b6eeca3e5ea0d

Independent tests include Schwinger integrals, explicit complex
Dirac matrices, fixed-scale mass derivatives, dimensional finite-part
extraction, both heavy radial regions and frontier mutations.
The proofs are written arguments with exact algebra, not
formalization or independent peer review. Exact 22-file staging
excludes later continuations and unrelated P4/P9 changes.

## Remaining work and nonclosure

The two quadratic nonlocal primitive remainders have passed direct
science and native replay; their fresh independent replays are
running. A private finite-MS-slope calculation has passed 220
science tests and awaits finalization. These later calculations
are not counted as published independent results here.

Finite mass references, vacuum rows, other parameter/counterterm
and canonical conversions, the full matched two-loop pole/error
and truncation control remain open. So do V contours, finite-gravity
G and common-parent B. No all-orders UV construction is added
to the adopted finite-EFT/necessary-positivity contract.
Scoped P8(a) and A.20-A.23 are unchanged.
Original P8(b) and original P8 remain OPEN.
