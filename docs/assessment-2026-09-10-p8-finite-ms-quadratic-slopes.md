# P8: finite MS slopes of both quadratic primitives

Original P8 remains OPEN. This follows the
[quadratic nonlocal forest checkpoint](assessment-2026-09-10-p8-quadratic-nonlocal-forests.md).
No user intervention is required for the current research.

## Finite local slopes, not just subtracted remainders

[S6.146](../problems/P8/s6/continuation/s6_146/FORMULATION.md)
computes the complete dimension-symbolic quadratic primitive
slope tensors. Both self-energy placements, the vertex word and
its cross derivative are retained. Gaussian massive/massive/massless
masters reduce every term exactly before d=4 is taken.

The proper fermion mass, kinetic and Yukawa MS counterterms are
included once. Fixed-scale differentiation precedes mu=m. After
overall MS pole subtraction, the finite Euclidean slopes are

    49/12 in NY^2/Q^2 units for scalar exchange,
    -37/6 in NYaC_F/Q^2 units for gauge exchange.

The Minkowski slopes have the opposite sign. Finite epsilon-times-pole
products cancel the pi-squared contributions exactly. The scalar
massless exchange is an auxiliary reference only. Its replacement
by the actual mass-one scalar is a fully convergent difference,
bounded jointly over both loop momenta and the chord diagonal.
The previous soft-tail disc bounds the shift to the on-shell slope.

The combined actual absolute on-shell slope is below 1e-410.
Its massless-reference contribution is bounded by approximately
2.82923410520425e-413; the actual-scalar correction and the
zero-to-on-shell allowance are both retained. Only the two
quadratic primitive slope statuses advance. Finite mass references
and full canonical matching are not inferred from this result.

## Independent verification

The report pins 18 source/proof/test files and 20 fields:
44 named identities, 44 scalar entries, 24 proof gates,
9 controls and 77 rejected inputs. Final private science passed
220 tests in 22.64 seconds; fresh repository science passed
220 in 22.88 seconds. Ordinary replay passed 245 tests in
2445.63 seconds. Independent native command-line replay passed.
The complete P8 snapshot passed 22455 tests in 3211.30 seconds.

All 545 captured test files were present and unchanged. Path-list SHA-256:

    c0abfd3649371bb22fe52173c0861e3f077b3094952b87e0bb7cc4bada679c7b

The full-run adapter passed 128 original tuple comparisons.
Counters were 34019 domain fallbacks, 7340 exact descents and
88 mixed fallbacks. This snapshot predates S6.147. Only full
regression uses the separately audited exact GCD adapter;
native, ordinary, CLI and direct science use unmodified SymPy
with interpreter-only allowances.

The 29188-character native report was transferred losslessly,
and all eighteen source hashes were independently verified.
Report SHA-256:

    8df7827728af867a091e4532574132d530f0e3815f4473cb9404e918c27ae202

Independent checks include explicit Dirac-matrix derivatives in
dimensions three, four and five, convergent Schwinger integrals,
complex-regulator finite-part extraction and exact fixed-scale
vacuum derivatives. Private numerical tests use a precision
tolerance for numerical zeros, while analytic zero identities
remain exact. All fixes preceded the final private run and freeze.
No frozen scientific source or report changed.
Written analytic proofs and exact tests are not formalization or
independent peer review. Exact 22-file staging excludes later
continuations and unrelated P4/P9 changes.

## Remaining work

The subsequent direct-primitive mass and older-insertion local
reference checkpoints have passed their private, repository and
native checks. Their independent replays are running and are not
counted as published results here. The vacuum-row calculation is
still private research.

Vacuum references, other parameter/counterterm and canonical
conversions, the complete matched two-loop pole/error, truncation
control, V contours, finite-gravity G and common-parent B remain
open. No all-orders UV construction is added to the adopted
finite-EFT/necessary-positivity scope. Scoped P8(a) and A.20-A.23
are unchanged. Original P8(b) and original P8 remain OPEN.
