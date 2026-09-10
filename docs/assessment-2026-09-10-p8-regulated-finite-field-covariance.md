# P8: regulated first finite-field covariance

Original P8 remains OPEN. This follows the
[scalar MS-slope checkpoint](assessment-2026-09-10-p8-complete-scalar-ms-slope.md).
No user intervention is required for the current research.

## Full first Phi normalization and its regulator coefficient

[S6.155](../problems/P8/s6/continuation/s6_155/FORMULATION.md)
derives the complete first normalization
k_D=r_D-fprime_MS,D through its positive-epsilon coefficient.
The full fermion trace reduces by a total parameter derivative,
with vanishing endpoints, before the finite pole products
are extracted. At the actual parameters,

    |k0| < 3.84413165810667438e-207,
    |k1| < 3.28629105179398148e-205.

The complete scalar and fermion proper MS poles are retained.
A fully regulated Phi4 graph has net field weight K^-2:
its vertex and internal-line factors cancel to that external
weight. The correlated kinetic and mass insertion preserves
pE^2+1, not just pE^2.

The positive-epsilon coefficient times the raw one-loop
pole cancels only against the paired counterterm product.
The inherited canonical first-order amplitude is not
multiplied by another first-order LSZ factor.

## Fixed-input coefficients versus complete bare matching

For a total vertex coefficient X of Phi weight w, this
checkpoint gives the fixed-input finite transformed
coefficient

    Xc2=[w(w+1)k0^2/2-w t0]X-w k1 p1.

This is not yet the full second bare-parameter map.
The first pole p1 must also be re-expressed at the full
regulated shifted coordinates. S6.156 performs that step;
its results are not counted as published here.

The first-field part of the raw-MS paired amplitude is
bounded by approximately 3.43692152926923e-213 relative
to tree. The unknown full second normalization remains
symbolic in this checkpoint. Fundamental cubic and
Yukawa coupling squares are retained.

## Independent verification

The report pins 18 source/proof/test files and 20 fields:
62 named identities, 62 scalar entries, 23 proof gates,
9 controls and 104 rejected inputs. Final private science
passed 210 tests in 21.63 seconds; fresh repository science
passed 210 in 21.84 seconds. Ordinary replay passed 235
tests in 2220.37 seconds. Independent native CLI replay
passed. The full P8 snapshot passed 24446 tests in
3317.81 seconds with final exit code 0.

All 563 captured test files were present and unchanged.
Path-list SHA-256:

    795eee7b16e6d18a9edf1b0a465a94b9db27016a57ec96b0e68392417611e92a

The full-run adapter passed 128 original tuple comparisons.
Counters were 34087 domain fallbacks, 7340 exact descents
and 88 mixed fallbacks. This snapshot predates S6.156.
Only full regression uses the separately audited exact
GCD adapter; native, ordinary, CLI and direct science
retain unmodified SymPy with interpreter-only allowances.

The 59166-character native report was transferred losslessly,
and all eighteen source hashes were independently verified.
Report SHA-256:

    76f6fac15d1f81af93b428f2d074ea12d12f5bbfcf08a6414a9dd2035822ed82

Independent checks include the dimension-dependent fermion
trace identity, regulated coupling products, finite-field
vertex/line covariance and negative controls against early
epsilon truncation or excess scope closure. All cleanup
preceded final private verification and freeze.

No frozen scientific source or report changed. Written
analytic arguments and exact replay are not formalization
or independent peer review. Exact 22-file staging excludes
later continuations and unrelated P4/P9 changes.

## Remaining work

The complete second Phi normalization/unit-disc pole and
matched through-two-loop forward-coefficient checkpoints
have passed direct science and native checks and are
completing fresh regressions. Neither is counted as
published here. The complete vacuum-reference assembly
and regulated first-source-square ownership are in private
development.

The complete source reference, physical finite-EFT
truncation, V contours/cuts, finite-gravity G and
common-parent B remain open. A finite-order positive
flat-vacuum coefficient does not discharge those
obligations. No all-orders UV construction is added
to the adopted finite-EFT/necessary-positivity contract.
Scoped P8(a) and A.20-A.23 are unchanged; original P8(b)
and original P8 remain OPEN.
