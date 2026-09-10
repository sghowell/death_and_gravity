# P8: complete double-bubble interaction-forest MS conversion

Original P8 remains OPEN. This follows the
[wineglass checkpoint](assessment-2026-09-10-p8-wineglass-ms-conversion.md).
No user intervention is required for the current research.

## All twenty-four double-bubble refinements

[S6.152](../problems/P8/s6/continuation/s6_152/FORMULATION.md)
converts the full twenty-four-refinement double-bubble family
to MS interaction references at fixed canonical coordinates.
It evaluates the actual restricted forests before removing
the regulator. Shared-cycle subtractions are not replaced
with a forbidden product of overlapping cores.

If I=P+f0+epsilon f1+... and P is the simple MS pole,
the pole projections of I^2 and PI retain their f0-dependent
simple poles. The six legal shared-core forest terms cancel
to (I-P)^2; the four disjoint-core terms do likewise. Mixed
triangle terms retain the full regulated finite triangle.

The resulting per-channel expression is

    [C^3 (IR+d)^2+C^2 (IR+d)(TA+TB)]/4,  d=ell/Q.

Both heavy triangles and all three channels are retained.
The linear scale-conversion b2 bound is approximately
6.67799149312032e-612. The exact C^3 forward coefficient
controls the quadratic scale term, including ordinary
repeated-heavy counterterm insertions. It is not a new
nonlocal counterterm. The complete converted family is
bounded by 6.68910760115438e-612, or
1.67227690028859e-12 relative to tree.

This discharges an assigned scalar interaction-family
conversion. It does not identify the fixed canonical
coordinates with the full boundary-MS field/coupling map.
All nine fermionic primitive rows remain unchanged.

## Independent verification

The report pins 18 source/proof/test files and 20 fields:
65 named identities, 65 scalar entries, 23 proof gates,
9 controls and 94 rejected inputs. Final private science
passed 205 tests in 20.51 seconds; fresh repository science
passed 205 in 20.34 seconds. Ordinary replay passed 230
tests in 2121.31 seconds. Independent native CLI replay
passed. The complete P8 snapshot passed 23802 tests in
3393.65 seconds.

All 557 captured test files were present and unchanged.
Path-list SHA-256:

    639fe8c1bbc5137c3bd1b9ca55b9fa40d1d6a1e28737e43f6872d1fdf4345d20

The full-run adapter passed 128 original tuple comparisons.
Counters were 34071 domain fallbacks, 7340 exact descents
and 88 mixed fallbacks. This snapshot predates S6.153.
Only full regression uses the separately audited exact
GCD adapter; native, ordinary, CLI and direct science
retain unmodified SymPy with interpreter-only allowances.

The 50339-character native report was transferred losslessly,
and its eighteen source hashes were independently verified.
Report SHA-256:

    df058b81fec53d4d6d2ee20c83935a59694bda97dcbabece17ad0ab1fc785505

Tests cover the actual twenty-four forest topologies,
regulator-first pole projections, complex finite-part
extraction, the full triangle factorization, independent
forward contour coefficients and all-radius integrals.
An initial private test exposed unreduced exact-rational
zeros in three factorization assertions; exact factorization
resolved those assertions without changing any formula.
All changes preceded the final private run and freeze.

No frozen scientific source or report changed. Written
analytic arguments and exact replay are not formalization
or independent peer review. Exact 22-file staging excludes
later continuations and unrelated P4/P9 changes.

## Remaining work

The scalar insertion and finite scalar MS slope checkpoints
have passed direct/native checks and are completing fresh
regressions. The first finite-Phi-field covariance and
projection-commutator calculation is in private development.
None is counted as published here.

Global field/coupling matching and cross terms, full GY14
vacuum/source and canonical amplitude/pole assembly,
physical higher-order truncation, V contours/cuts,
finite-gravity G and common-parent B remain open.
No all-orders UV construction is added to the adopted
finite-EFT/necessary-positivity scope. Scoped P8(a) and
A.20-A.23 are unchanged; original P8(b) and P8 remain OPEN.
