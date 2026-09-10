# P8: complete scalar insertion MS reference

Original P8 remains OPEN. This follows the
[double-bubble checkpoint](assessment-2026-09-10-p8-double-bubble-ms-conversion.md).
No user intervention is required for the current research.

## All sixty-four insertion refinements

[S6.153](../problems/P8/s6/continuation/s6_153/FORMULATION.md)
converts the remaining outer interaction reference of the
complete scalar on-shell insertion family to MS at mu=mF,
holding the physical inner OS condition and canonical
reference coordinates fixed.

The full dimensional inner slope fixes both alpha0 and
alpha1 before multiplication by the outer logarithmic pole.
Its finite outer reference is

    Falpha=(g/Q^2) integral b [2ell-log D],
    D=xM+(1-x)^2, b=x(1-x)/D.

Using only alpha0 ell/Q would miss alpha1/Q. Both insertion
positions and both convergent decaying kernels remain in
the same sixty-four-refinement family.

The complete local conversion is Falpha times the
three-channel sum C^2. Its exact heavy-tree variation gives
the relative forward coefficient
Falpha[-2L+3g/(M-2)]. At the actual parameters, the reference
is bounded by 4.41515887225116e-207. The resulting absolute
b2 conversion is below 8.08439734909270e-1011; adding the old
insertion family gives 6.24500451351651e-615, below 1e-614.

The equivalent linear second-order G, g, M and L shifts
are fixed by this same reference. They are an alternative
representation, not extra contributions or another LSZ factor.

## Four raw scalar interaction families assembled

The 88 UV-finite, 64 insertion, 24 double-bubble and
16 wineglass refinements are disjoint and exhaust all
192 raw scalar interaction refinements. With the last
three families converted to interaction-MS references,
their combined b2 bound is approximately
2.69375534479471e-607, or
6.73438836198678e-8 relative to tree, below 1e-7.

This is a fixed-coordinate interaction-family statement.
It does not yet complete the global field/coupling map or
the full GY14 canonical amplitude. The isolated finite
contact insertion was handled in S6.150 and is not added
again. All nine fermionic primitive rows remain unchanged.

## Independent verification

The report pins 18 source/proof/test files and 20 fields:
38 named identities, 38 scalar entries, 23 proof gates,
9 controls and 90 rejected inputs. Final private science
passed 165 tests in 23.92 seconds; fresh repository science
passed 165 in 22.38 seconds. Ordinary replay passed 190
tests in 2172.15 seconds. Independent native CLI replay
passed. The complete P8 snapshot passed 23992 tests in
3366.78 seconds.

All 559 captured test files were present and unchanged.
Path-list SHA-256:

    3b0e4559c450ca6397cf32293e817985094aea2f7d9e583976beda2c99fbb416

The full-run adapter passed 128 original tuple comparisons.
Counters were 34089 domain fallbacks, 7340 exact descents
and 88 mixed fallbacks. This snapshot predates S6.154.
Only full regression uses the separately audited exact
GCD adapter; native, ordinary, CLI and direct science
retain unmodified SymPy with interpreter-only allowances.

The 47959-character native report was transferred losslessly,
and its eighteen source hashes were independently verified.
Report SHA-256:

    46311af32da0d789a75132b6c9515c9e34f04966a4422deb794e31bd35e34a0d

Independent tests differentiate the nonzero-regulator
mixed bubble at fixed scale, extract the finite pole product
on a complex circle, check large-momentum OS asymptotics,
and independently recover the forward coefficient and
literal G/M coordinate shifts. Negative controls detect
early alpha truncation and inappropriate scope closure.
All cleanup preceded final private verification and freeze.

No frozen scientific source or report changed. Written
analytic arguments and exact tests are not formalization
or independent peer review. Exact 22-file staging excludes
later continuations and unrelated P4/P9 changes.

## Remaining work

The scalar finite-MS-slope and first finite-field covariance
checkpoints have passed direct/native checks and are completing
fresh regressions. Reconstruction of the complete two-loop
Phi normalization and unit-disc canonical pole is in private
development. None is counted as published here.

The remaining matched amplitude and vacuum/source assembly,
physical higher-order truncation, V contours/cuts,
finite-gravity G and common-parent B remain open. The field
and coupling correspondence must still be fully accounted
for; primitive and family bounds alone do not close it.
No all-orders UV construction is added to the adopted
finite-EFT/necessary-positivity contract. Scoped P8(a) and
A.20-A.23 are unchanged; original P8(b) and P8 remain OPEN.
