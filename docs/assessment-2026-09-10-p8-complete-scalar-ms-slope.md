# P8: complete scalar quadratic MS slope

Original P8 remains OPEN. This follows the
[scalar insertion checkpoint](assessment-2026-09-10-p8-scalar-insertion-ms-reference.md).
No user intervention is required for the current research.

## Complete scalar quadratic reference

[S6.154](../problems/P8/s6/continuation/s6_154/FORMULATION.md)
bounds the finite outer-MS slope of all thirty-two scalar
quadratic refinements while keeping the inner physical OS
conditions and interaction-MS references fixed.

The local sunset slope has the regulated simplex form

    L^2 e^(2 gamma epsilon) mu^(4 epsilon) Gamma(2 epsilon)
      integral w U^epsilon (1-sP)^(-2 epsilon) / (6 Q^2),
    U=xy+xz+yz, P=xyz/U, w=xyz/U^3.

Its residue is L^2/(24 Q^2 epsilon). The finite slope at
s=1 is L^2[ell/2+T1/2+C]/(6 Q^2), with T0=1/2,
|T1|<6 and 0<C<=1/16 proved by six explicit sectors.
No numerically suggested closed form is assumed.

The remaining mixed, repeated and distinct heavy terms,
nested insertions and proper interaction-reference changes
are retained. Their integrated analytic radius-two bounds
give the finite slope and the additional unit-disc OS
remainder. The complete scalar finite-MS-slope bound is
approximately 3.34630300149847e-19, below 1e-18.
The additional interaction-MS OS remainder coefficient is
below 1.14978095631541e-609. This is the complete scalar
reference, not yet the complete GY14 field normalization.

## Independent verification

The native report pins 18 source/proof/test files and
20 fields: 55 named identities, 55 scalar entries,
23 proof gates, 9 controls and 77 rejected inputs.
Final private science passed 194 tests in 12.48 seconds;
fresh repository science passed 194 in 12.59 seconds.
Ordinary replay passed 219 tests in 2209.10 seconds.
Independent native CLI replay passed. The full P8 snapshot
passed 24211 tests in 3374.33 seconds with final exit code 0.

All 561 captured test files were present and unchanged.
Path-list SHA-256:

    6fc4df21799fb83d09d7262cb46e25cfb9e833cb94500070018e213b0bd2b0b2

The full-run adapter passed 128 original tuple comparisons.
Counters were 34081 domain fallbacks, 7340 exact descents
and 88 mixed fallbacks. This snapshot predates S6.155.
Only full regression uses the separately audited exact
GCD adapter; native, ordinary, CLI and direct science
retain unmodified SymPy with interpreter-only allowances.

The 40816-character native report was transferred losslessly,
and its eighteen source hashes were independently verified.
Report SHA-256:

    67cd5c13683398935e985722d2018966d3bfcd25ca0f7657d949fdf619c08a0e

Independent checks include regulated sunset differentiation,
sector integration, Laurent extraction, explicit heavy-mass
derivatives and complex-disc OS remainder controls.
All cleanup preceded final private verification and freeze.
No frozen scientific source or report changed. Written
analytic arguments and exact replay are not formalization
or independent peer review. Exact 22-file staging excludes
later continuations and unrelated P4/P9 changes.

## Remaining work

The first finite-field covariance and complete two-loop Phi
normalization/unit-disc-pole checkpoints have passed direct
science and native checks and are completing fresh regressions.
The matched two-loop four-point assembly is in private
development. None of these successors is counted as published
here.

The vacuum/source assembly, physical higher-order truncation,
V contours/cuts, finite-gravity G and common-parent B remain
open. Finite-order amplitude and pole estimates alone cannot
close those obligations. No all-orders UV construction is added
to the adopted finite-EFT/necessary-positivity contract.
Scoped P8(a) and A.20-A.23 are unchanged; original P8(b) and
P8 remain OPEN.
