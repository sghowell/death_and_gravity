# P8: curved quadratic vector mass-insertion pole

Original P8 remains open. This follows the
[selected-state clock-profile audit](assessment-2026-09-08-p8-selected-state-clock-profile.md).
The selected action, vector state and original clock profiles are unchanged.

## New result

[S6.61](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/FORMULATION.md)
derives the four-dimensional quadratic retained-vector mass-insertion
UV residue, including curved two- and four-derivative terms.

The four-derivative mass dependence is obtained independently from
the scalar determinant on the auxiliary metric. Its full curvature
basis is retained. The component calculation uses independently
constructed coordinate curvatures and ordered frame derivatives.
All three flat coefficients replay the frozen Feynman calculation;
the curved scalar-mass conformal action, constant anisotropic mass
scaling and the independent two-derivative tensor basis also agree.

A literal fourth-order reference transcription fails its independent
curved scalar-mass check, with exact nonzero fixture 237/20. It is
retained as a withheld control, never used as the accepted pole.
The source audit does not assert an erratum or identify the ultimate
cause of that discrepancy.

The actual lapse mass profiles are inserted with all time derivatives
retained. Exact integration by parts gives explicit coefficients of
n''^2, n'^2 and n^2 on the original clock, in the stated Euclidean
local convention. Lorentzian conversion is a separate required step.

The existing full flat potential already includes the zero-derivative
quadratic term; it is checked, not added twice. The new derivative
operators have zero value and first variation at zero mass deviation,
so they do not change the selected background cancellation.

## Remaining work

The dimensional counterterm continuation and its evanescent part
are now separately constructed and undergoing their own regression.
A direct physical-signature continuation is also being checked.
Neither unfinished continuation is included in this frozen claim.

The finite bare/contact/retarded kernel, compatible varied-state
data and integrated feedback bounds remain to be calculated.
Coupled quantum stability, cones, interactions, cutoff and V/G/B
remain open. There is no user-intervention blocker.

## Verification

The report pins 17 sources and fully rebuilds S6.60. It verifies
63 named identities containing 153 scalar entries, 16 algebraic
audit checks and 170 rejected inputs. The focused science suite
passes 17 tests in 57.87 seconds. A test-construction type mismatch
was fixed before the final run and freeze.

The ordinary suite passes 48 tests in 512.66 seconds without the
broad GCD adapter. The independent read-only CLI passes. The full
regression through S6.61 passes 5246 tests in 1174.30 seconds,
with no frozen checkpoint excluded. The adapter passes 128 original
tuple comparisons and records 6472 domain fallbacks and 6509 exact
descents. This run collected before the S6.62 tests were added.

Report SHA-256:
6ac2f5d493737934590e3aa64786549419197ba8b74365660ae8bcbf09abd53a.
The evidence is exact symbolic/rational computation and written
local geometric arguments, not proof-assistant formalization.
