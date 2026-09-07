# P8: what the next formal order changes

This continues the [source-preserving reduction audit](assessment-2026-09-07-p8-source-preserving-reduction.md).
Two separate checkpoints sharpen that result. They do not close original
P8 or establish a controlled match to its original CD/M1 action.

## General ordinary interactions do not repair the retained order

[S6.31](../problems/P8/s6/matching/variable/reduction/general/FORMULATION.md)
extends the calculation to all five smooth ordinary two-metric potential
coefficients and unequal positive constant Einstein weights. The branch
must be a positive simple proportional root of the hidden metric's own
algebraic equation. It is not assumed to be a vacuum or a rolling solution.

The full ten-component Hessian remains invertible precisely when the
specified root polynomial is simple. Completing the hidden Einstein
source again produces the variable conformal curvature-square action.
After all coefficient derivatives and integrations by parts are retained,
its specified scalar-tensor normal form has A1=2F2_X and A3=0, with a
separate curvature-square operator. Choosing ordinary beta2 or beta3
therefore does not by itself repair the retained-order CD mismatch.

But this universality is not an all-order no-go. An explicit smooth
deformation preserves the root, its Hessian, the entire action through
degree four, and the first hidden-metric correction, while changing the
degree-six action. A compact conformal variation gives a nonzero bulk
effect. This is genuine action freedom, not merely a change of basis.
The deformation does not claim to preserve the original example's
full rolling background.

The full curvature-linear potential-cubic table is also derived. Its
higher-operator complement must be kept: a coefficient that vanishes
at the center can have nonzero derivatives after integration by parts.
An independent total-derivative example has zero full Euler expression
but a nonzero apparent projected F2_X. This rejects center-first
coefficient extraction as a test of the full higher-order action.

## The actual first omitted term is not simply small, and it cancels

[S6.32](../problems/P8/s6/matching/variable/reduction/sixth/FORMULATION.md)
returns to the unchanged parent and computes its complete first omitted
action on a specified off-shell tensor probe. The hidden Einstein
second variation, potential cubic and both boundary terms are retained.
The resulting quadratic action and full Euler equation are explicit.

At the actual coefficient center, with delta=c-2, all three action
coefficients tend to zero. Nevertheless the coefficient of the second
tensor derivative in the full sixth-order contribution tends to -2M^2.
An exact polynomial gives an individual magnitude greater than
(793/400)M^2 for 0<delta<=1/100. The finite term comes from derivatives
of the inverse-potential coefficient, not a nonzero limiting action
coefficient. Independent curvature, coframe, root-series and canonical-
clock calculations agree, including all physical time-scale factors.

The previous formal order contributes +2M^2 in the same limit. Their
combined center coefficient is exactly

    M^2 [7delta/2-173delta^2/32-143delta^3/32-349delta^4/128].

The cancellation is part of the result. The individual non-small term
does not establish a non-small combined defect, and the cancellation
does not establish convergence of the entire expansion.

This probe is not the actual FLRW tensor background or a full covariant
Xi calculation. The distinction persists at the bounce: the Hubble
rate vanishes, but its derivative does not. The original physical
metric and free matter coupling are not changed by either checkpoint.

## Remaining discriminating work

A controlled reduction must compare the complete higher-operator/source
map on a stated physical domain and bound its omitted effects. The
formal degree-six freedom, isolated coefficient projections and the
single tensor cancellation do not settle that comparison. A full
physical CD dictionary, selected inverse/state, remainder control and
the adopted positivity/V/G requirements remain open.

The next read-only exploration uses the hidden metric's exact equations
on the same off-shell probe to understand the cancellation without
assuming the formal series converges. It is not yet promoted to a
certificate or substituted for an actual complete healthy cosmology.
The scoped P8(a) photon result is unchanged. No new user authority is
needed for the current mathematical checks; original P8 stays open.

## Verification

S6.31 has 14 frozen sources, 69 exact identities and 22 rejected-input
controls. Its 70 ordinary tests passed in 131.31 seconds; the separate
ordinary report replay passed in 114.56 seconds. Its report SHA-256 is
`4ce6af3c87f527cac16ed3a79d893a315e4b80cbdcd986542a696feb28632cb4`.

S6.32 has 14 frozen sources, 70 exact identities and 26 rejected-input
or finite-derivative controls. Its 54 ordinary tests passed in 140.35
seconds; its separate ordinary replay passed in 119.84 seconds. Its
report SHA-256 is
`d392a734f46d8f4c1c39523c7dc3d41a88634ec0e3dfa645365c3b4d14eb344c`.

The first combined regression attempt stopped at collection because
two independent frozen checkpoints share a test-module basename.
No mathematical test ran in that attempt. Importlib collection avoids
the collision without changing frozen files. The opt-in regression
runner now selects it by default while preserving explicit overrides;
all 21 runner tests pass. Its exact algebra adapter is unchanged.

The full-P8 import-isolated regression passed all 3,181 tests in 646.66
seconds. Its unchanged exact-GCD adapter recorded 6,509 checked exact
descents and 4,552 original-domain fallbacks after 128 comparisons with
the untouched reference implementation. The broad run selected
importlib explicitly; the runner's equivalent default argument handling
was checked separately by its 21 passing tests in 0.11 seconds.

Both seeds were fixed at zero and the host's problematic faulthandler
timer plugin was disabled. These are CERTIFIED written-proof/exact-
algebra results, not FORMALIZED or controlled-EFT claims. Frozen
ancestors and unrelated P4/P9 changes remain untouched.
