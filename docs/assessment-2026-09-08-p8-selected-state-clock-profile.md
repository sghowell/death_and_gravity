# P8: selected-state clock-preserving scalar-profile candidate

Original P8 remains open. This follows the
[fifth-order vector regularity audit](assessment-2026-09-08-p8-fifth-order-vector-regularity.md).
This checkpoint introduces a new named lower-scalar action; it does not
revise the frozen predecessor action, vector state or finite prescription.

## New result

[S6.60](../problems/P8/s6/matching/affine/kinetic/aligned/margin/tadpole/FORMULATION.md)
defines two fixed smooth scalar coefficient profiles from the actual
S6.55 vector Cauchy problem, subtracted physical stress integrals and
S6.53 finite local terms on the original complete clock history.
A compact clock-norm window localizes the new scalar addition away
from an open zero-gradient region.

The profiles are fixed before variation. They cancel the selected
vector energy, pressure and Ward clock source exactly on that history.
Thus the unchanged complete clock geometry solves the reduced
retained-vector Gaussian semiclassical background equations for
this new, deliberately state-tuned action. This is a physical action
change, not a state-dependent renormalization rule. Another state
or metric does not cause the profiles to be recomputed.

The construction is motivated by an exact off-clock vector-block
check: the particular earlier fixed-source approximate initial
metric has a superluminal retained longitudinal block. This is
not a no-go for other initial data or a verdict on a full quantum cone.

Global smooth definition is proved on every compact real-time
interval. Quantitative smallness is restricted to the bounce window
|u|<=1/2. At M*tau=10^24,m0*tau=1000, every mixed (u,x) derivative
of the normalized scalar addition through total order five is
below 10^-18, for all real x. The literal lapse-square change is
bounded by 15eta0/8. The exact physical point-chart quadratic
terms are retained. The window changes no local field jet near
x=0; this does not prove that a healthy vacuum exists there.

## Remaining work

The selected first variation cancels, but the contact, retarded
and compatible initial-covariance second variation remains.
The working continuation is testing the curved quadratic
mass-insertion pole needed for that calculation against independent
flat and conformal-mass checks. A failing curved check is not
imported as a certified formula.

Coupled quantum stability and cones, controlled interactions and
cutoff, omitted sectors and higher loops, and the common-parent
V/G/B obligations remain open. There is no user-intervention blocker.

## Verification

The report pins 14 sources and fully rebuilds S6.59. It verifies
30 named identities with 30 scalar entries, 50 continuous proof
checks and 110 rejected inputs. Focused science passes 14 tests
in 139.36 seconds. The ordinary suite passes 50 tests in
514.74 seconds without the broad GCD adapter; the independent
read-only CLI passes.

The full regression through S6.60 passes 5198 tests in
1075.95 seconds without excluding a frozen checkpoint. The
adapter passes 128 original tuple comparisons and records
6458 domain fallbacks and 6509 exact descents.

Report SHA-256:
dba6c1f44957adaa54e15557bbf14eac6eeee5f911524caea301315d91425241.
The evidence is exact symbolic/rational computation and written
continuous proofs, not proof-assistant formalization.
