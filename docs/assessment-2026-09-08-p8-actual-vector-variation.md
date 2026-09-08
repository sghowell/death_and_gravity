# P8: actual vector response and initial-state compatibility

Original P8 remains open. This follows the
[fixed-vector-stress response audit](assessment-2026-09-08-p8-fixed-vector-stress-response.md).
The named action, fixed-background state, finite prescription and
adopted closure contract are unchanged.

## New result

[S6.58](../problems/P8/s6/matching/affine/kinetic/aligned/margin/response/variation/FORMULATION.md)
derives the actual vector mode and observable variations on
homogeneous physical metric histories. It eliminates the temporal
constraint before normalizing the two transverse and one longitudinal
modes. In physical canonical variables, the causal first-order source
depends on the metric variations themselves. The equivalent oscillator
equation retains the differentiated canonical map and initial data.

Physical lapse and scale variations reproduce the previously certified
energy and pressure weights. Their local contact matrices retain both
second lapse derivatives of the actual mass functions. The covariance
response includes its initial term; the contact and retarded momentum
integrals are not separately declared finite.

There is a concrete initial-state obstruction to the simplest feedback
prescription. At L=10^12,R=1000 the actual initial force exceeds
5*10^-15 and forces a nonzero negative initial lapse perturbation.
On its original physical chart lift, the longitudinal leading canonical
impedance ratio obeys 0<r_L^2<1. Copying the old physical coordinate and
momentum covariance therefore leaves a nonzero high-frequency
negative-frequency coefficient, with a quartically divergent excitation
tail. This fails the necessary leading adiabatic preparation condition.

That result does not rule out a compatible state family. The transverse
physical impedance has a leading conformal cancellation, which is not
an all-order state proof. Freezing normalized oscillator data is a
different prescription and already fails the transverse leading test.
Off clock, the longitudinal principal operator is not silently identified
with ordinary Proca on the physical metric.

## Remaining work

The next step is a compatible all-order covariance family, with explicit
regularity and state-change bounds, followed by the combined renormalized
contact/retarded second variation. Higher time-derivative estimates of the
existing fixed-background state are now in progress. No numerical
feedback norm, self-consistent quantum solution, corrected coupled cones,
interacting cutoff or V/G/B verdict is claimed. There is no
user-intervention blocker.

## Verification

The report pins 13 sources and fully rebuilds S6.57. It verifies
31 named identities comprising 58 scalar entries, 17 continuous proof
checks and 22 rejected inputs. Focused science passes 9 tests in
2.60 seconds. The ordinary suite passes 47 tests in 436.79 seconds
without the broad GCD adapter; the independent read-only CLI passes.

The full regression through S6.58 passes 5098 tests in 973.21 seconds
without excluding a frozen checkpoint. The adapter passes 128 original
tuple comparisons and records 5879 domain fallbacks and 6509 exact
descents. This run collected before any S6.59 tests existed.

Report SHA-256:
e4aa8aa397c2722168d7f475cd6b6fb73d742d3e143596d344955c5bd1feeb99.
The evidence is exact symbolic/rational verification with written
continuous proofs, not proof-assistant formalization.
