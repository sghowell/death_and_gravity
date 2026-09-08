# P8: homogeneous finite local mass-response subtraction

Original P8 remains open. This follows the
[bare retarded-response audit](assessment-2026-09-08-p8-retarded-mass-response.md).
The original clock, selected vector state and fixed scalar profiles
remain unchanged.

## New result

[S6.65](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/FORMULATION.md)
computes the local adiabatic subtraction for the homogeneous prepared
mass-response sector, including its finite dimensional counterterm part.

Independent physical D-dimensional Hamiltonians give the moving-mode
frequency, normalization and readout variations. The differentiated
Riccati expansion retains all mass-profile and source time derivatives
through adiabatic order four. Radial integration keeps D-1 transverse
polarizations and the first dimensional jet.

All three resulting homogeneous curved poles reproduce S6.63 exactly.
Combining the radial finite part at mu=m with the specified evanescent
counterterm gives a self-adjoint local response. The zero-derivative
part independently matches the already-frozen full finite potential;
it is not added twice.

The compact finite coefficients of n''^2,n'^2,n^2 are explicit.
Omitting the finite counterterm fails the self-adjointness controls,
with nonzero half-time fixtures -218103808/512578125 and
-67108864/512578125 at unit algebraic reference mass.

On [-1/2,1/2], the normalized local operator has a continuous
C4-to-C0 source norm below 10^-39 at M tau=10^24,m0 tau=1000.
This is a coefficient-wise rational interval proof, not point sampling
or a coupled feedback contraction.

## Remaining work

The exact nonlocal subtracted homogeneous response now has its own
derived bounds and is entering a separate checkpoint verification.
That result is not imported into this frozen local-only claim.

Spatially varying finite response, full metric and second mass-source
blocks, arbitrary initial covariance variations, no-loss coupled
feedback, quantum stability/cones, interactions, cutoff, finite Wilson
matching and V/G/B remain open. There is no user-intervention blocker.

## Verification

The report pins 16 sources and fully rebuilds S6.64. It verifies
43 named identities containing 43 scalar entries, 14 continuous audit
checks and 92 rejected inputs. The focused science suite passes
14 tests in 59.21 seconds.

The ordinary suite passes 44 tests in 553.80 seconds without the
broad GCD adapter. The independent read-only CLI passes. The full
regression through S6.65 passes 5427 tests in 1193.22 seconds, with
no frozen checkpoint excluded. The adapter passes 128 original tuple
comparisons and records 6833 domain fallbacks and 6509 exact descents.
This run collected before the S6.66 tests were added.

Report SHA-256:
1a6a16fcfd54827aefdfbe6e5a42f5a90fd001e611ae926c22bca01e29a66c70.
The evidence is exact symbolic/rational computation with written
dimensional and continuous-bound arguments, not proof-assistant
formalization.
