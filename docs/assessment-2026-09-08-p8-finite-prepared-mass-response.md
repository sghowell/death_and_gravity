# P8: finite homogeneous prepared mass response

Original P8 remains open. The completed scoped P8(a) objective is
unchanged; P8(b) still requires the remaining quantum, cutoff and
common-parent V/G/B work. This follows the
[finite local mass-response audit](assessment-2026-09-08-p8-finite-local-mass-response.md).
There is no user-intervention blocker.

## New result

[S6.66](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/FORMULATION.md)
combines the already-fixed local subtraction with a continuous bound
on the exact selected-state remainder for prepared homogeneous
mass sources. The original clock, all-order Gaussian preparation and
fixed scalar profiles are unchanged.

The source is smooth, spatial K=0, compactly supported in time after
u0=-1/2, and zero on an initial neighborhood. The source norm controls
time derivatives zero through ten on I=[-1/2,1/2].

An eighth-order WKB comparison frequency is differentiated using the
exact linearized Riccati equation. It is an auxiliary reference, not
a new state or new subtraction: the subtraction remains orders 0,2,4.
Continuous rational majorants control the varied residual and the
reference readout tail.

The exact mixing transport retains the selected state's nonzero
initial mixing. Zero initial response does not erase its interference
or varied phase; the squeezed control remains -15/4. The slowest
varied finite tail is inverse momentum to the fourth power and is
integrable in three spatial dimensions.

A written parameter-dependence and dominated-convergence argument
justifies the common finite dimensional limit in a sufficiently
small complex-D neighborhood of three. The diagonal mixing
coefficients stay bounded; the off-diagonal initial coefficients
are O(nu^-6). No uniform numerical claim on the whole older
dimensional strip is made.

At M tau=10^24,m0 tau=1000, the nonlocal subtracted C10-to-C0
operator bound is below 10^-43. Including the S6.65 finite local
piece gives a complete homogeneous mass-block bound below 10^-39.
It is retarded in time on this prepared-source domain.

## Remaining work

The full homogeneous metric local block has separately been derived
and is undergoing its independent regression. Its full state response,
including lapse/log-scale contacts and the existing fixed tadpole
term, is in progress. These are not imported into the frozen S6.66
mass-only claim.

Arbitrary spatial-momentum response, arbitrary initial covariance
variations, a no-loss coupled inverse, quantum stability/cones,
interactions, cutoff, finite Wilson matching and common-parent
V/G/B remain open. Small high-regularity homogeneous norms do not
settle those obligations.

## Verification

The report pins 16 sources and fully rebuilds S6.65. It verifies
40 named identities containing 43 scalar entries, 32 continuous
audit checks and 129 rejected inputs.

The focused science suite passes 14 tests in 344.50 seconds.
The ordinary suite passes 44 tests in 796.19 seconds without the
broad GCD adapter. The independent read-only CLI passes.

The full regression through S6.66 passes 5471 tests in 1506.27
seconds, with no frozen checkpoint excluded. The adapter passes
128 original tuple comparisons and records 9041 domain fallbacks
and 6509 exact descents. This run collected before S6.67 tests
were added.

Report SHA-256:
81a0419a304c689847d2005438080e584b7dba73d42db29780127da8c937d387.

The evidence is exact symbolic/rational computation with written
continuous and dimensional arguments, not proof-assistant formalization.
