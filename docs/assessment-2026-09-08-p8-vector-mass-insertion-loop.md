# P8: vector derivative pole and finite subtracted loop

Original P8 remains open. The scoped P8(a) photon objective, linear
CD/M1 classification, physical matter frame and adopted V/G/B contract
are unchanged. This follows the [local quantum checkpoint](assessment-2026-09-07-p8-local-vector-quantum.md).

## New result

[S6.48](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/bubble/FORMULATION.md)
calculates the retained vector's two-mass-insertion quantum component.
It keeps the full Proca numerator and the actual first mass variation,
with the physical metric held fixed. This external mass variation is
not the free physical lapse of the gravitational system.

The exact pole contains nonzero terms of momentum degree two and
four. The fourth-degree weight satisfies (k^2)^2/125<F4<(k^2)^2/25
for nonzero real Euclidean k on the compact clock interval. A direct
shifted tensor integral, an independent isotropic denominator reduction
and the previous potential's linear-mass second variation agree.
The second variation of the original mass tensor supplies a separate
tadpole; the bubble alone is not the full scalar two-point function.

The finite result is also derived, rather than inferred from the pole.
Actual radial Laurent coefficients fix the logarithmic form factor.
After subtracting its Taylor polynomial through momentum degree four,
the vector loop remainder on the explicit complex ball
kappa=sum_mu |k_mu|^2<=m0^2 obeys

    |R_loop(k)| <= kappa^3/(1920*pi^2*h^2*m0^2),
    h=(1+u^2)^3, |u|<=1/2.

An exact Feynman-parameter envelope gives 4852/229635<1/40;
convergent logarithm tails prove the bound over the whole domain.
An independent isotropic comparison verifies every Taylor order
above the subtraction order, not only selected momentum samples.

With Lp=M*tau, R=m0*tau and Q=kappa*tau^2<=R^2, the bound
relative to the reference density M^2/tau^2 is
Q^3/(17280*Lp^2*R^2). Lp=10^12, R=1000, Q=1 gives
1/(17280*10^30). This is a quadratic Fourier-kernel bound: a
band-limited source-functional estimate additionally carries its
squared Fourier/L2 norm. It is not a pointwise energy correction.

## What remains

This subtraction removes local scheme/scale terms through fourth
momentum order. It does not remove or bound unknown higher local
matching operators. No counterterm or frozen action is changed.
The complex ball supports stationary analytic continuation of this
component, not transfer to a curved, evolving, in-in state.

The next checkpoint derives the vector's on-clock local curvature
counterterms on the actual physical FLRW geometry, with an independent
spectral control. A finite state-dependent curved remainder, full
mass/metric variations, other field loops, quantum-corrected constraints
and cones, interacting cutoff and V/G/B estimates still require work.
The original exactly luminal matter-cone condition is not protected
by a small unsigned error alone. No user intervention is needed now.

## Verification

The report pins 13 local sources and fully rebuilds S6.47 and its
ancestry. It checks 27 named exact identities comprising 42 scalar
entries, 13 continuous/domain proof checks and 49 rejected inputs.
The scientific suite passes **21 tests in 2.35 seconds**; the ordinary
suite passes **52 tests in 225.32 seconds**, without the broad GCD
adapter. The separate seeded read-only CLI passes. Missing, extra
and mutated report fields are rejected.

The full P8 regression through S6.48 passes **4651 tests in 816.47
seconds**, with no frozen checkpoint excluded. The unchanged exact
GCD adapter passes 128 original tuple comparisons and records 5717
domain fallbacks and 6509 exact descents under the per-test seed recipe.
This is exact symbolic verification with written continuous integral
and analytic-domain proofs, not proof-assistant formalization or an
external peer review. Original P8 is not marked finished or closed.
