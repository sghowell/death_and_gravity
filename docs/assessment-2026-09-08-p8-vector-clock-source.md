# P8: differentiated vector matching and homogeneous clock source

Original P8 remains open. The scoped P8(a) photon objective,
original linear classification and adopted physical-frame V/G/B
contract are unchanged. This follows the
[finite clock-mass matching audit](assessment-2026-09-08-p8-finite-clock-mass-matching.md).

## New result

[S6.54](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/evolution/FORMULATION.md)
keeps the exact prepared Gaussian modes and S6.53's finite
prescription. It derives continuous bounds for the first time
derivatives of the actual vector energy and pressure, then
derives and bounds the homogeneous clock Euler source.

A direct rational calculation gives the residual and its time
derivative envelopes C/omega^4 and C1/omega^4, with
C=2000000,C1=100000000. Removing the diagonal evolution phase
by an exact transformation and integrating the remaining
oscillatory mixing once by parts, including both endpoints,
gives a stronger bound on actual mixing:

    abs(Bcal)<=K/nu^6, K=676000000,
    nu^2=m^2+k^2/(25/16)^2.

The physical quadratic readout is differentiated before momentum
integration. Its phase-independent product identities and this
stronger mixing estimate give an integrated derivative difference
below 40K/m. The full reference subtraction tail derivative adds
at most 100000000/(72m^2). The resulting uniform integrable
envelope establishes a C1 finite state energy and pressure term.

The ordinary-Proca comparison obeys conservation separately at
the exact mode, full adiabatic-order and matched local levels.
The actual clock-mass energy does not separately obey that
ordinary identity. Direct homogeneous clock variation of the
original mass functions and the same counterterm action gives

    J_clock=-[rho_mass'+3H rho_mass]
           =-[rho'+3H(rho+p)].

Compact variations keep initial Gaussian data fixed. S6.53's
uniform regulator limit also controls the distributional source
limit; the new C1 estimate identifies its continuous value.
The physical clock remains phi=tau*u.

Actual-clock local envelopes and the finite state estimates give
energy/pressure derivatives and vector clock-source magnitudes
below 10^-14 of M^2/tau^3 on u in [-1/2,1/2], for
M*tau=10^12,m0*tau=1000. All error terms are explicit.

## Remaining work

A bounded first source is not a solved quantum-corrected system.
Higher quantum functional variations, other fields/higher loops,
corrected constraints/cones, interacting cutoff and the adopted
vacuum/finite-gravity V/G/B conditions are still not established.
Small unsigned values do not protect the saturated classical cone.

The next separate state construction adds an all-order Cauchy
correction with controlled observable errors. The frozen
fourth-order state is not relabelled as Hadamard. That construction
is undergoing its own certificate checks. No old action or finite
matching prescription has changed and no user-intervention blocker
is present.

## Verification

The certificate pins 13 local sources and fully rebuilds S6.53
and its ancestry. It checks 54 named exact identities comprising
57 scalar entries, 28 proof checks and 70 rejected inputs.
Focused science passes **9 tests in 76.04 seconds**; the ordinary
suite passes **44 tests in 339.83 seconds**, without the broad
GCD adapter. The separate seeded read-only CLI passes.

The full P8 regression through S6.54 passes **4919 tests in
905.96 seconds**, without excluding a frozen checkpoint. The
exact GCD adapter passes 128 original tuple comparisons and
records 5850 domain fallbacks and 6509 exact descents under
the per-test seed recipe. Report mutation and source-manifest
checks pass. The report SHA-256 is
845f1ca80239f75a03d59d5a4cfcb05b3e2fdeefe084e49b71a57480e7087a47.

This is exact symbolic verification with written continuous
proofs, not proof-assistant formalization or original P8 closure.
