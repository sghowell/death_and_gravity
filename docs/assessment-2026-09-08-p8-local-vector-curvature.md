# P8: actual local vector curvature counterterms

Original P8 remains open. The scoped P8(a) photon objective, original
linear CD/M1 classification, physical matter frame and V/G/B contract
are unchanged. This follows the [mass-insertion loop result](assessment-2026-09-08-p8-vector-mass-insertion-loop.md).

## New result

[S6.49](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/curvature/FORMULATION.md)
derives the retained vector's local curvature pole on the actual
rolling clock. Its fixed-light on-clock Hessian is ordinary Proca
with constant mass, since the original dictionary gives a=b=1 all
along the clock. Off-clock variations do not vanish and are not
reconstructed from that on-clock value.

The Euclidean determinant prescription retains one Hodge-vector
trace minus one scalar trace. A finite Hodge complex independently
checks the operator/determinant factorization, its contact factor and
its scalar constant mode. It rejects dropping the finite contact
factor without specifying a regulator. Twenty-parameter curvature
tensor contractions independently fix the endomorphism and connection
curvature traces in the local Laplace-type coefficient formula.

With epsilon_DR=(4-d)/2 and the common 1/(32*pi^2*epsilon_DR)
normalization, the negative Euclidean pole weight is
3*m0^4/2+m0^2*R/2+a4. Here

    a4=-R^2/8+29*Ricci^2/60-Riemann^2/15-boxR/15.

Direct Christoffel/Riemann contractions of the original physical
FLRW metric, not a de Sitter replacement, give

    tau^4*a4=-8*(77*u^4+14*u^2-3)/(1+u^2)^4.

The exact square identity gives 0<R*tau^2<=49 on all real u;
a convex binomial-coefficient proof gives |a4*tau^4|<=308/3.
For Rm=m0*tau, the curvature correction relative to the flat
pole weight is bounded by 49/(3*Rm^2)+616/(9*Rm^4). This is
a local pole-coefficient ratio, not a finite quantum error.

The displayed divergence is kept. The boundary fluxes for boxR
and Gauss-Bonnet grow as u^9; their vanishing at infinite time
cannot be assumed. No finite spacetime-integrated action is claimed.

An independent unit-S4 spectral calculation agrees with the three
local coefficients. Keeping the scalar constant mode, its heat-trace
difference from 1/(2s^2)-1/s-11/30 is at most 1609s/1728 on
0<s<=1, proved with an explicit Euler-Maclaurin remainder. This is
a bounded diagnostic on S4, not a finite bound on the FLRW loop.

## Remaining work

Local counterterm functions can be evaluated covariantly without
thereby controlling the finite Lorentzian state-dependent part. The
proof preserves that distinction. Next work derives the actual vector
contribution to the lapse equation, including nonzero mass derivatives,
and bounds the finite difference between exact rolling modes and an
explicit fourth-order WKB reference across all momenta.

State admissibility, full subtraction/counterterm matching, a finite
renormalized curved quantum contribution, other field loops, corrected
constraints/cones, interacting cutoff and V/G/B remain separate work.
No frozen action or finite counterterm has been changed, and no user
intervention is currently needed. Original P8 is not finished or closed.

## Verification

The report pins 12 local sources and fully rebuilds S6.48 and its
ancestry. It checks 44 named exact identities comprising 81 scalar
entries, 21 continuous/spectral proof checks and 22 rejected inputs.
The scientific suite passes **20 tests in 2.07 seconds**. The ordinary
suite passes **50 tests in 233.89 seconds**, without the broad GCD
adapter; the separate seeded read-only CLI passes.

The full P8 regression through S6.49 passes **4701 tests in 811.65
seconds**, without excluding any frozen checkpoint. The unchanged
exact GCD adapter passes 128 original tuple comparisons and records
5724 domain fallbacks and 6509 exact descents under the per-test seed
recipe. Source manifests and missing/extra/mutated report fields are
checked. This is exact symbolic verification with written continuous
proofs and explicitly identified mathematical inputs, not proof-assistant
formalization or external peer review.
