# P8: the Ricci difference is a curved auxiliary deformation

Original P8 remains open. The completed scoped photon objective and
linear CD/M1 classification are unchanged. This follows the
[constant two-trace classification](assessment-2026-09-07-p8-constant-two-trace-family.md),
without modifying any frozen ancestor or enlarging the original
completion requirements.

## Result

[S6.40](../problems/P8/s6/matching/affine/kinetic/ricci/FORMULATION.md)
adds only the constant-coupling square of the antisymmetric difference
of the two specified Ricci contractions to the unchanged unrestricted
S6.37 parent. It is projectively invariant, vanishes on every
Levi-Civita connection and is independent of the two distortion-trace
curls. Its definitions follow
[Barker and Marzo, equations (2), (3) and (12)](https://arxiv.org/html/2402.07641v2);
no spectral-health conclusion is imported from that paper.

The original rolling solution is preserved exactly. At quadratic order
on it, the complete connection solve has a null flat differential Schur
symbol, but a nonzero curved commutator. For an antisymmetric two-form Z,

    C[nabla Z] = K_curv Z,
    K_curv = R_background/3 = 8(1+7u²)/(1+u²)².

Every one of the 64 connection equations and every original coefficient
variation is retained. The actual source is
Cstar_0i = -10u/(1+u²)^4 * partial_i n, with Cstar_ij=0.
All time derivatives of the lapse cancel. The six multiplier components
are algebraic on 1+lambda*K_curv!=0, leaving exactly

    Delta L/a³ = q*A*n²,
    A = 50lambda*u² /
        [(1+u²)^6*((1+u²)^2+8lambda*(1+7u²))].

Here lambda=lambda_physical/(M² tau²) and q=(tau*k_physical)².
This is an actual curved, constrained quadratic result, not a flat
pole calculation or a new heavy-particle spectrum.

## Physical sign cases and the crossing

For lambda>=0 the connection denominator is positive at every time.
After both scalar constraints, including the free matter perturbation,

    L_kin/a³ = (s_dot+w*v_dot/Theta)²/2
                +(J+q*A)*v_dot²/Theta²,

so both scalar kinetic directions are positive at finite u!=0, q>0.
The principal gradient matrix is also positive throughout that domain.
Its proof uses the positive Gss pivot and determinant: their exact
numerators have 37 and 180 nonnegative, even-time monomials with positive
constants. The time-dependent spatial boundary and gyroscopic mixing
are retained; an independent Legendre transform checks the resulting
positive principal Hamiltonian square.

At the center the regular first-order denominator is J+q*A, not Theta.
Since A=A_dot=0 there, the Hamiltonian and first coefficient jet agree
with the old crossing system. A direct momentum-Hessian calculation
reproduces the old positive center velocity chart for q>6. The
first-order system is regular for every q>0; q=6 is not inverted in
that particular velocity chart.

For every lambda<0 there is a regular tail with A<0. At q>J/(-A),
the fully constrained scalar kinetic matrix has one negative direction.
This conclusion does not identify an algebraic pole with a ghost.
Singular multiplier charts are separately rejected. Lambda=0 returns
the original auxiliary theory without division by lambda.

Positive principal matrices do not establish all-frequency stability,
nonlinear stability, a time-independent conserved energy or UV viability.
The original tensor sector and quadratic mode count are unchanged on
the regular rolling branch.

## Matching limit and the nearby-field warning

Every nonzero regular coupling changes the original physical quadratic
action away from the center: its lapse Hessian changes by 2q*A and its
physical v-velocity coefficient by q*A/Theta². This is not exact CD/M1
matching in the original metric and matter frame.

There is a limited uniform comparison. Exact polynomial positivity gives
r²/(2J)<10 for every real time, where r=-10u/(1+u²)^4. Hence, for
lambda>=0 on the punctured chart and a declared momentum range q<=Q,

    10lambda*Q<=epsilon implies
    K_original <= K_new <= (1+epsilon)*K_original.

The physical bound is 10lambda_physical*k_physical²/M². It is only a
quadratic kinetic-coefficient bound, not a full action/gradient/loop
remainder or a justified frequency cutoff.

An additional off-clock control prevents an unjustified open-tube claim.
With p left generic and k=(omega,0,0,k_z), the full connection Schur
response is k_z² diag(-a_p,-a_p,0,0,b_p,b_p), where

    a_p=(2p-1)²(3p+1)/[2p(8p²-1)],
    b_p=(2p-1)²/(8p²).

It is nonzero inside the original tube away from p=1/2, but has no
frequency dependence. It is therefore neither a nonlinear inverse
theorem nor, by itself, evidence of a new propagating ghost. The
complete metric constraints and variable coefficients are essential
to a separate nearby-background conclusion.

The next matching research must retain that distinction. A positive
quadratic deformation cannot be promoted to a healthy UV parent merely
by assigning it a mass scale. Changes to the auxiliary mass terms are
another possible branch, but require a new literal action, exact
matching, full constraints and quantitative control of any heavy sector.
The adopted V/G/B obligations remain open; no user choice or external
authorization is needed for the next local calculation.

## Verification

The frozen report covers 14 sources, 36 named exact identities comprising
1463 scalar entries, 13 proof/interface checks and 28 rejected-input
controls. Independent checks include dense unrestricted curvature
contractions, exact Lorentz covariance, a coordinate FLRW commutator
with arbitrary two-form functions, joint scalar constraint solves,
time-dependent boundaries, nonunit scales and singular-chart controls.

All **71 ordinary tests pass in 137.82 seconds**, without the broad
GCD adapter, and the separately seeded read-only CLI passes. The full
**4152-test P8 regression passes in 720.28 seconds**, with no checkpoint
excluded. The unchanged exact GCD adapter passes all 128 original
normalized-tuple comparisons; the run records 5168 domain fallbacks and
6509 exact descents. Its reproducible per-test seed and exact adapter
are unchanged from the
[previous verification recipe](assessment-2026-09-07-p8-constant-two-trace-family.md#verification).

This is exact symbolic verification with written proofs and independent
internal calculations, not proof-assistant formalization or external
peer review. No original P8 closure is asserted.
