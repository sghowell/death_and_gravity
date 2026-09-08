# P8: exact vector elimination and full-source remainder

Original P8 remains open. The scoped photon objective, linear CD/M1
classification, physical matter frame and adopted V/G/B contract are
unchanged. This follows the [local nonlinear constraint theorem](assessment-2026-09-07-p8-nonlinear-constraints.md).

## New result

[S6.46](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/FORMULATION.md)
uses the unchanged source-aligned action. At fixed light fields its
retained vector is exactly Gaussian. The ordered induced operator is
zeta D-zeta^2 D(M+zeta D)^-1 D, without commuting the mass and curl
operators. Completing the square retains the full source-centered term.

The stationary-action identity and retarded initial-value response are
kept separate. An explicit causal countercontrol shows that inserting
a retarded kernel into a one-copy action gives an advanced contribution
on variation. No general differential-inverse bound is inferred from
the exact algebraic identity or the old background-frequency floor.

The original Q coefficient ODE yields a continuous cubic Taylor
remainder and a separately derived quadratic bound for its derivative.
For |N-1|, |K_hat-3H| and their physical spatial-gradient norms
<=epsilon<=1/25, the full normal source differs from the actual
quadratic source by at most 104 epsilon^3; the coordinate source by
75 epsilon^3; its spatial derivative by 230 epsilon^3; and its physical
electric curl by 300 epsilon^3. The full source remains purely normal.

The first local Maxwell density evaluated on that source differs from
the same-physical-metric quadratic-source substitution by at most
14325 zeta epsilon^5 per physical volume. Physical normal errors
divide by tau, electric errors by tau^2, and normalized densities
multiply by M^2/tau^2. Nonunit tests retain every scale factor.

The exact nonlinear source preserves the prepared first three zero
time jets and the nonzero third jet for the prior profiles used as
specified fields. This does not turn the linear light solution into
a full nonlinear one. The rational lapse map also produces further
Fourier harmonics and variable vector coefficients; the old diagonal
fixed-output response estimate is not applied to those by assumption.

## Remaining work

These are genuine full-source and local-density substitution bounds,
not a bound on the nonlocal induced action, its causal variation,
full nonlinear evolution, quantum loops or an interacting cutoff.
The explicit source neighborhood is not a numerical radius for the
preceding nonlinear implicit-function theorem. V/G/B remains open.

The exact Gaussian formula identifies a nontrivial light-dependent
vector determinant even where the classical source vanishes. Work
continues on its first scheme-explicit local quantum coefficient,
with the temporal constraint and all three massive-vector polarizations
retained. Its sign conventions are being checked directly from the
dimensionally regulated integral. No frozen action or physical
counterterm is silently changed, and no user intervention is needed.

## Verification

The report pins 13 sources and fully rebuilds S6.45 and its ancestry.
It checks 15 named identities comprising 26 scalar entries, 31
continuous/ordering/interface proof checks and 42 rejected inputs.
The independent scientific suite passes **32 tests in 0.53 seconds**;
the ordinary suite passes **57 tests in 226.37 seconds**, without
the broad GCD adapter. The separate seeded CLI passes.

The full P8 regression through S6.46 passes **4554 tests in 853.96
seconds**, with no frozen checkpoint excluded. The unchanged exact
GCD runner passes 128 original tuple controls and records 5717 domain
fallbacks and 6509 exact descents under the per-test seeded recipe.
Verification is exact symbolic checking with written continuous proofs
and independent controls, not proof-assistant formalization or external review.
