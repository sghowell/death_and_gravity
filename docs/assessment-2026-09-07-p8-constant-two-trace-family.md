# P8: the constant two-trace family is now closed, not the UV problem

Original P8 remains open. The completed scoped photon objective and
linear CD/M1 classification are unchanged. This extends the
[three specified kinetic screens](assessment-2026-09-07-p8-affine-kinetic-screens.md)
while leaving every ancestor action and certificate frozen.

## Result and why the exceptional cases matter

[S6.39](../problems/P8/s6/matching/affine/kinetic/traces/FORMULATION.md)
keeps the exact unrestricted-connection CD/M1 parent and adds only
a constant real symmetric kinetic matrix Z for the two specified
projectively invariant distortion-trace curls. Off-diagonal mixing
and every matrix rank/sign are included. The original physical metric,
clock, free chi and complete rolling solution are unchanged.

The independently reconstructed eight-trace Schur matrix is

    D8 = D2 tensor diag(1,-1,-1,-1),
    D2 = (3/8)*[[1,-7],[-7,1]].

The all-64 connection equations are retained. The full-rank trace
problem has a nondegenerate 52-dimensional auxiliary complement.
For a rank-one term Z=zeta*(A,B)^T*(A,B), however, the relevant Schur
scalar is gamma=3(A²-14AB+B²)/8. Its zero values are separate
constraint charts, not smooth inverses of the ordinary Proca chart.

| Nonzero kinetic form | Decisive physical obstruction |
| --- | --- |
| Rank one, gamma nonzero, zeta negative | Negative unconstrained transverse velocity term |
| Rank one, gamma positive, zeta positive | Negative clock Schur pivot for q>2 gamma J/e² |
| Rank one, gamma negative, zeta positive | Negative longitudinal pivot for q>−1/(gamma zeta) |
| Rank one, gamma zero, either zeta sign | Nondegenerate higher-derivative Legendre map; Hamiltonian affine in an unconstrained momentum |
| Rank two, Z not positive definite | Negative transverse velocity direction |
| Rank two, Z positive definite | Indefinite physical longitudinal submatrix for q>(2/9)(1,1) Z^-1 (1,1)^T |

Here e=(A−B)/(1+u²)^3 and q is normalized squared spatial momentum.
The scalar reductions use finite u!=0, where the original Theta is
nonzero, and avoid each explicitly singular temporal chart. The zero
matrix is the unchanged auxiliary theory, not a propagating completion.

The null rays A/B=7±4sqrt(3) are the important new case. The full
connection equations force T=Tstar, but the curl of Tstar contains a
time derivative of the lapse. Only after the actual lapse/shift and
free-matter equations are retained can its effect be determined. The
joint highest-derivative Hessian in (v_ddot,s_dot) has determinant
a^6*zeta*q*e²/Theta². Its exact Hamiltonian is affine in an independent
P0, so no remaining constraint removes the higher-derivative mode.
The complete rolling Ustar source is not merely minus Vstar; its
background-coefficient contribution is included before this reduction.

The rank-two argument is also coupled, not an isolated mass test.
The full lapse/shift/both-temporal-vector auxiliary determinant is
−4q²Theta² det(D2^-1+qZ). The negative longitudinal submatrix is a
restriction of the physical four-scalar velocity form after those
constraints, which makes its negative value a physical obstruction.

## Scope and next research

This exhausts this particular constant two-trace family as a literal
everywhere-healthy parent of the original rolling solution. It does
not prove that an unhealthy frequency lies below a justified EFT
cutoff. It does not classify field-dependent curl matrices, other
curvature/torsion operators, additional mass tunings or new fields.
The zero-Schur result is quadratic on the actual rolling trajectory,
not a nonlinear open-tube inverse.

A primary-source check of the unrestricted Ricci-type spectrum in
[Barker and Marzo, sections III.C and IV](https://arxiv.org/html/2402.07641v2#S3.SS3)
does not supply a ready-made healthy massive parent here. Its explicit
healthy unrestricted example retains a massless vector, including the
already tested homothetic Maxwell case; healthy massive examples use
kinematic restrictions on the connection. The paper's Minkowski
quadratic analysis is explicitly not an exhaustive rolling/nonlinear
classification. Those restrictions must not be imported into our
unrestricted matched action.

A bounded next calculation is the antisymmetric difference of the
two Ricci contractions. It is independent of the two trace curls and
vanishes on a Levi-Civita connection. Preliminary flat-symbol algebra
indicates a null Schur response, so its full curved commutator and
actual rolling source must be derived before any new-mode or health
claim. It is not selected as a healthy UV parent on the strength of
the literature or a raw kinetic term.

The adopted vacuum, finite-gravity and controlled matching V/G/B
obligations remain research work. No new assumption, spending approval
or user intervention is needed for the next local calculation.

## Verification

The report covers 18 sources, 98 named exact identities comprising
7741 scalar entries, 20 continuous/rank/interface proof checks and
32 rejected-input controls. The scientific suite passes 118 tests,
including independent all-component inversions, generic full scalar
embedding, a time-dependent variational chain-rule check, both null
rays, both kinetic signs and nonunit normalization.

After the final documentation review made both CLI RNG seeds explicit,
the regenerated report and separate seeded CLI pass. The complete
ordinary suite passes **136 tests in 87.59 seconds**, without the broad
GCD adapter. The full **4081-test P8 regression passes in 731.20 seconds**,
with no checkpoint excluded. Its unchanged opt-in exact GCD adapter
passes 128 original normalized-tuple comparisons, with 5097 domain
fallbacks and 6509 exact descents during the tests.

The broad regression retains the per-test SymPy RNG reset that resolved
the previously recorded stationary-series arithmetic stall. It changes
no formula, assertion, cache policy or frozen source. From the repo root:

```sh
PYTHONHASHSEED=0 .venv/bin/python - <<'PY'
import sys
from sympy.core.random import seed
seed(0)
sys.path.insert(0, "scripts")
import p8_exact_regression as runner
import pytest
class ReproducibleSeed:
    def pytest_runtest_setup(self, item):
        seed(0)
print("Exact GCD adapter self-check:", runner.self_check())
with runner.exact_runner() as counts:
    try:
        result = pytest.main(runner.pytest_arguments(
            ["problems/P8", "-q", "--full-trace", "-p", "no:faulthandler"]),
            plugins=[ReproducibleSeed()])
    finally:
        print("Exact GCD runner counters:", dict(counts))
raise SystemExit(result)
PY
```

These are exact symbolic certificates with written proofs and an
independent internal audit, not Lean formalization or external peer review.
