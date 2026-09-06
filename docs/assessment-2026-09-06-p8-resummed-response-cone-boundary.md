# P8 after the resummed-response and all-mode cone screens

Recorded 2026-09-06; extends the
[A.12/S6.7 checkpoint](assessment-2026-09-06-p8-exact-qsei-vector-screen.md).
Original P8 remains open. Neither of the following scoped theorems is a
replacement for its realistic-field focusing or EFT matching obligations.

## New certified results

| Track | Result | Scope boundary |
|---|---|---|
| [A.13](../problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/continuation/FORMULATION.md) | Exact resummation of the stiff Einstein/logarithm block, with its complete causal poles/cut and an explicit weighted inverse bound | The full rolling/state remainder is retained but not bounded on a longer interval; the comparator pole is not an actual nonlinear runaway theorem |
| [S6.8.COMPOSITE](../problems/P8/s6/matching/composite/perturbations/cones/FORMULATION.md) | All-degeneracy no-bounce theorem under the explicit fundamental tensor/vector principal-cone requirements | It is not automatically a light-only EFT theorem below a heavy threshold, or a general composite/UV exclusion |

Certificate SHA256 values:

    A.13: 1569abacc14b4c53c951fc10e07eea7e5f17625812705bfb5709772a26261d97
    S6.8.COMPOSITE: 89145c7a91981f72941a30e902b7e83aa231f4f4c97440c8e6f454ba339c340a

A.13 hashes 18 files and verifies 16 exact identities plus independent
Fraction bounds. S6.8 hashes 13 files and verifies 17 exact identities,
six independent Fraction coefficient identities, and nine independent
covariant/API audits in its test suite. Their written contour, functional,
characteristic and connected-interval arguments are not Lean-formalized.
Old source-hashed subtrees and unrelated P4/P9 work remain unchanged.

Validation: **1,377 P8 tests passed** in 403.77 seconds with
`PYTHONHASHSEED=0` and diagnostic `faulthandler_timeout=60`, including all
40 A.13 and 36 S6.8 tests. Both read-only certificate CLI replays, full P8
Ruff, local artifact links and diff checks pass. An earlier unseeded run
was interrupted after 905 passing tests when an unchanged S5 symbolic
simplification took excessive time; no mathematical failure or old-source
edit was involved.

## A: the stiff response must be controlled, not removed

The full trace and integrated fixed-point map allow a constant
preconditioner with transform `2/[log(s)+beta-c/s^2]`. The variable Einstein
coefficient, anomaly, auxiliary Wick product, nonlinear actual-state
history and baseline residual remain explicit in the other side of the
equation. This is an exact rearrangement, not a full rolling linearization.

The principal logarithm has one positive pole and two damped poles. A
contour calculation retains all three and the positive cut; a separate
positive Volterra series proves the complete kernel's positivity. Exact
rational bounds place the growing pole between `10^6` and `2*10^6` in the
named dimensionless calibration. The unweighted comparator norm already
exceeds 1000 at duration `10^-5`. Deleting the growing term fails the
causal inverse identity and is not done.

An exponential weight instead gives a half-line norm below `240/769` at
weight rate `2*10^6`. That is useful only with complete remainder and
self-map estimates: removing the weight multiplies errors exponentially.
The next test keeps the original fixed preparation and data and targets
a longer, but still local, actual solution. It must preserve pointwise
geometric barriers, the actual state and propagated density constraint.
No no-runaway condition, state reset, order reduction or new physical
prescription has been selected silently.

## B: two physical tensor cones impose a strong restriction

For constant positive affine couplings, the two tensor speeds relative
to composite matter are `r/s` and `c*r/(y*s)`. Requiring both to be at most
one forces `c=y`. The complete pressure-dependent tensor stiffness then
becomes `mu=y*Q`. Positive regular vector inertia and strictly positive
vector principal speed give `Q>0`, so the undivided branch equation forces
`B=0`. The scale ratio becomes constant and the physical Hubble is
nonincreasing by the actual null equation. This covers arbitrary
degenerate transitions, all fixed HR coefficients and independent mass
scales within the stated class.

An actual local bounce with luminal tensor cones but zero vector speed
shows why strict vector positivity matters. A nonbouncing proportional de
Sitter control satisfies the stated TT/vector subset. Neither control
establishes a scalar stability result. The source-aware pressure and both
fundamental tensor roots are essential to the proof.

The original S6 matching question still permits studying a lower-energy
EFT in which a genuinely heavy mode has been integrated out. Fundamental
principal cones cannot simply be imposed on that light-only theory without
checking the heavy threshold and the reduction error. The next calculation
therefore retains the rolling tensor diagonalization, normalization
derivatives, retarded heavy response and heavy initial data. Setting the
two tensor fields equal is only a conditional leading approximation.

## Status

The frozen 32-row linear classification and these scoped gates are
complete. Original realistic-field cosmological focusing and controlled
matter/positivity matching are not. No missing credential, user choice or
external approval currently blocks the in-scope follow-on calculations.
P8 is not finished or closed.
