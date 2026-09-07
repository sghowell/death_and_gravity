# P8: physical tensor response and regular scalar constraints

Continuation of the [literal-vacuum and variable-parent checkpoint](assessment-2026-09-06-p8-vacuum-domain-and-variable-parent.md).
The fixed-slice tensor response is now proved with its actual physical
source and data maps. The scalar action and regular constraint reduction
have passed independent checks, but their physical stability interpretation
requires a uniform analysis that has not yet been established. Original P8
remains open; neither result replaces the adopted matching/UV contract.

## Fixed-slice tensor response

The [S6.21 response formulation](../problems/P8/s6/matching/variable/response/FORMULATION.md)
and [full proof](../problems/P8/s6/matching/variable/response/notes/proof.md)
retain the complete coupled tensor operator of the unchanged S6.20 action.
Write `delta=c-2>0`, `u=T/tau` and `kbar=tau*k_com`. The two physical slices
are fixed at `u=+-1/100`, and the spatial momentum band is `1<=kbar<=2`.
This improves the earlier result, which only controlled a shrinking
physical-time interval.

An explicit endpoint map includes the actual, delta-dependent conversion
from canonical fields to both metrics and their proper-time derivatives.
It also includes the full coupled punctured outer wave operators, defined
by a convergent Volterra series with a factorial truncation bound. In the
resulting physical-data norm the exact transfer differs from the matched
light/relative-mode transfer by at most

```
40000 delta^(1/3),       0<delta<=10^-9.
```

This is a convergence estimate, not a small-error claim on that whole box:
the bound is 40 at its upper endpoint. The explicitly smaller range
`delta<=10^-21` gives an error at most `1/250`. The endpoint maps and their
inverses are uniformly bounded; the norm is not a quotient that discards
physical data. The matching proof keeps all four components and both
moving-weight terms. Its bounds use a continuous exact interval cover,
not sampled propagation.

The inner relative mode is a Poeschl--Teller equation. The exact connection
has `|A|^2-|B|^2=1`, `B=i/sinh(pi*sqrt(39)/2)`, and a nonvanishing
log-delta transmission phase. The branch and frequency conventions are
checked against the [DLMF connection formula](https://dlmf.nist.gov/15.10.E21)
and [Gamma reflection identity](https://dlmf.nist.gov/5.5.E3).
The phase persists even if reflection is set to zero. This is a classical
response statement, not a particle-production calculation.

A conserved infinitesimal stress coupled only to the actual matter metric
g can prepare either a relative-mode input or a regular-light input on the
fixed preceding interval. The construction uses one fixed smooth compact
source for the entire positive-delta family, with its finite source norm
retained explicitly. Source-loading error is bounded by
`3,000,000 delta ||sigma||_L1`. Prepared-light leakage therefore has a
vanishing bound, whereas nonzero prepared relative data retain distinct
responses along two log-phase subsequences. After the source switches off,
the map from both metrics' Cauchy data to the g-only four-jet
`(g,g',g'',g''')` has positive determinant `(U_TT/K_g)^2`. Observability is
not assumed from an arbitrary two-component projection.

None of this selects the low-energy state or proves a low temporal band.
Fixed source duration and a fixed spatial momentum band do not do that.
The current leakage bound is also not sharp enough to compare an omission
error to the center's order-delta locked-cone excess. A merely bounded
pole-subtracted mass remainder is explicitly distinguished from the
analytic canonical correction: its two mass-only path limits are -32
and -141, so a smoothness argument at their joint origin would be invalid.

## Scalar/vector action and the nonuniform-symbol warning

The source-hashed [constraint implementation](../problems/P8/s6/matching/variable/perturbations/src/p8_variable_constraints/action.py)
starts from the literal quadratic ADM action, including both metric
lapses/shifts, the sourced scalar clock, the free spectator and the
interaction. The independently authored
[audit](../problems/P8/s6/matching/variable/perturbations/tests/test_constraints_independent_audit.py)
reconstructs its diagonal metric densities and interaction roots before
elimination, retains the Einstein time boundary, and separately performs
the Legendre/auxiliary reduction. It gives `H+N_f C_f`, with no spurious
quadratic f-lapse term.

A time-dependent canonical chart turns the remaining primary constraint
into a coordinate constraint without dividing by H. Its secondary Schur
coefficient is independent of spatial momentum; its center value is
`-5/24`, and a continuous polynomial bound gives `D<-1/8` on the frozen
local domain. This establishes regular linear constraint reduction, not
positive energy. In particular the sign of this elimination coefficient
is not a propagating scalar's kinetic sign.

The physical observable map retains both its first and second time
derivatives. In the checked c=4 example, at the center the chosen three
Bardeen observables commute, but they do not commute at nearby nonzero
times in the reduced Dirac
bracket. Its Cauchy determinant has a numerator of degree two in
`K=kbar^2` at the center and degree three at a fixed nearby time. Thus
freezing its apparent high-momentum equation and then taking the bounce
limit is not a uniform operation. Apparent center sound speeds and an
instantaneous Hamiltonian momentum block are not promoted to a scalar
health or characteristic-cone theorem.

The vector calculation independently uses the unit-Jacobian transverse
spatial pullback and both shifts. It retains the full kinetic coefficient,
relative stiffness, canonical normalization pump and boundary term.
Its positive-branch formal vector kinetic sign is useful, but it does not
establish scalar health or a verified subcutoff frequency window.

## All-time background exploration is separate

The unpromoted [G1 feasibility note](../problems/P8/s6/matching/variable/global/notes/feasibility.md)
constructs a different, variable-lapse action with the CD-shaped physical
scale factor over all real time. Exact reconstruction, continuous
positivity bounds and tail estimates give a globally analytic canonical
clock and a complete physical g metric. This is not a continuation of the
constant-c action: the lapse and coupling second derivatives differ.

Two limitations are recorded explicitly. The f metric has finite null
affine length in both tails despite bounded curvature, and its relative
algebraic tensor mass tends to zero as `4/u^30`. That coefficient cannot
supply a uniform positive heavy threshold. Neither fact alone proves a
physical instability or rules out every coupling-suppressed reduction.
G1 has exploratory tests, but no independently reviewed source-hashed
certificate or certified ledger row.

## What must happen next

The immediate research tasks are a uniform scalar reduction with a genuine
physical energy/propagation interpretation, and sharper prepared-light
response estimates on an explicitly justified source/state and frequency
window. A successful candidate would still need the common-action domain,
canonical operator and matter-frame match to the original C/D witness,
with gap/cutoff and omitted-term errors, and the adopted vacuum plus
finite-gravity positivity gates. A result about this named parent cannot
decide an entire DHOST row's UV status.

No new user choice or external permission is currently required for these
mathematical checks. The completed scoped photon result is unchanged.
The project must remain open until its outstanding original requirements,
not just these candidate-level checkpoints, have been met.

## Verification record

The frozen tensor child passes 81 ordinary tests in 149.80 seconds and its
standalone `--check` replay. Its separately authored physical-action audit
contains 38 tests. The report retains 33 exact residuals, 18 continuous
coefficient margins, 23 omission/domain controls and 17 verified source
hashes. Report SHA-256:
`7ac41e1ca23c9fe9649e6fc703c5d3d236fb599481f921db319416e10865f118`.

The frozen scalar/vector report retains 440 exact residuals, 1,264
independent output comparisons, 140 independently derived positive
Bernstein coefficients, 18 warmed-cache invalid-input controls and 18
verified source hashes. Its separately authored ADM/physical-map audit
passes 30 ordinary tests. Report SHA-256:
`9cb56291ce32499293a65ee1355e9b13d3abc4d429643f76c0dcf63e9b71a9b1`.
The full frozen suite passes 76 ordinary tests in 178.90 seconds and its
standalone `--check` replay; no arithmetic adapter was needed.

G1 passes 19 ordinary exploratory tests in 0.72 seconds. This is not a
source-hashed or independently audited theorem package.

The full checkpoint regression passes **2,452 tests in 504.91 seconds**
with the existing separately tested exact Gaussian GCD adapter. Its 128
startup comparisons pass; the run records 6,509 exact descents and 3,154
domain fallbacks. This is an exact-arithmetic-adapted full-tree run, not
an unmodified SymPy claim. The unchanged M1 nonlinear check emits its
60-second diagnostic stack before passing. Active prepared-light and
microlocal research children are excluded from that checkpoint run and
from its commit manifest. The two new frozen children separately pass
their ordinary suites as recorded above.

All 35 new frozen source hashes and 72 local Markdown links are verified.
Ruff, default whitespace checks and exactly-one-final-newline checks pass
for the 47-file checkpoint. No frozen ancestor was edited; no unrelated
P4/P9 work is included.
