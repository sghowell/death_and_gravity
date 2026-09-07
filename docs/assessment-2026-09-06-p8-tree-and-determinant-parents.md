# P8: tree and determinant parent-matching obstructions

Recorded 2026-09-06 after the
[auxiliary/star checkpoint](assessment-2026-09-06-p8-star-parent-no-bounce.md).
These are extensions of named parent exclusions, not a universal
UV-completion theorem or a change to the completed linear classification.

## The physical source placement is essential

The new work tests whether extra interacting metrics and extra separately
conserved NEC matter sectors evade the preceding no-bounce results.
It also checks a recent genuinely non-pairwise interaction, without
confusing its physical matter metric with its summed coframe.

| Result | Action/source setting | Essential boundary |
| --- | --- | --- |
| S6.15 auxiliary | Two Einstein leaves; actual matter on the nondynamical auxiliary metric | At least one nonzero link and a regular positive physical chart |
| S6.16 star | General constant pairwise star; sole central NEC matter | Positive central kinetic term or a genuine link; algebraic leaves cannot be added to a fixed dynamic K |
| S6.17 tree | Constant pairwise interaction tree; separately conserved NEC matter on any vertices | Nonnegative Einstein coefficients and a strictly positive coefficient at the chosen physical vertex |
| S6.18 determinant | Rank-one determinant of a coframe sum; separate matter on the actual Einstein metrics | Aligned positive flat-FLRW metrics and both summed lapse and scale nonzero |

The new formulations and written proofs live in
[TREE](../problems/P8/s6/matching/star/tree/FORMULATION.md) and
[DETERMINANT](../problems/P8/s6/matching/determinant/FORMULATION.md).
Neither conclusion requires a tensor-cone screen or a stationary massive
spectrum. No such spectrum or perturbative-health result is inferred.

## Tree: an instantaneous component, not a global kinetic sum

For an oriented pairwise edge i to j, use `y=a_j/a_i`, `c=N_j/N_i` and
`P=2(beta1+2 beta2 y+beta3 y^2)` in the actual `-2 beta` action convention.
Literal lapse/scale variations give the interaction balance C_i and its
reciprocal C_j:

```
C_i=3P c(y H_j-H_i),     C_j=-C_i/(c^2 y^3).
```

The edge flux has weight `N_i^2 a_i^3 C_i`; the opposite endpoint gets
its negative. This is distinct from the `N_i a_i^3` weight that cancels
interaction null stresses. Separately conserved matter and the full
Einstein equations, including algebraic zero-Einstein-coefficient
equations, give zero divergence of these fluxes at every vertex.
Leaf elimination on a finite tree forces every edge flux to vanish.
Cycles have a nonzero incidence kernel, so that argument cannot simply
be reused for a cycle. The cycle control is not a claimed bounce solution.

At each time take the component containing the physical vertex using only
edges with nonzero current P. On its internal edges velocity locking and
its derivative hold on a neighborhood; every boundary edge has zero null
stress at the time in question. Sum only this component's null equations.
With physical time `dT=N_r dt`, `y_i=a_i/a_r`, `c_i=N_i/N_r`,

```
A_C=sum_(i in C) G_i y_i^2 >= G_r > 0,
-2A_C H_r'+H_r (A_C')_fixed=sum_(i in C) c_i y_i^3 nu_i >= 0.
```

The derivative holds the vertex subset fixed. It does not differentiate
a changing branch label. On any compact regular interval the coefficient
`(A_C')_fixed/(2A_C)` is bounded by the largest absolute logarithmic
scale-ratio rate, independently of which component is selected. The
positive part of H_r therefore obeys a Gronwall inequality. This excludes
every nonpositive-to-positive transition without assuming finitely many
branch changes or omitting root points.

An actual four-vertex Einstein--scalar solution includes a zero-kinetic
intermediate vertex. A separate actual algebraic-branch attachment shows
why including all vertices in one exact K formula would be wrong. The
strictly positive kinetic coefficient at the physical vertex is necessary
for this particular proof: an undetermined disconnected zero-kinetic
center supplies an explicit outside-hypothesis control.

## Determinant: same regular obstruction, different matter metric

The non-pairwise model uses `L_int=-lambda S_N S_a^3`, with
`S_N=sum beta_i N_i`, `S_a=sum beta_i a_i`, and constant `lambda>0`.
Its physical matter fields couple separately to the Einstein metrics,
not to `U=sum beta_i e_i`. Direct coframe variation gives

```
rho_int,i=lambda beta_i S_a^3/a_i^3,
p_int,i=-lambda beta_i S_N S_a^2/(N_i a_i^2),
C_i=3lambda beta_i S_a^2/(N_i a_i^3)
    *[dot S_a-S_N dot a_i/N_i].
```

For active nonzero beta and both sums nonzero, C_i=0 locks all proper
scale-factor velocities, including their zeros. The weighted interaction
nulls cancel exactly. With the physical clock of an active Einstein metric,
the positive active-only K consequently satisfies
`(H/sqrt K)'=-sum c_i y_i^4 nu_i/(2K^(3/2))<=0`, where this model's
`c_i=N_i/(N_r y_i)` is a speed ratio, unlike TREE's lapse ratio.
Disconnected zero-beta metrics independently obey the flat GR null equation.

The paper's same-sign beta restriction guarantees regular sums for positive
coframes. The direct mathematical proof also admits mixed signs whenever
both sums stay nonzero; it does not establish those parameters' health.
Actual proportional rolling canonical-scalar and vacuum solutions check
the full equations, with all reference-clock factors retained.

There is an exact stationary auxiliary rewriting, but matter must stay on
the original Einstein leaves. Moving it onto the auxiliary field changes
the auxiliary equation and the model. Nor may the proof divide by a
vanishing summed lapse: a checked point with `S_N=0,S_a!=0` satisfies
every interaction Bianchi equation while failing velocity locking. That
point also fails the specified Einstein lapse equations. It is a guard
against an invalid inference, not evidence for a viable singular branch.

## Consequence for the remaining P8 work

These exclusions rule out actual contraction-to-expansion in their stated
parent domains. They cannot be repaired merely by computing a different
vacuum positivity coefficient for the same action. An actual-background
matching claim must respect the physical metric/clock dictionary; formal
operator matching alone does not supply the missing solution.

The original local-UV question remains open. Variable interaction
coefficients, derivative operators, nonaligned or singular-sum sectors,
different matter couplings and actual quantum stress require new analyses.
Some change the allowed model rather than answer the adopted question.
No omitted-operator or loop error is assigned an invented numerical bound,
and no untested alternative is counted as a UV-compatible C/D witness.

## Verification

The frozen TREE report has SHA256
`d3d4c612fd0942bc66f0005ed21b76d23ca86285a461d87b32e1a2d79f242bb2`.
Its 15 source files include 34 separately authored audit cases. The report
checks 60 primary identities, six independent coefficientwise identities,
fourteen Fraction/first-jet/Gaussian fixtures and 87 primary comparisons.
All 65 ordinary child tests passed in 117.75 seconds, as did the standalone
read-only certificate replay.

The frozen DETERMINANT report has SHA256
`97ee3e0d70864dc6c1b8d4a5088dd3b245b2f94ab57d166a88f9912c5fb74906`.
Its 14 source files include 22 separately authored audit cases. The report
checks 41 primary identities, ten independent coefficientwise identities,
192 literal full-coframe directions, sixteen full auxiliary stationary
directions and 35 primary comparisons. Independent actual scalar and
vacuum Einstein fixtures are retained. All 49 ordinary child tests passed
in 130.27 seconds, as did the standalone read-only certificate replay.

All 29 source hashes match and Ruff passes for both children. Neither
ancestor certificates nor unrelated P4/P9 work were changed for these gates.
The 56 separately authored audit cases also passed together in 0.75 seconds;
all 64 local Markdown links in the checkpoint manifest resolve.

The combined P8 suite passed **2,154 tests in 428.41 seconds** with the
explicitly opt-in [exact regression runner](p8-exact-regression-runner.md).
Its startup checked 128 original Gaussian-domain GCD tuples; the full run
recorded 6,509 exact descents and 2,509 domain fallbacks. This is an adapted
exact full-suite pass, not an ordinary unmodified SymPy run. All ten ordinary
runner tests passed in 0.17 seconds. The known unchanged long symbolic test
emitted the configured 60-second diagnostic traceback and subsequently
passed; the run was not stopped or counted from partial output.
