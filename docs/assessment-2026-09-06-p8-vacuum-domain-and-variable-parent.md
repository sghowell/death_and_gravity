# P8: literal vacuum domain and a variable-coupling parent

Recorded 2026-09-06 after the
[tree/determinant checkpoint](assessment-2026-09-06-p8-tree-and-determinant-parents.md).
The new results identify a genuine obstruction for an unchanged published
action and a genuine local bounce outside the constant-interaction
no-bounce theorems. Neither supplies the missing common-parent UV and
finite-window matching certificate. P8 remains open.

## Two different questions, kept separate

| Checkpoint | Established | Not established |
| --- | --- | --- |
| S6.19.A26 | The literal June-2026 DHOST action cannot acquire a smooth constant-clock vacuum chart while remaining unchanged near nonzero X on an open finite-field neighborhood | Exclusion of a separately specified off-tube repair, cosmological instability, or a UV no-completion theorem |
| S6.20.VARIABLE | An actual local scalar-dependent-interaction bimetric bounce with a separately conserved free spectator; full canonical tensor operator and a nonadiabatic inner limit | Scalar/vector health on a verified frequency window, a uniform heavy gap/cutoff, original C/D operator matching, or vacuum/finite-gravity positivity |

The [A26 formulation](../problems/P8/s6/matching/a26/FORMULATION.md) and
[written proof](../problems/P8/s6/matching/a26/notes/proof.md) concern a
literal source action. The
[VARIABLE formulation](../problems/P8/s6/matching/variable/FORMULATION.md)
and [proof](../problems/P8/s6/matching/variable/notes/proof.md) concern a
separately named reconstructed parent. The latter is not a completion of
the former, and neither is silently substituted for the frozen free-M1
C/D witness.

## A26: test the full directional action, not isolated poles

The source uses signature `(-,+,+,+)` and `X=(partial phi)^2` without a
half. Its rolling clock has `X=-1`; a constant-field Lorentz-invariant
vacuum has `X=0`. At a defined finite coefficient point set
`f=1/2-g1(phi0)` and `r=2g1(phi0)-a1(phi0)`. For `f!=0`, the residues of
`A3,A4,A5` are `r,-r,-r^2/(2f)`. On fixed Minkowski space,
`phi=phi0+epsilon psi` has leading higher-derivative action

```
epsilon^2 r L[psi],       L=(L3[psi]-L4[psi])/Xpsi.
```

The A5 term first enters at fourth, not second, small-field order.
Coefficient poles alone would not prove an obstruction: degeneracy and
boundary cancellations must be checked in the full action.

Here the complete symmetric-Hessian Euler variation on an embedded
2+1-dimensional quadratic jet gives `E=-2d0*d1*d2/X`. The timelike-gradient
examples have `E(psi1)=2`, `E(psi2)=1/2`, but `E(psi1+psi2)=16/9`.
Their nonadditivity `-13/18` is incompatible with a bilinear action
Hessian; an integration-by-parts boundary term cannot change an Euler
expression. Setting all but the temporal Hessian entry to zero instead
gives zero, demonstrating why homogeneous-only tests would miss this.

A second calculation uses
`psi=t+lambda T(t)(cos x+cos y)` with smooth compact-time T and small
lambda. The quadratic spatial average is the time boundary `(T T')'`.
After every time boundary is removed, the quartic coefficient is
`-integral(T^2 T'^2)/2<0`. A fourth finite difference cancels every
quadratic functional and every common-boundary contribution but retains
this bulk term. The independent replay uses exact Fourier coefficients
and a different logarithmic integration-by-parts identity, not sampled
quadrature or the primary symbolic expression.

For the specified source background, the analytic numerator controlling
`r` has derivative `-1/20` at zero. Its zeros cannot contain an interval;
neither can the exceptional denominator or Einstein-coefficient zeros.
Every open finite-field interval therefore contains regular nonzero-r
points with the obstruction. This excludes an unchanged smooth vacuum
chart, not just one unfortunate expansion point. An isolated `r=0` point
is not classified by this leading-order test.

The source's printed A4 parentheses differ from its earlier Ia expression.
Both versions are explicitly transcribed and their difference is regular
at the nondegenerate X=0 points; the leading result is the same. No
cosmological degeneracy calculation is silently repaired. These source
inputs are documented against the
[June action](https://arxiv.org/html/2606.03302v1) and
[earlier Ia formula](https://arxiv.org/html/2501.09985v2).

The optional fractional-power audit chooses `b=9/500` itself; this is not
attributed as a printed numerical parameter. Source equation (21) maps
that choice to `n_s=24/25` for the stated background epsilon=10. In these
literal coordinates Q is C2 but not C3, and the named nonconstant entropy
profile makes W2 C1 but not C2. W2 is finite. Crucially, all reconstructed
Q-dependent terms cancel off shell on `chi=B(phi)`. After the smooth
entropy shift, the linear source also vanishes on the rolling clock.
The result is a coefficient-regularity statement, not a claimed physical
instability or a replacement for reducing the perturbations.

## VARIABLE: the source term really evades the old locking identity

The new action has positive Einstein terms for g and f, a canonical clock
phi and the free canonical spectator chi on the actual g metric, and
`-2 sqrt|g| sum beta_n(phi) e_n(sqrt(g^-1 f))`. Its scalar-dependent
interaction is motivated by the general Einstein-frame construction in
[Wood's scalar-tensor multigravity action](https://arxiv.org/html/2501.16442v2).
This fixture chooses unit matter conformal factors; it is not that paper's
specific dimensional-deconstruction model or a proven UV theory.

Write `y=a_f/a_g`, `c=N_f/N_g`, `P=2(beta1+2beta2*y+beta3*y^2)`, and let
U and V be the usual two lapse polynomials. With both possible canonical
clocks retained in the general derivation, the combined balance is

```
B=3P(y H_f-H_g)+2(U_phi_f q_f-V_phi_g q_g)=0.
```

The reciprocal edge weights remain correct, but the last term is not
zero for variable coefficients. Thus vanishing B no longer implies the
constant-tree velocity lock. The clock is sourced by the interaction;
the physical free chi remains separately conserved. Conflating those
two conservation laws would manufacture a false no-bounce conclusion.

In g proper time, put `u=T/tau`, `d=1+u^2`, `a_g=d^2`, `a_f=2/d^2`,
`N_g=1`, `N_f=c`, and equal Einstein coefficients M^2. Explicit analytic
beta functions are reconstructed using a strictly monotone canonical
clock. They are functions of phi through its local analytic inverse,
not externally prescribed functions of time. All four metric equations,
the clock equation, the free-chi equation and the sourced Bianchi relation
are checked exactly on `abs(u)<=1/10`, for `c=1` or `2<c<=4`.

The free charge is `a_g^3 dot(chi)=M/(10tau)`, matching the normalized
free-spectator background of the frozen witness. The physical g scale
factor also has the CD form. These are background agreements only:
they do not establish the required field/operator/canonical dictionary.
This construction is local, not a geodesically complete global solution.

## Retain the actual rolling tensor operator

Literal relative-shift variation gives the coefficient
`a_g^5 y^2 P/[2N_g(c+y)]`. Exact two-shift elimination makes its sign the
formal high-frequency relative-vector kinetic sign. The c=1 bounce has
a negative coefficient and is a negative control, not a healthy candidate.
This formal principal warning is not promoted to a subcutoff instability
without a cutoff estimate. At c=4 the coefficient is positive, but the
f principal tensor speed squared relative to g matter is four. Neither
control proves scalar health or excludes every light-only reduction.

The full TT action retains both Einstein kinetic terms, both spatial
gradients and the literal interaction stiffness. The moving common and
relative canonical variables l,Q have normalization rates thetaS,thetaR
and mixing omega. After its explicit normalization boundary term, the
equations have the form

```
l''+VLL*l+B_adjoint*Q=0,
Q''+VHH*Q+B*l=0,
B=2omega*(partial_u-thetaS)+q^2*D.
```

VHH includes the relative mass, the second normalization derivative and
`-4omega^2`; the adjoint includes `-2omega'`. The physical source couples
to the g-metric projection of l and Q, not automatically to l alone.
At c=4,u=0 the algebraic relative mass squared is `24/tau^2`, whereas
the zero-momentum canonical heavy diagonal is `22/tau^2`. Neither number
by itself is a rolling spectral-gap theorem.

Let `delta=c-2>0` and use the shrinking variable `u=sqrt(delta)*x`.
For fixed finite spatial momentum and compact fixed x intervals, the
complete canonical operator converges after time rescaling to

```
l_xx=0,
Q_xx+80/(1+8x^2)*Q=0.
```

This is coefficient and normalized-solution convergence on shrinking
physical intervals, not a fixed-duration EFT matching theorem. Although
the center frequency grows as `1/(tau*sqrt(delta))`, its variation occurs
on the same shrinking time scale. The exact center limits are
`(m^2)''/(m^2)^2 -> -1/5` and `m_frequency''/m_frequency^3 -> -1/10`;
the canonical diagonal has the same limits. Increasing this center
frequency does not make the usual time-derivative expansion parametric.

For the limiting forced equation `R''+80R/(1+8x^2)=1`, the exact
particular solution is `(1+8x^2)/96`. Algebraic inversion instead gives
`(1+8x^2)/80` and leaves residual 1/5. This checks a real omitted derivative,
but it is only a specified limiting-ODE control: it is not yet an actual
physical-metric retarded pulse, prepared state or finite-window error
bound. The locked light-cone excess at the center (and on fixed inner
intervals) is itself O(delta), not uniformly so on a fixed physical
interval. The relative order of the actual response correction must be
computed before making a light-only cone verdict.

## What remains on the route to P8 closure

The scoped P8(a) photon/cosmological-strength objective and the frozen
linear M0/M1 classification are unchanged. The original P8(b) question
still needs an action and physical matter-frame dictionary joining the
rolling EFT to an applicable vacuum/finite-gravity positivity calculation,
or a theorem excluding the remaining allowed class.

For the variable parent, the next concrete checks are a genuinely
physical sourced tensor response and the reduced scalar/vector system
on the same regular background. Any proposed light-only reduction must
carry its state dependence, frequency band and omitted-response error.
A smooth vacuum completion would be a separately specified action,
requiring its own spectrum, interactions, decoupling and matching bounds.
The finite-gravity Regge/IR and loop remainders are not assigned guessed
numbers. Failure of one reconstructed parent would not close a DHOST row.

Two older primary sources help route this calculation but do not discharge
it. The scalar/metric mass Hessians in
[Cusin, Khosravi and Noller, section 4.2 and appendices B/C](https://arxiv.org/pdf/1608.06643)
provide comparison formulas; that section explicitly leaves the full scalar
stability analysis to separate work. The no-ghost inequality in
[De Felice, Mukohyama and Uzan, equation (31)](https://arxiv.org/pdf/1702.04490)
is derived there for a constant-scalar de Sitter normal branch. Neither
that stationary result nor the paper's environment-dependent mass-scale
argument is a rolling-bounce stability or S6 common-parent certificate.

## Verification

The frozen A26 report has SHA256
`5c32597acb4f62ea4c2ad566480fd1be34fcaa29591c0ad771eab9d338067beb`.
Its 16 direct sources include 32 separately authored audit cases. It
checks 34 exact primary residuals and ten independent Fraction Laurent
fixtures, plus the separate Fourier/IBP replay. All 67 ordinary child
tests passed in 6.66 seconds; the standalone read-only replay, Ruff and
scoped diff check passed.

The frozen VARIABLE report has SHA256
`335cd52028baf56b30c75377aa7e6edd7ec86db2799f49b13f467242674fc4d0`.
Its 16 direct sources include 30 separately authored audit cases. It checks
64 exact primary residuals, six independent coefficient identities, seven
full-background and five canonical-second-jet Fraction fixtures, and 148
primary/independent output comparisons. All 55 ordinary child tests passed
in 121.71 seconds; its standalone read-only replay and Ruff passed.

All 32 direct source hashes match and all 65 local Markdown links in the
37-file checkpoint manifest resolve. The A26 tests and separate VARIABLE
audit also passed together: 97 tests in 8.31 seconds. All ten ordinary
regression-runner tests passed in 0.15 seconds. No ancestor certificate or
unrelated P4/P9 work is changed by these additions.

The combined P8 suite passed **2,276 tests in 418.20 seconds** with the
explicitly opt-in [exact regression runner](p8-exact-regression-runner.md).
Startup checked 128 original Gaussian-domain GCD tuples; the full run
recorded 6,509 exact descents and 2,562 domain fallbacks. This is an adapted
exact full-suite pass, not an ordinary unmodified SymPy run. The unchanged
long M1 nonlinear symbolic test emitted its configured 60-second diagnostic
traceback and subsequently passed. No partial output was counted as a pass.
