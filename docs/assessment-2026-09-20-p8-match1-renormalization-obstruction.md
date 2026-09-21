# MATCH-1: finite-normalization obstruction and four-value matching target

Date: 2026-09-20. Scope: the QG2-H8A420 branch of the
[closure-driven plan](p8-closure-plan.md), following the
[first matching decision audit](assessment-2026-09-20-p8-matching-decision-audit.md).

## Decision

The retained data alone give **CONDITIONAL ONLY**, not a physical exclusion.
They do not determine the finite-gravity matching combination

```text
C_pos = alpha/(8*pi^2*kappa^2)
      + 4*g*(n+2*mu)*c_RH/[kappa*(n-2*mu)^3].
```

Here `mu` and `n` are the light and heavy squared masses. This is the
specified two-direction projection, not the whole physical positivity
coefficient. The M gate and physical MATCH-1 enclosure remain OPEN.

Two concrete results sharpen the next task:

1. A covariant, order-hbar finite counterterm preserves the retained
   vacuum value normalization and reference-clock jets, but changes
   `C_pos`. Even a stringent clock-tube smallness condition permits both
   signs of a selected vacuum reference coefficient.
2. Four nonforward amplitude values suffice to extract this projection
   while eliminating the independent constant and Newton-shape anchors.
   Their uniform absolute error amplification is less than `15/mu^2` for
   `n >= 32*mu`. This remains bounded at the original enormous heavy mass.

Neither result supplies physical matching values, an omitted-order bound,
a gravitational positivity verdict, or a quantum bounce. No frozen S6
source or certificate is changed, and no new numbered S6 result is claimed.

## 1. What the existing prescription fixes

The relevant distinction is between a specified calculation and a fully
matched quantum parent:

| Retained input | What it establishes | What it does not establish |
|---|---|---|
| [S177](../problems/P8/s6/continuation/s6_177/FORMULATION.md) | Classical fixed-canonical gravity-decoupling family | Quantum gravitational decoupling or finite-gravity matching |
| [S182](../problems/P8/s6/continuation/s6_182/FORMULATION.md), [S238](../problems/P8/s6/continuation/s6_238/FORMULATION.md) | Specified Gaussian Proca retuning and full classical heavy parent | Complete interacting quantum counterfunctional |
| [S239 H8A420-VAC-OS4](../problems/P8/s6/continuation/s6_239/notes/renormalization.md) | First scalar loop of the classical limiting action; generated off-shell UV subtraction, selected derivative couplings, mass/residue and symmetric four-point value conditions | Finite-gravity derivative coefficients absent from that limiting calculation |
| [S240](../problems/P8/s6/continuation/s6_240/FORMULATION.md) | Specified heavy Gaussian state/prescription and QG2 reference profile | Full light-scalar curvature prescription or interacting quantum mean |
| [S285](../problems/P8/s6/continuation/s6_285/FORMULATION.md), [S290](../problems/P8/s6/continuation/s6_290/notes/endpoint.md) | Explicit light-curvature and endpoint matching coordinates | A physical value for the independent light-curvature or `c_RH` coordinates |
| [S302](../problems/P8/s6/continuation/s6_302/FORMULATION.md) | Known finite gravity reference plus explicit `alpha`, `beta`, `delta_kappa` | Values for those three independent finite anchors |
| [S336](../problems/P8/s6/continuation/s6_336/FORMULATION.md), [S347](../problems/P8/s6/continuation/s6_347/FORMULATION.md) | Selected known curvature contribution and its subsequent calculations | An assignment of the independent parent matching data |

Consequently, setting an unassigned finite coefficient to zero would add a
renormalization condition to the candidate; it would not be a deduction
from the existing loop calculation. Conversely, changing subtraction
coordinates while retaining the same physical matching data does not
change an observable. The comparison below varies an unspecified finite
datum, not merely the name of a fixed physical theory.

This is the usual EFT distinction between calculable nonlocal loop effects
and local coefficients requiring matching or normalization input; see
[Donoghue and Holstein, *Low Energy Theorems of Quantum Gravity from Effective Field Theory*](https://arxiv.org/html/1506.00946).
The specific obstruction below is derived from this repository's parent
and conditions, not inferred solely from that general observation.

## 2. An explicit surviving finite direction

Let `kappa0 = 10^800`, `Y = (nabla Phi)^2`, `X = Y/kappa0` and
`u = Phi/sqrt(kappa0)`. Use the already adopted S239 profile

```text
V(X) = (1-X)^1024 / [X^1024 + (1-X)^1024].
```

Consider the generally covariant action-density addition, with the volume
factor understood,

```text
Delta L_eta = hbar * eta/kappa^2 * V(Y/kappa0)
              * [Y^2 - mu^2*Phi^4/9].
```

`eta` is a real comparison parameter, not an assigned value in the original
candidate. These are alternative formal finite quantum completions of the
retained data; no UV realization of this family is asserted. Keep
`kappa0` fixed when taking the classical family's fixed-canonical
`kappa -> infinity` limit. The addition then vanishes, and it also vanishes
in the formal `hbar -> 0` classical action. Thus it does not alter the
derivative-coupling conditions of the already taken scalar limiting action.

### Vacuum normalization and on-shell effect

At the vacuum, `V = 1 + O(X^1024)`. Its localizing factor therefore does
not change the four-light vertex. There are no quadratic or heavy
one-point terms. At the displayed order in hbar, it changes no light
mass/residue, linear Newton pole, or already computed one-loop
discontinuity: putting this counterterm inside a loop is a higher order.

The literal 24-permutation vertex of `Y^2` is
`2*(s^2+t^2+u_M^2-4*mu^2)`, where `u_M` is the Mandelstam variable,
distinct from the clock coordinate `u`. The vertex of `mu^2*Phi^4` is
`24*mu^2`. After stripping the common hbar marker, the amplitude shift is

```text
Delta A_eta = 2*eta/kappa^2
              * [s^2+t^2+u_M^2 - 16*mu^2/3].
```

It vanishes at `s=t=u_M=4*mu/3`, so even imposing the retained symmetric
value condition at finite gravity does not remove this direction. In the
S302 coordinates,

```text
delta alpha = 32*pi^2*eta,
delta beta = -16*delta alpha/3,
delta(delta_kappa) = 0.
```

With `s=2*mu-t/2+v` and `u_M=2*mu-t/2-v`, its forward `v^2` coefficient is

```text
delta b20 = delta C_pos = 4*eta/kappa^2.
```

The nonzero on-shell polynomial is not removable by declaring this
direction redundant under an equation of motion. A field redefinition
must carry all induced couplings, source and measure changes; it cannot
remove a physical amplitude difference between these comparisons. For the
higher-order qualifications to equation-of-motion substitutions, see
[Criado and Perez-Victoria, *Field redefinitions in effective theories at higher orders*](https://arxiv.org/html/1811.09413).

### Clock jets and a quantitative tube bound

At `X=1`, the denominator of `V` is nonzero and its numerator has a zero
of order 1024. Thus every Taylor variation of the added density through
total order 1023 in the light and metric perturbations vanishes on the
exact reference clock. This is a local factorization statement, not a
statement about arbitrary fluctuating states or nonlinear evolution.

There is also a useful finite-neighborhood bound. At the original
`mu=1`, `kappa=kappa0`, the normalized scalar addition is

```text
Delta F = eta/kappa0 * V(X) * (X^2-u^4/9).
```

For real `|u|<=1`, `7/8<=X<=9/8`, use a joint complex radius `1/64`.
Then `|z|<=73/64`, `|u_complex|<=65/64`, and
`|(1-z)/z|<=9/55`. The denominator of `V` is nonzero there,
`|V|<2*(9/55)^1024`, and `|z^2-u_complex^4/9|<3`. Cauchy's inequalities
give, for `i+j<=4`,

```text
|partial_u^i partial_X^j Delta F|
  <= 6*|eta|/kappa0 * i!*j!*64^(i+j) * (9/55)^1024
  < 10^-590                         if |eta| <= 2*10^1000.
```

All 15 rational majorants are checked exactly. This is a bound on the
added coefficient function, not a bound on the full quantum stress,
response inverse, or distance between nonlinear bounce solutions.

### The small tube bound does not fix the selected vacuum sign

At the original parameters,

```text
lambda = 10^-600, D = 10^200/512, n = D+2, g = 1/8192,
b20_tree = 2*g^2/D^3 = 4*lambda.
```

S239's selected limiting-matter first-loop correction satisfies
`0 < delta b20_matter/(4*lambda) < 10^-203`. The two comparison values
`eta = +/- 2*10^1000` produce `delta b20 = +/- 8*lambda`. Hence
tree plus that selected matter correction plus this counterterm has
opposite signs in the two comparisons, while both obey the above tube
bound and the retained value/clock-jet conditions.

This does not sacrifice the simple contact-only potential control either.
With S238's `q = g^2*(D-1)/[6*D^2*(D+2)]`, the retained quartic lower
margin is greater than `q/2`, whereas
`|eta|/(9*kappa0^2) < q/4`. At `Y=0` the potential shift is
`eta*Phi^4/(9*kappa0^2)`, so this particular contact-only lower margin
remains positive. This is not a bound on the full effective potential.

The sign comparison is deliberately limited. It is not the sign of the
full physical finite-gravity coefficient; other contributions and the
gravitational dispersion allowance remain unknown. Large `eta` can also
reorganize loop power counting. No higher-loop smallness, positivity pass,
UV health, or quantum-bounce existence follows from either comparison.

The exact normalization conditions and reference-clock jets alone leave
`eta` arbitrary. Imposing the displayed finite tube budget restricts its
size but still permits both comparison signs. Thus neither the exact
conditions nor that smallness budget supplies the decision-relevant
physical enclosure sought here.

## 3. Four matched values isolate the needed combination

The next matching calculation need not separately fit two badly scaled
coefficients, nor fix the nuisance constant and Newton normalization to
zero. Set `x=n/mu >= 32`. Choose four physical, nonforward, 90-degree
massive two-body kinematic points, in units of `mu`:

```text
q0=(8,-2,-2), q1=(10,-3,-3), q2=(12,-4,-4), q3=(16,-6,-6).
```

They obey `s+t+u_M=4*mu` and avoid both forward poles and the heavy
resonance. For the original `mu=1` they lie in the S302 compact window.
This kinematic observation does not by itself turn S302's analytic soft
reference into a physical infrared-safe observable.

For a dimensionless triple `q=(a,b,c)`, define the cyclic sums

```text
S(q)   = a^2+b^2+c^2,
H(q;x) = sum_cyclic (a+2)/(x-a),
T(q)   = sum_cyclic (2-2*a-b*c)/a.
```

`T` is the dimensionless S288/S302 Newton shape. The fully crossed S290
endpoint, with `f2(a)=h/(n-a)` and `h=-2*g*c_RH`, gives the amplitude
term `(2*g*c_RH/kappa)*H`. This retains all three crossed channels.

In a common physical normalization and infrared prescription, let `r_i`
be the real matched residual at `q_i` after the specified known terms
are subtracted. The four-shape ansatz and its remaining error are

```text
r_i = A*S(q_i) + B*H(q_i;x) + C + D_N*T(q_i) + e_i,

A   = alpha*mu^2/(16*pi^2*kappa^2),
B   = 2*g*c_RH/kappa,
D_N = delta_kappa*mu/kappa^2.
```

`C` includes the independent constant anchor and constant improvement
contributions. The terms `e_i` must include every unmodeled sector,
higher-derivative structure, subtraction error and relevant physical
infrared conversion error. They are not presumed small.

Let the four rows of `M(x)` be `[S(q_i),H(q_i;x),1,T(q_i)]`, and define

```text
L(x) = [2, 2*(x+2)/(x-2)^3, 0, 0],
w^T  = L(x)*M(x)^(-1).
```

Then the desired projection is exactly

```text
C_pos = sum_i w_i*(r_i-e_i)/mu^2.
```

In particular, `sum w_i = sum w_i*T(q_i) = 0`. Neither nuisance anchor
has been assigned a physical value. The matrix is nonsingular throughout
the stated range:

```text
P(x) = 3659*x^4 + 309616*x^3 + 815324*x^2
       + 30852016*x - 93987840,

det M = -2*P(x) /
        [5*(x-16)*(x-12)*(x-10)*(x-8)*(x+3)*(x+4)*(x+6)].
```

The rational weights have signs `(+,-,+,-)` and obey the uniform bound
`sum_i |w_i| < 15` for every `x>=32`. These are exact half-line results:
substituting `x=32+z` makes the numerator and denominator coefficients of
the relevant sign-adjusted rational expressions strictly positive.
They are not conclusions drawn from a finite numerical grid.

Their limits as `x -> infinity` are

```text
(42447/14636, -24900/3659, 65923/14636, -4385/7318).
```

The projection therefore remains well conditioned even when a separate
fit of the heavy-pole coefficient would be poorly scaled. Use the exact
derived weights or controlled arithmetic, not a floating-point inverse
of the nearly degenerate raw matrix at the original huge `x`.

### Error threshold

If each matched residual, including the entire error just described, is
known to absolute accuracy `epsilon>0`, then

```text
|delta C_pos| < 15*epsilon/mu^2.
```

At `mu=1`, allocating `lambda` (one quarter of the tree coefficient)
to this matching error alone requires `epsilon <= lambda/15`. Equivalently,
residuals expressed in units of `lambda` require an absolute error at
most `1/15`. This is a proposed error allocation, not a result asserting
that such precision or the central values are available. Other terms in
the necessary test and its gravitational allowance still need their own
budgets.

## 4. Next input and stop rule

The next subtarget is **MATCH-1.PARENT-4**: obtain these four matched
residuals and a joint error bound from genuinely specified quantum-parent
input, or prove why a proposed parent cannot provide them. This is an
input contract, not a requirement for an all-orders UV construction.

The available routes differ materially:

| Route | Effect on MATCH-1 | Decision |
|---|---|---|
| Refine already known loop terms | Changes the known subtraction but leaves the demonstrated finite direction | Not a matching-resolution task |
| Set finite gravity coefficients to zero by convention | Defines additional candidate conditions without deriving physical values | Only as a separately named conditional candidate, never retroactively |
| Compute or constrain the four residuals from an explicitly specified parent/physical normalization | Can supply `C_pos` with the error bound above | The next physical matching target |
| Retain a parameterized renormalized family | Can state conditional tests and compare assumptions | Useful classification evidence; does not close the physical M gate |

The retained files do not provide the third route's missing input. Running
the existing EFT amplitude with free coefficients through this projector
would merely recover those same free coefficients. No arbitrary value is
requested from the user, and absence of input is not an all-parent no-go.
Before another numbered scientific sequence, a proposed additional parent
input must be named and justified, or a different route must be compared.

For the bounce gate, an on-shell four-value match is not sufficient. The
counterterm comparison already shows that identical reference-clock jets
can coexist with different vacuum data; conversely, four on-shell values
cannot determine the full curved/state-dependent counterfunctional.
Any selected parent must also carry its finite terms, matter/source
contacts and omitted-order control into the B calculation in the same
physical frame. M/V/G/B and original P8 remain OPEN.

## Reproduction and evidence level

Run the read-only diagnostic from the repository root:

```sh
.venv/bin/python scripts/p8_match1_audit.py
.venv/bin/ruff check scripts/p8_match1_audit.py
.venv/bin/ruff format --check scripts/p8_match1_audit.py
```

The [diagnostic](../scripts/p8_match1_audit.py) pins 11 parent reports and
their source manifests, 209 files total, before and after its calculation.
It imports no scientific parent package and writes no files. Original
SymPy checks the literal vertex, symmetric-value and forward identities,
the projector identity, determinant, weight signs and uniform error bound.
An independent `Fraction` polynomial ring reconstructs all 24 vertex
permutations; an independent rational Gauss-Jordan solve cross-checks all
four weights at `x=32,64,1000` and the original heavy mass, for 16 exact
weight comparisons. Nonzero nuisance coordinates and all 16 unit error
corners are also exercised. Clock and potential majorants use exact
rational arithmetic. The order-1024 jet statement is the factorization
proof above, not a claim of evaluating 1024 symbolic derivatives.

The new ledger entry is **VERIFIED_N** for this finite algebra and its
written scoped argument. This diagnostic is not a frozen physical matching
certificate, a formal proof-assistant development, or a replacement for
the separate long-running validation of earlier S6 results.
