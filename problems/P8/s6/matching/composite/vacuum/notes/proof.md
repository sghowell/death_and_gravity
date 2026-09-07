# S6.12: an exact physical-source-preserving vacuum branch

The proportional GR and massless physical-source results have explicit
[prior literature](https://arxiv.org/pdf/1409.3146), audited in
[sources.md](sources.md). The theorem below concerns a specified exact
branch, not generic nonlinear equivalence to GR.

## 1. One physical metric, a relative coordinate, and the domain

The physical metric is `mathcalG=g_eff=g(I+S)²`, `S=sqrt(g^-1 f)`, with
the positive real branch continuously connected to S=I. The equal unit
composite weights are fixed; no matter-frame change is made. Define

`Delta=(I-S)(I+S)^-1`.

The matrices S and Delta commute, since each is a rational function of
the other. Solving gives

`S=(I+Delta)^-1(I-Delta)`,

`g=mathcalG(I+Delta)²/4`, `f=mathcalG(I-Delta)²/4`.

These products have mixed-index Delta on the right of the covariant
metric. Require Delta to be mathcalG-self-adjoint, both I+-Delta
invertible, and g,f Lorentzian in the same real positive-root chart.
Near Delta=0 this is an ordinary invertible algebraic field map. More
generally the claimed domain is precisely the connected chart where
these properties and the selected square root hold. A Lorentz-self-adjoint
matrix need not have arbitrary globally admissible roots; no such claim
is used. A nonzero det mathcalG alone is insufficient.

The inverse formulas verify `g S²=f` and `g(I+S)²=mathcalG`. Exchange
g and f sends S to S^-1 and Delta to -Delta, leaving mathcalG unchanged.
The proof is not restricted to simultaneous diagonal metrics; the root
audit includes non-diagonal Lorentz congruences.

## 2. Full determinant, sources and exchange parity

Write `E_n=e_n(Delta)`, E0=1. The exact generating identity is

`sqrt|g| det(I+tS)/sqrt|mathcalG|`

`=det[(I+Delta)+t(I-Delta)]/16`

`=sum_(n=0)^4 E_n(1+t)^(4-n)(1-t)^n/16`.

This characteristic-polynomial identity holds for general matrices.
Its coefficient of t^n is the original potential density
`sqrt|g| e_n(S)/sqrt|mathcalG|`, and summing its coefficients gives 1.
Thus every prescribed light/source action depending on mathcalG alone
is exactly independent of Delta, including its physical volume factor.

Equal Einstein coefficients make the two complete EH actions invariant
under exchange. The potential is invariant when constant
`beta_n=beta_(4-n)`. Appropriate identical gravitational boundary terms
are exchanged as well. Alternatively use compactly supported variations.
For either choice the full action obeys

`S[mathcalG,Delta,psi;J]=S[mathcalG,-Delta,psi;J]`.

There is no assumption that light fields or physical sources are static.
All derivatives of the algebraic metric map are retained inside this
identity; terms linear in Delta or its derivatives cancel between the
exchanged actions, even on a time-dependent light metric.

For bare beta2=1 with all other beta_n zero, the potential is exactly

`-M²m² sqrt|mathcalG| [6-2E2+6E4]/16`

`=-M²m² sqrt|mathcalG| [6+tr(Delta²)-(tr Delta)²+6det Delta]/16`.

It is an exact quartic expression, not a TT-only mass fit. In particular
the constant physical-volume term has no hidden relative tadpole.

An external source remains the literal prescribed `Sprobe[mathcalG;J]`.
At a physical Minkowski vacuum and for one unit TT polarization,
`delta mathcalG_ij=-gamma eij`, `eij eij=1`. With convention
`delta Sprobe=(1/2) integral T^{ij} delta mathcalG_ij`, its term is
`-Pi_TT gamma/2`. Our normalized `j gamma/4` therefore means
`j=-2Pi_TT`, in physical Minkowski time and volume. It is not an
identification j=Pi_TT. No source couples directly to Delta in these
coordinates. The metric and all light actions, not just a TT coupling,
are preserved.

## 3. Exact branch, classical elimination and its limits

Evenness implies the relative Euler functional satisfies

`E_Delta[mathcalG,0,psi;J]=0`

for every light history in the chart, whether or not its light equations
hold. Thus Delta=0 is an exact off-shell stationary branch. At Delta=0,
g=f=mathcalG/4. In four dimensions constant scaling gives
`sqrt|g|R[g]=sqrt|mathcalG|R[mathcalG]/4`; the two EH terms therefore
reduce to `-M² sqrt|mathcalG|R[mathcalG]/4`. Consequently

`S_tree[mathcalG,psi;J]=S[mathcalG,0,psi;J]`

`=integral sqrt|mathcalG|[-M²R[mathcalG]/4-3M²m²/8+L_light]+Sprobe`.

Differentiating the stationary substitution yields the full light
equations because the chain-rule term multiplies E_Delta=0. This is
not the operation of solving the g equation for f and treating that as
massive-mode integration. Nor was a retarded inverse inserted into a
symmetric single-copy action: the selected branch is explicitly zero,
so its substitution is local without a Green-function manipulation.

For the canonical-light action defined in FORMULATION, no C/D DHOST or
curvature-squared operator is induced at tree level, to any derivative
order. If such an operator is already inserted into L_light, the
identity retains it; it does not prove its absence. Pure light GR/matter
tree interactions also remain. The assertion concerns new operators
generated by the relative sector on this branch.

A diagrammatic check gives the same tree statement. Expand about
Delta=0. Every vertex with relative fields has positive even relative
degree. With no external relative legs, any nonempty connected heavy
subgraph has E>=V, hence cycle number E-V+1>=1. It cannot contribute
a tree graph. A loop is allowed: a light-dependent quadratic heavy
operator can produce a nonconstant log determinant. No loop cancellation
or radiative closure is implied by parity.

The identity establishes a branch, not automatically the unique
physical elimination prescription. Nonzero relative Cauchy data,
boundary sources, coherent states or other stationary saddles are not
set to zero by symmetry. In particular the pinned asymmetric composite
bounces have Delta!=0; they are not contradicted by the exact branch.

## 4. The precise functional uniqueness and continuation theorem

There are two distinct valid formulations.

**Causal formulation.** Fix the common gauge, solve the relative
constraints with their rank preserved, and specify an admissible
well-posed relative initial-value problem. Zero relative fields and
physical canonical data, including the necessary auxiliary constraints,
select Delta=0 by uniqueness, since it solves the complete homogeneous
relative problem. This remains true for dynamical light fields and
sources. Neither constraint regularity nor nonlinear well-posedness on
an arbitrary rolling CD domain is asserted in this gate.

**Stationary functional formulation.** Let lambda vary on a connected
parameter domain of prescribed light configurations and source data.
After fixing gauge/constraints and a common zero-relative boundary or
state prescription, let `E(lambda,h):X->Y` be a continuously differentiable
map between specified Banach spaces of relative histories/residuals.
Require `E(lambda,0)=0`, a continuous solution branch h_*(lambda), and
`D_h E(lambda,h_*)` to be an isomorphism wherever continuation is
claimed. Suppose h_*(lambda0)=0 at a reference vacuum point. Then the
set where h_*=0 is closed by continuity. It is open by the implicit
function theorem: near each such point the nearby solution is unique,
and h=0 is already a solution for the nearby light configurations.
Connectedness therefore implies h_*=0 throughout that domain.

This is a statement about a selected solution functional in
configuration space. It does not assume or require a physical-time
history connecting the vacuum and bounce backgrounds. A loss of chart
regularity, operator invertibility, branch continuity, common boundary
conditions or the zero-data prescription invalidates the premise, not
the original DHOST witness.

For quantitative control, an algebraic mass is not enough. One needs
the full differential operator with its temporal, spatial, constraint
and boundary structure. For example `h''+h=0` has strictly positive
algebraic mass squared 1, yet with h(0)=h(pi)=0 the nonzero sin(t)
solution is in the kernel. The same equation with h(0)=h'(0)=0 has a
unique zero causal solution. A source or nonzero initial derivative
selects a different solution. The finite-dimensional branch equation
`h(h²-lambda)=0` similarly permits branch joining precisely where its
linearization at h=0 ceases to be invertible.

If a *supplied* constrained inverse has norm at most K, `N(0)=0`, and in a stated
ball `||L^-1(N(h)-N(k))||<=eta||h-k||`, eta<1, then a residual problem
`Lh+N(h)=r` obeys

`||h||<=K||r||/(1-eta)`.

The contraction additionally maps the ball into itself when
`K||r||<=(1-eta)*radius`. This includes appropriately normalized
boundary/data residuals in Y, not just a local potential derivative.
The exact API replays this arithmetic and rejects inexact inputs. Its
fixture K=2, eta=1/4, r=3/100 yields 2/25; these are demonstration
constants, **not** measured bounds for the bimetric CD operator. No
rolling heavy-gap, inverse-norm or nonlinear remainder estimate is
silently supplied by this theorem.

## 5. Conditional proportional Minkowski vacuum and source residues

For one fixed light potential to realize a proportional Minkowski
vacuum with constant scalar fields, its value must be

`V_v=-3M²m²/8`, `V_phi(phi_v)=0`.

Canonical clock masses are its canonically normalized Hessian. Require
them nonnegative (strictly positive for the massive-clock vacuum test);
the exact free chi is massless. The constant potential contributes
`-V_v sqrt|mathcalG|=-V_v sqrt|g| sum e_n(S)`. Thus it shifts **all**
five beta_n, giving

`beta_eff=(-3/8,-3/8,5/8,-3/8,-3/8)`.

Both proportional Minkowski lapse equations vanish with these shifted
coefficients. Applying the standard relative mass formula in the g
proper clock first gives

`m_FP,g²=2m²(beta_eff,1+2beta_eff,2+beta_eff,3)=m²`.

But mathcalG=4g, so physical proper time is twice the g proper time.
The physical mass is consequently `m_FP²=m²/4`, not m². Using the bare
beta2 coefficients also gives the incorrect physical value m²: it
misses the vacuum source shifts and is not an independent mass check.

A second route uses the exact Cayley FP invariant from section 2.
At quadratic TT order `h_g=gamma+2Delta`, `h_f=gamma-2Delta`. The
physical-time kinetic coefficients are M²/32 for each original
dimensionless TT field. Their diagonal form has

`M_eff²=M²/2`, `M_relative²=2M²`, `m_FP²=m²/4`,

with action after integration by parts

`-[M_eff² D gamma² + M_relative²(D+m_FP²)Delta²]/8 + j gamma/4`,

`D=partial_T²+k_phys²`.

The exact quadratic mass invariant is
`-M²m²[tr(Delta²)-(tr Delta)²]/16`, the full Fierz–Pauli form, not
only its TT restriction. The usual divergence constraint, then its
trace, impose transverse tracelessness on source-free massive modes;
the timelike rest-frame polarization space has five positive-norm
spatial traceless components. The root audit derives these full
constraints independently. Both kinetic coefficients are positive.
The zero heavy physical-source residue must not be confused with a
zero mass, a negative kinetic residue, or absence of the five modes.

Stationarity gives `gamma=j/(M_eff²D)`, Delta=0 for this physical source,
and its stationary quadratic functional is `j²/(8M_eff²D)`. For a
general conserved source, the massless trace projector is
`(P2-P0/2)/(M_eff²D)`, not the massive `P2/(D+m_FP²)` projector. A
hypothetical independent relative diagnostic source would see positive
residue `1/(2M²)` and the physical massive pole. It is not an additional
physical matter coupling in the model. The full scalar/constraint
projector and source absence are checked separately from the TT mass.

This is conditional flat quadratic spectrum information. It is not
full-parent nonlinear or quantum health, and it does not establish that
the already reconstructed local rolling potential has the required
vacuum elsewhere. One could choose a separate potential with this
minimum as a vacuum example, but that would not prove matching to the
pinned bounce.

Indeed, the pinned one-scalar family at its initial centre has

`V_b/(M²m²)=9epsilon²/2-epsilon/(1+epsilon)²`.

For `0<epsilon<=1/10000`, `V_b/(M²m²)+3/8>=3749/10000>0`.
The constant-vacuum choice cannot reproduce even this initial value.
A smooth same-potential extension may or may not satisfy both the
rolling profile and vacuum conditions; constructing and controlling
that one action is an unfulfilled matching obligation. The symmetry
identity above holds for any one such fixed potential, if supplied.

## 6. Distinct matter content and direct old CD/M1 comparison

The pinned composite parent has one internal scalar. The old CD/M1
system has a DHOST clock and an additional exactly free, rolling
canonical chi. They are not the same light content or action.

For a direct M1 test define a separate canonical-light extension:

`L_light=Z_AB(psi) partial psi^A.partial psi^B/2 -V(psi)`

`+ (partial chi)²/2`,

with positive Z and no chi dependence in Z or V and no mixed chi
kinetic term. One clock psi suffices for two light scalar fields. The
parity and tree-branch arguments apply unchanged, but a rolling
solution of this extension is a new background problem. This gate
does not import the old one-scalar composite solution into it.

On Delta=0 the physical Einstein coefficient is constant and there
are no second derivatives of a clock in the canonical input matter.
In the frozen physical metric and a stated nondifferential clock
dictionary, F2_X=A3=0. The old CD/M1 certificate instead gives, in
its exact M=tau=1 normalization,

`F2=-1/2+(1-X)/(2(1+phi²)³)`, `A3=1/[(1+phi²)³ X]`.

At `(phi,X)=(0,1)`, these jets are -1/2 and 1. Matching the target
as this branch plus a remainder therefore requires

`|R_(F2X)(0,1)|>=1/2`, `|R_(A3)(0,1)|>=1`.

These are budgets in the prescribed operator/clock normalization, not
field-redefinition-invariant observables. Point scalar reparametrizations
preserve the zero branch values; arbitrary derivative redefinitions,
metric changes or EOM trades require a separate dictionary and error
accounting. The physical-metric obstruction below does not depend on
this operator-coordinate comparison.

Match `M_target²=M_eff²=M²/2`. For flat homogeneous backgrounds of
the zero-relative branch with the canonical-light extension,

`-2M_target² Hdot = Z_AB psidot^A psidot^B + chidot² >=0`.

The potential and constant vacuum energy drop out of this null
equation. There are no additional NEC-violating background sources.
Hence H is nonincreasing. The original physical CD target has

`H_CD(T)=4T/(tau²+T²)`,

`H_CD(-tau/2)=-8/(5tau)`, `H_CD(tau/2)=8/(5tau)`.

If both endpoint errors were strictly less than 8/(5tau), the candidate
H would be negative at the first endpoint and positive at the second,
contradicting monotonicity. Equivalently the maximum endpoint error
must be at least 8/(5tau). No Hdot error assumption is needed. Static
Minkowski saturates the boundary for nonrolling matter, furnishing a
sharp scope control; it does not satisfy the old nonzero M1 charge.

For an actual equation-level correction define

`-2M_target² Hdot = n_clock + chidot² + R_null`, `n_clock>=0`.

The exact old charge restores as
`chidot=M_target/[10tau(1+(T/tau)²)^6]`. At the exact CD bounce,

`R_null(0) = -(801/100)M_target²/tau² -n_clock(0)`.

Thus its magnitude is at least `(801/100)M_target²/tau²`. This residual
must include all omitted heavy-state, operator, loop or external
background effects before using the bound. Dropping the free chi
understates the required normalized magnitude by 1/100.

If instead `|tau² Hdot(0)-4|<=epsilon_H<4` and
`|tau chidot(0)/M_target-1/10|<=epsilon_chi<1/10`, then

`tau² |R_null(0)|/M_target²`

`>=2(4-epsilon_H)+(1/10-epsilon_chi)²`.

These are necessary error budgets, not bounds computed for quantum
corrections. They exclude a small-error match of this branch under the
stated hypotheses, not every state or branch of the parent.

## 7. What this gate settles, and what remains research

The source-preserving zero-relative tree branch is exactly GR plus
its input canonical light matter. It cannot reproduce the old CD/M1
action or physical bounce with the stated small errors. A continuous
stationary branch attached to that vacuum cannot evade the result
while all the explicitly stated functional uniqueness hypotheses hold.
This is a direct named-branch matching obstruction, not another
particular-centre tensor-frequency screen.

Nonzero relative backgrounds, other boundary/state choices, loss of
functional invertibility, symmetry-breaking actions, loops, different
parents and explicit higher-operator light actions are not excluded.
The old asymmetric composite solutions remain valid within their
certified scope. Their possible matching requires its own complete
field/source/state dictionary and errors; a different background is
not automatically the zero-relative vacuum EFT branch.

Nothing here establishes the full parent cutoff, an all-scale absence
of the BD mode, a common potential realizing vacuum and bounce, a
rolling heavy functional inverse, Regge/IR dispersion remainders, or
UV completion. S6 does not require a physical-time vacuum-to-bounce
history, and none was used. The original DHOST row verdicts are not
overwritten, and neither S6 nor P8 is declared closed.
