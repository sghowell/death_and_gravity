# S6.13: full auxiliary-parent equations before matching

## 1. Frozen action and exact variational convention

The action and units are fixed in FORMULATION. Write
`d=det u`, `U=u^-1`, `Q=p_g e+p_f v`, `s=B+tr(UQ)`. All inverse
formulas below refer to regular real oriented vierbeine. The stronger
positive-root and symmetrization restrictions of the parent metric
formulation are retained wherever that interpretation is used.

Define the Frobenius gradients by

`delta S=integral tr(E_e^T delta e+E_v^T delta v+E_u^T delta u)`

plus the matter-field variation and gravitational boundary terms.
Use compact variations or appropriate EH boundary terms. For general
spacetime-dependent fields the complete gradients are

`E_e=G det(e) eta e Einstein[g]^{up,up} -2p_g d U^T`,

`E_v=F det(v) eta v Einstein[f]^{up,up} -2p_f d U^T`,

`E_u=-2d[s U^T-U^T Q^T U^T]+epsilon J_u`,

`J_u=delta L_m/delta u`.

The scalar Euler functional is `epsilon delta S_m/delta psi`. In the
chosen model there is no u EH term and no direct e/v matter term.
The EH formulas follow from the metric variation
`delta S_EH=(G/2) integral sqrt|g| Einstein[g]^{mu nu} delta g_mu nu`
and `delta g=delta e^T eta e+e^T eta delta e`, in the repository
curvature convention. They are not an omission of the kinetic action.

The interaction gradients follow from
`delta det u=d tr(U delta u)` and
`delta tr(UQ)=-tr(U delta u UQ)+tr(U delta Q)`.
No simultaneous diagonalization or Euclidean commutativity is assumed.
The independent Fraction engine differentiates all 48 e/v/u components
at noncommuting ambient matrix fixtures. The root audit separately
varies the literal action in noncommuting directions and verifies the
EH coframe dictionary. These fixtures test algebra, not solutions of
the parent Lorentz constraints.

For comparison with the source's inverse-vierbein convention define

`T=u J_u^T u/(2d)`, `t=tr(UT)`, `calT=(t/3)u-T`.

Then let `C=Q+(B/3)u-epsilon calT`. The full identity is

`E_u=2d U^T[C^T-tr(UC)u^T]U^T`.

Consequently `E_u=0` is equivalent to `C=0`: multiply by u^T on each
side, transpose, and trace with U. The trace of
`C-tr(UC)u` is `-3tr(UC)`, so the trace vanishes and then C=0.
This equivalence divides neither B nor a source amplitude. T is the
explicitly defined coframe-source matrix, not an unexplained reuse of
a physical stress convention.

## 2. Exact no-flat prerequisite theorem

For constant e,v,u the metrics g,f,h are constant Lorentzian matrices
and both Einstein tensors vanish. The equations reduce in particular to

`E_e=-2p_g d U^T`, `E_v=-2p_f d U^T`.

Since u is invertible, each is an invertible nonzero matrix whenever
its corresponding link is nonzero. Equivalently,

`u^T E_e=-2p_g d I`, `det E_e=16p_g^4 d^3`,

and likewise for f. Thus if either link is nonzero there is no regular
constant flat triple satisfying the full equations. This includes the
regular proportional constant vacuum needed for common Lorentz
invariance, but the algebra actually needs no proportionality assumption.
The signs of the links and B do not affect this contradiction.

Nothing was divided by B, Q, a Hubble rate or the auxiliary equation.
Any matter gradient depending only on u leaves these two equations
unchanged. Therefore neither a choice of constant scalar potential nor
a nonlinear source treatment can fix them within the stated action.

For a literal constant scalar potential,
`L_m=-V det u`, `J_u=-V d U^T`. The full action and equations simply
replace `B` by `B_eff=B+epsilon V/2`. The e/f equations remain exactly
the same. The convenient fixture e=v=u=I, p_g=p_f=1, B=-6 already has
E_u=0 but E_e=E_v=-2I. Checking only the eliminated auxiliary equation
would incorrectly accept it as a flat vacuum.

There are two essential scope boundaries. A singular u is not a regular
counterexample: it loses the metric/inverse used in the action. If both
links vanish, the theorem's premise fails; the dynamical metrics can
be free Einstein theories, and tuning B_eff=0 leaves u undetermined
rather than supplying the required invertible auxiliary construction.

This is a failure of a selected S6 prerequisite for this particular
parent. It is not an exclusion of curved vacua, of every theory discussed
by the source, of other trimetric potentials, or of all consistent
cosmologies. No physical-time vacuum-to-bounce trajectory is required
or used. The source's ghost-free degree-of-freedom construction is not
being contradicted by the additional Minkowski requirement.

## 3. Exact elimination, geometric coefficients and a false truncation

When B!=0 and Q is invertible, the source-free equation gives
`u0=-3Q/B`. Substitution into the full density gives

`L_eff,int=54 det Q/B^3`,

or potential `V_eff=-54 det Q/B^3` with action convention `-V_eff`.
The leading physical metric is `h0=u0^T eta u0`, with volume
`sqrt|h0|=81 det Q/B^4` on the oriented chart. Its weights are
`a_actual=-3p_g/B`, `b_actual=-3p_f/B`. The source paper defines
a,b with the opposite common sign for its metric-square expression;
the actual u0 is kept when differentiating matter.

For p_g!=0 the effective potential is

`V_eff=det(e) sum beta_n e_n(e^-1 v)`,

`beta_n=C r^n`, `C=-54p_g^4/B^3`, `r=p_f/p_g=b/a`.

The equivalent polynomial formula
`beta_n=-54p_g^(4-n)p_f^n/B^3` remains meaningful at a zero link.
These coefficients are not five independent HR parameters. In
particular `beta1 beta3-beta2^2=0`; the old bare beta2-only parent
instead gives -1, and its S6.12 shifted-vacuum coefficients give -1/4.
Thus this chosen auxiliary model cannot simply be imported as that
parent. Adding only separate beta0/beta4 endpoint terms also does not
alter this middle-coefficient identity.

For positive ratio r and positive proportionality y, the geometric
flat g/f lapse polynomials are `C(1+ry)^3` and `Cr(1+ry)^3`.
They cannot vanish for C!=0. This is a secondary proportional check;
the undivided full-action theorem also covers excluded denominators
where this eliminated formula is unavailable.

Constant matter vacuum energy can be resummed exactly, not guessed
from a first-order action. Whenever B_eff!=0 the exact solution is
`u(V)=-3Q/B_eff`, and

`L_eff(V)=54 det Q/B_eff^3`.

Relative to the fixed leading h0 volume this is

`L_eff(V)/sqrt|h0|=(2B/3)[B/(B+epsilon V/2)]^3`.

Its first-order expansion is `2B/3-epsilon V`. Choosing
`epsilon V=2B/3` would cancel that truncated expression, including
the entire geometric potential; treating this cancelled action as a
two-Einstein flat theory would give no relative FP interaction. But
it is not the exact parent at that point. The full coefficient equals
`9B/32`, not zero, and the exact u equals `3u0/4`. Moreover
`epsilon V/B=2/3` does not tend to zero in that cancellation prescription.
A small-source expansion cannot promote it to exact matching.

This distinction prevents both wrong conclusions: neither a positive
FP mass computed around a nonsolution nor a zero FP mass in an
incorrectly cancelled truncation establishes the chosen full parent's
Minkowski spectrum. The full no-flat equations must be checked first.

## 4. Actual source metric and Lorentz constraint

For the canonical scalar let k_mu=partial_mu psi, `Y=k^T h^-1 k` and
`w=h^-1 k`. Direct matter variation gives

`J_u=d[(Y/2-V)U^T-eta u w w^T]`,

`t=Y/2-2V`,

`calT=u[-(Y/12+V/6)I+(w k^T)/2]`.

At fixed e,v,k and B!=0, with regular u0, the trace-reduced equation
has u-Jacobian `(B/3)I` at epsilon=0. The finite-dimensional analytic
implicit-function theorem therefore gives a local auxiliary branch
for analytic canonical matter on this chart. No explicit finite source
radius, cutoff, or rolling error bound is claimed. Its first derivative
is fixed exactly:

`u=u0+epsilon u1+O(epsilon^2)`, `u1=(3/B)calT(u0)`.

The actual metric is h=u^T eta u, not h0. Its first coefficient is

`h1=u0^T eta u1+u1^T eta u0`

`=[3kk^T-(Y0/2+V)h0]/B`.

This expression is derived from the definition and full auxiliary
equation, not by importing a printed metric-correction prefactor.
With constant V and k=0 it agrees with the exact rescaling
`h=[B/B_eff]^2 h0` to first order. Generally it has a source-dependent
rank-one part as well as a conformal part. Identifying h with h0 while
claiming the exact prescribed matter frame would omit that correction.

For the following pair of Lorentz constraints assume both links are
nonzero. The e/v Lorentz equations imply, on the original parent chart,

`e^T eta u-u^T eta e=0`, `v^T eta u-u^T eta v=0`.

The separate EH and beta4 densities are individually Lorentz invariant,
and matter depends on u, so the nonzero link is what supplies these
constraints in either the chosen model or the named extension.
Substituting the auxiliary expansion into the first constraint gives

`p_f(e^T eta v-v^T eta e)`

`=epsilon[e^T eta calT(u0)-calT(u0)^T eta e]+O(epsilon^2)`.

In vacuum the e-v square-root symmetry is recovered. In a sourced
solution it generally receives a first-order correction. For the
off-shell dictionary fixture e=diag(1,2,3,4), v=I, B=-6,
p_g=p_f=1, V=0, k=(1,1,0,0), the leading e-v antisymmetric part
vanishes and Y0=5/9>0, but the required epsilon coefficient in its
01 component is -1/6. Also h1_01=-1/2. This does not purport to be a
complete background solution; it directly falsifies keeping the old
unsourced Lorentz condition and frozen h0 under a generic source.

The paper's small-source vierbein construction therefore is not
silently identified with the old metric-only composite action. On
special commuting backgrounds these differences can vanish, but
their perturbations and constraints still require their own audit.
No first-order source calculation here supplies quantitative original
free-M1 or CD operator matching.

## 5. A real separately named beta4 countercontrol

Now explicitly **change** the action by adding

`-2 beta4g det e-2 beta4f det v`.

These are the source's omitted n=4 link terms, independent of u,
and therefore add `-2beta4g det(e)e^-T` and its f counterpart to
the dynamical metric equations without changing the auxiliary equation.
Choose q>0,

`p_g=p_f=q`, `B=-6q`, `beta4g=beta4f=-q`, `V=0`.

At e=v=u=I all three Euler matrices and the total vacuum density
vanish. A constant V may also be included by choosing
`B=-6q-epsilon V/2`. Thus the original theorem cannot exclude this
extension. Omitting either compensating term restores a nonzero
corresponding flat equation.

For V=0 its exact auxiliary solution is `u0=(e+v)/2`. The reduced
interaction density is

`-q det(e+v)/4+2q(det e+det v)`.

To verify the full FP invariant, write e=r(I+Delta),
v=r(I-Delta), where `r=(e+v)/2` is the physical auxiliary coframe
in the source-free problem. On the symmetric local chart,

`L_int=4q det(r)[e2(Delta)+e4(Delta)]`.

Its quadratic part is exactly

`-2q det(r)[tr(Delta^2)-(tr Delta)^2]`.

The quartic determinant is retained and the computation is not only
a TT fit. At the flat vacuum the two TT fields h_g,h_f have quadratic
action after integration by parts

`-[G D h_g^2+F D h_f^2+q(h_g-h_f)^2]/8`,

`D=partial_T^2+k_phys^2`, in the actual h=eta clock.

The kinetic and spring matrices are diag(G,F) and
`q[[1,-1],[-1,1]]`. The generalized eigenvalues are 0 and
`m_FP^2=q(G+F)/(GF)>0` for G,F,q>0. The positive relative kinetic
coefficient in the coordinate h_g-h_f is `GF/(G+F)`.
For G=F the physical common mode gamma=(h_g+h_f)/2 and the relative
coordinate Delta_TT=(h_g-h_f)/4 have coefficients 2G and 8G;
the mass is `2q/G`. For G!=F the physical common source coordinate
is not itself the massless eigenmode; source rediagonalization would
be required. No source-decoupling result is borrowed across that case.

This flat quadratic countercontrol does not establish a controlled
rolling branch, scalar/matter stability, a cutoff, quantum positivity,
or original DHOST matching. The full constraint and EFT analyses for
such an extension are separate work. Its purpose here is to keep the
chosen-parent exclusion logically sharp while identifying an actual
surviving vacuum variation of the action.

## 6. Verification and remaining boundary

The replay checks literal diagonal Euler derivatives, the full arbitrary
source trace identity, canonical source and physical-metric coefficients,
exact elimination and geometric minors, the resummed cancellation
control, all extended flat equations, and the full FP invariant/mass.
The independent engine uses exact Fraction first jets in all 64
e/v/u/canonical-source matrix directions and a separate sparse
polynomial ring. Root-owned covariant audits use additional routes.
All sources/tests/proofs are hash-pinned; ancestors are replayed without
writes. The analytic proof and finite symbolic checks are not labelled
proof-assistant formalization.

What is settled is a concrete prerequisite failure of the stated
unextended auxiliary parent, together with a genuine extension outside
its hypotheses. No claim closes every parent class, old row, S6 or P8.
