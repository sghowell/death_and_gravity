# SOURCE-1: source completion feasibility and the unchanged-action EFT route

Date: 2026-09-28. The user approved **SOURCE-1**, a bounded feasibility
study. This is not approval to replace the physical action. The
[study specification](candidates/p8-source-1-study.json) records that
distinction; all frozen parents, prepared states and prior evidence remain
unchanged. The already approved COVZERO first-order boundary rule remains
approved. It does not need to be selected again.

**Decision:** recommend a separately named
**QG2-H8A420-RATE4-COVZERO-S0** candidate, replacing the entire retained
source square `(W-S)^2/2` by `W^2/2`, with everything else held fixed.
It passes this bounded classical admission study. It has NOT been adopted,
and its full common quantum reference has NOT been constructed.

The study establishes three useful distinctions:

- The original source obstruction survives an independent constraint-chain
  audit. At nonzero null gradients there is a second mechanism: a surviving
  scalar primary fails to commute with the temporal-vector primary.
- S0 and a first-derivative-only covariant source family avoid both failures.
  They preserve the exact classical clock and its quadratic action. A
  restricted all-tilt result rules out repairing the old Hessian-linear,
  gradient-parallel source merely by changing its two coefficients.
- The old source has no source-tagged flat-vacuum four-scalar diagrams
  through one loop. Removing it therefore preserves that **formal
  perturbative projection**, but changes the curved nonlinear response.
  An unchanged-action order-reduced EFT is still a possible alternative,
  not an admitted common prescription or a no-go theorem.

Original M/V/B/R/P8 remain OPEN and G remains UNTESTED. Completing this
feasibility study does not close those gates.

## 1. Admission question and fixed data

The missing item is the classical constraint prerequisite to the proposed
common vacuum/bounce reference, not a smaller known spectator loop. Success
means a concrete source candidate with a valid local constraint chain,
preserved classical data, and an explicit list of changed quantum work.
A scoped obstruction or a controlled EFT alternative would also answer the
study. Absence of a global secondary-pivot bound must remain visible.

Retain the physical metric/signature +---, X=du.du, the entire S238/S240
functions and fixed profiles, the original regular lower dictionary and
both canonical boundary momentum shifts, zeta=10^-6 and vector mass 1000,
the scalar/heavy masses, the four original RATE4 targets and detector 1/256.
R below is the tensor coefficient, not the Ricci scalar. Its curvature
density is `-R*Ricci_P8/2`, as corrected in the preceding packet.

The literal source is the one in the
[previous obstruction](assessment-2026-09-24-p8-covzero-source-constraint-obstruction.md):

```text
S_mu=u_mu (R-1)[-3H_clock+3R_u/(4R)+Box(u)/X
                         +(-1/X^2+3R_X/(2RX))*uHu],
uHu=u^mu u_munu u^nu.
```

No timelike-only restriction, numerical cutoff, new state or rule deleting
the extra solutions is assumed. These would be additional physical data.

## 2. Independent constraint-chain audit

### 2.1 Regular tilted and spacelike branch

The earlier covariant calculation establishes a rank-seven scalar/metric
velocity block where its determinant is nonzero, while the source-free
block has rank six. Its temporal-vector Hamiltonian pivot is `b^2/X`,
nonzero on that rank-lifted branch. Thus the usual vector primary has a
secondary that solves its temporal component; it does not impose an
extra scalar constraint. Minimal matter adds regular velocity blocks.

For clarity, after eliminating the spatial auxiliary-gradient pairs,
the first-order configuration variables are

```text
metric h_ij:6; lapse and shift:4; u and A_*=n.du:2;
W_mu:4; M1 and heavy scalar:2. Total:18.
```

The four diffeomorphism generators and their lapse/shift primaries give
eight first-class constraints. On the regular rank-lifted branch the
vector pair supplies two second-class constraints, with no scalar primary
left. The generic local Dirac count is therefore `18-8-2/2=9`, rather than
the intended eight. Preservation fixes the vector multiplier; it does not
start another chain. This is a local regular-rank count, not a theorem
about a nonlocal boundary prescription, a frequency of an on-shell ghost,
or an instability on the prepared bounce.

The original timelike, scalar-aligned clock calculation remains valid.
It cannot be used across the vacuum by treating scalar-unitary gauge as
regular at du=0. Nor does a result for a specific timelike shadow-mode
model supply the missing spacelike/null prescription here. The published
shadow-mode analysis explicitly distinguishes that domain issue.
See [De Felice, Mukohyama and Takahashi](https://arxiv.org/html/2110.03194).

### 2.2 New null-gradient check: a primary-bracket failure

Take the regular limit X=0 **before** constructing the Hessian. Define
`A(u)=[(R-1)/X^2]_(X=0)`. The actual retained parent's formula gives

```text
A(u)=1024[1-(1+u^2)^(-3)],
R=1, R_X=0, A3=2A, A4=-2A, A5=0 at X=0.
```

The X=0 values of A3 and A4 are not zero just because R_X is zero.
Set `u_mu=(a,a,0,0)`, a!=0, `Z=a V+a^2 K11`. The principal pure-light
density and source scalar are independently obtained as

```text
L0=(Kij Kij-K^2)/2+2A[(V+a K)*a Z-Z^2],
q=-A a Z,                 S_mu=u_mu q + velocity-independent terms.
```

A null velocity, normalized to Z=1, is

```text
V=(1-2A a^4)/a, K11=2A a^2, all other Kij=0.
```

Suppress the common positive kappa*volume density here. The metric
momentum conjugate to h_ij enters the null contraction with the factor
two relating dot(h_ij) to Kij; this does not change its bracket with p_T.
Its contraction with the Hessian is identically zero. The metric minor
has determinant `-16(1-2A a^4)^2`; use its regular branch. The source
projection is `q_null=-A a`. Since `S^2=X q^2=0`, the total principal
Hessian still has rank six. But the **linear** term `-W.S` changes the
scalar primary to

```text
Psi=n_null.(p-p_light,lower)+(a T-a W1) q_null,
{Psi,p_T}=a q_null=-A a^2.
```

The two primaries have an invertible antisymmetric bracket matrix when
A!=0. Their preservation fixes their two multipliers instead of producing
the two expected secondaries. On this regular stratum they provide only
two second-class constraints, again giving the generic local count nine.
At the actual coefficient point `u=1,a=1/128`, A=896 and the normalized
bracket is `-7/128`. This is not a generic-R surrogate. At u=0 the bracket vanishes;
no extension of this count onto that special stratum is asserted.

This check distinguishes loss of a kinetic null direction from loss of
its secondary constraint. That distinction is emphasized in
[Garcia-Saenz](https://arxiv.org/pdf/2106.14960). His generalized-Proca
no-go assumes coupling only through the metric, so it is not imported as
a no-go for the explicitly scalar-dependent P8 source.

## 3. A bounded source-only classification

Hold the pure-light and Maxwell terms fixed and examine

```text
S_mu=u_mu[c(u,X) Box(u)+d(u,X) uHu+f(u,X)].
```

On a non-null slice, normalize the light null by Z=1 as in the preceding
audit. Its contractions are

```text
(Box u)_null=D/a, (uHu)_null=a,
D=1-3X R_X/(2R),          q_null=c D/a+d a.
```

The rank-one Lorentz mass update lifts the scalar null whenever q_null is
nonzero and X!=0. On a tilted slice the temporal-vector pivot then remains
`b^2/X`, so it does not restore the pair. To preserve the intended local
count for **all tilts at fixed u,X**, two distinct a^2 values require

```text
c D+d a_1^2=0, c D+d a_2^2=0  ==>  d=0, c=0, if D!=0.
```

This applies on the open regular set containing the vacuum neighborhood
and the actual spacelike obstruction. Continuous coefficients extend the
vanishing through X=0 there. It makes no claim on an exceptional open
D=0 branch, nor for a joint scalar-vector-metric completion outside this
ansatz. In particular it is **not** a classification of all possible
scalar-vector theories.

Consequently a source-only repair within this class cannot keep the old
acceleration dependence and all its clock cubic vertices. With
`N=1+epsilon*n`, `K_hat=3H+epsilon*k1`, `h=(1+u^2)^3`, the old source is

```text
S_normal=epsilon^2[-2n(k1+3Hn)/h]+O(epsilon^3).
```

A source f(u,X)u_mu, with no scalar second derivatives, has no k1. Its
removal or replacement keeps the classical and quadratic clock action,
but not these interactions.

## 4. Viable classical candidates and their secondary constraints

### 4.1 S0: the recommended minimal candidate

S0 sets S_mu=0 in the **entire** source square. This is not the inconsistent
operation of dropping S^2 while retaining -W.S. It leaves minimally
coupled massive Proca alongside the unchanged degenerate light sector
and two minimally coupled scalar matter fields.

In the auxiliary-gradient Hamiltonian, Proca contains neither A_* velocity
nor metric velocities, so the light primary Psi is unchanged. It is
independent of T and the vector momenta: `{Psi,p_T}=0`. Its self bracket
and its brackets with the other primaries are the original light ones.
Thus its preservation gives a secondary, not an equation for the vector
multiplier. The vector temporal primary gives the usual Proca secondary.
This is the Hamiltonian mechanism, not just a Hessian-rank claim. The
general pure-light primary-to-secondary construction is reviewed in
[Deffayet and Garcia-Saenz, section II.1](https://arxiv.org/html/2004.11619).

There must still be a nonzero second-class pivot. At the actual quadratic
vacuum, after normalizing the positive kappa*volume factor to one, the
relevant Hamiltonian is

```text
H=p_u A_*-A_*^2/2+u^2/2-T^2/2-T G
  +|W_i|^2/2+|pi_W|^2/(2 zeta)+spatial terms.
```

For primaries `(p_*,p_T)`, the secondaries are `(A_*-p_u,T+G)`.
Their four-by-four bracket determinant is one. Continuity gives a regular
local vacuum neighborhood, including small timelike, spacelike and null
gradients. The eight-mode count is `18-8-4/2=8` on such regular branches.

At each finite point of the original timelike clock, source removal leaves
the full quadratic action and the original joint lapse/vector pivot
`diag(-2J,-1)` unchanged, with J>0 from the pinned clock proof. Therefore
its local eight-mode count also persists. **No uniform secondary-pivot
bound throughout a connected vacuum-to-clock domain, or through infinite
clock tails, follows from these two local statements.** The common reference
must justify the further domains and branch transport it actually uses.
It need not prove an unnecessarily global bound if a controlled local or
patched construction suffices for the original closure contract.

### 4.2 SALG: a nonzero source algebraic in the scalar gradient

Take `f=lambda*(R-1)*(X-1)`, with lambda an unspecified fixed finite
coefficient, and `S_mu=f u_mu`. No value is chosen. This family illustrates
that S0 is a minimal choice, not the only source with safe primary algebra.

For unit positive density, its complete temporal dependence after the
Maxwell Legendre transform is

```text
H_vec=-(T-f a)^2/2+|W_i-f u_i|^2/2-T G+E_electric+E_magnetic,
T*=f a-G,
H_vec,red=|W_i-f u_i|^2/2+G^2/2-f a G+E_electric+E_magnetic.
```

The temporal pivot is exactly -1. Psi remains vector-independent and
commutes with p_T, so the light secondary is generated as well. Its
subsequent bracket must still be checked on the desired branch. At the
vacuum and clock it has the same regular quadratic data as S0; this gives
local neighborhoods, not a global positivity or well-posedness theorem.
Its clock source starts `4 lambda epsilon^2 n^2/h`, which differs from
the old n*k1 term. It adds an unnecessary free coefficient and additional
response work, so S0 is preferred for the next candidate.

### 4.3 Exact-gradient equivalence control

A source `S_mu=partial_mu J(u,X)` can be removed by the explicitly invertible
substitution `W'=W-dJ`: `F(W)=F(W')`, and its entire mass term is W'^2/2.
This gives a constrained benchmark outside the gradient-parallel ansatz
when J_X!=0. For example J=(R-1)^2 has zero first clock jet. It is not a
distinct interacting physical cure once all sources, observables, boundary
terms and initial data are transported. For derivative-dependent J, one
must use that full transformation, not pretend it is a metric-only point
map. No quantum Jacobian or state conversion is set to one by assertion.

The old source is not generally exact: on the unitary slice S_i=0 but
S_0 varies spatially with K and lapse, so `(dS)_(i0)=partial_i S_0` need
not vanish. Renaming W-S as a free vector would then change its curl.

## 5. Preservation and the exact local affine lift

The alternatives can be represented before eliminating the 56 connection
complement directions by adding

```text
Delta L/(kappa sqrt(-g))
 =W.(S_old-S_new)+(S_new^2-S_old^2)/2,
W=T_trace-B(u,X)du.
```

It depends on the retained trace, not the complementary connection
coordinates. The completed-square identity gives precisely the new mass
square; the complementary nondegenerate block and four projective gauge
directions are unchanged. Keep the full B and regular q dictionary; do not
reset q on the clock. This is an explicit local covariant **proposed
action change**, not a redefinition of the unchanged original theory.

Both S_old and the proposed sources have zero value and first variation
on the clock. W=0 there. Their action difference begins at cubic order,
so the same classical clock, stationary connection, physical matter
metric, canonical quadratic dynamics and quadratic prepared-state data
are preserved. Vacuum mass/residue data are likewise preserved. Higher
interaction vertices are not. Classification/operator-support coverage
must be reassessed for the separately named parent; no automatic all-row
P8 statement is made.

### Formal flat-vacuum four-scalar projection

There is a stronger source filtration for the actual S238 parent than for
the older S174 parent. With `u=epsilon*phi`, `X=epsilon^2 Y`, the literal
retained R has

```text
R-1=1024 epsilon^6 Y^2(3phi^2-7Y)+O(epsilon^8).
```

The S238 addition cancels the old degree-four coefficient. The localizer
A=10^420 is fixed and finite before taking this formal jet; no convergence
radius or finite field-domain approximation is inferred. Substitution
into the **regular** source gives minimum scalar degree six:

```text
S_old,leading=1024 du_mu(3u^2-7X)[X Box(u)-uHu].
```

Every source-free vector vertex has even W parity. With no external W,
a source-tagged graph needs at least two W.S vertices of total degree
seven each, or one S^2 vertex of degree twelve. For any connected graph,
`E=2-2L+sum_v(d_v-2)`. At one loop this gives E>=10; extra vertices cannot
lower it. Hence no source-tagged four-scalar graph exists at this order.
The SALG source starts still later, at degree seven. The more conservative
old-S174 degree-four bound would already suffice for the four-leg result.

This preserves the formal action-vertex vacuum calculation and the chosen
four rate **targets**. It does NOT prove equality of already completed
physical amplitudes: no original global constrained reference was
completed. Gauge/measure conversion, complete hard/radiative bounds and
finite matching still need their own common reference.

At the clock the source starts at degree two, not six. Its known cubic,
lapse-generated quartic and one-loop response therefore change. In S0
the source-tagged vertices of the
[previous response packet](assessment-2026-09-24-p8-source-response-pole-match.md)
are absent. Free Proca gravitational stress/response, all light/metric and
matter sectors, their interactions and known scale conversions remain.
Neither vector tadpoles nor curved/state-dependent counterterms can be
carried over by claiming that all quantum response is unchanged.

## 6. Unchanged-action EFT comparison: possible, not yet controlled

The S6 contract already treats the action as an EFT. The exact-rank
obstruction is therefore **not** an exclusion of that contract. A
perturbative treatment must specify its expansion, constraints, observables
and domain, and bound its omitted terms. Formal reduction of order alone
does not furnish such a bound. These distinctions and the role of field
redefinitions are discussed by
[Glavan](https://arxiv.org/html/1710.01562); the separate adiabatic condition
on time-dependent backgrounds is explained by
[Burgess and Williams](https://arxiv.org/html/1404.2236).

Here is an explicit sufficient diagnostic, without choosing a cutoff.
After gauge/constraint reduction on a specified common domain and with
fixed causal/state boundary data, let G0 be a bounded source-free Green
operator and DeltaD the retained-source correction. IF

```text
rho=||G0 DeltaD||<1,
G_p=sum_(j=0)^p (-G0 DeltaD)^j G0,
```

then a Neumann identity gives

```text
||G-G_p|| <= rho^(p+1)/(1-rho) * ||G0||.
```

This is a conditional linear-response error bound, not a completed
nonlinear or quantum EFT. The actual domain, uniform rho, state/constraint
transport, interaction bounds and loop remainder would all be needed.
One cannot use an inverse on a real on-shell pole without specifying the
causal prescription and functional spaces.

The tiny actual coefficient on the earlier spacelike interval does not
by itself supply rho. As a diagnostic only, freezing the acceleration
coupling produces the two-coordinate form

```text
L=K_s dot(v)^2/2+zeta dot(w)^2/2-(w+alpha ddot(v))^2/2,
det D(z)=z[-K_s+K_s zeta z-alpha^2 zeta z^2], z=omega^2.
```

Its relative higher-derivative entry grows as alpha^2 z/K_s. Small
alpha does not control unrestricted frequencies. K_s and these roots
are **not** the physical P8 reduced residue or pole masses. No numerical
cutoff, ghost mass or instability timescale is inferred from this toy.
The previous finite-regulator identities and small spectator bounds do
not establish the missing uniform operator/error estimate. The EFT route
remains CONDITIONAL, not excluded; it is less ready than S0 for the next
bounded common-reference calculation.

## 7. Decision and next work

The bounded SOURCE-1 study is complete. Recommend adoption of **S0 as a
new candidate**, keeping the old candidate and its obstruction intact.
This requires explicit approval because it changes the physical classical
source. It does not require inventing a source coefficient or a cutoff.

If adopted, the next milestone is **S0-REFERENCE**, with three admission
tasks: establish the required regular constraint domains and any branch
transport (without assuming global control from the two local neighborhoods),
complete its regulated
gauge/measure and Ward/finite-conversion prescription, and recompute/bound
the changed curved response using the original preparations. Preserve the
approved COVZERO rule relative to that newly completed reference; do not
reset known legacy terms. Then assemble the existing RATE4 hard/radiative
projection and B error budget. V/G still require their actual observable,
spectral/contour and remainder inputs; S0 supplies none by fiat.

The independent audit, source-class restriction and preservation checks
are in [the read-only diagnostic](../scripts/p8_source1.py). The
[validation receipt](validation/p8-source1-2026-09-28.json) records replay
and preservation evidence. Finite symbolic controls accompany the written
proofs; they are not a formal quantum-field-theory certificate.
