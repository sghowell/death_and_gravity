# S0 adoption and the canonical/vector reference bridge

Date: 2026-09-28. The user explicitly approved adopting the recommended
S0 candidate and proceeding with S0-REFERENCE. This supersedes the
**adoption-pending status**, not the scientific evidence, in the frozen
[SOURCE-1 study](assessment-2026-09-28-p8-source1-feasibility.md).

**Result:** adopt **QG2-H8A420-RATE4-COVZERO-S0 v1** as a separately named
candidate. Establish its compact-clock local auxiliary admission, exact
source-free Gaussian split, retained-state vector probe vertices and
ordinary Proca reference bridge. Derive the full finite-regulator change
from the old source response, including its local contact. The common
light/metric quantum reference and continuum response bounds are still
in progress. There is no further user-approval blocker for that work.

The [new specification](candidates/p8-rate4-covzero-s0.json) is physical
candidate data, not a modification of the old parent or its reports.
The [diagnostic](../scripts/p8_s0_reference.py) pins 507 inherited inputs;
the [validation receipt](validation/p8-s0-reference-2026-09-28.json) records
the scoped replay. Original M/V/B/R/P8 remain OPEN and G remains UNTESTED.
This packet does not complete PRESCRIPTION-2 or original P8.

## 1. Exactly what has been adopted

In the physical +--- metric, retain the massive vector but replace its
**entire** source square:

```text
kappa sqrt(-g)[-zeta F(W)^2/4+(W-S_old)^2/2]
    -> kappa sqrt(-g)[-zeta F(W)^2/4+W^2/2],
zeta=10^-6,  A_mu=sqrt(kappa*zeta) W_mu,  m_V=1000.
```

The exact affine lift is

```text
Delta L/(kappa sqrt(-g))=W.S_old-S_old^2/2,
W=T_trace-B(u,X) du.
```

It depends on no complementary connection coordinate. All 56 algebraic
complement directions, the projective gauge and the full original B/q
dictionary remain as before. This does not authorize dropping either
canonical primitive or either scalar boundary-momentum shift.

Keep the entire S238/S240 light/matter functions and fixed profiles,
the physical metric, all masses, original S176/S240/S251 preparations,
the four RATE4 targets and detector1/256. Keep the already approved
COVZERO first-order finite-functional rule and selected dimensional
light/metric continuation. Set the vector source to zero at every d in
that continuation. The curvature density is `-R*Ricci_P8/2`, with R the
tensor coefficient and X=+g.du.du. No cutoff, state reset, field-domain
restriction or additional prescription for an extra mode is adopted.

The old source-free flat four-scalar filtration and classical/quadratic
preservation proofs are inherited from SOURCE-1. They are not a computed
full physical amplitude or a claim that the curved response is unchanged.
The earlier SOURCE-1 specification and receipt remain historical snapshots
of a study which, at that time, did not adopt the candidate.

## 2. Full trace/temporal reduction and compact-clock admission

Use the original regular scalar-unitary chart only where du is timelike.
With gamma the hat spatial metric, V its volume, and normalized momenta
as in [S257](../problems/P8/s6/continuation/s6_257/notes/canonical.md), put

```text
M=R^(1/4), U=R^(-3/4), Cchi=R^(-1/4), C3=R^(3/4)/2,
a_trace=-M/3, p=2 tr(pi_metric)/(3V), G=partial_i pi_W^i/V,
B=-U R_u/(2N)-I,
Fhat=U[F+9R_u^2/(16R N^2)]-I_u/N,
I_N=3U R_u R_N/(4RN), I(u,1)=0.
```

The complete S0 trace density and remaining energy are

```text
Ltrace=a_trace*K^2+B*K+Fhat+U*T^2/2,
Eother=2|pi_TF/V|^2/M+(p_M1^2+p_H^2)/(2U)
       +Cchi(|grad M1|^2+|grad H|^2)/2
       +U(m_H^2*H^2/2-j_H*H)-C3 R3
       +|pi_W/V|^2/(2 zeta Cchi)+zeta M |Fij|^2/4
       +Cchi |Wi|^2/2.
```

All gradients, curvature, electric/magnetic terms, potentials and the
entire retained heavy source are present. The normal Hamiltonian is
`Hraw=N[(pK-Ltrace)|Kstar-TG+Eother]`. Thus

```text
Kstar=(p-B)/(2a_trace), Tstar=-G/U,
Hred=N[(p-B)^2/(4a_trace)-Fhat+G^2/(2U)+Eother].
```

For D=d_(N,T)^2 Hraw, direct differentiation at fixed canonical data
and finite spatial jets gives

```text
D_TT=-N U,
D_NT|Tstar=-D_TT*partial_N Tstar,
(D_NN-D_NT^2/D_TT)|Tstar=partial_N^2 Hred.
```

In particular the last derivative acts on Eother; a trace-only lapse
test would be insufficient. The spatial constraint generator and its
primary-momentum transport terms are inherited unchanged.

At formal tree grade on the actual clock, W=G=0, N=R=1. The auxiliary
block is diag(-2J,-1) in the normalized variables. The literal S174
clock polynomial is

```text
J(u)=P(u^2)/[800(1+u^2)^18],
degree(P)=17, every coefficient positive, P(0)=1215.
```

Consequently on every finite interval |u|<=T,

```text
J(u) >= 243/[160(1+T^2)^18] > 0.
```

The time-dependent background and coefficients are smooth. Positive
N,R,zeta and this compact bound give nonzero local kinetic/auxiliary
pivots. The parameter-dependent implicit function theorem, with the
canonical fields and the finitely many spatial jets as parameters,
therefore supplies local branches along the clock; on a compact segment
one can choose a common sufficiently small bounded finite-jet
neighborhood by a finite cover. On overlaps the nearby roots agree by
local uniqueness. This is auxiliary algebraic admission, not nonlinear
PDE well-posedness, an evaluated continuum inverse norm or a global
field-space gauge theorem. It imposes no new physical field restriction.

The fixed stress profiles retain formal quantum grade one, and all their
background dependence must still be varied. The displayed bound is the
tree pivot used in the formal first-loop construction; it is not a bound
on an unexpanded interacting mean solution. Nor is it uniform as T tends
to infinity. At du=0 use SOURCE-1's regular covariant auxiliary-gradient
vacuum chart and its nonzero secondary pair, not scalar-unitary gauge.
No connection between arbitrary field-space branches is inferred merely
from having these local patches.

## 3. Complete homogeneous Gaussian reduction and its measure

The [S253 algebra](../problems/P8/s6/continuation/s6_253/src/p8_vacuum_affine_physical_background_vertices/coupled.py)
defines a general homogeneous scalar/longitudinal-vector quadratic
Lagrangian. Independently start from its full light/M1 `base_lagrangian`
and add **only** `K*(A_dot-P*A0)^2/2+Z*A0^2/2-Yv*A_L^2/2`.
Legendre-transform all three velocities, before eliminating the lapse,
shift divergence and temporal vector. The kinetic and auxiliary Hessians
are exactly

```text
Kinetic=diag(-6D,Z,K),
Daux=diag(-2J,-2D/3,-Z),
det(Daux)=-4DZJ/3,
det(Kinetic)*det(Daux)=8D^2 Z^2 J K.
```

Let chi_i=partial_(n,b,A0) Hraw. Their mutual bracket E is generally
nonzero. With primary momenta followed by chi, the full second-class
matrix and its inverse are

```text
C = [ 0    -Daux  ],   C^-1 = [ Daux^-T E Daux^-1   Daux^-T ],
    [ Daux^T  E   ]           [ -Daux^-1                 0  ].
det(C)=det(Daux)^2.
```

Since the retained phase coordinates commute with the auxiliary primary
momenta, the lower-right zero of C^-1 makes their Dirac brackets exactly
canonical. Integrating
`delta(p_aux) delta(chi) sqrt(det C)` over the auxiliaries cancels the
Jacobian |det Daux| of chi on the regular local branch. This establishes
the reduced Liouville measure at the stated common finite regulator;
it does not omit E or infer the full covariant gravity gauge measure.

For every admitted homogeneous background jet the reduced phase matrix
splits into a coupled **four-dimensional light/M1 block** and a
**two-dimensional longitudinal Proca block**, with the latter
`diag(Yv,1/K+P^2/Z)`. The two light modes do not split from each other.
Together with two transverse vectors, two TT modes and H this remains
16 phase dimensions/eight modes. The heavy source's existing formal
clock-jet filtration does not delete its free metric determinant.

Setting the module's source slots r,rN to zero is checked against the
independent transform. It is **not** setting the actual R, R_N, Theta,
J or their physical-probe derivatives to constants. Both old scalar
boundary shifts and all time-dependent state transports remain.
At Wbar=0, evenness of the vector action likewise makes its unreduced
Hessian block diagonal from the other fields; extending that observation
to a complete covariant quantum measure still needs the remaining
light/metric construction.

### 3.1 The exact change from the old finite Gaussian response

Use the original S253 notation and define

```text
L0=Theta*(p_v+3c*sigma)/D-w*p_sigma/Z-Lnv*v,
bS=r*Theta/D+rN*dH,  r=R-1, dH=Hhat-H_clock.
```

The old-source minus S0 reduced Hamiltonian is exactly

```text
Delta H=P*p_A*[r*(p_v+3c*sigma)/(2D)-3*L0*bS/(2J)]
        +P^2*p_A^2*[9*bS^2/(4J)-3*r^2/(4D)].
```

This retains the lapse-generated quadratic-source contact. At the clock
r=dH=0. Its first probe Hessian is light/vector off diagonal, while its
second includes a nonzero vector diagonal contact. The unchanged initial
product between the coupled light block and the vector block therefore
gives zero mean for the first difference at a common finite regulator.
The metric/light first Gaussian mean jets agree there; the second response
does not. Cross contractions of a block-diagonal source-independent
vertex with an off-diagonal first difference also vanish.

Thus the old-minus-S0 second response is the two-first-source-vertex
connected kernel **plus its complete second-vertex contact**, in the
same state and regulator. The previous source-squared calculation can
be used for this comparison; one must not retain its contact after
deleting its two-vertex contribution, or infer a zero vector stress.
This is a finite-regulator identity and a bookkeeping reduction for the
future common subtraction, not a new bound on a renormalized response.
It does not assert that all old tadpoles or all covariant finite terms
are unchanged.

## 4. Physical vector probes, original state and Ward check

Work in the physical coordinate amplitude A and conjugate Pi, holding
comoving k fixed when varying the metric. The hat scale is
ahat=R^(1/4)*a. Substituting the actual D/Z/K factors into the preceding
longitudinal block cancels **all** R dependence. For `H=z^T Hmat z/2`,

```text
Hmat_L=N*diag(a/zeta, 1/a+zeta*k^2/a^3),
Hmat_T=N*diag(a/zeta+k^2/a, 1/a).
```

There are two T polarizations and one L. For N=1+eta and a=a0*exp(v),
the physical probe matrices at zero probe obey

```text
H_eta=H0, H_eta_eta=0,
H_v=a*partial_a H0, H_eta_v=H_v,
H_vv=(a*partial_a)^2 H0,
H_L,vv=diag(a/zeta,1/a+9*zeta*k^2/a^3), H_T,vv=H_T,0.
```

The nonzero H_eta,v and H_vv contacts are essential. A fixed-physical-
metric clock variation has no direct vector-action dependence. That
statement is not obtained by freezing the metric dependence of the
hat-variable map.

Use precisely the original state at t0=-1/2. If g denotes the reduced
kinetic coefficient (called g^2 in S176's bridge note), its endpoint map is

```text
g_L=a/(1+zeta*k^2/a^2), g_T=a, d=dot(g)/(2g),
(A,Pi)=(v_osc/sqrt(g),sqrt(g)*(dot(v_osc)-d*v_osc)),
d_L=Hubble*(1/2+z), d_T=Hubble/2,
z=zeta*k^2/(a^2+zeta*k^2).
```

The map is symplectic and supplies the same S55/S176 mode operators.
Do not apply the shear twice to a covariance which already includes it,
or replace the all-order prepared state by an instantaneous vacuum.
For compact probes away from the preparation surface the initial state
is held fixed. The finite CTP trace follows the old continuous
metaplectic phase, not separately selected principal square roots.

For the original Wightman phase covariance W(t,s), the vector contribution
to two quadratic insertions is

```text
C_AB(t,s)=1/2 sum_polarizations Tr[H_A(t) W(t,s) H_B(s) W(t,s)^T],
chi_AB=2 Im C_AB,
retarded term=-theta(t-s)*chi_AB,
contact=-Tr[H_AB(t) V(t,t)]/2.
```

The transpose is ordinary, not Hermitian. These follow directly from
[S253's common-regulator formula](../problems/P8/s6/continuation/s6_253/notes/ctp.md).
For each unperturbed mode, let E=H0 and pressure=-partial_a E/(3a^2).
Hamilton's equations give
`dot(E)=Hubble*a*partial_a E=-3*Hubble*a^3*pressure`.
The diagnostic checks the vanishing canonical evolution contribution
independently. This is the homogeneous metric Ward identity before the
momentum sum. Covariant counteractions preserve its renormalized
interior version; a noncovariant cutoff needs its own transport.

S0 also has exact W->-W parity, including its canonical momenta and
temporal pair. The original even vector state, metric-only gauge and
parity-preserving canonical measure preserve this symmetry. Odd vector
expectations and mixed linear W/light response at Wbar=0 vanish in such
a reference. This does **not** set Gamma_WW, vector stress loops, or the
response at a nonzero mean vector to zero.

## 5. Ordinary Proca quantum bridge: what it proves and what it inherits

For the constant physical mass, the ordinary Proca sector uses the
one-form-minus-scalar determinant, with a common regulator and compatible
boundary prescription. A transparent finite-complex identity is useful.
Let d1*d0=0 and put

```text
P=d1^T*d1+m^2 I1,
K1=P+d0*d0^T, K0=d0^T*d0+m^2 I0.
```

Then P*d0=m^2*d0, so the determinant lemma gives

```text
det(P)=det(K1)*(m^2)^n0/det(K0).
```

The constant normalization of the constrained canonical measure cancels
the displayed nondynamical mass factor. The continuum local comparison
is the ordinary Proca `Tr_1 log(Delta1+m^2)-Tr_0 log(Delta0+m^2)`, with
the half factor for a real field and the chosen Lorentzian/CTP branch.
See [Ruf and Steinwachs, section V](https://arxiv.org/pdf/1806.00485).
Three exact finite-complex controls check the algebra; they do not alone
prove continuum determinant identities or boundary matching. The latter
ordinary-vector/state convention is the retained
[S176 bridge](../problems/P8/s6/continuation/s6_176/notes/bridge.md) and
[renormalization prescription](../problems/P8/s6/continuation/s6_176/notes/renormalization.md).
No global Euclidean continuation of the bounce is assumed: the elliptic
comparison determines the local UV coefficients, while the original
Lorentzian state fixes the nonlocal finite parts.

At fixed physical g the source-free vector has no direct u dependence.
A covariant regulator transforms its form operators by conjugation under
compactly supported diffeomorphisms. The variation of the regulated trace
of a commutator vanishes with compatible boundary conditions. Its local
poles are covariant and its interior metric stress is conserved. This is
an ordinary nonchiral Proca-sector result, not the missing proof of an
anomaly-free constrained light/gravity measure.

The following finite conversion was **already derived** in the COVZERO
seed and included in S176. This packet rechecks it, not a new loop result.
With R_o=-Ricci_P8 and ell=log(m^2/mu^2),

```text
a4s=R_o^2/72+(Riemann^2-Ricci^2)/180,
h0=d-1, h1=(d-7)*R_o/6,
h2=(d-1)*a4s+Ricci^2/2-R_o^2/6-Riemann^2/12,
P(d)=m^4*h0-2*m^2*h1+2*h2.
```

Keeping d=4-2epsilon through subtraction produces the finite evanescent
term `-2 partial_d P=-2[m^4-m^2 R_o/3+2a4s]`. Together with the ordinary
heat finite part this is exactly the literal S176 polynomial

```text
(5/2-3ell)*m^4+(5/3-ell)*m^2*R_o-2ell*a4V-4a4s,
a4V=-R_o^2/8+29Ricci^2/60-Riemann^2/15.
```

At the common mu=1 retain `+log(10^6)*P(4)/(64pi^2)` to recover the old
mu=1000 normalization. Do not add the evanescent finite term a second
time or erase it under the complementary-zero rule. All surface terms
remain; displayed bulk comparisons use compact variations away from
preparation/final surfaces. The uncomputed light/mixed residues cannot
be set to zero by analogy with this completed ordinary-vector check.

## 6. A restricted frame shortcut is inadmissible across the clock

In the regular hat chart the ADM quartic coefficients are
`A4hat=-R^(1/4)/2`, `B4hat=R^(3/4)/2`, in the convention
`A4hat*(K^2-Kij Kij)+B4hat*R3`.
Try keeping the hat spatial metric fixed and replacing only the lapse
by F(u,N)>0. The transformed coefficients are

```text
A4new=(F/N)*A4hat, B4new=(N/F)*B4hat.
```

The quartic Horndeski condition is
`A4new+B4new+F*partial_F B4new=0`; the ADM relation follows from
[Gleyzes et al., Eq. 16](https://arxiv.org/html/1404.6495).
On a regular lapse map it yields

```text
F*F_N=-N*(B4hat+N*B4hat_N)/A4hat
     =N*sqrt(R)*(1-3X*R_X/(2R)).
```

On the actual tree clock R=1, X=1, R_X=(1+u^2)^-3. The right side
vanishes at `u^2=(3/2)^(1/3)-1`. Positive F would require F_N=0, losing
local invertibility. Thus this **fixed-hat-spatial-metric lapse-only**
route cannot transport a Horndeski measure over the full clock segment
through those points. This is not a theorem about every field
redefinition or a pathology of S0: the original chart and J stay regular
there. Do not replace the remaining measure work with this singular
shortcut or divide the canonical evolution by this new factor.

## 7. Remaining executable obligations and stopping boundary

The next reference task is the light/metric sector, not a new source
choice or another finite-boundary approval:

1. Construct its gauge and second-class measure with the chosen
   dimensional continuation, retaining canonical state/boundary and
   ordering transport. Show compatibility of the needed vacuum and
   clock patches; global coverage of irrelevant field space is not a
   new requirement.
2. Verify the needed Ward and UV-locality hypotheses, and compute the
   pole/finite conversions for the light/metric projections entering
   the four rates and curved response. A finite-dimensional Liouville
   identity is not this proof. General BRST renormalization results
   explicitly assume an invariant/anomaly-free measure and local
   divergences; they cannot supply those hypotheses by citation.
   See [Barvinsky et al., sections 2 and 4](https://arxiv.org/pdf/1705.03480).
3. Assemble the fixed-state causal finite parts and actual response/error
   bounds using the S0 vertices, then the complete hard and radiative
   four-rate inputs. Preserve all other sectors and the independent
   finite-gravity same-observable error requirement.

No missing light/gravity pole, finite part, norm, Regge envelope or
omitted-order allowance is assigned a default. The full-J gravity route
remains optional. The explicit vector bridge removes the former
source obstruction from this candidate but does not by itself close
physical matching, an interacting/global bounce or classification coverage.

## 8. Validation scope

Run `.venv/bin/python -B scripts/p8_s0_reference.py`. The diagnostic checks
adoption and retained data, 507 protected input hashes, the full auxiliary
Schur identity and literal compact J bound, an independent Legendre
reduction, nonzero lower constraint bracket and physical Dirac brackets,
exact old-source difference, physical first/second vector probes, state
symplecticity, homogeneous Ward algebra, parity, finite Proca determinant
controls, inherited heat/scale conversion and the restricted frame-map
zero. Eleven invalid-input controls reject source/state/cutoff/finite-rule
changes and premature gate promotion.

The receipt records all 16 focused root diagnostic replays, lint/format
and whitespace checks, source hashes and preservation of all 45
preexisting non-status files, including unrelated P4/P9 work. These are
exact algebra and written scoped arguments at **VERIFIED_N**, not a Lean
formalization, independent external review or rerun of the old 115,487-test
snapshot. Historical packets and the old sourced parent are not rewritten.
