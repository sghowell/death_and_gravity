# RATE4-COVZERO-1: authorized boundary and a regular reference seed

Date: 2026-09-24. The user approved the explicit finite-boundary candidate
recommended by the [closure audit](assessment-2026-09-24-p8-response-closure-blockers.md).
The approval removes that audit's user-decision blocker. It does not
certify a completed quantum reference, gravity allowance or original P8.

**Result:** adopt the first-order RATE4-COVZERO-1 boundary rule as new
physical candidate data. Construct a dimension-dependent light/metric
action seed that preserves its kinetic degeneracy, instead of silently
continuing the four-dimensional coefficient unchanged. Extend the
regular chart's cotangent bookkeeping to independent scalar gradients.
Derive the explicit finite scale transport needed to retain the old
Proca prescription at a common subtraction scale. Full quantum-reference
construction remains in progress; no new user decision is needed for
that work. Gravity still requires its independent physical certificate.

The [candidate specification](candidates/p8-rate4-covzero-1.json) records
the adopted rule, selected reference seed and unresolved obligations.
It deliberately distinguishes an adopted physical boundary rule from
a complete executable quantum model. No frozen scientific packet is
rewritten. M/V/G/B/R and P8 remain OPEN, G UNTESTED.

## 1. Authority, conventions and admission

Retain S238/S240's whole classical parent, fixed profiles, physical
matter metric, original prepared states, masses, four rate targets and
detector. At first quantum order the additional finite functional is

```text
DeltaGamma = integral sqrt(-g_phys) V(X) sum_(a=0,S,U,SS) c_a O_a,
V(X)=(1-X)^1024/[X^1024+(1-X)^1024].
```

Use the original proposed fully symmetrized covariant four-slot lift
and no extra independent complementary finite functional in the
completed reference. The c_a must solve the unchanged four physical
rate equations using complete reference inputs. No missing coefficient,
loop value, radiative input or error bound acquires a numerical default.
The old S239 OS4 value is reference data, not a fifth condition on the
full new amplitude. Reference changes transport the full boundary
functional oppositely; zeros are not reset in a changed reference.

**Convention correction:** section 3 of the preceding closure note
wrote X=-g^(mu nu)u_mu u_nu. That minus sign is a prose typo. The retained
P8 signature is +--- and the correct definition is
X=+g_phys^(mu nu)u_mu u_nu, hence X=N^-2 in scalar-unitary gauge.
Both previous witness calculations already used X=N^-2, so their jets,
response determinant and obstruction are unchanged. The old note and
its manifest remain historical evidence; this explicit correction is
the current convention. R below denotes the positive tensor coefficient,
NOT spacetime curvature.

Admission: the addressed M/B obligation is one shared reference for
flat normalization and curved response. A successful seed must recover
the unchanged four-dimensional action, preserve the kinetic constraint
under its regulator continuation, and preserve already fixed finite
Gaussian prescriptions. A nonzero Schur complement or unmatched finite
stress would reject that particular implementation, not the physical
parent. The naive continuation fails the first check; the selected one
and the explicit scale conversion pass the checks below.

## 2. A dimension-dependent covariant light/metric seed

Let d=4-2epsilon, m=d-1, and retain the exact coefficient functions
F(u,X), R(u,X), their fixed profiles and parameters. In the physical
metric scalar-tensor density R*Ricci_P8/2+F+sum A_i L_i select

```text
A1=A2=0,
A3=R_X/X,
A4=-R_X/X-(3m-2)*R_X^2/[2(m-1)R],
A5=R_X^2/(R X).
```

Here L3=(Box u)u^mu u_munu u^nu,
L4=(u_munu u^nu)(u^murho u_rho),
L5=(u^mu u_munu u^nu)^2, with the repository's inherited conventions.
At m=3 the coefficient in A4 is 7/4, exactly the
[S238 action](../problems/P8/s6/continuation/s6_238/notes/family.md).
This is a chosen evanescent continuation, not a derived physical theory
in noninteger dimension. No other bare evanescent operators are added
to this seed. Explicit regulator conventions must subsequently be used
through the complete loop and counterfunctional calculations.

### 2.1 Arbitrary-tilt kinetic check

This calculation is not confined to scalar-unitary gauge. On a general
ADM slice let A_*=n.du, A_i be the spatial gradient, b^2=A_i A_i,
X=A_*^2-b^2, K_ij the metric velocities, V_* the scalar normal
second-derivative velocity, and

```text
Z=A_* V_*+A_i K_ij A_j.
```

The principal velocity density in the scalar/metric sector is

```text
R/2*(K_ij K_ij-K^2)-2 R_X K Z
 +A3*A_*(V_*+A_* K)Z+(A4+A_*^2 A5)Z^2.
```

For A_* nonzero use Z as the scalar velocity coordinate. The mixed
tensor and scalar-square coefficient are

```text
T_ij=(-2R_X+A_*^2 A3)delta_ij-A3 A_i A_j,
beta=A3+A4+A_*^2 A5.
```

The DeWitt metric on K_ij is invertible for m>1, R>0. Its Schur
complement is beta-[T_ij T_ij-(tr T)^2/(m-1)]/(2R).
Rotate A_i onto the first axis and write l=b^2/X. Then

```text
T_11=-R_X,   T_aa=R_X(l-1) for a>1,
(T_ij T_ij-(tr T)^2/(m-1))/(2R)
    =R_X^2/R * [l-m/(2(m-1))]
    =beta.
```

Thus the full scalar/metric kinetic matrix has precisely the intended
null direction wherever this regular metric block and coordinate
description apply. The polynomial original expression also covers
A_*=0 without dividing by A_*; at a purely spatial gradient the V_*
row vanishes. The diagnostic checks nine exact finite-dimensional
ranks, including this separate branch. These fixtures accompany the
invariant Schur proof; they do not replace it.

If instead one keeps A4=-R_X/X-7R_X^2/(4R) at all d, the Schur
complement is

```text
-(m-3)R_X^2/[4(m-1)R].
```

It is generically nonzero off d=4. This rejects the naive continuation
as an implementation that purports to preserve the original kinetic
constraint. It does not prove that every constraint-breaking regulator
is unusable with suitable repairs, or that a physical four-dimensional
ghost has been found.

The independent proof is for the pure scalar/metric kinetic sector.
It is not a full arbitrary-background secondary-constraint count for
all matter/source fields, a BRST anomaly calculation or quantum health
theorem. The geometric-frame literature supplies useful context but
not this dimension-dependent completion by citation.
[Langlois, Noui and Roussille](https://arxiv.org/pdf/2012.10218).

### 2.2 Exact clock-chart cancellations and the primitive

Generalize the regular [S174 chart](../problems/P8/s6/continuation/s6_174/notes/chart.md)
by choosing

```text
C=R^(-1/(m-1)), D=(1-C)/X,
g_phys=C*g_hat+D*du*du, X_phys=X_hat,
omega=(log C)/2, s=sqrt(X).
```

In scalar-unitary gauge h_phys=C*gamma, with the same lapse and shift.
The trace and squared extrinsic curvature become

```text
K_phys=K_hat+m*b,
(K_ij K_ij)_phys=(K_ij K_ij)_hat+2b K_hat+m b^2,
b=s*omega_u+2s*omega_X*V,  V=n.partial(s).
```

The inherited ADM coefficients generalize to

```text
B=-R/2, C_ADM=-s R_X,
D_ADM=-m X R_X^2/[2(m-1)R],
E_ADM=-X R_X+(3m-2)X^2 R_X^2/[2(m-1)R].
```

Substitution cancels both K_hat*V and V^2. With
delta=X R_X/[(m-1)R], the full spatial-curvature integration by parts
gives the separately checked cancellation

```text
E_ADM+(m-1)(m-2)(R/2)delta^2
  +2(m-1)(R/2-X R_X)delta=0.
```

No spatial lapse-gradient term is being dropped. The remaining linear
V coefficient before volume multiplication is
m X R_u R_X/[(m-1)R]. Define

```text
U=R^(-m/[2(m-1)]), M=U R,
I_s=U*m*s^2*R_u*R_X/[(m-1)R], I(u,1)=0,
B_hat=-U*R_u/(2N)-I,
F_hat=U*[F+3m R_u^2/(8(m-1)R N^2)]-I_u/N.
```

The entire time and spatial primitive boundary remains. At m=3 these
are exactly [S257's conventions](../problems/P8/s6/continuation/s6_257/notes/source.md),
including the differentiated primitive, not just its value on the clock.
The trace kinetic coefficient is -M(m-1)/(2m), matter/electric spatial
factor U/C, and intrinsic-curvature coefficient U R/(2C).

### 2.3 Source continuation, regularity and evanescence

Write R-1=X^2 A(u,X), using the inherited regular factor. Select the
following covariant one-form in the hat variables:

```text
S_mu=A(u,X)*[X*Box_hat(u)-u^alpha*u^beta*nabla_hat_alpha(u_beta)
                         -m*H_clock(u)*X^2]*u_mu.
```

It is regular at X=0, including nonzero null gradients. Its unitary
normal component is (R-1)(K_hat-m H_clock/N), exactly the retained
four-dimensional source at m=3. The entire Proca source square, local
contact and physical-metric matter couplings are retained, not aligned
away before variation. This defines their continuation after the
four-dimensional algebraic connection reduction; no unproved
d-dimensional reconstruction of all affine projectors is asserted.

The joint trace/source velocity factor is

```text
Gamma_m=1-m(R-1)^2/[(m-1)R].
```

For 5/2<=m<=7/2 and 1/2<R<6/5, it exceeds 1/6, since
(R-1)^2/R<1/2 and m/(m-1)<=5/3. This is a useful regular trace-block
neighborhood of the regulator dimension. It is not the full lapse
pivot, complete spectrum or continuum inverse estimate.

The chart remains regular: C>0, C+DX=1, and the metric-position
Jacobian is C^[d(d+1)/2-1]. Its exceptional rank-one eigenvalue is
C-X C_X-X^2 D_X=1. Because R-1=X^2 A,
D=A X/(m-1)+O(X^2) near zero. The apparent divided coefficients in
the action likewise extend smoothly.

Do not discard the epsilon derivatives. In particular

```text
partial_epsilon A4|0=-R_X^2/(4R),
partial_epsilon log(C)|0=-log(R)/2,
partial_epsilon log(U)|0=partial_epsilon log(M)|0=-log(R)/4,
partial_epsilon log(U R/(2C))|0=+log(R)/4.
```

They can multiply ultraviolet poles to leave finite functionals. The
actual light/mixed pole residues still need calculation; this packet
does not assign their conversion terms zero.

## 3. Extend the chart bookkeeping without freezing the scalar gradient

S257's exact cotangent lift is already available in the scalar-unitary
chart. To address the apparent derivative dependence away from that
chart, first introduce an independent covector A_mu and enforce
A_mu=partial_mu u with its delta functional/Lagrange multiplier. This
is an auxiliary first-order presentation, not a new propagating field
or state. In the enlarged position variables (g_hat,A,u) the disformal
map is an ordinary point map, regular also on the null locus.

Its position Jacobian is block triangular. The metric block has the
rank-one form above, while the A,u blocks are identity. For its complete
point Jacobian J, transport momenta by p_old=J^(-T) P_new. The canonical
one-form identity p_old.dq_old=P_new.dq_new then gives unit phase
Jacobian and symplecticity. The momenta conjugate to A_mu and u acquire
the induced shifts; keeping only the metric momentum rescaling would
not be this transformation.

In differential form the phase Jacobian is

```text
[[J,0],[-J^(-T) H,J^(-T)]],
H_ab=sum_i p_old_i partial_a partial_b f_i = H_ba.
```

The diagnostic independently builds the actual 15-position Jacobian
and checks the 30-phase identity in timelike, spacelike and nonzero-null
fixtures with X=1/16,-1/16384,0. The generic-R fixture is not substituted
for the original parent. The written identity holds for every regular
point map. Regular second-class constraint redefinitions transform
their delta functions and square-root bracket determinant oppositely
on the constraint surface; this coordinate fact is not a computation
of the model's still missing full regulated gauge determinant.

Retain the constraint that introduced A, source transformations, all
primitive boundaries and the transported initial density. A nonlinear
change need not preserve a Gaussian representation of that density.
Unit classical Liouville volume is not an all-order quantum Weyl-ordering
theorem and does not allow the state to be re-minimized. The need to
transform sources and matching data as well as the action is emphasized
by [Criado and Perez-Victoria, section 2](https://arxiv.org/html/1811.09413v2#S2).

## 4. Explicit Proca conversion to the common scale

The original [S176 prescription](../problems/P8/s6/continuation/s6_176/notes/renormalization.md)
uses mu=m_V=1000, whereas
[S240](../problems/P8/s6/continuation/s6_240/notes/renormalization.md)
uses mu=1. A common reference at mu=1 must transport the first finite
prescription, not silently replace it.

Use the heat-kernel curvature convention R_o=-Ricci_P8; on the FLRW
clock R_o=6(H'+2H^2). This agrees with S240's explicitly named R_old
in its bounds note. Squared curvature invariants are unchanged under
the signature bridge. To avoid confusing the heavy mass with the
dimension m, write the vector mass m_V explicitly.

For the source-independent ordinary vector determinant, the minimal
one-form-minus-scalar bulk heat coefficients at general d are

```text
a0=d-1, a2=(d-7)R_o/6,
a4=(d-1)*a4s+Ricci^2/2-R_o^2/6-Riemann^2/12,
a4s=R_o^2/72+(Riemann^2-Ricci^2)/180.
```

These follow by tracing the vector endomorphism and connection
curvature, then subtracting the scalar. At d=4,

```text
a4V=-R_o^2/8+29 Ricci^2/60-Riemann^2/15,
P_V=3 m_V^4+m_V^2 R_o+2 a4V.
```

The resulting bulk pole agrees with the ordinary Proca calculation in
[Ruf and Steinwachs, section V, Eq. 39](https://arxiv.org/pdf/1806.00485),
after conventions are matched. The diagnostic binds the actual finite
and pole polynomials directly to the frozen S176 report, not to an
unrelated vector model.

Expanding the full d-dependent coefficients before subtraction gives
the finite evanescent contribution

```text
-2*partial_d[m_V^4 a0-2m_V^2 a2+2a4]|d=4
 =-2m_V^4+(2/3)m_V^2 R_o-4a4s.
```

With ell=log(m_V^2/mu^2), the complete local matching density is

```text
L_V(mu)=[(5/2-3ell)m_V^4+(5/3-ell)m_V^2 R_o
                           -2ell a4V-4a4s]/(64*pi^2).
```

This reproduces S176, including its finite nonlogarithmic terms. It is
the local part matched to mode subtraction, not the full state-dependent
determinant or source/mixed response. The chosen subtraction keeps the
fixed four-dimensional pole coefficients multiplying d-dimensional
invariants, whose variations precede the limit. Subtracting their
evanescent finite coefficients as well would be a different prescription.

The required known conversion is therefore

```text
DeltaGamma_legacy=+log(10^6)/(64*pi^2)*integral sqrt(-g) P_V,
Gamma_V(mu=1)+DeltaGamma_legacy=Gamma_V(mu=m_V).
```

The same original finite state-dependent part is retained. This
conversion belongs to the reference, as required by the adopted
boundary rule's retention of existing finite conditions. It is not
an independently adjustable complementary Wilson coefficient, and
it is not multiplied by the clock-null switch.

There is an observable negative control. The bulk energy variation of
P_V on flat FLRW is

```text
rho[P_V]=-3m_V^4-6m_V^2 H^2+(H')^2-6H^2 H'-2H H''.
```

The independent scale-factor Euler variation, with q=log(a), gives
P=f-(partial_t+3H)f_H/3+(partial_t+3H)^2 f_H'/3 for the restricted
bulk density f(H,H'). This polynomial pressure agrees with the continuity
expression and verifies rho'+3H(rho+P)=0, regular also at H=0.
At the bounce H=0,H'=4,
omitting the conversion would shift the existing vector energy by
log(10^6)(3*10^12-16)/(64*pi^2). The adopted four-contact lift has zero
clock first variation and cannot cancel that error. Including the
conversion preserves the old vector mean and all its compact metric
variations; it does not retune the fixed scalar profile.

Keep the scope to bulk compact variations away from the preparation
and final surfaces. A nonminimal Proca determinant can carry
total-derivative/double-pole and factorization subtleties; those are not
erased by a finite-dimensional determinant identity.
[Barvinsky and Kalugin](https://arxiv.org/html/2408.16174v2).
No infinite-history boundary flux or initial-state boundary term is
declared zero here.

## 5. Reference completion: concrete remaining calculation

COVZERO-REF-DIM-1 now fixes a regular dimensional action seed, the
zero-complement boundary rule and the known Proca scale conversion.
Use physical-metric background de Donder gauge with unit gauge
parameter as the covariant target, keeping the projective trace gauge
and the canonical constrained measure. The metric/gradient point-map
calculation above handles one coordinate issue; it does not evaluate
the first-class determinant or complete continuum second-class measure.

The next admission target is the common regularized one-loop functional:
derive its gauge/constraint contribution, verify the Ward identities
and local counterfunctional, and only then combine source-squared,
source-linear mixed and source-independent response pieces. Flat rates
use the vacuum Feynman reference. Clock variations use the original
fixed-state closed-time contour, with no state reset; nonlocal finite
pieces are different even when local UV subtraction is shared.

General background-field renormalization results do not bypass this
verification. Their hypotheses include an invariant regulated measure
and locality of the remaining divergences; these have to be checked
for the actual constrained operator.
[Barvinsky et al., sections 2.3-2.4](https://arxiv.org/pdf/1705.03480).
The present packet does not import those hypotheses as conclusions.

After the common reference is validated, solve the four contact
normalizations with its complete hard/radiative inputs, evaluate and
bound the same-state response, and retain higher-order errors. Gravity
remains independently conditional on a quantitative same-observable
certificate. User authorization of COVZERO is not authorization to
assume a Regge envelope, spectral budget, contour error or UV completion.

## 6. Verification and publication boundary

The [new diagnostic](../scripts/p8_covzero_reference.py) pins 496 inherited
files and records the candidate specification's hash. It checks the
arbitrary-tilt Schur identity, nine finite kinetic ranks, full ADM
cancellations, regular source/chart, four chart epsilon derivatives,
three extended cotangent fixtures, pinned Proca heat/finite polynomials,
scale transport and local stress Ward identity. A synthetic four-rate
fixture checks only the retained normalization interface; it supplies
no physical missing rates. Ten unsupported inputs/statuses are rejected.

The [receipt](validation/p8-covzero-reference-2026-09-24.json) records
all thirteen focused diagnostic replays, lint/format and preservation
checks. The frozen 115,487-test suite is not rerun. Initial tests used
two inappropriate structural symbolic equalities and assumed a specific
printed denominator for a polynomial; exact cancellation/polynomial
checks corrected those test implementations without changing formulas.

All 34 prior non-status edited/untracked files, including unrelated
P4/P9 work, remain byte-identical. Frozen P8 and the original three
tracked diagnostics are unchanged. The work is local, uncommitted and
unpushed. Adoption is now recorded; full reference construction and
original P8 closure are not claimed.
