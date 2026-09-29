# PRESCRIPTION-2 continuation: reference reduction and finite-window gravity

Date: 2026-09-24. Continues the authorized model and gravity tracks from the
[September 22 packet](assessment-2026-09-22-p8-prescription2-gravity.md).
Neither track, nor original P8, is complete. No replacement parent,
finite-zero boundary data, regulator, state, detector or rate is adopted.

**Progress:** the reference problem now separates 56 algebraic connection
directions from the retained vector and constrained light/metric sector.
The original heavy *source* cannot change the reference background or
quadratic response at formal one-loop order, whereas the original Proca
source can. The next source-response calculation is therefore identified,
not hidden inside a generic instruction to compute the whole effective
action. Its actual value and remainder bound are still missing.

On gravity, a full positive forward moment is shown to be an unnecessarily
strong input for the previous comparison method: an explicit positive
even-spin spectral fixture has divergent J but a convergent finite-angle
arc. A finite-energy comparison plus a finite-angle Regge tail avoids
that full-J requirement, with a proved tail coefficient below 60 times
the residue/slope ratio. No physical Regge bound is inferred from this
example or from the small Newton coupling.

## 1. Which connection modes really are auxiliary?

The inherited [S174 reduction](../problems/P8/s6/continuation/s6_174/notes/affine.md)
and [S238 domain](../problems/P8/s6/continuation/s6_238/notes/domain.md)
apply to the actual retained classical parent. Of the 64 connection
components, four are projective gauge, 56 are algebraic, and four form
the retained, dynamical Proca variable. It is incorrect to discard a
60-component determinant as if every quotient mode were auxiliary.

For the already proved full quotient matrix M, trace map N and source
centering, let

```text
D=N M^-1 N^T=diag(tau,sigma,sigma,sigma),
tau=3(2p^3-1)/p,  sigma=(8p+5)/(8p^2),  p=sqrt(R)/2,
L=M^-1 N^T D^-1,  NL=1,  P=1-LN.
```

The source-centered update Xi=eta-D^-1 gives
`P^T M_new L=0`, `L^T M_new L=eta`. Thus a local basis Q of ker(N)
and y=Q a+L w separate a nonsingular 56-dimensional algebraic block
`C=Q^T M Q` from the four-vector. The inherited domain has
`1/3<p<3/5`, hence tau<0, sigma>0 and determinant ratio
`det(M_new)/det(M)=-tau sigma^3>0`. This determinant sign is not a
physical propagating-ghost test. The new diagnostic checks the scalar
update identity; the full 60-row identities remain inherited frozen
evidence, not claimed as independently rerun here.

For an arbitrary retained-field quadratic operator A and mixed block B,
the exact block identity is

```text
det [[C,B],[B^T,A]] = det(C) det(A-B^T C^-1 B).
```

Completing the square has unit triangular Jacobian even if B acts by
derivatives. In a local component measure, integrating the algebraic
Gaussian leaves an ultralocal factor `(det C)^(-1/2)`. Field-dependent
basis changes and the regular metric-chart Jacobian are also ultralocal.
The projective trace gauge has variation `delta Gamma^a_ac=4 U_c`, so
its Faddeev-Popov determinant is field independent in four dimensions.
If continued in this form it becomes d per gauge component, also constant.

In the **specified perturbative dimensional component-measure convention**,
an ultralocal determinant is

```text
Tr log C = integral d^d x tr(log C(x)) integral d^d k/(2pi)^d 1 = 0.
```

The last equality is the scaleless-integral prescription, not a
regulator-independent statement about the measure. A cutoff can instead
produce local power terms, which require a finite-reference conversion.
This reduction removes no nonlocal determinant or physical matching data.

In particular it does **not** derive the lapse/normal-vector/scalar
second-class measure. A canonical algebraic toy with constraints
`pi_a=0, C a+J=0` has `sqrt(det{chi,chi})=|det C|`; its constraint
delta function cancels that factor after a-integration. A naive Gaussian
configuration measure instead produces `(det C)^(-1/2)`. The two are
related by a measure choice; this toy cannot determine the actual DHOST
constraint measure. S174 proves a local classical constraint rank, not
the common covariant quantum measure. Nor does this argument continue
four-dimensional degeneracy/epsilon-tensor identities to d dimensions.

**Decision:** isolate the 56-mode ultralocal factor and retain the complete
propagating/constrained Schur operator. Do not identify dimensional
scalelessness with completion of PRESCRIPTION-2.

## 2. Heavy Gaussian reduction, with the mixed loops retained

For the actual added scalar the physical-metric action has the form

```text
S[q,H]=S0[q] - (1/2) H K[q] H + H J[q],
K=Box_g+n,  J=(g/2) Phi^2 (1-X)^8 exp(-10^420 X^2).
```

Here q denotes all retained light/metric variables. The displayed K uses
the covariant inner product; for variations, K and J include their metric
densities so that the pairing is fixed. With a specified Green function G=K^-1,
Gaussian elimination gives

```text
H_c=G J,
S_eff[q]=S0[q]+(1/2) J G J,
Gamma_H^(1)=(i/2) Tr log K
```

up to field-independent phases in Lorentzian conventions. This is a
functional reduction, not replacement of H by a classical saddle in
the complete quantum problem. The light-field determinant must be taken
from the **full** reduced Hessian. For one repeated variation delta,

```text
delta^2[(1/2)J G J]
 = H_c delta^2 J - (1/2) H_c (delta^2 K) H_c
   + (delta J - delta K H_c) G (delta J - delta K H_c).
```

It is the Schur complement of the coupled light/heavy Hessian. In
particular the last term retains mixed light/heavy loops when the
remaining determinant is expanded; `Tr log K` alone omits them. This
organization is consistent with the primary functional-matching account
of [Henning, Lu and Murayama, sections 2.3 and 3.1](https://arxiv.org/html/1604.01019v1#S3.SS1).
The displayed variation identity is derived here and independently
checked against a field-dependent finite-dimensional Hessian.

The inverse and determinant must have the appropriate contour/state:
Feynman boundary data for the vacuum amplitude; a closed-time-path
construction with the original prepared covariance for causal clock
response. Writing G does not make these two objects equal or reset the
[S240 heavy state](../problems/P8/s6/continuation/s6_240/FORMULATION.md).
The same local ultraviolet subtraction can be compared between them;
their nonlocal, state-dependent finite parts cannot be glued together
without the actual causal construction.

### 2.1 A useful source-sector zero at the original clock

At X=1, the original J has field degree at least eight in arbitrary
compactly supported light/metric variations. The smooth exponential and
Phi^2 do not lower this degree. Assuming the retained Green function and
gauge/constraint reduction are regular on the finite domain considered,
`(1/2)JGJ` therefore starts at degree 16, even when G is varied.
Its Hessian starts at degree 14. Expanding the remaining one-loop trace
log about the source-free Hessian cannot lower that degree.

Thus the source-dependent tree functional has zero variations through
degree 15, and its formal one-loop contribution has zero variations
through degree 13. In particular this **source sector** gives no direct
reference-mean or quadratic-response change at that order. This uses the
original order-eight source zero, not the proposed order-1024 finite lift.
It concerns retained light/metric external variations after H integration.
It bounds the field degree of these source-induced diagrams and their
local subtractions, not an arbitrary independent finite boundary action.

The free heavy determinant still varies with the metric and is not zero.
Its mean stress and finite local density in
[S240](../problems/P8/s6/continuation/s6_240/notes/renormalization.md) remain:

```text
[(3/2-log n)n^2+(log n-1)n R/3-2log(n) a2]/(64pi^2),
a2=R^2/72+(Riemann^2-Ricci^2)/180.
```

Its causal stress response remains necessary. The zero does not bound
higher loops, a nonlinear neighborhood, all-time tails or omitted terms.
No uniform inverse bound is inferred from finite-time regularity.

### 2.2 The Proca source does not have the same protection

The actual [S174 source](../problems/P8/s6/continuation/s6_174/src/p8_vacuum_analytic_affine_parent/source.py)
has zero value and first variation on the clock, but its second jet is

```text
S_n^(2) = -2 nu (K_1+3 H_ref nu)/(1+u^2)^3.
```

nu is the lapse perturbation, not the heavy mass n. At the bounce this
is `-2 nu K_1`, not zero. The source-centered action retains BOTH
`-kappa W.S` and `kappa S^2/2`. Eliminating W with its full inverse
therefore gives a quartic source functional, whose Hessian is quadratic
and can contribute to the one-loop quadratic light/metric response.

A simple frozen-coefficient Euclidean symbol detects why a homogeneous
test would miss it. With p=(omega,k), zeta>0, the normalized Proca operator
and inverse on this two-dimensional subspace are

```text
P=(1+zeta p^2)I-zeta p p^T,
P^-1=(I+zeta p p^T)/(1+zeta p^2).
```

For a temporal source the local-plus-exchange coefficient is

```text
(I-P^-1)_00 = zeta k^2/[1+zeta(omega^2+k^2)].
```

It vanishes at k=0 but is nonzero at nonzero spatial momentum. This is
a symbol-level control of the source cancellation, **not** a Euclidean
definition of the curved clock state, a physical constraint reduction,
or a sign/bound on the actual response. The original canonical boundary
transport must still be retained.

**Next model calculation:** the source-induced Proca/light/metric Hessian
with the actual constraint measure and prepared causal propagators, plus
the free heavy/vector metric response and needed common finite conversions.
This is smaller and more informative than recomputing the negligible
flat heavy remainder or declaring every Gaussian sector free.

## 3. Finite reference changes must transport the boundary functional

The elementary Laurent identity

```text
(A/epsilon+B)(1+epsilon E) = A/epsilon + (B+A E) + O(epsilon)
```

exhibits the finite term lost by taking tensor/dimensional factors to
four dimensions too early. S240 already retains this effect in its
scalar counterstress; it must not be discarded when reconciling sectors.
If a common reference changes by local Delta, the same proposed physical
action requires `F_boundary -> F_boundary-Delta`. Resetting the new
complement to zero instead defines a different proposal. The identity
does not assign the actual missing finite conversion coefficients.

Still needed are the decision-relevant d-dimensional continuation,
gauge/constraint determinants, source contacts, causal response, complete
hard/rate conversion and omitted-order bounds. The RATE4-COVZERO-1
boundary proposal remains unadopted and is not yet an executable common
quantum parent.

## 4. Why full J can fail while the finite-angle arc exists

Retain mass one, E=1/256, q_max=1, s0=2 and the original weight
`psi(q)=q phi(q)`, `phi(q)=1-3q/2+q^3/2`. For a split T>=32 and
C,alpha>0, construct the spectral **fixture** at s>=T:

```text
b=alpha(s-4)log(s/T)/2,
rho(s,x)=C s^2 exp(-b) cosh(b x).
```

It is even in x. Rodrigues' formula and integration by parts give

```text
integral_-1^1 exp(bx) P_l(x) dx
 = b^l/(2^l l!) integral_-1^1 exp(bx)(1-x^2)^l dx >=0.
```

Consequently every even-spin coefficient of the fixture is nonnegative
and every odd coefficient vanishes. This is an all-spin proof. At x=1,
`rho>=C s^2/2`, so

```text
J >= (C/pi) integral_T^infinity s^2/(s-3)^3 ds
  >= (C/pi) integral_T^infinity ds/s = infinity.
```

For the actual angle x=1-2q^2/(s-4), however,

```text
0 <= rho(s,x) <= C s^2 (s/T)^(-alpha q^2),  E<=q<=1.
```

The reflected exponential is no larger because `s-4-q^2>=q^2` for
s>=32. The bound in the next section proves a finite arc. This fixture
has positive partial-wave spectral coefficients; it is **not** asserted
to obey the full unitarity upper bound, energy crossing, a common analytic
amplitude, or to be a UV completion. It establishes only that the previous
positive-reference assumptions plus finite-angle convergence do not
force full-J finiteness. It does not establish that the actual parent's
J diverges or refute the earlier conditional J/80 estimate.

## 5. A finite-window comparison and a finite-angle tail bound

The weight and infrared arc are unchanged, so the previous four-rate
projection and its conditional D8-tail estimate still apply. No additional
normalization condition is introduced by splitting the spectral integral.

Use the positive all-spin comparison only on the finite interval [32,T).
Its proof is pointwise in s and therefore gives the same error constant:

```text
J_32,T=(2/pi) integral_[32,T) rho_positive(s,1)/(s-3)^3 ds,
a_E >= -C_E J_32,T - R_tail - delta_contour - delta_observable,
C_E<1/80.
```

The finite moment is still required physical input; finite energy alone
does not prove infrared finiteness of a forward gravitational observable.
The high-energy interval is [T,infinity), avoiding double counting any
spectral atom at the split. An empty finite window has zero moment.
The positive-reference mismatch and low-cut/IR conversion errors must
be included in delta_observable. Regulated limits require uniform bounds.
All errors refer to the normalized arc, with no double counting.

For the high-energy tail, suppose BOTH actual cut discontinuities have
the uniform finite-angle envelope

```text
|rho(s,-q^2)| <= C s^2 (s/T)^(-alpha q^2) + remainder,
s>=T>=32, E<=q<=1, alpha>0.
```

Let delta_Regge bound the complete normalized integral of the remainder.
This is a specified input to derive, not a universal Regge theorem. The
two unequal crossed denominators obey

```text
(s-2)^(-3)+(s-2-q^2)^(-3) <= 2/(s-3)^3,
[s/(s-3)]^3 <= [T/(T-3)]^3,
integral_T^infinity (ds/s)(s/T)^(-alpha q^2)=1/(alpha q^2).
```

Therefore Tonelli's theorem for the nonnegative envelope gives

```text
R_tail <= [2/(pi W_E)] [T/(T-3)]^3 (C/alpha) P_E + delta_Regge,
P_E=integral_E^1 phi(q)/q dq
   =log(1/E)-4/3+3E/2-E^3/6.
```

No q->0 limit or infinite forward moment was used. At the unchanged
E=1/256, `P_E<6`, `W_E>1/11`, `2/pi<2/3` and T>=32 give

```text
R_tail < 60 C/alpha + delta_Regge  (C>0).
```

For C=0 the leading term is zero. The executable bound uses a sharper
rational dyadic logarithm majorant. The complex contour remains a
separate error: a real-cut envelope alone does not bound that contour.
Likewise this tail estimate does not by itself prove a positive spectral
reference at the physical finite coupling.

### 5.1 What the hierarchy does and does not buy

IF a calculation established `C/alpha<=gamma/kappa`, the leading tail
would be below `60 gamma/kappa`. With the original kappa=10^800 and
lambda=10^-600, even gamma<=10^196 would put this piece below lambda/100.
This is a sensitivity calculation and a sufficient target, **not** an
assumed order-one coefficient, a physical gamma bound, or the whole
gravity allowance. The tail, graviton pole, low cuts and contour must
still enter the same physical dispersion relation before taking limits.

Raising T alone does not suppress the bound when C is defined at that T:
the energy integral is independent of T apart from its kernel prefactor.
Changing the spectral split must transport the envelope normalization.
T is an integration split, not a new physical cutoff or detector choice.

The primary literature motivates, but does not supply, these inputs:
[Bellazzini et al., sections 4-5](https://arxiv.org/html/2512.13780v2#S5)
distinguish the finite-resolution observable from its limiting unitary
reference and retain finite-coupling errors. The published
[gravitational Regge bounds](https://arxiv.org/pdf/2202.08280)
have additional semiclassical assumptions and a d>4 scattering scope;
they are not numerical four-dimensional bounds for this parent. No
Froissart bound or massless example is imported into the present massive
calculation.

## 6. Closure decisions and reproduction

| Question | New answer | Still required |
|---|---|---|
| Must every affine determinant be computed as a propagating loop? | No: 56 directions reduce to an ultralocal factor in the stated dimensional measure convention. | Actual constrained light/metric/vector measure, dimensional continuation and finite conversions. |
| Is the added heavy source the next one-loop bounce-response obstacle? | Not at the reference mean/quadratic jets: its source-dependent contribution starts too high in field degree. | Free heavy metric response, regularity and higher-order/domain bounds. |
| Can the retained Proca source be dropped on the clock? | No: its quadratic source generates a quartic exchange with a nonzero spatial-momentum control. | The constrained, state-preserving one-loop response and its error bound. |
| Is a finite full forward J necessary for this research route? | No: finite-window comparison plus a finite-angle tail can replace it. | Actual finite-window moment, uniform tail/envelope remainder, contour and observable errors. |

**M/V/G/B/R and original P8 remain open; G remains UNTESTED.** The original
scoped A and linear classification are unchanged. No user choice of an
unknown numerical coefficient is needed to continue the identified
calculations. Adoption of a new physical boundary proposal is separate.

Run the [read-only diagnostic](../scripts/p8_reference_regge.py):

```sh
.venv/bin/python -B scripts/p8_reference_regge.py
.venv/bin/ruff check scripts/p8_reference_regge.py
.venv/bin/ruff format --check scripts/p8_reference_regge.py
```

The diagnostic pins inherited inputs, checks block Gaussian and actual
source-jet identities, a nonzero Proca symbol, reference transport,
Regge antiderivatives, dyadic logarithm bounds, finite Rodrigues controls,
budget sensitivity and rejection of missing physical inputs. The written
proofs carry the functional-degree and all-spin/all-energy statements.
Evidence level is VERIFIED_N, not FORMALIZED or physical MATCHED status.
All ten focused matching/prescription diagnostics and their lint/format
checks pass. The new diagnostic protects 342 inherited inputs; all 25
pre-existing edited/untracked files outside the three living status
documents were retained byte-for-byte. The frozen P8 science is unchanged.
The historical 115,487-test release snapshot was not rerun for this
standalone diagnostic. This work remains local, uncommitted and unpushed.
See the [validation receipt](validation/p8-reference-regge-2026-09-24.json).
