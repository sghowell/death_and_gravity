# S0 reference: gauge operator, ghost ultraviolet input and response correction

Date: 2026-09-29. This continues the already approved S0 and COVZERO
construction. It does not change the physical candidate, state, cutoff,
field domain, finite boundary rule or finish line.

**Result:** the dimensional clock-gauge algebra, actual de Donder ghost
operator, an all-momentum retarded inverse bound on the bounce, and the
ordinary ghost's local ultraviolet polynomial are now explicit. The full
coupled light/metric quantum reference and physical response bounds remain
**IN PROGRESS**. This is a required measure component, not a substitute
known-sector smallness estimate or a claim of P8 closure.

**Correction:** section 3.1 of the
[previous S0 packet](assessment-2026-09-28-p8-s0-reference.md)
has an exact **source-aligned** Hamiltonian difference. Its unqualified
extension to the complete held-physical-vector second response is not
valid: the old normal-vector one-point embedding contact must also be
retained. Section 1 below supplies it from the existing S253 calculation.
The old packet, script, specification and receipt remain frozen; the
living ledger carries the correction. S0 adoption, its constraints,
source-free vector probe matrices and parity are unaffected.

Evidence: [read-only diagnostic](../scripts/p8_s0_gauge_measure.py) and
[validation receipt](validation/p8-s0-gauge-measure-2026-09-29.json).
The 45 accumulated P8 files were published first as
`bb4632f7f9dbd0f4a7e99768ab37f8147de04d9d`; unrelated P4/P9 work was excluded.
That publication does not certify an overbroad physical interpretation.

## 1. Restore the held-vector embedding contact

The pinned [S253 held-vector algebra](../problems/P8/s6/continuation/s6_253/src/p8_vacuum_affine_physical_background_vertices/coupled.py)
already distinguishes the background paths. For the old parent,

```text
Wbar_aligned = S = 3*(R-1)*(Hhat-H_clock).
```

At the reference, S and its first derivative vanish. For a homogeneous
physical lapse/scale probe (eta,v_phys), with time u=t in clock units,

```text
r1 = R_N|clock = -2/(1+t^2)^3,
S_second = 6*r1*eta*d_t(v_phys+r1*eta/4).
```

This expression gives the same-direction second derivative; polarization
gives mixed probes. It follows by differentiating the old background
embedding, not by setting the actual R or its derivatives to constants.
At the bounce, smooth compact probes can have eta=1, eta_dot=0 and
v_phys_dot=1, giving S_second=-12.

For any finite-regulator functional Gamma, the chain rule is

```text
Gamma_old,aligned,AB = Gamma_old,held,AB + Gamma_old,W*S_AB,
Gamma_old,held,AB - Gamma_S0,held,AB
  = (Gamma_old,aligned,AB - Gamma_S0,aligned,AB) - Gamma_old,W*S_AB.
```

S0 alignment is Wbar=0, so its aligned and held paths coincide. The
two-first-vertex kernel in the earlier comparison is unaffected because
S_first=0. Its source-aligned second contact, however, is **not the entire
held response contact**. For the Gaussian Hamiltonian contribution,
Gamma_old,W=-expectation(h_W), so the missing term has sign
**+expectation(h_W)*S_AB**. Complete measure, other-sector and local
counterfunctional contributions to Gamma_old,W must also be carried;
they have not been computed or assigned zero here.

An independent envelope derivative of the full held temporal-vector square
agrees entry by entry with S253. At D=Z=1 it gives

```text
nstar = [Theta*(p_v+3c*sigma)-w*p_sigma-Lnv*v]/(2J),
h_W = rN*nstar*(3Theta*nstar-3c*sigma/2-p_v/2)
      +p_v*A_L/(2P)+(Z_N*nstar+3v)*P*p_A.
```

The first term is a purely light quadratic operator. The original centered
light/vector product kills the cross-block means, not this light mean.
In particular, Lnv=2E*q+3T must not be replaced by zero. At the actual
grade-zero bounce with q=1, J=243/160, Theta=0, c=1/10, w=1/20,
Lnv=-1 and rN=-2, that term is

```text
[(80v-4p_sigma)/243]*(p_v+3sigma/10).
```

Its v*p_v coefficient is 80/243. This is a nonzero-operator control,
**not** a calculation of its expectation in the original state. Neither a
nonzero original tadpole nor its cancellation by another sector is claimed.
The old profiles retain their original formal grades; the tree control
does not set the fixed higher-grade profiles to zero in the full problem.

The diagnostic also checks the full S253 physical phase-density map,

```text
Cphase = sqrt(kappa)*diag(1,1,sqrt(zeta),a^3,a^3,a^3/sqrt(zeta)),
H_W,physical = kappa*a^3*Cphase^(-T)*H_W,normalized*Cphase^(-1).
```

Both original scalar boundary momentum shifts remain in state transport.
The displayed 1/P term is on P>0 with the inherited physical-shift infrared
domain, not an unrestricted zero-mode formula. Exact S0 vector parity
still gives zero odd vector mean in its even reference; it does not erase
the old sourced-parent one-point term or the S0 gravitational response.

## 2. Continue the actual clock gauge in dimension

Put d=4-2epsilon and m=d-1. On the timelike scalar-unitary patch, use

```text
Q_m(gamma) = det(gamma)^(1/m)*gamma^(-1),
chi^i = partial_j Q_m^ij,       density weight w=2/m.
```

The already selected light/metric continuation has
h_physical=C*gamma, C=R^(-1/(m-1)). Therefore Q_m(h_physical)=Q_m(gamma)
for every integer m and by the stated analytic continuation. This retains
the physical/auxiliary conformal cancellation rather than assuming it
survives an unchanged four-dimensional gauge expression.

Differentiating the contravariant density Lie derivative gives the full
operator, before imposing chi=0,

```text
M xi = xi.grad chi - chi.grad xi + w*chi*div(xi)
       -Q^jk partial_j partial_k xi^i
       +(w-1)Q^ij partial_j div(xi).
```

The lower-order chi terms are necessary for off-gauge vertices. For Q_ell,
incoming ghost momentum k and outgoing p=k+ell, its full Fourier matrix is

```text
M(p,k)=(p.Q_ell.k)I+(Q_ell k)p^T-(Q_ell p)ell^T-w*(Q_ell p)k^T.
```

The incoming k=0 column is -(Q_ell ell)ell^T and need not vanish off the
slice. The outgoing p=0 row vanishes. This is consistent with the existing
[S260 full algebra](../problems/P8/s6/continuation/s6_260/notes/algebra.md),
not permission to remove an incoming column before differentiation.

For Q positive and k nonzero, r=k.Q.k,

```text
sigma(M)=rI+(m-2)/m*(Qk)k^T,
det sigma=2*(m-1)/m*r^m,
M_flat^(-1)=[I-(m-2)/(2*(m-1))*k*k^T/k^2]/k^2.
```

The determinant is the rank-one determinant lemma; Q^(-1) symmetrizes the
principal matrix. On the exact slice and a periodic/decaying domain,
integration by parts gives the flat form norm(grad xi)^2 plus
(1-2/m)*norm(div xi)^2. For integer m>=2 and
delta=norm(Q-I)_op, the perturbation costs at most (m-1)*delta times
norm(grad xi)^2. In physical m=3, delta<=1/8 gives the inherited lower
bound 3/4. On a torus whose first nonzero momentum is 1/L, the mean-zero
inverse estimate is at most 4L^2/3. This is local gauge admission, not a
global Gribov theorem or a new physical field restriction. No Hilbert
space/coercivity assertion in noninteger dimension is made.

### 2.1 Evanescent terms cannot be discarded

If Q_fixed=det(gamma)^(1/3)*gamma^(-1) is copied unchanged instead, then

```text
Q_fixed(C*gamma)=C^(m/3-1)*Q_fixed(gamma),
d_epsilon log[C^(m/3-1)]|epsilon=0 = log(R)/3.
```

Even though R=1 on the clock, its first lapse probe is nonzero:
partial_N(log R/3)|clock=-2/[3(1+t^2)^3]. Also, at gamma=a^2 I,
Q_fixed=a^(2m/3-2) I, with epsilon derivative -4log(a)/3 of its logarithm.
These O(epsilon) probe and volume terms can multiply ultraviolet poles.
The continued gauge has

```text
d_epsilon(2/m)|0=4/9,
d_epsilon[(m-2)/(2*(m-1))]|0=-1/4,
d_epsilon log[2*(m-1)/m]|0=-1/3.
```

An alternate fixed-exponent gauge is not intrinsically forbidden; it must
retain its evanescent vertices and full finite reference transport.
What fails is identifying it with the same conformally cancelling gauge
while dropping those terms. This packet adopts no extra finite physical
boundary coefficient through a gauge change.

## 3. Temporal clock Jacobian and its limit

Geometrically, delta u=A_* xi_perp+A_i xi^i. On u=t, A_i=0 and
A_*=sqrt(X)>0, so the time/spatial Faddeev-Popov block is triangular:

```text
FP = [[A_*,0],[v,M_spatial]],      det FP=A_* det M_spatial.
```

The lower-left leg v remains in the inverse and response vertices; only
the determinant factors. In a local canonical clock chart, the Jacobian
from integrating delta(H_perp) over its conjugate clock momentum cancels
the clock factor. The bracket must be the correctly reduced/Dirac bracket,
not a naive bracket before auxiliary constraints are imposed.

A separate finite auxiliary control makes this distinction explicit:
H=p_u*A-A^2/2+E(u,q,p), with second-class pair (p_A,A-p_u).
The Dirac clock bracket is p_u; the full second-class-plus-clock constraint
matrix has determinant p_u^2. On the positive branch,
abs(p_u)*delta(p_u^2/2+E) dp_u cancels the clock solve Jacobian. This toy
is a check of the mechanism, not a replacement P8 Hamiltonian.

The covariant auxiliary-gradient Hamiltonian in
[Deffayet and Garcia-Saenz, section II.1](https://arxiv.org/html/2004.11619)
provides the relevant normal clock generator; its signature conventions
are not copied blindly. The result here does not construct the complete
regularized BFV measure or discard a second-class determinant by declaring
it ultralocal. At du=0 the unitary clock gauge is inadmissible; the previously
admitted covariant vacuum chart is still required.

## 4. Actual de Donder operator and causal inverse

Use the selected physical-metric background de Donder gauge, parameter one.
For the linear physical metric perturbation h, define
F_mu=nabla^nu h_numu-(1/2)nabla_mu h. Directly substituting
h=Lie_xi g gives the vector ghost operator

```text
M xi = Box_P8 xi + Ricci_P8 xi.
```

This fixes the convention by metric variation rather than a curvature-sign
guess. On ds^2=dt^2-a^2 dx^2, H=a_dot/a, q=k^2/a^2,

```text
(M xi)^0 = xi0_ddot+mH xi0_dot+[q-m(H_dot+2H^2)]xi0
           -2iH k_i xi^i,
(M xi)^i = xii_ddot+(m+2)H xii_dot+q xii-2iH k_i xi0/a^2.
```

All components and the orthonormal transport are independently checked by
Christoffel/metric variation in d=3,4,5. The geometric identity gives the
general dimensional form. There is no division by Theta, J or A_*;
at the flat vacuum it is the ordinary vector wave operator.

Set y=(xi0,a*xii), Forth=(f0,a*fi). The equation is

```text
y_ddot+mH*y_dot+(qI+V+K)y=Forth,
V00=-m(H_dot+2H^2),       Vii=-H_dot-(m+1)H^2,
K0i=Ki0=-2iH*k_i/a,      other K entries zero.
```

For the actual four-dimensional clock a=(1+t^2)^2, exact bounds are

```text
abs(H)<=2,
V00=-12*(1+7t^2)/(1+t^2)^2,       abs(V00)<=49/2,
Vii=-4*(1+15t^2)/(1+t^2)^2,       abs(Vii)<=225/14<49/2.
```

The two potential bounds follow respectively from the nonnegative squares
(7t^2-5)^2/[2(1+t^2)^2] and (15t^2-13)^2/[14(1+t^2)^2]
for the upper bound plus the negative potential. With
E=abs(y_dot)^2+(1+q)*abs(y)^2, differentiating gives four homogeneous
contributions: damping <=12E; I-V <=(51/2)E; the K cross term <=4E;
and q_dot=-2Hq <=4E in absolute value. Hence

```text
E_dot <= (91/2)E + 2*sqrt(E)*abs(Forth),
sqrt(E(t)) <= exp[91*(t-s0)/4]*sqrt(E(s0))
              + integral_s0^t exp[91*(t-s)/4]*abs(Forth(s)) ds.
```

This proves a momentum-uniform finite-interval retarded estimate, including
k=0. The stated weighted energy loses no spatial derivatives for this
ghost inverse. Twelve exact Hermitian-slack controls supplement the
inequalities; they are not used to infer a continuum bound from a grid.
At the bounce the time component of M acting on a constant time vector is
-12, explicitly rejecting an unchanged flat/scalar ghost operator.

For compact/retarded sources, zero Cauchy data selects a unique inverse.
Choosing xi=-M_ret^(-1)F(h) sets the linear de Donder condition and preserves
the initial preparation neighborhood when the source vanishes there.
It does **not** generally vanish at the final CTP endpoint. Final sewing,
source observables and residual gauge-volume normalization must be
transported as well; no second endpoint is silently reset. The bound is
not uniform in infinite elapsed time, a positive ghost state, or a norm
bound for the coupled physical response.

In particular, the [S251 reduced physical light block](../problems/P8/s6/continuation/s6_251/notes/ultraviolet.md)
has two distinct scalar principal cones, one strictly subluminal and one
luminal. It is not thereby the same normally hyperbolic operator as this
ghost. No ordinary Einstein-graviton kernel is transplanted into S0.

## 5. Ordinary ghost local ultraviolet polynomial

This section computes only the ordinary complex Grassmann de Donder pair
in the target gauge, not the coupled bosonic/constraint measure. Use the
local Euclidean comparison operator P_gh=-(nabla^2 I+Ricci), so E=+Ricci.
The inherited heat convention is R_o=-R_P8. This comparison supplies local
UV coefficients, not a global Euclidean rotation or a replacement state.

The universal Laplace-type coefficient in
[Vassilevich, equation (4.28)](https://arxiv.org/pdf/hep-th/0306138)
includes the divergence terms. Here tr I=d, tr E=R_o,
tr E^2=Ricci^2 and tr Omega^2=-Riemann^2. Therefore

```text
a4_scalar = R_o^2/72+(Riemann^2-Ricci^2)/180+Box_o R_o/30,
a4_gh(d) = d*a4_scalar+R_o^2/6+Ricci^2/2-Riemann^2/12+Box_o R_o/6.
```

The diagnostic checks the bundle trace on independent non-Einstein
algebraic curvature tensors in d=3,4,5. Restricting only for comparison
to Ricci^2=R_o^2/d and Box_o R_o=0 gives

```text
a4_gh,Einstein = (5d^2+58d+180)*R_o^2/(360d)
                 +(d-15)*Riemann^2/180.
```

This agrees with [Bastianelli et al., equation (2.41)](https://doi.org/10.1007/JHEP10(2023)152);
their ghost action and relative weight are in (3.19) and (2.38).
The actual bounce is not assumed Einstein, and their pure-gravity bosonic
coefficients are not used for the S0 coupled block.

Relative to a real boson's half determinant, the ghost pair has weight -2.
In the inherited pole-polynomial normalization before 1/(64*pi^2), this is

```text
P_gh(d)=-4*a4_gh(d),
P_gh(4)=-8R_o^2/9-86Ricci^2/45+11Riemann^2/45-6Box_o R_o/5.
```

Keeping d=4-2epsilon in the numerator gives the explicit component-count
evanescent finite contribution

```text
-2*partial_d P_gh = 8*a4_scalar
  = R_o^2/9+2*(Riemann^2-Ricci^2)/45+4Box_o R_o/15.
```

This is **one identified finite contribution**, not the complete finite
ghost determinant, covariant measure conversion or boundary functional.
It is already present if the full dimensional expression is expanded;
do not add it twice. Invariants and regulator transport can supply other
terms. Keep the surface terms and vary before any four-dimensional
Gauss-Bonnet or on-shell simplification. A temporary massless infrared
regulator used to identify UV poles would require its own cancellation;
no physical cutoff or mass is introduced here. Ghost polynomials are not
standalone physical energy or positivity observables.

## 6. On-clock determinant agreement does not fix response

Differentiating the exact gauge identity E_i R^i_alpha=0 gives

```text
K_ij R^j_alpha + E_j partial_i R^j_alpha = 0.
```

Thus K R=0 is an on-stationary-background simplification. It cannot be
differentiated as if the E term vanished on nearby physical probes.
A finite gauge-invariant control, with no continuum/regularization issue,
is S=(x^2+y^2-1)^2/2, R=(-y,x), evaluated at (r,0). Gauge F_b=(b,1)
and parameter one have ghost determinant r. If
D_b=det(K+F_b^T F_b), the one-loop difference is
DeltaGamma1=(1/2)log(D_b/D_0). Exactly at r=1,

```text
DeltaGamma1=0,
partial_r DeltaGamma1=b^2/2,
partial_r^2 DeltaGamma1=-(13b^2+b^4)/2.
```

Equality of the Gaussian determinant on the clock is therefore
insufficient to identify its first or second physical derivatives. This
control is not a P8 anomaly or a proof that its physical observable is
gauge dependent. It explains why the existing
[S260 source/state Nielsen transport](../problems/P8/s6/continuation/s6_260/notes/ward.md)
and [endpoint preparation](../problems/P8/s6/continuation/s6_260/notes/state.md)
must be implemented, not bypassed with on-shell determinant counting.
The BRST locality theorem of
[Barvinsky et al., section 2.3](https://arxiv.org/pdf/1705.03480)
assumes an appropriate invariant/anomaly-free regulated measure and local
divergences; citing it does not verify those hypotheses for this system.

## 7. Next closure-directed calculation and validation scope

The next required calculation is the **coupled S0 light/metric/constraint
measure bridge** in the selected dimensional continuation, with its
off-shell source and endpoint transport. Combine its actual bosonic and
second-class contributions with the explicit ghost block; establish the
needed local pole/finite conversion relative to the approved COVZERO
boundary; then compute or bound the same-state causal projections used by
the four rates and B. The S253 held-vector correction remains required
when comparing to the old sourced parent, not a new source term in S0.

Only required projections need bounds. No global field-space gauge theorem
or all-orders UV completion is added to the contract. The separately needed
same-observable finite-gravity error certificate is still absent. S0 and
COVZERO approvals are resolved; no further user decision is needed for
this ongoing reference construction. Physical M/V/G/B/R/P8 are not closed.

The focused validation covers 17 root diagnostics, including exact new
algebra and analytic inequalities, with 517 pinned inputs. It checks 49
preserved preexisting files, including the old S0 packet and unrelated
P4/P9 work, plus lint/format, source hashes and local links. The historical
115,487-test covering snapshot is **not rerun**. A successful replay of a
frozen report is not endorsement of a superseded physical interpretation.
No background research is claimed after handoff.
