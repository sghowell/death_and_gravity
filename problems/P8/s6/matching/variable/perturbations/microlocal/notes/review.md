# Uniform scalar reduction: proof of P8-S6.24.MICROLOCAL

This note supplies the written proof accompanying the exact replay. Its
action and physical metric are exactly those of frozen S6.20/S6.22. In their
dimensionless units, `u=T/tau`, `K=(tau*k_com)^2`, `k=sqrt(K)`, and the
physical spatial wave number is `k/(tau*a)`. The positive branch has
`2<c<=4`, `|u|<=1/10`, and a strictly positive canonical clock kinetic
coefficient. The free chi source remains free; beta(phi) sources the
other canonical scalar. No derivative operator is added.

## 1. Exact differential field map

Starting from the frozen unitary-clock, common-spatial gauge, define

```
pi = -e/K,
xg = a² B_g,
xf = (b²/c²)(B_f-pi'),
q_g = psi_g-h*xg,     q_f = psi_f+h*xf,
z = chi_B-w*xg,
n_g = Phi_g-xg',      n_f = Phi_f-xf'.
```

The map is invertible as a differential change of eight perturbation
fields. `pi=E_f-E_g` is the relative spatial scalar, `xg` is the physical
clock displacement, and `chi_B` is the physical g-metric Bardeen matter
perturbation. Unlike a flat-curvature clock gauge, no formula divides h.
The two lapse variations and both original shift equations are retained.

The literal substituted action has **zero K² coefficient identically**.
Write it as `K L1+L0`. Exact polynomial primitives F1,F0 obey

```
L1bar = L1-d_u F1,
L0bar = L0-d_u F0.
```

L1bar contains no velocity except pi'; L0bar contains no xf'. The latter
identity is essential before substituting a value of xf containing pi'.
Otherwise an artificial acceleration is introduced by the substitution.
All background derivatives in these boundaries are taken before u=0.

## 2. The leading stationary fields

Set `U=a² y P/2`, `V=a²P/(2y)`, where `P=2 beta1` in the normalized action.
The five leading stationary fields are

```
psi_g=-U*pi,        psi_f=V*pi,
Phi_g=-U*pi,
Phi_f=-(1-2y/c)*V*pi,
xf=(2/c)*pi'+xg+J(u,c)*pi.
```

The source-derived J and all other expressions are exposed without a
center substitution by `leading.derive()`. This is a stationary solution
of L1bar only, not the exact solution of all full constraints.

After this substitution, L1bar still contains `f(u,c) pi*pi'`. Its boundary
is `Fpi=f pi²/2`. The final leading gradient is obtained from
`L1red=L1bar|stat-d_u Fpi`, not from the value of L1bar|stat at the bounce.
At c=4, f(0)=0 but f'(0)=352. Omitting this derivative changes the helicity
gradient from 800 to 448. The latter number is an explicit failed control.

## 3. Exact completion and auxiliary weights

Subtract the leading stationary fields, naming the curvature remainders
`r_g,r_f`, the lapse remainders `lambda_g,lambda_f`, and the xf remainder
`t_f`. The exact leading completion is

```
L1bar = L1bar|stat + a*r_g²+c*b*r_f²
       +2a*r_g*lambda_g+2c*b*r_f*lambda_f+(C_x/2)*t_f²,
C_x = a*c²*P/(c+y) > 0.
```

There is no physical-field/remainder cross term in this completed K
piece. After the substitution, L0bar is still first order: only the three
physical velocities and the two curvature-remainder velocities occur.
It has smooth rational background coefficients for each fixed c>2.

The useful semiclassical weights are

```
r_g,r_f,t_f = k^-1 times order-one amplitudes,
lambda_g,lambda_f = order-one amplitudes.
```

Divide the physical Euler equations by k². Divide the two curvature
equations by k², and the two lapse and t_f equations by k. In columns
`(lambda_g,lambda_f,k*r_g,k*r_f,k*t_f)` the leading auxiliary block is

```
diag(2a,2cb,2a,2cb,C_x).
```

It has a uniform inverse on every compact parameter rectangle with
`c>=2+delta_min>2`, independently of bounded normalized time frequency
omega/k. The remaining auxiliary differential terms have at least one
power of k^-1. Derivatives of background coefficients also have that
suppression in this semiclassical scaling. This explains why substituting
the stationary fields into the order-one velocity Hessian is the candidate
physical leading kinetic calculation: errors in the curvature, temporal
and lapse responses contribute below the k² physical equation order.

This finite-cone parametrix argument is not, alone, a proof excluding an
additional nonuniform fast branch. Sections 5–6 complete it using the exact
six-dimensional reduced Cauchy system and a uniform rational-symbol
inverse. They do not assume that freezing the unreduced differential-
algebraic equations preserves its time-dependent secondary constraint.

## 4. Principal matrices and independent exact bridge

In q=(pi,xg,chi_B), the stationary action gives

```
Lprincipal = q'^T K_s q'/2 - K*q^T G_s*q/2.
```

At u=0,

```
K_s=diag(30720/[c²(c-2)²], (6400-801c)/(100c), 1),
G_s=diag(256(9c²-22c+144)/[c²(c-2)²], (6400-801c)/(100c), 1).
```

Thus the center helicity speed squared is
`(9c²-22c+144)/120`, equal to 5/3 at c=4 and tending to 17/15 as c->2+.
These are not the discarded 14/15 boundary-omission result or the old
nonuniform Psi_f-chart value 5.

An independent exact route takes the frozen canonical phase Z and defines
q=OZ with pi=-e/K, xg=3(a³wz-P_E)/(2aK), chi_B=z+w*xg. It retains O', O''
and the full Hamiltonian derivative in the exact q'' equation. At c=4,u=0,

```
det C = 5(K+60)(K²+26K-24)/(19176 K³),
Mq/K -> diag(-5/3,-1,-1),    Mv=O(1),
Sigma_qv -> diag(480,799/100,1),    Sigma_vv=O(K^-2).
```

This chart is not valid at every finite K; its center zero near K=0.89
must not be discarded. It is regular at the center for K>=1. Point checks
at signed punctures have the same degree orders and agree with the action
kinetic matrix. They are independent fixtures, not the proof of uniformity.

## 5. Exact generic symbol and a constructive inverse threshold

The native rational calculation keeps u,c,K indeterminate throughout the
pinned original Hamiltonian, time-dependent canonical boundary, secondary
Schur elimination and physical Cauchy map. It proves the rational identities
and degree bounds

```
Mq = -K K_s^-1 G_s + Rq,        deg_K Rq <= 0,
deg_K Mv <= 0,
Sigma_qv = K_s + O(K^-1),
Sigma_vv = O(K^-2),             Sigma_qq = O(1),
deg_K detC = 0,                leading(detC)*det(K_s)=1.
```

These are polynomial identities over Q(u,c,K), not interpolation or
assertions inferred from numeric samples. The backend is exact rational
arithmetic using native multivariate FLINT polynomials. Separate checks
compare its arithmetic and matrix inversion against SymPy and the root
audit reconstructs the literal scalar action/Cauchy system independently.

Fix `0<delta_min<=2`, `2+delta_min<=c<=4`, `|u|<=1/10`. Write
`detC=N(u,c,K)/(K³ D0(u,c))` in the exact reduced rational normalization
used by bounds.py. Its cubic numerator coefficient is

```
n3 = c^4 d^18 (c*d^4-2)^2 (u²-1)/344064.
```

Consequently `|n3| >= 99 delta_min²/2150400`. Direct weighted coefficient
norms of the three lower numerator polynomials on `|u|<=1/10, |c|<=4`
give `|n0|<16`, `|n1|<14`, `|n2|<1`. For

```
K >= K_star(delta_min) = 2000000/delta_min²
```

the relative correction to the cubic leading term is smaller than
`6944/20625 < 1/2`. Thus detC has its nonzero leading sign throughout the
whole parameter box. The original phase map is already regular there:
its only old constraint denominators are the pinned positive clock
coefficient, nonzero P, finite K and the strictly negative D of S6.22.
No pole introduced by writing a reduced rational fraction is promoted to
a physical singularity.

Here is why the generic degree cancellations give uniform removability,
not merely pointwise limits. Before the physical Cauchy inversion, the
original phase generator, O and C have finite Laurent dependence on K;
their coefficient denominators are the nonvanishing background factors
just described. Use `C^-1=adj(C)/det(C)` in this original representation,
rather than declaring every factor of a separately reduced fraction to be
a physical denominator. After `zeta=1/K`, its only additional denominator
is a power of

```
Nhat = n3+n2*zeta+n1*zeta²+n0*zeta³,
|Nhat| >= |n3|/2 > 0.
```

Consequently every exact equation and conserved-form entry is a finite
Laurent numerator divided by smooth, nonvanishing background factors and
powers of Nhat. The generic degree identities state that the coefficients
of every forbidden negative power of zeta vanish identically, throughout
the parameter box. These are exact polynomial identities, so the
cancellations remain valid at h=0 and all other points of the box. For
K_s^-1 the required background determinant is bounded away from zero on
each fixed delta_min box by Section 6. Thus Rq, Mv and the appropriately
rescaled symplectic remainders extend smoothly to `zeta=0`, with all
needed time derivatives bounded on the compact set
`(u,c,1/K) in [-.1,.1] x [2+delta_min,4] x [0,1/K_star]`.
This supplies uniform O bounds, not just pointwise degree statements.
The exact six-state Cauchy evolution contains no omitted additional fast
physical scalar branch: all of its coefficients have the displayed
uniform first-order wave scaling after the following energy reduction.

The explicit K_star is deliberately loose. It is a mathematical symbol
threshold, not a computed physical cutoff, nor a claim that this range is
inside an EFT regime.

## 6. Continuous positivity and uniform exact evolution

Use the nonzero time-dependent field weight `Pi=U*pi` and leave xg,chi_B
unchanged. Let `T=diag(U,1,1)` and

```
Kbar = T^-T K_s T^-1,     Gbar = T^-T G_s T^-1.
```

These rational principal matrices extend smoothly to c=2 on the closed
coefficient box `|u|<=.1,2<=c<=4`. This does not extend the parent action
or the time connections to that singular value. The clock diagonal is
exactly `a³*kphi` in Kbar and `a*kphi` in Gbar; it is positive by the
pinned clock proof. The chi diagonals are a³,a, with no principal cross
terms involving chi. Four strict polynomial Bernstein certificates give

```
trace(Kbar_2), trace(Gbar_2) < 64,
det(Kbar_2), det(Gbar_2) > 16.
```

There are respectively 87,190,112,222 coefficients, all positive, after
`u²=x/100`, `c=2+2z`, `(x,z) in [0,1]²`. All clearing denominators are
explicitly positive. The positive clock diagonal and the positive two-
minor imply positive definiteness. The trace/minor inequalities bound
both eigenvalues between 1/4 and 64. The uncoupled chi eigenvalue is also
in this interval. Therefore

```
(1/4) I < Kbar,Gbar < 64 I
```

on the complete coefficient box, with no sampling in time or c.

For the exact normalized physical equation, retain the time connections
of T. Its form remains

```
Q'' = [-K Kbar^-1 Gbar + R] Q + V Q',
Q=(Pi,xg,chi_B),
```

where R,V and their required time derivatives are uniformly bounded at
fixed delta_min for K>=K_star. They include T'/T and T''/T; dropping those
terms would not be the exact equation. Set `W=(k Q,Q')`,
`S=diag(Gbar,Kbar)`, and `E=W^T S W/2`. Its exact first-order generator is

```
W' = k [[0,I],[-Kbar^-1 Gbar,0]] W + [[0,0],[R/k,V]] W.
```

The order-k part is skew with respect to S. Hence

```
|E'| <= C_delta E,
C_delta = 4 sup ||S'|| + 512 sup (||R||+||V||) < infinity,
E(u) <= exp(C_delta |u-u0|) E(u0).
```

The suprema are over the stated compact coefficient/1/K box; the proof
establishes finiteness and explicit rational formulas, not a numerically
optimized expanded value of C_delta. This is a genuine energy bound
uniform in all K above the stated threshold, through h=0. It permits
finite lower-order growth; it does not assert an instantaneous ground
state or positive exact Hamiltonian at every finite K.

This is the constrained, spatially weighted Goldstone energy. Reconstructing
the original metric components can multiply by fixed finite powers of k;
the result is not a uniform bound on every unweighted raw metric amplitude.
Those fixed spatial derivative weights are distinct from the discarded
H=0-dependent degree jump. A prescribed physical-source response would
still require its own source normalization and constraint-compatible data.

The action and conserved-form identities also fix the high-frequency
kinetic sign, not just a positive artificial norm. In the energy phase
W, the conserved two-form multiplied by k tends to
`[[0,Kbar],[-Kbar,0]]`; its other corrections vanish uniformly. The
positive Kbar is the principal canonical kinetic weight, with the same
orientation as the retained canonical matter. This proves positive
scalar principal kinetic and gradient forms in this precise sense, not
nonlinear or all-band health.

At the center the physical g speed equals the displayed eigenvalue
because a(0)=1. In general physical squared speeds are the eigenvalues
of `a² K_s^-1 G_s`. The center helicity excess satisfies

```
cs_helicity² - 17/15 = (c-2)(9c-4)/120 > 0  for 2<c<=4.
```

Thus the positive scalar principal sector is not fully subluminal with
respect to the prescribed matter metric. The two center matter speeds
are unity. This is a full-system principal statement, not a statement
about the result of integrating out the relative mode on a finite band.

## 7. Limits of the result

The proved formal full-parent principal cone is not
a finite-band light-only EFT matching verdict. The limit k->infinity at
fixed c>2 differs from c->2 with a shrinking inner window. Normalizing
Pi=U*pi removes the pole from the candidate principal matrices, but the
time connection U'/U is not bounded uniformly in the joint c->2 limit.
In particular K_star grows as delta_min^-2, whereas the tensor algebraic
mass squared at the center grows only as (c-2)^-1. The uniform principal
gate cannot be used as a low-energy heavy-mode matching computation.
No cutoff, physical mass gap, finite-band response, loop health or original
DHOST/free-matter operator match is supplied by this scalar calculation.
