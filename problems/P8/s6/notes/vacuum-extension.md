# S6.1 — Vacuum data are not fixed by the smooth clock tube

## 1. Exact, tube-preserving D family

Use P8's (+---), X=(partial phi)^2 convention. For the algebraic replay set
M=tau=1 and phi=u. Write F_D and a3_D for the pinned witness, and choose
Z>0, m>0, any finite u0 and any real lambda. Set

```
r(y) = 0 for y<=0, exp(-1/y) for y>0
B(X) = r(X-1/4)/(r(X-1/4)+r(3/4-X))
F_v = Z*X/2 - Z*m^2*(u-u0)^2/2 + lambda*Z^2*X^2
F_ext = B*F_D + (1-B)*F_v
F2_ext = -1/2; K_ext=A1_ext=A2_ext=0
A3_ext = B*a3_D
A4_ext = -A3_ext + X^2*A3_ext^2/4
A5_ext = -X*A3_ext^2.
```

Dependent Ia coefficients are recomputed, **not linearly spliced**. The
denominator F2-X*A1=-1/2 never vanishes. Since a3_D>0 at every finite u,
the smallest supporting optional group is still D, not the baseline.

The denominator defining B is everywhere positive: both r arguments cannot
be nonpositive. B=0 for X<=1/4 and B=1 for X>=3/4. Thus the entire original
tube 9/10<=X<=11/10 lies strictly inside an open set where **every function
and every derivative equals the old action**, for all u. All existing
background, principal, tail and finite-order interaction results there are
identical. This is local equality, not a nonlinear/global stability claim.

For completeness, r is C-infinity at zero. Its nth derivative for y>0 has
the form exp(-z)*P_n(z), z=1/y, with

```
P_0=1; P_(n+1)=z^2*(P_n-P_n').
```

Induction makes P_n a polynomial for every n. For each monomial z^k,
exp(z)>=z^(k+1)/(k+1)! implies z^k exp(-z)<=(k+1)!/z ->0 as z->infinity.
Every right derivative therefore matches the zero left derivative. The
positive denominator then makes B smooth. The replay checks eight derivative
recurrences and exact branch/degeneracy identities; the induction and bound
above, not those eight examples, prove the all-orders smoothness statement.

In physical units the same construction uses
F_D_phys=M^2/tau^2 F_D(phi/tau,X), A3_D_phys=M^2 a3_D(phi/tau),
F2=-M^2/2 and
F_v_phys=Z X/2-Z m^2(phi-phi0)^2/2+lambda Z^2 X^2.
Here [phi]=-1, [Z]=4, [m]=1 and [lambda]=-4 in mass units.
The dependent terms become A4=-A3+X^2 A3^2/(4M^2) and
A5=-X A3^2/M^2. This gives the same nonzero Ia denominator for M>0.

## 2. Vacuum and independent scalar test

Near X=0 the action is exactly Einstein gravity plus F_v. At phi=phi0:
F_v=F_v,phi=0, F_v,X=Z/2>0, F_v,phiphi=-Z*m^2. No second-clock-derivative
operator or varying Planck coupling remains there. With
psi=sqrt(Z)*(phi-phi0), its scalar Lagrangian is

```
L_v = (partial psi)^2/2 - m^2*psi^2/2
      + lambda*((partial psi)^2)^2.
```

It has a healthy massive quadratic scalar and the usual massless tensor for
either sign of lambda. This is perturbative vacuum health, not global
Hamiltonian boundedness or nonlinear health at arbitrary gradient. Taking
M->infinity while holding Z,m,lambda fixed is a scalar decoupling limit **in
this vacuum patch**; it is not a decoupling theorem along the bounce.

For all incoming k_i, sum k_i=0, k_i^2=m^2, write s=(k1+k2)^2,
t_transfer=(k1+k3)^2 and w=(k1+k4)^2, so s+t_transfer+w=4m^2.
Expanding four labelled plane waves directly in lambda*((partial psi)^2)^2
gives the contact amplitude (S=1+iT convention)

```
A = 8*lambda*((k1.k2)*(k3.k4)+(k1.k3)*(k2.k4)+(k1.k4)*(k2.k3))
  = 2*lambda*((s-2m^2)^2+(t_transfer-2m^2)^2+(w-2m^2)^2).
b2 = (1/2)*partial_v^2 A(v,t_transfer=0)|v=0 = 4*lambda,
v = s+t_transfer/2-2m^2.
```

There are no scalar cubic interactions in this vacuum patch. After gravity
decouples there is no light exchange pole in this tree calculation. Under
the *vacuum* dispersive hypotheses and controlled perturbative errors, b2
is nonnegative (strict with nonzero spectral weight). Positive/negative
lambda therefore give opposite leading scalar positivity diagnostics while
leaving the **same whole clock tube** untouched. A finite positive coefficient
is only one necessary test. At lambda=0, tree saturation is not a standalone
exclusion. Finite M and loop/higher-order errors have not been set to zero
in the physical P8 problem.

Independent known answer: add a healthy heavy scalar h of mass M_h with
interaction (g/2)*h*(partial psi)^2. Tree exchange gives
g^2/4 sum_z (z-2m^2)^2/(M_h^2-z), z=s,t_transfer,w. Eliminating h at leading
derivative order gives lambda_eff=g^2/(8M_h^2). The exact forward coefficient
is g^2/[2(M_h^2-2m^2)], agreeing at leading order and retaining its finite-mass
correction. Assume M_h>2m. This derivative-coupled mediator is a tree-level
matching/sign regression, **not a proposed complete UV theory**.

## 3. A stronger analytic branch has a different answer

Suppose instead F_ext is real analytic on a connected open domain containing
both the whole specified open tube and a proposed (phi0,0) vacuum, and equals
F_D there exactly. F_D itself is real analytic for all real u,X. The real
analytic identity theorem forces F_ext=F_D throughout that connected domain.
The pinned exact check has

```
F_D(u,0) = -(10u^14+62u^12+162u^10+263u^8
             +288u^6+186u^4+34u^2+19)/(1+u^2)^8 < 0.
```

It cannot satisfy the necessary flat-vacuum metric equation. Hence **exact
tube equality plus connected global real analyticity plus a finite
constant-clock Minkowski vacuum is impossible for this D witness**. This
does not exclude analytic models agreeing only to a finite jet or within an
error tolerance, extra-field completions/thresholds, other vacua, other D
witnesses, or the entire D row. A splice is nonanalytic at its switching
surfaces, exactly as that theorem requires. Vacuum momentum analyticity
does not, by itself, supply the global field-space hypothesis used here.

## 4. Consequence and unfinished gate

The smooth tube alone determines neither a finite constant-clock vacuum nor
the sign of its leading scalar forward coefficient. A separate healthy
vacuum can be attached algebraically, but this is **not a matching proof**.
In particular the transition 1/4<X<3/4 has not been shown to have a uniformly
healthy spectrum, a controlled heavy gap, or to arise from any common UV
parent. The vacuum's positive test does not certify the bounce's UV class.

The selected next ansatz is the positive-lambda member, keeping Z,m,lambda
symbolic until matched scales can be justified. S6.2 must relate this vacuum
patch to the tube in an actual controlled EFT construction or reject that
ansatz. Gravitational dispersion residuals, all-order power counting,
radiative stability, M1 and P8(a) remain open. No new user choice is needed
to begin those specified investigations.
