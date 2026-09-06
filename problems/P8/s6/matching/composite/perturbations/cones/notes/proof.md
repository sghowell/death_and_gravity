# Matter-cone ordering and the all-degeneracy obstruction

All action coefficients G,F,m4,alpha,beta and beta_n are constants in the
specified class. A field-dependent conformal coupling or Einstein
coefficient would change the equations and is not included.

## 1. Literal physical tensor characteristics

S6.7 derives the tensor kinetic and spatial-gradient terms directly from
the two Einstein actions, with their full positive lapse factors. In
physical time, their kinetic coefficients per tensor contraction are

    Kg,T=G*a^3*s/8,  Kf,T=F*a^3*y^3*s/(8c).

The spatial terms are respectively `G*a*k^2/(8s)` and
`F*c*a*y*k^2/(8s)`. The relative potential mu is algebraic. Derivatives
of the rolling normalization and this algebraic mixing do not alter the
principal symbol. Comparing with `k_phys=k/(a*r)` therefore gives

    cg²=(r/s)²,  cf²=(c*r/(y*s))².

These are formal tensor principal speeds of the two-field system, not
finite-momentum mass eigenvalues. Both polarizations have the same pair
of characteristic speeds. There is no high-k statement about a separate
low-energy truncation or its unknown validity range.

All ratios are positive. Hence `cg<=1` is equivalent to y<=c, while
`cf<=1` is equivalent to c<=y. Together they force c=y. An alternative
exact proof uses the positive weights

    wg=alpha/r, wf=beta*y/r, wg+wf=1,
    wg*cg+wf*cf=1.

Two speeds both at most one can have this weighted average only when
both equal one. This is a velocity identity, not an average of squared
speeds. The independent polynomial check verifies both order factorizations.

## 2. Retain the pressure in the relative stiffness

The pinned, source-aware expression is

    mu=y{m4[beta1+beta2(y+c)+beta3*c*y]-alpha*beta*r*s*p}.

For arbitrary constant potential coefficients, exact algebra yields

    mu=y*Q+y*(c-y)*[m4(beta2+beta3*y)-alpha*beta²*r*p].

Thus c=y gives mu=y*Q. In the regular Xi>0 vector chart,
`cV,eff²=(r/s)²*mu/Xi`. Strict positivity of that speed and c=y require
mu>0 and hence Q>0. In particular a pressure-branch slice Q=0 cannot
simultaneously meet both tensor-cone inequalities and the strict vector
condition. Omitting the homogeneous matter pressure would invalidate
this conclusion; its contribution was derived before imposing the branch.

Xi=0 is not a hidden escape inside the stated contract. In the undivided
two-shift action, its relative-shift coefficient is then zero and shift
elimination leaves zero relative-vector inertia. Xi<0 gives negative
inertia at sufficiently large formal k, beyond a singular constraint
denominator. These cases are outside the regular positive all-k chart.
No physical subcutoff ghost scale is inferred from this observation.

## 3. Connected-interval proof, including degenerate bounces

The assumptions apply at every time. Step 1 therefore gives c=y on the
whole interval, and step 2 gives Q>0 there. The full undivided Bianchi
identity Q*B=0 implies B=0 everywhere. This never divides by a Hubble,
the matter velocity or a polynomial near a root: Q has already been
proved strictly positive on the required domain.

At c=y the actual coordinate identity is

    ydot=B/(Ng*a).

Consequently y, and hence r, is constant. B=0 also gives Hg=y*Hf and
H_eff=Hg/r. The physical lapse is Ng*r. The full g null equation reads

    -2G*Dg Hg=alpha*r^3*(rho+p)+(y-c)*Q.

The last term vanishes, so conversion to physical time gives

    dH_eff/dT=(Dg Hg)/r²=-alpha*r*(rho+p)/(2G)<=0.

This proves monotonicity, not only a sign for the acceleration at an
isolated stationary point. A smooth nonincreasing function cannot change
from negative to positive even through a degenerate zero, zero interval,
or a nonisolated zero set. No assumptions on the multiplicity of Q are
needed. The result applies to any prescribed potential whose admitted
solution satisfies the canonical null inequality.

## 4. Controls on the assumptions

If y=2,c=1 and alpha=beta=1, the two squared tensor speeds are 9/4 and
9/16; if y=1,c=2 they are 4/9 and 16/9. These are algebraic cone fixtures,
not asserted background solutions. They show why checking only one tensor
cone does not suffice. If beta=0, the physical metric is proportional to
g, and its cone alone does not force c=y; strict positivity of both affine
couplings is therefore an explicit hypothesis.

A stronger control is an actual local analytic pressure-branch solution
when zero vector speed is allowed. Use the S6.6 beta2 reconstruction ODE
but prescribe the state-dependent physical Hubble

    h(y,rho)=[X(y,rho)+y²Y(y,rho)]/(1+y)².

Its reconstructed lapses obey Nf=y*Ng identically. Starting at y=1,
rho=1/2 gives positive lapses 1/2, positive null source one, h=0,
and the exact actual-clock derivative h'=1/6>0. The analytic ODE and the
positive scalar velocity provide an analytic locally reconstructed
potential just as in the pinned reconstruction argument. In a small
open neighborhood it is a nondegenerate composite bounce with both
tensor cones equal to the matter cone. However Q=0 and mu=0 identically:
the relative vector inertia stays positive near the initial point, while
its principal speed is identically zero. It deliberately fails the
strict vector premise. It is not the frozen free or CD solution, nor
claimed stable, globally complete, or controlled as an EFT.

The theorem cannot be reinterpreted as a general UV exclusion by declaring
the massive principal cone a light-field cone. A physical heavy threshold,
validity window and controlled reduction must first be established for
that different assertion. Nor does mu=0 alone determine a rolling
canonical gap, as the pinned normalization identity already demonstrates.
