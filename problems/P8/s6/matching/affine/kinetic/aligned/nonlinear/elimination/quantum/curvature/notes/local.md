# Local curvature calculation and actual rolling metric

## Literal vector Hessian, not a different background

The original coefficient dictionary gives p(phi,x=-1)=1/2. Substitution
into the unchanged mass response gives a=b=1 for every real clock
time, and their derivatives along that clock vanish. Metric compatibility
therefore makes the on-clock mass tensor m0^2*g covariantly constant.
The aligned source vanishes on the background; it does not enter the
fixed-light Gaussian Hessian. That Hessian is ordinary Proca on the
actual physical metric. This does not make the off-clock mass variations
vanish: S6.47's clock coefficient jets explicitly reject that inference.

In Euclidean signature define delta as the adjoint of the exterior
derivative d. On one-forms, K=delta*d+m0^2. The Hodge operator is
D1=d*delta+delta*d=-nabla^2+Ricci; D0=delta*d=-nabla^2 on scalars.
Since d^2=0,

    K*(I+d*delta/m0^2)=D1+m0^2.

The determinant of I+d*delta/m0^2 equals the corresponding scalar
determinant, with contact and zero-mode factors handled consistently.
An independent finite cochain complex verifies this factorization and
the determinant identity *including* its nontrivial contact factor.
Its scalar constant mode is kept. Discarding that factor in this finite
control fails at nonunit mass; the continuum prescription matters.

Equivalently introduce the redundant variables A,psi in
W=A-d*psi/m0, add the gauge-fixing square (delta*A+m0*psi)^2/2,
and the associated determinant of D0+m0^2. The mixed quadratic
terms cancel by adjointness. The A and psi Gaussian determinants,
with that gauge determinant, leave one vector minus one scalar trace.
With scaleless contact traces zero in dimensional regularization,

    Gamma=1/2 Tr log(D1+m0^2)-1/2 Tr log(D0+m0^2).

This is not massless Maxwell (vector minus two scalars) and is not
three scalar determinants. Boundary conditions and global zero-mode
normalizations must respect this prescription. The independent S4
check below uses this ratio literally, with no scalar mode omitted.

## Universal input and independent tensor contractions

For D=-(nabla^2+E), the local coefficients excluding (4*pi*s)^-2 are

    a0=tr I,
    a2=tr(E+R*I/6),
    a4=tr[ (5R^2-2Ricci^2+2Riemann^2+12boxR)*I
            +60R*E+180E^2+60boxE+30Omega_mu_nu*Omega^mu_nu ]/360.

This is the standard Laplace-type input, with its smooth-manifold
and local-interior scope; see [Vassilevich, Eqs. 4.26-4.28](https://arxiv.org/abs/hep-th/0306138).
For one-forms E=-Ricci. The implementation constructs a general
four-dimensional algebraic curvature tensor with twenty independent
parameters and checks tr E=-R, tr E^2=Ricci^2 and
tr(Omega_mu_nu*Omega^mu_nu)=-Riemann^2 by explicit matrix contractions.
In particular the last sign is negative. The scalar has E=Omega=0.
Subtracting its coefficients gives

    a0=3, a2=-R/2,
    a4=-R^2/8+29*Ricci^2/60-Riemann^2/15-boxR/15
       =-GB/15+13*Ricci^2/60-7R^2/120-boxR/15,

where GB=Riemann^2-4Ricci^2+R^2. Both forms are checked exactly.

The proper-time representation is -1/2 integral ds/s exp(-m0^2*s)
times the heat trace. Expanding the mass exponential through s^2,
the logarithmic UV coefficient is a4-m0^2*a2+m0^4*a0/2. Hence in
d=4-2epsilon_DR,

    Gamma_pole=-1/(32*pi^2*epsilon_DR) integral sqrt(g)
                 [3*m0^4/2+m0^2*R/2+a4].

The flat coefficient is exactly S6.47's -3*m0^4/(64*pi^2*epsilon_DR).
This fixes the convention directly, without importing a pole sign from
a source that prints the opposite definition of epsilon.

## Direct physical FLRW geometry

The code constructs Christoffel symbols and Riemann/Ricci contractions
from g=diag(-1,a^2,a^2,a^2), with a_dot=a*H, not from a de Sitter
replacement. It verifies R=6(H_dot+2H^2),
Ricci^2=12(H_dot^2+3H_dot*H^2+3H^4) and
Riemann^2=12[(H_dot+H^2)^2+H^4], including zero Weyl squared.
For the unchanged a=(1+u^2)^2, H=4u/(1+u^2) in normalized time,

    R=24*(1+7u^2)/(1+u^2)^2,
    Ricci^2=192*(1+10u^2+37u^4)/(1+u^2)^4,
    Riemann^2=192*(1+6u^2+25u^4)/(1+u^2)^4,
    boxR=48*(63u^4-14u^2-5)/(1+u^2)^4,
    GB=1536*u^2*(1+3u^2)/(1+u^2)^4.

Here boxR=-R''-3H*R' is the physical Lorentzian scalar d'Alembertian.
The resulting local invariant coefficient is

    a4=-8*(77u^4+14u^2-3)/(1+u^2)^4.

These are evaluations of local counterterm functions. Their covariant
Lorentzian use for renormalization does not establish a finite
Lorentzian determinant or control its nonlocal/state dependence;
that distinction is also explicit in [Vassilevich, Sec. 2.3](https://arxiv.org/abs/hep-th/0306138).

Keep the divergence. Direct differentiation proves
a^3*boxR=-(a^3*R')' and a^3*GB=(8a^3*H^3)'. Neither boundary
flux tends to zero at infinite cosmological time: they grow like u^9.
Dropping boxR would instead give
-8*(259u^4+98u^2-5)/[5*(1+u^2)^4]. There is no assertion of an
integrable full spacetime action on this noncompact background.
Local variational arguments require their own boundary conditions.

## All-time coefficient bounds, with units

The exact square identity
49-R=(7u^2-5)^2/(1+u^2)^2 gives the sharp 0<R<=49 on all real u.
For other coefficients use P(u^2)/(1+u^2)^n. Express it as the
weighted mean of c_j/binomial(n,j), with nonnegative weights
binomial(n,j)*u^(2j)/(1+u^2)^n summing to one. The largest absolute
coefficient is therefore a continuous, not sampled, bound. Exact
reconstruction is checked for every expression. In particular
a4's coefficients at n=4 are (24,-28,-308/3,0,0), so |a4|<=308/3.

Restore R_physical=R/time_scale^2 and a4_physical=a4/time_scale^4.
For Rm=m0*time_scale>0, the absolute curvature addition to the
pole weight divided by its flat 3*m0^4/2 is bounded by

    49/(3*Rm^2)+616/(9*Rm^4).

For Rm=1000 this is 18375077/1125000000000. It is a ratio of
local pole coefficients, not a finite quantum uncertainty. In particular
the pole must be renormalized; its residue is not itself an observable.
No finite matching coefficient, full quantum stress, scalar constraint,
cone, cutoff, vacuum/Regge condition or original P8 verdict follows.
