# Complete covariant curvature variations

Use the -+++ convention of the fixed finite action. All spatial
coordinate partial derivatives vanish on the homogeneous family.
Covariant spatial derivatives generally do not vanish.

The implementation constructs the full ordered inverse time jets,
all Christoffel jets0..3, Riemann(1,3) jets0..2, Ricci jets and
scalar jets. For rank-r covariant T,

    (nabla_a T)_I=delta_a0 T'_I
                  -sum_over_r_slots Gamma^j_(a i) T_(I:i->j).

Its second derivative differentiates this entire expression and
subtracts connections on all r+1 indices, including the first
derivative index. For a scalar,
Hess(R)_ba=delta_a0 delta_b0 R''-Gamma^0_ba R'.
In particular the spatial Hessian is not dropped.

For compact metric variations the lower-index covariant Euler
tensors, whose raised versions pair with delta g_mu_nu, are

    E_R=-(Ric-g R/2),
    E_R2=2 Hess(R)-2 R Ric+g(-2 box(R)+R^2/2),
    E_Ric2=Hess(R)-box(Ric)-2 Ric^(pq) R_mu_p_nu_q
             +g(-box(R)/2+Ricci^2/2).

They follow from the full curvature variation and integration
by parts. Independent tests also derive their spatial pairing
by differentiating the scalar density with respect to independent
metric jets0,1,2 and applying the alternating Euler time derivatives.
This second calculation uses a generic positive non-diagonal
spatial metric and does not repeat the Christoffel implementation.

Here are conservative bounds in the factorial-weighted amplitude
norm, with enough remaining time jets for each displayed derivative.
From N(dt^r g),N(dt^r g^-1)<G,

    Gamma_j <= (3/2) 2^j G^2, j=0..3,
    Riem_j <= 2 Gamma_(j+1)
               +2 sum(l=0..j) binomial(j,l) Gamma_l Gamma_(j-l)
             <100 G^4, j=0..2.
    Ric_j <100 G^4, R_j <400 G^5, j=0..2.

The scalar bound includes all inverse-metric product derivatives.
For any rank-r covariant tensor with time bounds T0,T1,T2,
a complete second covariant derivative is bounded by

    T2+r Gamma1 T0+r Gamma0 T1
       +(r+1) Gamma0 (T1+r Gamma0 T0).

Using Gamma0<2G^2, Gamma1<3G^2 gives
Hess(R)<1200G^7, box(R)<1200G^8 and box(Ric)<5000G^9.
Raising Ricci twice gives100G^6; lowering Riemann once gives100G^5.
Ricci^2<10000G^10 and R^2<160000G^10. Every term in each Euler
tensor is then bounded before contraction:

    E_R <300G^6,
    E_R2 <1e6 G^11,
    E_Ric2 <1e6 G^11.

tensors.constants() records the unrounded sums and checks every
strict display. These are bounds for the full tensors, not only
their tracefree projections or their highest derivative terms.
