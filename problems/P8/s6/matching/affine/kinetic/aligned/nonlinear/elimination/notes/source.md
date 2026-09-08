# Full source Taylor and spatial-curl remainder

Use the actual frozen source, not a freely chosen forcing. On each
scalar-clock slice u is constant, h=(1+u^2)^3>=1, |H|<=2 and
h_phi/h=3H/2. Write y=x+1, where x=-N^-2. The ODE is

    Q_x+A Q=F,   Q(u,-1)=0,
    A=1/(2x)+3/[4(h-y)],
    F=-3h_phi*y/[4h^2(h-y)].

Its exact linear forcing coefficient is f1=F_x(u,-1)
=-3h_phi/(4h^3)=-9H/(8h^2). Consequently

    F-f1*y = f1*y^2/(h-y).

On |y|<=1/10, the frozen S6.42 bounds are |A|<=25/18 and
|Q|<=(45/31)y^2/h^2. Also h-y>=9h/10 and |f1|<=9/(4h^2).
Thus |F-f1*y|<=(5/2)y^2/h^2. For
R_Q=Q-f1*y^2/2, integration from y=0 on either side gives

    |R_Q| <= (140/93)|y|^3/h^2,
    |R_Q,x| <= (140/31)|y|^2/h^2.

The second inequality follows directly from R_Q,x=F-f1*y-A Q;
it is not inferred by differentiating an upper bound. This is a
continuous Taylor remainder of the original ODE solution.

## Normal and coordinate sources

Let r=N-1, k=K_hat-3H, N=1+r, s=N^-1 and
y=r(2+r)/N^2. The exact source and its quadratic expression are

    S=-y*k/h-3H*y*(1-s)/h+(3/2)s Q,
    S2=-2r*k/h-[6H/h+27H/(8h^2)]r^2,
    S0=N S.

For rho=1/25 and |r|<=rho, X=N^-2 lies strictly in [9/10,11/10],
N>=24/25 and |y|<(9/4)|r|. Exact identities in source.py give

    S-S2 = r^2(3+2r)k/(hN^2)
          +3H r^3(5+6r+2r^2)/(hN^3)
          -(3/4)f1 r^3(16+39r+40r^2+20r^3+4r^4)/N^5
          +(3/2)R_Q/N,

    S0-S2 = r^2 k/(hN)+3H r^3(3+2r)/(hN^2)
           -(3/4)f1 r^3(12+23r+16r^2+4r^3)/N^4
           +(3/2)R_Q.

The four normal coefficients are bounded respectively by 4,36,37,27;
the coordinate coefficients by 2,21,26,26. Each is checked using an
exact rational numerator envelope at rho and the positive lower N
denominator, with the ODE remainders above. Dropping h^-1,h^-2<=1,
for |r|,|k|<=epsilon<=rho this proves

    |S-S2|<=104 epsilon^3,    |S0-S2|<=75 epsilon^3.

Also |S2|<=(83/4)epsilon^2 and |S|<=25 epsilon^2.
The upper bounds are conservative, not estimates from samples.

## Actual spatial derivative and local curl density

Write R0=S0-S2. In the clock chart no spatial derivative acts on
h,H or f1. Exact differentiation, including R_Q,x y_r, gives

    R0_r = y*k/h+3H r^2(9+11r+4r^2)/(hN^3)
          -(3/2)f1 r^2(18+40r+40r^2+20r^3+4r^4)/N^5
          +3 R_Q,x/N^3,
    R0_k = r^2/(hN).

The three r^2 coefficients in R0_r are <65,82,78, totaling 225.
Its k coefficient is <3|r|/h; |R0_k|<=2r^2/h. In any one common
positive physical spatial norm, if also ||D r||,||D k||<=epsilon,

    ||D R0|| <=230 epsilon^3,
    ||D S2|| <=(83/2)epsilon^2.

The full electric source covector is E=-s D S0; define E2=-D S2
using the same actual physical spatial metric for both norms.
Since |s-1|<=25|r|/24 and s<=25/24,

    ||E-E2|| <= (4525/16)epsilon^3 <300 epsilon^3,
    ||E|| <=54 epsilon^2.

There is no magnetic component because the full source is normal.
The first local term suggested by the exact ordered elimination is
-zeta F(S)^2/4=zeta ||E||^2/2 per physical volume. Its difference
from zeta ||E2||^2/2, in the same actual metric, is at most

    (zeta/2)||E-E2|| (||E||+||E2||)
       <=14325 zeta epsilon^5.

This bounds substitution of the full source into a specified local
curl density. It is not a bound on the remaining nonlocal induced
action, on its variation, or on the difference between two metrics.
It also does not bound a time derivative merely from a spatial norm.

Physical normal-source errors divide by tau, electric errors by tau^2,
and normalized action-density errors multiply by M^2/tau^2. With
zeta=zeta_physical/(M^2 tau^2), the last physical density bound is
14325 zeta_physical epsilon^5/tau^4. Nonunit tests keep these factors.

## Preparation and Fourier boundary

If r=r'=r''=0 initially, the full S and S0 have those same three
zero initial time jets, and S'''=-2r''' k_initial/h. Smooth u-dependent
coefficient variations cannot enter this third jet. The S6.44
nontrivial prepared profiles therefore keep the exact source's
preparation when used as specified N,K_hat fields. This does not turn
their linear light solution into a full nonlinear parent solution.

Even the full source's coefficient of k is rational in r. For
r=epsilon*cos(theta), it is not a finite Fourier polynomial. Higher
harmonics and light-dependent vector coefficients must be controlled
before extending the old fixed-output, diagonal-mode response theorem.
The explicit source neighborhood is not a quantified neighborhood
for nonlinear stability or for S6.45's auxiliary implicit function.
