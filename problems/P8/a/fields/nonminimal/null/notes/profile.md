# Continuous profile and full nonminimal null-stress estimate

The estimate covers every real |s|<=1, not sampled times.
Write v=B v0 for the functions in state.md. Exactly

    v0(0)=2, v0'(0)=1.

The coefficient triangle bounds are |C1|<=4+1/epsilon,
|C2|<=2+1/epsilon, hence

    |v0''(s)| <=12epsilon^2+5epsilon
       =1253/250000 =D2.

Taylor's integral formula on the real segment yields

    |v0(s)-(2+s)|<=D2/2,
    |v0'(s)-1|<=D2,
    |v0(s)|<=3+D2/2,
    |v0'(s)|<=1+D2.

Because b is even and nonnegative of mass one, B is real,
|B|<=1, and on the same segment

    |B-1|<=width^2/2, |B'|<=width^2, |B''|<=width^2.

These are the integrated elementary cosine/sine bounds and the
second moment of b. The three actual mollified profile estimates are

    |v-(2+s)|<=D2/2+(width^2/2)(3+D2/2),
    |v'-1|<=D2+(width^2/2)(1+D2)+width^2(3+D2/2),
    |v''|<=D2+2width^2(1+D2)+width^2(3+D2/2).

Each rational right-hand side is STRICTLY less than 1/100.
No unknown small-epsilon remainder or numerical quadrature is used.

Let v=x+iy and set e0=1/100. These imply

    1-e0<=x<=3+e0, |y|<=e0,
    x'>=1-e0, |y'|<=e0, |x''|,|y''|<=e0.

From the actual squeezed covariance, divided by A_R^2,

    w=(2/3)(-x^2+3y^2),
    <:(partial_ell Phi)^2:>=(2/3)(-x'^2+3y'^2).

The full null contraction of A.20's nonminimal covariant stress
is checked before restriction:

    T_ell,ell=(partial_ell Phi)^2-xi partial_ell^2(Phi^2).

Thus its expectation, divided by A_R^2, is exactly

    (2/3)(1-2xi)(-x'^2+3y'^2)
       +(4xi/3)(x x''-3y y'').                           (1)

In particular the improvement is not omitted or replaced by
the minimal tensor when xi=1/6.

Define the positive rational numbers

    A=(1-e0)^2-3e0^2,
    D=(3+e0)e0+3e0^2.

For 0<=xi<=1/4, (1) is bounded above by
-(2/3)(1-2xi)A+(4xi/3)D, increasing in xi. At its upper endpoint
this is (-A+D)/3. The exact resulting bounds are

    w/A_R^2 <=-1633/2500 < -3/5,
    T_ell,ell/A_R^2 <=-4747/15000 < -3/10,
    at xi=1/6: T_ell,ell/A_R^2 <=-4823/11250 < -2/5.

Together with A_R^2=hbar R^2 gamma0^2 these prove uniform
segment negativity and the divergence of every fixed nonzero
compact g^2-weighted null average. The same R family works for
all such samplers. Its strict negative w satisfies ANY fixed
nonnegative one-sided upper cap. It does not satisfy a
uniform two-sided bound on |w|.
