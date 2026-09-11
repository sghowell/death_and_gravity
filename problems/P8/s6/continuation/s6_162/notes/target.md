# The actual finite-kappa target and its normalization

Use the S6.109 rank-regular target, not S6.108's older polynomial
switch. With n=1024, kappa=10^800, lambda=10^-600, write

    T=X^n/[X^n+(1-X)^n],
    w=n X^2 exp(-n X^2),
    B=T+(1-T)w, S=(1-T)(1-w),
    R=1+B(X-1)/(1+u^2)^3.

The exact retuned tree function is imported literally from the
frozen parent. Its denominator is 200(1+u^2)^12 and its value
at the origin is F0=-28/25. Set

    Fv=(X-u^2)/2+(lambda kappa-n F0)X^2+(-n/3-F0)u^4,
    F=Ftree+exp(-u^4) S (Fv-Ftree),
    A3=R_X/X,
    A4=-R_X/X-7 R_X^2/(4R),
    A5=R_X^2/(R X).

The apparent divisions by X are removable: R_X is divisible
by X at the vacuum. The canonical conversion is
u=Phi/sqrt(kappa), X=Y/kappa, where Y=(partial Phi)^2.
The dimensionless rolling clock X=1 therefore has canonical
gradient invariant kappa, not one.

At fixed flat physical metric the canonical scalar density is

    L=kappa F+(A3 L3+A4 L4)/kappa+A5 L5/kappa^2,
    L3=(Box Phi) Z,
    L4=(raised gradient)^T Hessian eta Hessian (raised gradient),
    L5=Z^2,
    Z=(partial^mu Phi)(partial^nu Phi) partial_mu partial_nu Phi.

The curvature term vanishes on this prescribed fixed metric;
its metric variation does not vanish. No gravity or stress
estimate is obtained by discarding that variation.

The zero and linear vacuum terms vanish. The quadratic term
is (Y-Phi^2)/2. All coefficient functions and contractions
are even under simultaneous reversal of the scalar jets.
