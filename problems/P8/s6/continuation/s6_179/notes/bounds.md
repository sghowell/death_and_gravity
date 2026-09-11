# Explicit current norm, with all contacts retained

The physical source norm is the maximum over both sources and
time derivatives0..10 on I. With the analogous clock-chart norm,
the exact Leibniz formula for v=vhat+w1 n gives

    ||(n,v)||_C10 <= M ||(n,vhat)||_C10,
    M=46090764897/8.

Each derivative of w1 is rational in t with a positive power of
1+t^2 in its denominator. The continuous coefficient majorant
is recomputed on [-1/2,1/2], and all eleven polynomial
reconstructions vanish. No sampling grid defines this bound.

Pointwise write r=1/h. Then 64/125<=r<=1,
0<w1=r/2<=1/2 and

    w2=r^2-3r/2=(r-3/4)^2-9/16,
    |w2|<=9/16, |w1+w2|=r(1-r)<=1/4.

The stress output map J^T diag(-1,3) has row norm at most3.
For unit clock-chart source norm, |v|<=3/2. The contact rows
therefore obey, using the exact unperturbed physical stress,

    |contact_N| <= (9/2)|rho|+(15/2)|P|,
    |contact_vhat| <= (33/2)|P|.

Let B_rho,B_P be the S6.176 zero-derivative physical stress
bounds divided by the fixed kappa, and Rphys the maximum of
the two actual S6.85 physical response bounds. Then

    Rclock = 3 M Rphys
             + max((9/2)B_rho+(15/2)B_P, (33/2)B_P)
             < 1e-784.

The code uses exact source-pinned rationals for every term;
both nonlocal pieces, matched local pieces and both background
contacts remain in the report. No old profile response or
background-cancelled output is included. Rphys<5e-795.

These bounds refer to physical stress and reference-volume
mean-current variations, respectively. Both are normalized by
kappa=1e800, never by the vanishing bounce density or by H.
The result is a derivative-losing C10-to-C0 operator bound,
not a bounded inverse on the same function space.
