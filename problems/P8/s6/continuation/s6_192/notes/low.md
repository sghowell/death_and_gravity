# Actual finite-band covariance response

The exact real balanced evolution matrix, derived from the constrained
three-mode dynamics, is

    A=[[R+S,omega I],[-omega I,R-S]].

R is real skew and S is real symmetric. Hence
(A+A^t)/2=diag(S,-S). S6.190 gives ||S||<=11 on the fixed unit slab,
so the exact propagator satisfies ||U(t,s)||<=exp(11|t-s|) at every
momentum. The fast omega and R part is norm preserving; this is not
an adiabatic gap or a restriction to a finite momentum grid.

The unchanged S55 preparation has |f|^2<=9/omega and |p|^2<=27omega.
At the common initial surface the balanced covariance consequently has
trace at most 3(9+27)=108, and positive covariance has operator norm
at most its trace. The whole initial shear neighborhood is common,
so its first two amplitude derivatives are zero in this same frame.

Let nu=nu_minus>=1000. The S6.191 mixed jets imply

    ||A_e|| <= 3nu+212.04+216.545 < 4nu,
    ||A_ee|| <= 6nu+10617.64+10640.665 < 28nu.

For Sigma'=A Sigma+Sigma A^t, differentiation gives

    Sigma_e'=A Sigma_e+Sigma_e A^t+A_e Sigma+Sigma A_e^t,
    Sigma_ee'=A Sigma_ee+Sigma_ee A^t+A_ee Sigma+Sigma A_ee^t
               +2 A_e Sigma_e+2 Sigma_e A_e^t.

Combine time exponents before integrating: exp(22(t-s))exp(22s)
is exp(22t), not two independent full-slab factors. Zero initial
parameter data and Duhamel give, for slab length T<=1,

    ||Sigma|| <=108 exp(22T)<4e12,
    ||Sigma_e|| <=108 exp(22T) 2(4nu)T <4e13 nu,
    ||Sigma_ee|| <=108 exp(22T)[2(28nu)T+4(4nu)^2 T^2]
                  <3e14 nu^2.

The last term is the iterated first response and cannot be discarded.

Use |tr(XY)|<=6||X||||Y||, omega<=3nu, and all three metric-vertex
bounds from vertices.md, including the frequency derivatives.
Complete product differentiation gives

    |J_D|<1e14 nu, |(J_D)_e|<1e15 nu^2,
    |(J_D)_ee|<1e16 nu^3.

These are exact-state finite-momentum bounds, not ultraviolet
majorants by themselves. With b=.99/(25/16)^2 and
nu^2=m^2+b|k|^2, b^(-3/2)<4 and pi^2>9, radial measure is bounded by
(4/18)nu^2 dnu. Thus the complete band nu<K has integral bounds

    (4/18) CURRENT[a] K^(4+a)/(4+a), a=0,1,2.

Extending the lower integration endpoint from m to zero only enlarges
the positive integral. No physical mode or longitudinal piece is removed.
