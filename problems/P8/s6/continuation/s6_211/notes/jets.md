# Actual four-degree mode and source-derivative jets

Put x=1/r, k_hat=n and ell_hat=-n+xP; use P=p e3 when evaluating the exact azimuthal coefficients. All intermediate coefficient expressions keep symbolic t,m,p and u=n.e3.

Normalize Obar=a*x*omega_physical and Wbar=a*x*W_physical. Then

    Obar_k^2=1+m^2 a^2 x^2,
    Obar_l^2=1-2pux+(p^2+m^2a^2)x^2,
    z=1-m^2a^2x^2/Obar^2.

Both positive analytic frequency branches have constant1. The actual reference gives, through degree4,

    Wbar=Obar+a^2 x^2 P1(t,z)/Obar
                 +a^4 x^4 P2(t,1)/Obar^3.

Only P1(t,1), its first z derivative, and P2(t,1) enter; the second z change in P1 and all P3/P4 terms are degree6 or higher. This filtration is a proof about the UV coefficients, not replacement of finite modes. Every independent finite test retains all four W8 orders and their time/z derivatives.

Define pbar=-i Wbar-x a b, where

    a*b_T=a*partial_t Wbar/(2Wbar),
    a*b_L=a*partial_t Wbar/(2Wbar)+a*z*H.

These follow directly from the original physical momenta and their transverse/longitudinal dilution rates. In particular the mixed amplitude has a nonzero imaginary second-degree term; no scalarization or Maxwell substitution is made.

The five normalized amplitudes Abar,Bbar,TLbar,LTbar,LLbar are the exact original mass/electric/magnetic/constraint amplitudes multiplied by a. Their denominators are2sqrt(Wbar_1 Wbar_2); the TT numerator is m^2a^2x^2-pbar_1 pbar_2, Bbar numerator1, the TL/LT numerator includes ma x times its actual longitudinal frequency and momentum quotient, and LL retains Obar_1 Obar_2-m^2a^2x^2 pbar_1 pbar_2/(Obar_1 Obar_2).

For each of seven geometric contractions, keep detector Abar_sharp(t) fixed. Begin with source Abar(s)Gamma(s)/a(s) and iterate

    L_source F=partial_s[a(s)gbar(s)F],
    gbar=1/(Wbar_1+Wbar_2).

The normalized current coefficient is -Im[-i*i^j*gbar(t)*Abar_sharp(t)*(L_source^j F)(t)]. The detector is never differentiated by this operation. For0<=j<=4 every source derivative r<=j is retained:105 amplitude rows in total. Finite truncated-series arithmetic is exact because multiplication, reciprocal and positive square-root recurrences preserve all lower coefficients; differentiation in t does not lower inverse-radius order.

All odd endpoint labels remain. Their actual homogeneous contribution can be nonzero even when their spatial difference vanishes.
