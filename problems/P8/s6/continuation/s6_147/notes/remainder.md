# Entire finite mass-ratio correction

The exact finite correction from the complete S6.144 outer reference
at b=1 scales to general b by r=b/T, T=4m^2. In CY/Q units,

    delta F(r)=-12 Integral_0^1 (2z-1)/(z sqrt(1-z))
       [h(rz)(-log(rz)-1)+rz A(rz)] dz,
    A(k)=(1-k)^(-2), h(k)=A(k)-1.

The difference is convergent before the regulator limit. This
scaling introduces no extra log b finite term from a pole:
the difference itself has no pole. On 0<r<=1/16, the frozen
pointwise and log-moment estimates give

    |delta F(b/T)| <= (b/T)[300+100 log(T/b)].

The physical domain m>=720 is contained in this domain. Integrating
the complete correction with the contraction and coupling factors
from forests.md yields

    |delta Delta f| <= NY^2/Q^2 * (175+50 log T)/T.

Here Integral_0^1 b db=1/2 and Integral_0^1(-b log b)db=1/4.
There is no logarithmic endpoint omission at b=0.

For exact numerical certification use the least dyadic exponents
n_m,n_T with m^2<=2^n_m and T<=2^n_T. Since log 2<1,
log(m^2)<=n_m and log T<=n_T. Strict rational inputs are
required, with m>=720,Y>=0 and 0<Q_lower<=144.
The bound uses N=6 and a lower bound Q_lower on physical 16pi^2.
At the actual candidate n_m=1329 and n_T=1331.

The complete scalar-mass difference is approximately bounded by
5.11704625009938e-410. The finite-ratio portion alone is below
2.00344382855632e-810. These are rigorous rational upper bounds
printed approximately, not floating-point fit coefficients.
