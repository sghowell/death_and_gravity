# Exact zero-soft kernel and finite outer momentum differences

A constant Phi shifts the two active fermion masses by +/-y Phi.
Two background derivatives of the two-point determinant therefore give

    Gamma0_MS(q^2)=Y d^2 f_MS(-q^2)/dm^2,

where mu is held fixed during the mass derivatives, and only then set
to m. The raw identity includes all six boxes and their same proper
quartic counterterm. Its value at q=0 is -64 N Y^2/Q, not an independently
adjustable quartic matching coefficient.

Let C=2NY/Q, T=4m^2. Differentiating the known two-point density gives
the signed density nu(v)=12CY(2T/v-1)/sqrt(1-T/v). The density and its
first mass derivative vanish at threshold before the second derivative;
there are no delta-function boundary terms. The once-subtracted identity is

    Gamma0(q^2)=Gamma0(0)+q^2 integral_T^infinity nu(v)/[v(v+q^2)] dv.

This density is not positive and is not a normalized exact KL measure.
Since |nu|<=12CY/sqrt(1-T/v), put z=T/v and t=q^2/T.
On z<=1/2 bound 1/sqrt(1-z)<2. On z>=1/2 bound t/(1+tz)<=2.
The latter integral is below three. This proves

    |Gamma0(q^2)| <= NY^2/Q [136+48 log(1+q^2)], T>=1.

For the outer bubble subtract the entire zero-momentum F_D reference
before taking epsilon to zero, as defined in the next proof. At four
dimensions the difference of the physical light product and 1/S^2 has
numerator z(S/2-q_axis^2)-z^2/16 and denominator P S^2.
Here |S/2-q_axis^2|<=S/2. The previous gap proves its norm <=33/S^3.
The radial moments integral x/S^3 dx=1/2 and integral x log(S)/S^3 dx=3/4
give |B_Gamma(z)-B_Gamma(0)|<=3432 NY^2/Q^2.

Each internal-heavy term is absolutely convergent. The light and heavy
bounds give 32 NY^2/Q^2 times
136 J0(M)+48 Jlog(M), where J0 and Jlog integrate
x/[(1+x)^2(M+x)] and the same times log(1+x).
Split at M>=1. Below M use 1/(M+x)<=1/M and x/(1+x)^2<=1/(1+x).
Above M use the integrand <=1/x^2 and log(1+x)<=log(2x).
With ell_M=log(2M),

    J0 <=(ell_M+1)/M,
    Jlog <=(ell_M^2/2+ell_M+1)/M.

Thus each triangle is bounded by
32 NY^2/(Q^2 M) [24 ell_M^2+184 ell_M+184].
There are two such terms per channel, weighted by g. The bubble is
weighted by |C|<=L+g/(M-3). All three channels are included.
