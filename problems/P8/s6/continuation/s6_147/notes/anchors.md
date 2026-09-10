# Massless-exchange anchors before taking d=4

Write e=epsilon, d=4-2e, Q=16 pi^2 and
A(e)=exp(gamma_E e) Gamma(1+e). The massive/massive/massless
Gaussian master and trace conventions are exactly those of S6.146.

In units N Y m^4/Q^2 for the scalar chord and N a C_F m^4/Q^2
for the gauge chord, the raw vacuum rational factors multiplying
Gamma(e)^2 exp(2 gamma_E e)(mu/m)^(4e) are

    V_s = 4/[(1-e)(2e-1)] + 1/(1-e)^2,
    V_g = -4/[(1-e)(2e-1)] + (2-2e)/(1-e)^2.

These follow directly from the two-propagator traces and the
I(1,1,1) and factorized tadpole masters, including the one-half
contraction weight and the two gauge vertex phases. Differentiating
twice with respect to the background fermion mass, at fixed mu,
multiplies each rational factor by (4-4e)(3-4e). Multiply by Y
for the two external scalar insertions.

Independently, the complete zero-momentum quadratic scalar tensor
is 2 tr(S_k^3 S_l)+tr(S_k^2 S_l^2). The gauge tensor is
-[2 tr(S_k gamma_mu S_l gamma_mu S_k^2)
  +tr(gamma_mu S_k^2 gamma_mu S_l^2)].
The equal self-energy placements are combined only using symmetry
of the integrated k,l measure. Exact dimension-symbolic Gaussian
reduction agrees with both fixed-mu vacuum derivatives.

The proper mass/kinetic/Yukawa counterterms in forests.md contribute
c_s=2(6-3e)[4-2/(e-1)] and c_g=-4(6-3e)[4-2/(e-1)].
The complete normalized forests are

    [A(e)^2 R_s(e)+A(e)c_s(e)]/e^2
        =36/e^2-48/e-56+O(e),
    [A(e)^2 R_g(e)+A(e)c_g(e)]/e^2
        =-72/e^2+24/e+40+O(e).

Here R=(4-4e)(3-4e)V. In both finite terms the pi^2 parts
cancel between the raw double pole and the proper-counterterm
pole products. Set mu=m only after the derivative and then
subtract the remaining overall MS poles, not the finite terms.

Thus f_MS(0)=N m^2[-56Y^2+40Ya C_F]/Q^2 for massless exchange.
The gauge chord really is massless. The scalar chord is not;
the next notes restore its physical mass. Along Y=19a/45 the
combined coefficient is N Y a m^2(1336/45)/Q^2. Its sign is
not a sign claim for the entire matched two-loop amplitude.
