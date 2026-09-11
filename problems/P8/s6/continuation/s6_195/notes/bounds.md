# Finite-band bounds uniform in external momentum

Everything here is evaluated at the original isotropic clock.
In the actual balanced variables Sigma_k=T_k^-1 C_k T_k^-t,
the unchanged initial state has trace at most108. The free
generator has a skew fast-frequency part and symmetric part
of norm at most11 on the unit slab, from the S190/S192
actual mode estimates. Hence ||U_k(t,s)||<=exp(11|t-s|)
uniformly in k.

An elementary exponential-series bound proves e<68/25 and
(68/25)^22<4e9. Thus the product of the two propagator norms
is below4e9, while ||Sigma_k||<108*4e9<4e12.
These are value bounds, not second-amplitude derivative bounds.

Set nu_k=sqrt(m^2+|k|^2/a(t0)^2). On this slab
nu_k<=omega_k(t)<2nu_k. The balanced forcing vertex has norm
at most sqrt(omega_k omega_q)||Gamma_(k-q)||. Both source
terms and both propagators therefore give

    ||delta Sigma_kq(t)||
      < [2*4e12*4e9*2] sqrt(nu_k nu_q)
            integral_(t0)^t ||Gamma_(k-q)(s)|| ds
      <1e23 sqrt(nu_k nu_q) integral ||Gamma_(k-q)||.

For a unit operator-norm detector, the six-dimensional trace
inequality and the factor1/2 in the ordered current give
3*2*1e23 nu_k nu_q<1e24 nu_k nu_q times the source integral.
The local contact is bounded by
3*2*4e12 nu_k<3e13 nu_k times the two direction norms.
The constants retain both frequency and canonical factors.

Now restrict the field phase space with the common Galerkin
projection P_K, |k|<=K, and differentiate P_K H_gamma P_K
BEFORE any ultraviolet limit. The propagated term at external
p integrates over |k|<=K,|q=k-p|<=K. For K>=m=1000,
nu_k,nu_q<2K. Its overlap measure is no larger than
K^3/(6pi^2)<K^3/54, for every p. This proves

    memory <1e24 K^5 integral ||Gamma_p(s)|| ds,
    contact <1e13 K^4 ||Gamma_p(t)||.

The contact integral has only |k|<=K. It is not restricted by
the two-mode propagation overlap, since D Gamma carries zero
total momentum for the output/input Fourier pair.

For K=1e16 and slab length1 the complete projected current
is below2e104 sup_I||Gamma_p||. Dividing by a^3>=1 and
retaining both canonical factors2/sqrt(kappa) gives
8e-696 sup_I||Gamma_p|| for the canonical linear response.
For smooth spatially compact inputs one may integrate this
uniform multiplier estimate in p using the norm
integral d^3p/(2pi)^3 sup_I||Gamma_p||. Finite Fourier
fixtures are supplementary, not the proof of this bound.

These estimates deliberately do not integrate the growing
internal-momentum majorant over an infinite domain.
