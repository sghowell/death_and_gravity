# The other five modes and the unchanged vector preparation

## Two tensors

Use each real TT polarization E with E:E=2 and dual
E/2. Its reduced action is a^3*(gamma'^2-q*gamma^2)/4.
The canonical oscillator y=a^(3/2)*gamma/sqrt(2) obeys

y''+(q-pump)y=0,
pump=3H'/2+9H^2/4=(6+30u^2)/(1+u^2)^2.

On I, 0<pump<=15 and |pump'|<=65. The latter follows
from |H|<=2, |H'|<=5 and |H''|<=13. The value bound
uses H'<=4 and H^2<=4. For q>=10^4, Omega^2=q-pump>=q/2,
and |(Omega^2)'|<=4q+65<=5q. Thus the positive oscillator
energy E=(|y'|^2+Omega^2|y|^2)/2 satisfies |E'|<=10E.

Choose y_i=1/sqrt(2*kappa), y'_i=-i*sqrt(kappa/2),
where kappa=sqrt(q_i). Then E_i<=kappa. Even using the
whole interval I, E<3^10*kappa<10^5*kappa. The global
scale ratio gives q>=kappa^2/4, hence the generous bounds

|y|<2000/sqrt(kappa), |y'|<500*sqrt(kappa).

The actual local canonical phase components are

gamma=sqrt(2)*a^(-3/2)*y,
p_gamma/a^3=a^(-3/2)*(y'-3H*y/2)/sqrt(2).

Keep the H*y term, the TT tensor/dual factors and the
fixed-center Fourier change. Every component is then
less than 10^5*U0, far below S. Both polarizations are
kept. This does not establish tensor loop control or
an all-momentum tensor Hadamard-state theorem.

## Three Proca modes

Keep the ORIGINAL selected Borel-prepared state, fixed
dimensionless m=1000 and zeta=10^-6. The parent replay
retains the S6.59 preparation and S6.65/S6.68 evolution
bounds relative to the eighth-order reference:

v=A_mix*f+B_mix*f*, p_v=A_mix*p_f+B_mix*p_f*,
|B_mix|<=C_mix/nu^6<1, |A_mix|<2,
C_mix=1813229+5347035781757616/m^4,
nu^2=m^2+|k_com|^2/(25/16)^2.

The initial B_mix is not set to zero. The reference is
not reselected as the exact state. The frozen reference
bounds give |f|<=omega^(-1/2), |p_f|<=2*omega^(1/2),
where omega^2=q+m^2 and p_f=f'-(g'/g)f. Therefore

|v|<=3/sqrt(omega), |p_v|<=6*sqrt(omega).

For either transverse polarization g_T^2=a*zeta.
For the longitudinal scalar g_L^2=a^3*q/(q+m^2).
The local physical coordinate and momentum density are

W_T=a^(-3/2)*m*v,
Pi_T=a^(-3/2)*p_v/m,
W_L=a^(-3/2)*omega*v,
Pi_L=a^(-3/2)*p_v/omega.

These follow from W_local=W_com/a, Pi_local=Pi_com/a^2
and the original longitudinal W_com=grad(sigma)
decomposition. Their product of weights is a^-3 in both
sectors, as required by the local phase density. The
symbolic checks reproduce both canonical weights, not
only transverse Maxwell normalization.

Since 1<=omega<=5U0 and a>=1, all phase columns, after
the fixed-center Fourier factor below two and unit
polarization components, are less than 10^5*U0<S.
The common seed thus covers two scalars, two tensors
and all three Proca polarizations.

The band is above m, not a light-only effective theory
below its vector threshold. No heavy-mass elimination
expansion is applied at these momenta.
