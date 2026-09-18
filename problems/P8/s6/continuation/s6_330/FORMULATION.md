# Positive increments and normalized infrared cutoff

Retain the original action and fixed parameters kappa=10^800,
n=10^200/512+2,g=1/8192,mass1,nu=mu=1.
Fix E in[5/4,2] and a unit rest-frame outgoing direction u.
For a positive Borel angular-energy measure sigma of mass<=1/8,
C_sigma is S328's complete continuous logarithmic TT coefficient,
with its original recoil and both complex components.
Its physical normalization is C_sigma/kappa^(3/2).

If sigma=tau+rho with tau,rho positive and t=massrho, then
||C_sigma-C_tau||sup<=40000t. This applies to a positive added
tail at fixed E,u, not an arbitrary signed change or changing hard data.

Use the fixed-Born Poisson intensity Lambda(dn)*dw/w on0<w<=1,
a=int Lambda,0<=a<=1, including the original a<4/(5kappa).
Set R=masssigma, sigma_eta=sigma restricted to w>=eta,
R_eta=masssigma_eta and T_eta=R-R_eta. For0<eta<=x<=1/8,
P(R<=x)=exp(-EulerGamma*a)*x^a/Gamma(1+a) is kept unexpanded.
The exact conditioned missing-energy mean is
a*x/(a+1)*[1-(1-eta/x)^(a+1)]<=a*eta.
Its second moment keeps both the one-emission diagonal and ordered
two-emission contribution; see notes/cutoff.md.

Writing P_eta(x)=P(R_eta<=x), the relative conditioning loss obeys
0<=1-P(R<=x)/P_eta(x)<=a*eta/x.
Let B=23667. The actual finite-cutoff and limiting conditioned marks obey
||E[C_sigma_eta|R_eta<=x]-E[C_sigma|R<=x]||sup
 <=(40000+2B)*a*eta=87334*a*eta.
The difference of E[||C-C_Born||sup^2] under their respective
conditioning events is at most
(2B*40000+B^2)*a*x*eta=2453486889*a*x*eta.
These statements also bound each angular/polarization second moment.
At original a<4/(5kappa), physical errors are respectively
<70000*eta/kappa^(5/2) and<2000000000*x*eta/kappa^4.

At a=0 all emissions and changes vanish. At t=0 the coefficient
difference is exactly zero. The result does not discard a phase,
expand x^a, assume conditioned Poisson independence, or identify
the leading reference with the interacting physical state.
