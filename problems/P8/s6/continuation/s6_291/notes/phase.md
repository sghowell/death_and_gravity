# Entire D phase, independent B0 normalization and evanescence

Use D=4+2e and the frozen loop normalization nu^(-2e).
The on-shell radial delta integral gives
Phi2_D=(2pi)^(2-D)Omega_(D-2)p^(D-3)/(4sqrt(s)nu^(2e)),
Omega=2pi^((D-1)/2)/Gamma((D-1)/2),p=(s-n)/(2sqrt(s)).
Gamma duplication yields

Phi2_D=beta_Hh/(8pi)*(4pi nu^2)^(-e)s^e beta_Hh^(2e)*
Gamma(1+e)/Gamma(2+2e).

Independently, frozen B0_D=Gamma(-e)(4pi nu^2)^(-e)int[Delta-i0]^e dx.
For s>n its negative interval is0<x<beta_Hh with
Delta=s*x*(x-beta_Hh). The substitution x=beta_Hh*y gives the whole
Beta integral s^e beta_Hh^(1+2e)Gamma(1+e)^2/Gamma(2+2e).
Reflection -Gamma(-e)sin(pi e)=pi/Gamma(1+e) fixes
ImB0_D=8pi^2 Phi2_D. Hence optical half Phi2_D=ImB0_D/(16pi^2).
The independent derivations fix the phase, sign and factor2; no D0 guess.

The polar angle measure in D spatial kinematics has weight(1-x^2)^e.
Let I_e=sqrt(pi)Gamma(1+e)/[2Gamma(3/2+e)] and
N_e(B)=int0^1(1-x^2)^(e+2)/(1-Bx^2)^2 dx.
The normalized average isF_e=N_e/I_e. Euler's beta integral gives
F_e=4(e+1)(e+2)/[(2e+3)(2e+5)]*
2F1(2,1/2;e+7/2;B).
Pointwise comparison proves its lower bound atB0 and upper bound1.

Combining the full projector and all phase/angle factors gives EXACTLY
rho_e(s)=(s-n)^(-1+2e)K_e(s), with
K_e(s)=g^2(s-4mu)^2/(4pi kappa s)*
(1+2e)/(2+2e)*(16pi nu^2 s)^(-e)/Gamma(1+e)*N_e(1-4mu/s).
The phase divided byI_e is1/[(16pi nu^2 s)^e Gamma(1+e)], without
an extra factor2. At e0 this reproduces the full D4 cut.

Write L_B=int0^1[(1-x^2)^2/(1-Bx^2)^2]ln(1-x^2)dx.
Differentiating every factor BEFORE discarding e terms gives
K1(s)=g^2(s-4mu)^2/(8pi kappa s)*
[(EulerGamma+1-ln(16pi nu^2 s))*F0(B)+L_B].
The dimension-dependent TT factor, gamma phase and angular log all matter.
For fixed massive compact domains, the denominator has a positive lower
bound and logarithmic angular moments dominate the derivatives. This
justifies differentiation under the integral. Numerical full phase/B0,
hypergeometric and derivative tests are independent checks only.
