# Regulator-uniform infrared cutoff and joint error

Use the SAME full Poisson space, A,A_eta,p,p_eta as S332, but with
Z_e in both densities. Its generic nested-event L1 split remains valid.

The q_e increment plus its two normalization terms give
(2*18000+2*15000)*a*eta*L_eta/kappa
=66000*a*eta*L_eta/kappa.

The delta a*h_e part is dominated term by term by S332's logarithmic
argument: same-event14600*a*eta*L_eta/kappa and event-change
6300*a*eta*L_eta^2/kappa.

Only the q_e(tau)[g_e(y_eta)-g_e(y)] term is new. For eta<=x/4,
monotonicity R_eta*L_Reta<=x*L_x and S332's small-tail logratio bound give
<=2e*15000*x*L_x*a*eta/x*[4+2a*(2+ln(x/eta))]/kappa
<=240000e*a*eta*L_eta^2/kappa
<=30000*a*eta*L_eta^2/kappa.

For eta>x/4, use R_eta*L_Reta<=R*L_R and logratio<=|ln(x-R)|.
The needed mixed conditional moment is
E[R*L_R*|ln(x-R)||A]
 =a*x/(a+1)*[L_x*(ln(1/x)+H_(a+1))
             +ln(1/x)/(a+1)+H_(a+1)/(a+1)-psi1(a+2)].
This follows from two derivatives of the beta integral.
For0<a<=1, H_(a+1)<=3/2,1/(a+1)<=1,psi1(a+2)>0,
so it is <=3*a*x*L_x^2.
Consequently the extra exponent term is
<=2e*15000*3*a*x*L_x^2/kappa
<=45000*a*eta*L_eta^2/kappa, using2e<=1/4 and x<4eta.
At a=0 it vanishes. The mixed moment bound is integrable at both endpoints.

Combining66000+45000+14600+6300=131900<132000 gives
sup_(0<=e<=1/8) E|f_e,eta-f_e|
 <=132000*a*eta*L_eta^2/kappa
 <106000*eta*L_eta^2/kappa^2 at original a<4/(5kappa),
since132000*4/5=105600.


## Joint error and ordering of limits

Combine this with the earlier uniform cutoff bound:
||f_e,eta-f_0||_1
 <=[132000*a*eta*L_eta^2+400000e*a*x*L_x^2]/kappa
 <[106000eta*L_eta^2+320000e*x*L_x^2]/kappa^2
at original a<4/(5kappa).

The reference probability space, exact cuts, changing remaining energy,
Born subtraction and dimensional ADDITIONAL-soft convention are fixed.
This is an explicit two-parameter common-space signed-density bound.
It gives arbitrary joint and both iterated limits for fixed x, with
no ordering assumption between eta and e. It is not TV convergence
of configuration laws or a continuation of the outer hard amplitudes.
Hard/evanescent matching, interacting detector probabilities, quantum
state and all-N summability, absolute Regge, bounce and UV remain open.
