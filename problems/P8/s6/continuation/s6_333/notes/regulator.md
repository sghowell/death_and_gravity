# Energy-dependent explicit regulator error

Keep exactly S301/S296/S325's fixed-ball and phase prescription:
c(e)=e/[2(1+e)], A(e)=Gamma(3/2+e)/[Gamma(3/2)Gamma(1+e)],
p(e)=(4pi)^(-e)Gamma(3/2)/Gamma(3/2+e).
All differences below subtract the SAME Born state, andR<=1/8.
Write K0,U0 for the boundary differences and
I_H(e),I_U(e) for radial differences with the extra tau^e weight.
Then

q_e-Delta =1/(8pi^2*kappa) times
 {[(p-1)/e-p'(0)]*K0
  +[p*c/e-1/2]*U0
  +(p*A-1)*I_H(e)+(I_H(e)-I_H(0))
  +p*A*c*I_U(e)}.

This is an exact algebraic decomposition, not an assignment of hard
finite matching. Frozen bounds give
|K0|<=50000R, |U0|<=30000R,
|I_H(e)|<=400000R*L_R, |I_U(e)|<=240000R*L_R,
p<=1, |p-1|<=4e, |p-1-p'(0)e|<=17e^2/2,
A<=2, |A-1|<=2e, c<=e/2.
In particular |p*A-1|<=|p|*|A-1|+|p-1|<=6e and p*A<=2.

The only new radial majorant is
J1=integral_0^1 min(A0*R,B0*tau^(1/4))*|ln tau|/tau dtau.
Set z=ln(B0/(A0*R))>0 and split at tau0=(A0*R/B0)^4.
Exact integration gives
J1=A0*R*(8z^2+16z+16)
   =8A0*R*[(1+z)^2+1].
For A0=100000,B0=64000 and R<=1/8, z<=ln(1/R), so
J1<=1600000R*L_R^2.
Since |1-tau^e|<=e*|ln tau|,
|I_H(e)-I_H(0)|<=1600000e*R*L_R^2.

Also |p*c/e-1/2|=|p/(2(1+e))-1/2|<=5e/2.
The five numerator budgets are respectively
425000R,75000R,2400000R*L_R,1600000R*L_R^2,240000R*L_R,
each multiplied by e. As L_R>=1 their sum is <=4740000e*R*L_R^2.
Division by8pi^2>72 gives
|q_e-Delta|<=66000e*R*L_R^2/kappa.
The zero-radiation and e=0 cases have exact zero differences.
S331's positive Borel extension preserves the angular and radial
majorants, so this estimate extends to the limiting cloud on its cut.


## Conditional regulator error

For y=x-R, l=|ln y|,
|g_e(y)-1|<=2e*l and |h_e(y)-ln y|<=e*l^2.
The second inequality is exp(-t)-1+t<=t^2/2 for t>=0.
Therefore, using the actual Born-specific bounds,

|Z_e-Z_0|
 <=e/kappa*[66000R*L_R^2+20000R*L_R*l+1400R*l^2].

Conditional expectation is bounded by
(66000*5+20000*3+1400*5)*e*a*x*L_x^2/kappa
=397000e*a*x*L_x^2/kappa
<=400000e*a*x*L_x^2/kappa.


At epsilon=0 or zero radiation use the exact zero remainder. The strict
arithmetic margins for positive parameters also imply these non-strict
bounds throughout the closed regulator domain.
