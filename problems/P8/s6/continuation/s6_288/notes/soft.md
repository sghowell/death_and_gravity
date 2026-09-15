# Exact massive physical soft kernel and forward Coulomb phase

Use the complete S278 kernel
F(z)=V(z)M(z), V(z)=z^2-4mu*z+2mu^2,
M(z)=integral_0^1 dx/[4mu-z(1-x^2)],
B(s,t,u)=2[F(s)+F(t)+F(u)]-mu,
with its retained Feynman sheet.

For physical s>4mu put Q=s-4mu, beta=sqrt(1-4mu/s),
t=-tau, u=-Q+tau,0<=tau<=Q. The exact boundary values are

 M(s+i0)=[-atanh(beta)+i*pi/2]/(s*beta),
 M(-Q)=atanh(beta)/(s*beta), V(-Q)=V(s).

Hence B(s,0,-Q)=i*pi*V(s)/(s*beta), not zero.
Since the other two channels remain below threshold,
ImB=pi*V(s)/(s*beta) at EVERY physical angle.

## Positive attenuation density on the entire physical interval

Differentiate the parameter integrand twice:
F''(z)=4mu^2 integral_0^1
(1+6x^2+x^4)/[4mu-z(1-x^2)]^3 dx>0 for real z<4mu.

The real part is a convex finite difference. With A=4mu,w=1-x^2,
P=1+6x^2+x^4, its complete density is

 -ReB=4mu^2*tau*(Q-tau) integral_0^1
 P*(2A+wQ)/[A*(A+wQ)*(A+w*tau)*(A+w*(Q-tau))] dx.

Every factor is positive in the interior. There is no scan or
selected-angle assumption. Both endpoint values vanish, and
tau<->Q-tau leaves the result unchanged.

Equivalently the finite difference is twice the integral of F''
over the rectangle[0,tau]x[0,Q-tau], evaluated at minus the sum
of its coordinates. Since integral_0^1 P dx=16/5,

 128mu^2*tau*(Q-tau)/[5(4mu+Q)^3]
 <=-ReB<=2*tau*(Q-tau)/(5mu).

The forward attenuation slope is2[F'(0)-F'(-Q)]>0.
The exact zero jets are F(0)=mu/2,F'(0)=-11/12,F''(0)=1/(5mu).

## What soft division does and does not remove

Define eta(s)=V(s)/(8pi*kappa*s*beta).
The stipulated raw factor has
log W_E=B*(E/nu)^(2EP)/(8pi^2*kappa*EP).
At forward transfer it is i*eta*(E/nu)^(2EP)/EP.
For L=log(E1/E0), S278's conditional resolution ratio has modulus
exp[-ReB*L/(4pi^2*kappa)] and phase exp[-2i*eta*L].
For a coarser resolution L>0 the modulus increases, as expected
from that precise ratio convention. A single endpoint's charge
Ward identity does not cancel this four-leg Coulomb phase.

The two ordered boxes I4(t,s),I4(t,u), with their fixed coefficients,
have leading D4 term

 i*V(s)^2/[8pi*kappa^2*s*beta*t] *
 [-1/EP+log(4pi*nu^2/(-t))-Gamma_E].

Their logarithmic transfer coefficient, after the conditional soft
division, is A_tree*i*eta*log((-t)/E^2). It is not an ordinary
simple pole. The finite1/t constant is deliberately not inferred
from this D4 formula: evanescent box, tree and soft factors
multiplying IR poles affect it. The new rational pole comparison
subtracts the entire KNOWN scalar boxes and restores them in the
representative; it does not erase their cuts or choose a finite
physical phase convention.

The detector interpretation still assumes the S278 factorization
premise. This agrees with the need for regulated/smeared observables
emphasized in Bellazzini et al., Positivity with Long-Range
Interactions, https://arxiv.org/html/2512.13780v2 .
Their displayed gravitational eikonal example has massless
external states; the massive coefficients here are derived
directly above, not imported from that example. No finite-G error
bound, nonperturbative no-limit theorem or Regge closure follows.
