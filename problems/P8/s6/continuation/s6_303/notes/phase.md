# Full finite Coulomb phase, classical term and real quantum logarithm

Putq=s-4>0,D=sqrt(sq),V=s^2-4s+2 and
h=atanh(sqrt(q/s)). The physical massive master values are

J_s=(-4h+2i*pi)/D, J_-q=4h/D.

Thus J_s+J_-q=2i*pi/D and the complete forward soft kernel,
including the diagonal, isB_soft=i*pi*V/D.

## Leading imaginary part

The two boxes whose first channel ist give
[V_s^2 J_s+V_u^2 J_u]*(ell-ln tau)/t
-4[V_s J_s+V_u J_u]/t.
At t=0 their endpoint is
2i*pi*V^2*(ell-ln tau)/(D*t)-8i*pi*V/(D*t).

The finite D-tree soft division contributes4B_soft/t at its endpoint,
or+4i*pi*V/(D*t). The combined known principal phase is therefore

2i*pi*V^2/(D*t)[ell-ln tau-2/V].

The term-2/V includes both the dimensional box derivative and finite
D-tree division. Dropping either changes this finite reference.
All nonbox integer poles cancel by the exact regularized block in
notes/regularized.md. Against the leading Born V/(kappa*tau), the
known amplitude phase is

i*eta[ln tau-ell+2/V], eta=V/(8pi*kappa*D).

It is imaginary on the physical edge, not absent from the amplitude.

## Real small-channel coefficients

The complete S284 row gives at t=0

c00=-3(5V+6)/2, b00=(87V+118)/60.

For negative transfer the massless triangle is
C(-tau)=-pi^2/[sqrt(tau)*sqrt(tau+4)]
       -(1/2)int_0^1 ln[tau*x(1-x)]/[1+tau*x(1-x)]dx.

For the first term, x=(1-cos theta)/2 followed bytan theta on
the half interval gives exactly
int dx/[sqrt(x(1-x))(1+tau*x(1-x))]
=pi/sqrt(1+tau/4). This fixes the entire classical radical.

For0<tau<=1,

C(-tau)=-pi^2/(2sqrt(tau))-ln(tau)/2+1+R_C,
|R_C|<5sqrt(tau)/8+tau*|ln tau|/12+5tau/36.

Use1-(1+x)^(-1/2)<=x/2,pi^2<10,
int x(1-x)dx=1/6 and
int-x(1-x)ln[x(1-x)]dx=5/18.
The logarithmic remainder uses1-1/(1+tau*x(1-x))<=tau*x(1-x).
All endpoint terms are integrable.

Forf(a)=V_a^2 J_a on negativea,

f'(-q)=4V(s-2)(8-3V)h/D^3-2V^2/D^2.

The real endpoint of f(s)+f(-q) vanishes. Taylor expansion of the
small-channel boxes consequently suppliesf'(-q)*ln tau.
Combining it with-c00*ln tau/2 and-b00*ln tau gives

C_log(s)=4V(s-2)(8-3V)h/D^3-2V^2/D^2+(69V+76)/30.

The classical coefficient is-c00*pi^2/2=3pi^2(5V+6)/4.

## Full complex bounded remainder and analytic boundary

The remaining large-channel master values are smooth boundary
functions on the compact s interval, and the crossed channel stays
strictly spacelike near the endpoint. Their coefficients are regular.
The massive small-channel block extends continuously after exact
integer-pole cancellation. Small-channel coefficient differences
multiply only O(tau^-1/2) or O(ln tau); their residuals are bounded.
The box Taylor remainder divided byt isO(tau ln tau). The soft
divisor remainder after its endpoint is bounded as well.
These facts give the stated complexO(1), uniformly on the compact
energy interval, not merely an expansion at sampled energies.

Taking ln tau around its branch changes the principal amplitude by
-4pi^2 V^2/(D*tau). It is nonzero. Unknown local polynomial and Newton
simple-pole coordinates do not cancel this logarithmic monodromy.
A real interference bound therefore creates no holomorphic transfer
disk and supplies no complex Regge contour estimate.
