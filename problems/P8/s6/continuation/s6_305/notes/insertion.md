# Conserved sources and the full inverse-kernel insertion

Use eta=(+---), g=eta+2h/sqrt(kappa) and the original linear de Donder
propagator. For an all-incoming scalar pair p,q of mass squared mu, let

T^{mn}=p^m q^n+q^m p^n-eta^{mn}(p.q+mu).

This tensor is conserved against p+q. With channel invariants a+b+c=4mu,
the two conserved stress contractions are

N0=(a+2mu)^2/3,
N=2mu^2-2mu*a-b*c,
N2=N+N0/2=[3(b-c)^2-(a-4mu)^2]/12.

Their Einstein contraction is N2-N0/2=N, and the whole canonical
gravitational Born amplitude is AG=-sum_channels N/(kappa*a).
Independent exact component contractions check all three partitions
for multiple rational massive shells, including masses different from1.
The symbolic result agrees with the independently rebuilt S297 Born form.

Write p=-a on the Feynman sheet. S285's full fixed H/Proca kernels are

O2=p(1+dk/kappa)+p^2 A(p)/(16pi^2 kappa),
O0=-2p(1+dk/kappa)+p^2 H(p)/(192pi^2 kappa).

Insert a formal loop marker before expanding the inverse. Its derivatives
at zero are

d(1/O2)=-dk/(kappa*p)-A/(16pi^2 kappa),
d(1/O0)=dk/(2kappa*p)-H/(768pi^2 kappa).

Contract with (N2,N0)/kappa and sum all three channels. The COMPLETE
fixed spectator correction is consequently

MHV=-dk*AG/kappa
    -sum[N2*A(-a)/(16pi^2 kappa^2)+N0*H(-a)/(768pi^2 kappa^2)].

In particular the trace inverse sign is negative in its finite term,
and the fixed Newton insertion has the same sign in both physical
source sectors. The spin0 source sector is not an extra massless scalar
graviton. The original volume cancellation has already removed the
cosmological double pole before this inverse is formed.

For L=log(n)+2, the complete radial functions are

A(p)=L/60+sum_X integral_0^1 p*w_X2(v)/(4nu_X+p(1-v^2)) dv,
H(p)=2L+sum_X integral_0^1 p*w_X0(v)/(4nu_X+p(1-v^2)) dv,

with X=H,Proca and their actual masses. The weights are

w_H2=v^6/30, w_H0=v^2(3-v^2)^2,
w_V2=v^2(30-20v^2+3v^4)/30,
w_V0=v^2(3-2v^2+3v^4).

Every crossed integral and fixed local constant is retained. This is a
full convergent nonlocal Gaussian insertion, not a finite low-momentum
Taylor truncation. Its use here is still first loop in the scattering
amplitude, not an exact interacting gravitational theory.
