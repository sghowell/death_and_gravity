# Channel Cauchy lift and required pole running

Fix t in[0,4mu), and use the indicated original normal density
rho(sigma,t) on positive sigma. For t=0, sigma*rho is O(sigma^(-1/2));
for fixed0<t<4mu it is O(1)+O(sqrt(sigma)). It is therefore locally
absolutely integrable in either case. No estimate uniform in t down
to zero is asserted.

For any finite L>0 and z outside[0,L], define

    F_L(z;t)=c_L(t)/z+(1/pi) integral_0^L
                rho(sigma,t)[1/(sigma-z)+1/z]dsigma.

The identity in brackets equals sigma/[z(sigma-z)] proves absolute
convergence and uniform compact domination. Thus zF_L is an ordinary
Cauchy transform plus c_L; F_L is holomorphic off the integration cut.
The continuous interior normal density, including its removable external
threshold, gives jump2i*rho by the Cauchy boundary formula. This is a
well-defined CHANNEL construction. It is not claimed to be the full
crossed amplitude or a physical choice of finite counterterm.

The cut does not fix c_L(t). Multiplication by z has not removed the
massless pole from the amplitude: it has exposed its matching datum.
For a change in the upper cap to add only the ordinary shell,

    d_L F_L(z;t)=rho(L,t)/[pi(L-z)],

direct differentiation proves the necessary and sufficient condition

    d_L c_L(t)=-rho(L,t)/pi.

The initial c_L0(t) remains undetermined. Keeping it fixed adds the
spurious shell rho(L,t)/(pi*z). This is a theorem about matching the
subtraction convention, not an arbitrary choice of the original anchor.

Even at fixed forward transfer, a constant crossing-even anchor error is
visible: delta_c(1/s+1/u), s,u=2mu+/-v, has v^2 coefficient
delta_c/(4mu^3). Likewise half the second derivative of the compensated
kernel differs from1/(sigma-z)^3 by1/z^3. Thus using a compensated low
integral and then ignoring its pole term does not produce the correct
positivity coefficient.

Full crossed/double-spectral assembly can contain the same graph in
different channel cuts. Adding three independently reconstructed channel
functions need not avoid double counting. The transfer contribution,
original finite pole normalization and local polynomial must all be
matched before claiming the full subtracted b20. None is set to zero here.
