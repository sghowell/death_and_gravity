# Full massless-pair coefficient and angular-energy measure

Let pi be the four future massive momenta and eta=(-1,-1,+1,+1).
For arbitrary positive outgoing null rays qr=omega_r*(1,nr), define
a_i(n)=pi^0-vec pi.n and b(n,n')=1-n.n'. The massive pair function is

F_mu(d)=(d^2-mu^2/2)*acosh(d/mu)/sqrt(d^2-mu^2),
F_mu(mu)=mu/2.

Its equivalent Feynman-parameter integral is
(d^2-mu^2/2)*integral_0^1 dx/[mu+2x(1-x)(d-mu)].
The constant equal-velocity value is used directly, so no removable
0/0 is evaluated at forward/backward zero-radiation endpoints.

With sigma=sum_r omega_r*delta(nr), the COMPLETE finite coefficient is

K_sigma = 2mu+2sum_(i<j) eta_i*eta_j*F_mu(pi.pj)
 +2sum_i eta_i*integral a_i(n)*ln[2a_i(n)/sqrt(mu)] dsigma(n)
 +integral integral b(n,n')*ln[2b(n,n')] dsigma(n)dsigma(n').

Set0ln0=0. The double integral counts each unordered radiation pair
twice; self-pairs and exactly collinear pairs vanish. Dropping it is
incorrect when more than one noncollinear graviton is present.

To derive this expression, temporarily give each radiation leg a
positive auxiliary mass m_r=epsilon*omega_r. The massive pair asymptotic
is F_(mi,mr)(d)=d*ln[2d/(mi*mr)]+o(1).
One physical regulating family keeps each energy fixed, replaces its
spatial velocity by sqrt(1-epsilon^2)*n_r, and adjusts the massive pair
with the smooth recoil map. Dot products change by O(epsilon^2), so
the extra logarithmic terms vanish. Exactly collinear radiation pairs
have the common-velocity value epsilon^2*omega_r*omega_s/2 and vanish;
the resulting finite formula extends continuously to those pairs.
The coefficient of each ln m_r is proportional to

qr.(sum_i eta_i*pi+sum_s qs)=0.

This cancels every collinear mass logarithm before the masses vanish.
After cancellation, the massive/null terms give ln(2a_i/sqrt(mu)) and
the null/null terms give ln(2b). For one null leg, its common ln2 term
also vanishes because sum_i eta_i*pi.q=0. For many legs the
radiation/radiation contribution is essential to the cancellation.

A boost changes the reference energies omega_r in the logarithms.
The coefficient of each ln(omega'_r/omega_r) is the same conserved zero.
The displayed COM formula is therefore Lorentz invariant. Arbitrary
collinear splitting of one ray preserves sigma and thus K exactly.

The full current is J=sum_A P_A P_A/(P_A.l), with all momenta outgoing
and l=(1,n) a null test direction. Conservation gives l.J=0.
Its physical TT norm is nonnegative and integrates to K. The recoiled
four-massive-leg subset alone has l.J=-Q and is not conserved.

The angular limit needed for the dimensional pole is controlled, but
not its first dimensional derivative. Each massive leg contributes
transverse matrix norm at most Ei+|vec pi|<=4. Each massless leg
contributes at most2omega_r, including its collinear limit. Hence the
full transverse matrix norm and absolute trace are at most16+2R.
The D=4+2e projector therefore has absolute contraction bounded by
3*(16+2R)^2/2 for e>=0. The normalized angular measure tends to the
ordinary D4 sphere measure. For finite N the only potentially
discontinuous directions are the finitely many massless rays, a
zero-area set. Dominated convergence fixes the pole coefficient K;
it does NOT give the finite angular conversion Delta_sigma.

Finally, the mixed angular kernels are continuous because
a_i>=1/4, and b ln(2b) is continuous on the compact double sphere
after its zero extension. Finite positive angular energy measures
can be approximated weakly by atomic ones. The recoil map is smooth
in their total four-momentum, so K has a unique continuous extension.
Positivity passes to that limit. This is a statement about kinematic
coefficients, not the existence of a full quantum state or amplitude.
