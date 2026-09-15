# S6.282: full original massive graviton normal-channel continuation

Use the unchanged S281/S280 formal vacuum source, mass squared mu=1,
kappa=10^800, all three vacuum constants, heavy source, physical frame
and literal Einstein four-helicity normalization. This does not identify
the formal vacuum parameters with exact interacting LSZ poles.

Put Q=s-4mu, u=4mu-s-t, and define the analytic kernel

    M(z)=integral_0^1 dx/[4mu-z(1-x^2)], z outside[4mu,infinity).
    V(z)=z^2-4mu*z+2mu^2,
    A(s)=s^3-8mu*s^2+20mu^2*s-12mu^3,
    B(s)=s^3-10mu*s^2+30mu^2*s-20mu^3,
    R(s)=s^2-18mu*s+146mu^2,
    D(s)=13s^2-24mu*s-232mu^2.

The COMPLETE original two-TT-graviton normal cut is

    rho_gg(s,t)=1/(4pi kappa^2) {
      [V(t)^2 M(t)+V(u)^2 M(u)]/s
      -A(s)M(4mu-s)+R(s)/120
      +3tu [B(s)M(4mu-s)+D(s)/120]/Q^2 }.

The last quotient is regular at Q=0, with value29/40. A uniformly regular
parameter integral gives the extension, not numerical cancellation.
The formula is holomorphic for s outside(-infinity,0] and t,u outside
[4mu,infinity), on the component continuing the physical normal cut.
In particular it continues to real s>0,t<4mu,u<4mu. At fixed negative t,
u crosses its massive threshold when s falls below -t; use the specified
complex boundary value, not a real-angle shortcut.

The physical source angular integral is reduced in both axes and retains
its phase-sensitive odd double resolvent. Its entire forward restriction
matches S280 exactly. The invariant dictionary uses the same M as S278's
massive soft kernel, but does not replace the full graviton amplitude by a
soft approximation.

At the forward massless endpoint,

    rho_gg(s,0) ~ mu^(7/2)/(4kappa^2 s^(3/2)).

At fixed strict real angle z, the generic leading behavior is instead

    rho_gg(s,z) =
      mu^3(2z^2-1)^2/[4kappa^2 sqrt(1-z^2) s]
      +3mu^(5/2)(5z^2-1)/(16kappa^2 sqrt(s))+O(1).

These expansions are not uniform as z approaches+/-1. There is no ordinary
uncompensated full low-cut integral from0. For the stated fixed-transfer
channel density, multiplication by the channel energy gives a convergent
Cauchy construction

    F_L(z)=c_L/z+(1/pi)integral_0^L rho(sigma)
                         [1/(sigma-z)+1/z]dsigma.

If a change in L is to add only the ordinary upper shell,
c'_L=-rho(L)/pi is necessary and sufficient. Its initial value remains
an original matching obligation. A crossing-even shift
delta_c(1/s+1/u) changes the v^2 coefficient at s,u=2mu+/-v by
delta_c/(4mu^3). Dropping the pole is not innocuous.

The normal cut is NOT the full imaginary part in a region where other
channel cuts overlap. A channel Cauchy construction does not by itself
assemble the full crossed amplitude or determine local terms, a transfer
subtraction, b20, detector/Regge errors or original P8 closure.
