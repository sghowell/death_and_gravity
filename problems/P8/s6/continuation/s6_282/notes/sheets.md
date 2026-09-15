# Invariant sheet, uniqueness and removable external threshold

Set Q=s-4mu, h=s/Q, z=1+2t/Q. BEFORE choosing a square root or log,
the full Feynman denominator is

    h-a^2-(1-a)^2-2za(1-a)=[4mu-4t a(1-a)]/Q.

Its symmetric parameter integral is therefore I(h,z)=Q M(t).
The other seed gives Q M(u), while L0=Q M(4mu-s). Substitution into
all four angular coefficients yields the complete invariant formula in
FORMULATION.md. This is not a low-momentum or near-threshold expansion.

For any compact subset of C minus[4mu,infinity), the denominator
4mu-z(1-x^2) is bounded away from zero uniformly on x in[0,1].
Integration and all complex derivatives therefore commute. M is
holomorphic there, with M(0)=1/(4mu). For z<0 it equals

    atanh sqrt[-z/(4mu-z)]/sqrt[-z(4mu-z)].

For0<z<4mu use atan sqrt[z/(4mu-z)]/sqrt[z(4mu-z)].
For z>4mu approached from above, write beta=sqrt(1-4mu/z).
The denominator z(x^2-beta^2)-i0 has a single root in the interval.
The checked logarithmic principal-value primitive and delta Jacobian give

    M(z+i0)=[-atanh(beta)+i*pi/2]/(z*beta).

The lower boundary is its conjugate. Thus no discontinuity is inferred
by silently taking the real part of an overlapping cut.

The nominal Q^-2 term is also holomorphic at Q=0. Let q=Q and B=B(4mu+q).
The exact remainder identity

    1/[4mu+q(1-x^2)]-1/(4mu)+q(1-x^2)/(16mu^2)
      =q^2(1-x^2)^2/[16mu^2(4mu+q(1-x^2))]

and the moments2/3,8/15 give

    [B M(-q)+D/120]/q^2
      =(83mu^2+20mu*q-5q^2)/(120mu^2)
       +integral_0^1 B(1-x^2)^2/[16mu^2(4mu+q(1-x^2))]dx.

The integrand is uniformly regular near q=0. Its exact value there
is29/40. This proves removability even at fixed nonphysical transfer.

The invariant expression is holomorphic where s avoids(-infinity,0]
and t,u avoid[4mu,infinity), with the external threshold filled in.
It agrees on an open real physical region with the original normal cut.
The one-variable identity theorem applied successively in two variables
gives its unique continuation on that connected physical component.
The real set s>0,t<4mu,u<4mu is connected to that region, and the displayed
formula supplies the claimed continuation there. Equivalently, at fixed
real angle r=beta^2 has a uniformly holomorphic original angular integrand
for r outside[1,infinity), including continuation through r=0.

At fixed t<0, however, u crosses4mu as s passes -t. Its M(u) boundary
must be continued explicitly. Moreover massless crossed cuts of the full
amplitude can already overlap elsewhere. A normal Cutkosky channel is a
particular analytic discontinuity component, not the total imaginary
part of the crossed amplitude in every real region.
