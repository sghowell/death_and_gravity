# Complete optical cuts and the analytic continuation

S285's whole massive densities, with beta=sqrt(1-4nu/a), are

rho_H2=beta*(a-4nu)^2/(3840pi^2),
rho_H0=beta*(a+2nu)^2/(384pi^2),
rho_V2=beta*(13a^2+56nu*a+48nu^2)/(3840pi^2),
rho_V0=beta*(a^2-4nu*a+12nu^2)/(384pi^2).

The radial substitution sigma=4nu/(1-v^2) identifies each subtracted
integral with f_i*pi^2*p integral rho_i(sigma)/
[sigma^3(sigma+p)] dsigma, where f_2=64,f_0=768.
For p=-a-i0 and a>4nu, its imaginary part is equivalently
-pi*w_i(beta)/(2beta)=-f_i*pi^3*rho_i(a)/a^2.

The inverse signs in the insertion therefore give

Im M_X=pi[4rho_X2(a)N2+rho_X0(a)N0]/(kappa^2*a^2).

For scalar production set b=-(a-4mu)(1-z)/2. This is exactly
beta/(32pi)[a0^2+a2^2 P2(z)/5], where
a0=(a+2mu)(a+2nu)/(6kappa*a),
a2=-(a-4mu)(a-4nu)/(6kappa*a).
The same computation matches S281's full nine-polarization Proca cut.
Both include the identical-pair and optical halves. Neither replaces
the complete angular cut by its forward value.

For a closed radial representation put B=4nu+p, r=sqrt(p/B) and

J0=atanh(r)/(B*r),
J_j=(B/p)J_(j-1)-1/[p(2j-1)].

Then integral p*v^(2j)/(B-pv^2) dv=p J_j. Continue from p>0.
The denominator B*r is essential: replacing it with the independently
principal sqrt(p*B) flips the continuation above the massive cut.
Numerical tests use separate radial and spectral integrations at real
and complex arguments, plus lower-bank Feynman boundary values for both
spins and both massive species. A wrong-principal-product negative control
has the opposite imaginary sign. These finite calibrations supplement
the exact residue and dispersion identities, not a proof by sampling.

As nu tends to zero for a DIFFERENT minimally coupled scalar, the density
gives the original M1 cut. The complete crossed known logarithmic part is

MM1=-sum(4N2+10N0)log(-a-i0)/(3840pi^2 kappa^2)
   =-sum[a^2-bc+2mu*a+6mu^2]log(-a-i0)/(960pi^2 kappa^2).

This agrees exactly with S278 at reference scale squared1. For a>0,
log(-a-i0)=log(a)-i*pi. Formal symbolic logarithms in the report denote
this boundary prescription; ordinary principal evaluation at a negative
real argument must not replace it. The nonzero transfer logarithm
remains after universal gravitational soft division. No real local
polynomial is determined from this cut or set to zero physically.
