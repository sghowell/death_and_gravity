# Entire signed-series derivatives and connector regulator limit

Use a>=0 for the physical index and b=(a_e-a)/(2e) for the finite
regulator coordinate. This b is not a radiation energy. Set
 t_e=Gamma(1+2e)*y^(2e), lambda=t_e*(a/(2e)+b),
 Phi_e(lambda)=sum_(N>=0) lambda^N/[N!Gamma(1+2eN)],
 S_e(a,b;y)=exp[-a/(2e)]Phi_e(lambda).
The domain is0<e<=1/8,0<y<=1,|b|<5500/kappa.
Log convexity gives0<t_e<=1. The original hierarchy gives
exp(|b|)<2, and
 exp(|lambda|-a/(2e))<=exp(|b|).
The entire series and all its parameter derivatives converge
absolutely on compact sets. Negative a_e is allowed.

## A global reciprocal-Gamma derivative bound

For1<=z<=2, monotonicity of psi gives |psi(z)|<1. Also
psi'(z)=sum_(j>=0)(j+z)^(-2)>0 on the positive axis; see the
[NIST DLMF trigamma identity](https://dlmf.nist.gov/5.15.E1).
Gamma(z)>=int_1^infinity exp(-u)du=exp(-1).
For z>=2, psi(z)>0 and Jensen on a Gamma(z,1) random variable gives
psi(z)=E[ln U]<=ln E[U]=ln z. The recurrence yields
Gamma(z)=(z-1)Gamma(z-1)>=(z-1)/exp(1), and ln z<=z-1.
Both ranges therefore give
 |d(1/Gamma(z))/dz|<exp(1)<3.
In particular neighboring coefficients
g_N=1/Gamma(1+2eN) satisfy |g_(N+1)-g_N|<6e.

On[1,1+2e], Gamma<=1 and |psi|<1, so |Gamma'|<1.
This implies
 (1-t_e)/(2e)<=1+|ln y|,
using1-y^(2e)<=2e|ln y|.

## Uniform parameter derivatives

The derivative in b is
partial_b S_e=t_e exp[-a/(2e)]Phi'_e, whose absolute value is<6.
For a, the apparent1/e cancels after writing
 t_e Phi'_e-Phi_e=(t_e-1)Phi_e+t_e(Phi'_e-Phi_e).
The first term costs<6(1+|lny|); the second costs<6 by the
neighboring reciprocal-Gamma coefficient bound. Hence
 |partial_a S_e|<12(1+|lny|), |partial_b S_e|<6.
The argument bounds the absolute series, not a supposed positive
probability distribution in fractional dimension.

The segment joining the elastic and radiative (a,b_e) is in the same
convex domain. Combine the derivatives with |delta a|<40R/kappa
and |delta b_e|<8000R(1-lnR)/kappa to obtain
 |S_e(sigma;y)-S_e(0;y)|
 <50000R[1-lnR+|lny|]/kappa.

For the connector, y=x-R. The nonnegative elastic Born seed has
angular integral a0*dR/R. Thus the bound is integrable at BOTH
endpoints, with
 int_0^x [1-lnR-ln(x-R)]dR=x(3-2lnx).
The fixed-x entire-series limits from S309 can therefore be passed
under the seed integral by dominated convergence:
 E1(x)=int dB1[P_sigma_R(x-R)-P0(x-R)].
Neither unsubtracted Born-seed integral has been assigned a separate
finite value. The proof does not take a fixed-N limit before summing,
and makes no assertion for arbitrary simultaneous e*|lnx| scaling.

## Independent calibration and explicit logarithm domain

Six independent series/derivative calibrations include zero physical
index and negative regulated index.240- and360-term sums agree within
10^-60; numerical differentiation agrees with both explicit derivative
formulas within10^-60 at80-digit precision. Four positive/negative
finite-conversion connector examples agree with independent Beta
moment sums to the same tolerance.

The real integral int_0^1 -ln(1-z)dz is evaluated using the explicit
positive-domain substitution u=1-z, giving int_0^1 -ln(u)du=1.
An initial scratch CAS primitive returned1-2I*pi on the shifted form,
which the minimal original interpreter reproduced. No imaginary term
was discarded from a physical amplitude; the real integration domain
and its substitution were made explicit before evaluating the integral.
The backend, identity, physical parameters and tolerances were unchanged.
