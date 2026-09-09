# Polynomial transfer growth and spacetime test functions

All times belong to the same closed inner interval I of length
T=10^-7. S6.90/S6.91 give, uniformly for ordered endpoints and
k>=K0=4*10^31, the exact packet transfer factorization
T_Y=S_t C_t V(t,s) C_s^-1 S_s^-1.
Its infinity norm is at most 1000*(3/2)*2*2*1000=6*10^6.
In particular this is a bound on exact evolution, not only its
leading oscillatory approximation.

Using the actual coefficient bounds R in (1/2,2), abs(alpha)<1,
abs(r_N)>1/10, the explicit Laurent maps obey
norm(E)<=24 k^2 and norm(E^-1)<=21 k for k>=1.
The native majorants check every entry of E/k^2 and E^-1/k;
all remaining k powers are nonpositive. Their row-sum bounds
are 12*sqrt(2)<24 and 21*sqrt(2)/2<21.
Thus norm(U)<=24*21*6*10^6 k^3<10^10 k^3.

On 0<=k<=K0, use S6.91's complete original generator bound B,
including its finite-q polynomial terms. The two D_R endpoint
norms are at most eight, so norm(U)<=64 exp(TB), a finite
constant. B and the unevaluated exponential are retained in the
report; no numerical representability or modest size is needed.
Reverse-time transfers are controlled either by the same finite
interval argument or U^-1=-Omega U^T Omega. Infinity and transpose
norms differ by at most the fixed dimension factor four.

Together these give norm(U(t,s,k-vector))<=C(1+|k-vector|)^3
with a finite uniform C. Parameter-dependent linear ODE theory
gives smoothness in the Cartesian momentum variables, including
the origin, using the original polynomial generator.

Here is an all-orders bound rather than a finite jet check. For
a multi-index alpha of positive order, differentiate the original
ODE and use variation of constants:
partial^alpha U(t,s)=integral_s^t U(t,r)
 sum_{0<beta<=alpha} binom(alpha,beta)
 (partial^beta M_rho(r)) (partial^(alpha-beta) U(r,s)) dr.
Derivatives of M_rho vanish for |beta|>4; otherwise they are
polynomial of degree at most 4-|beta| with uniformly bounded
time coefficients. Induction gives the sufficient degree
d_alpha=3+7|alpha|. A nonzero summand's degree is at most
3+(4-|beta|)+(3+7(|alpha|-|beta|))
=10+7|alpha|-8|beta|<=d_alpha.
All finitely many coefficient and time integrals are bounded on I.
Repeated time derivatives add finite polynomial factors because
the actual background and generator coefficients are analytic.

Leibniz's rule now shows U and U^T map Schwartz Cauchy tests
continuously into Schwartz tests, uniformly on I; the same holds
for all test derivatives. For any smooth compact spacetime test f,
the Cauchy test obtained by integrating U(t,t_*)^T times its
Fourier transform is Schwartz. Smooth multiplication by V(t)
does not change this property. These estimates establish
tempered linear fields and well-defined compact spacetime
smearings. They are not Hadamard, energy or stress-tensor bounds.
