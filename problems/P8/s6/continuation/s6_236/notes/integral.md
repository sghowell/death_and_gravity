# Exact mass-ordered primitive and its physical branch

The full S235 on-shell triangle denominator is n z+(1-z)²L, with L=1-s xi(1-xi). Its ordered box denominator is n z+(1-z)²L-bz², with b=t eta(1-eta). Therefore define

T(n,L,b)=integral_0^1 (1-z)/[n z+(1-z)²L-bz²] dz,
Cbar(s)=integral_xi T(n,L,0),
Dbar(s,t)=-partial_n integral_xi,eta T(n,L,b).

The derivative holds with L and b fixed. It gives the entire z(1-z)/Q² box measure, not a mass-order interchange.

Let delta=sqrt(n²-4nL+4Lb), alpha=(n-2L+delta)/2 and beta=(L-b)/alpha. The exact factorization is

Q=(L+alpha z)(1+beta z), alpha-L beta=delta.

Partial fractions give the primitive

[(alpha+L)/alpha ln(L+alpha z)-(1+beta)/beta ln(1+beta z)]/delta.

After evaluating BOTH endpoints,

T=[(1+L/alpha)(ln(alpha+L)-Log(L-i0))-H(beta)]/delta,
H(beta)=(1+beta)ln(1+beta)/beta, H(0)=1.

The removable value is essential, including b=L fixtures. The parameters in the stated physical and symmetric domains keep alpha+L and1+beta positive. The only imaginary part there comes from Log(L-i0). For L<0, the upper physical boundary has Log(L-i0)=ln|L|-i pi and the triangle imaginary part is positive. An independent Feynman delta-function calculation at the unique root zstar=-L/alpha reproduces pi(1+L/alpha)/delta. The mass derivative keeps the complete box discontinuity.

The formula is first derived with positive denominators and then continued with the Feynman prescription. The sign and endpoint continuation are also fixed directly by the unique real light root. A double-pole integral through a zero is not treated as an ordinary absolutely convergent integral.

Write w=1/n, d=sqrt(1-4Lw+4Lb w²), a=(1-2Lw+d)/2, v=Lw/a and beta=(L-b)w/a. Then

T=w[P(w)(-ln w-Log(L-i0))+Q(w)],
P=(1+v)/d,
Q=[(1+v)ln((1+d)/2)-H(beta)]/d.

P and Q are analytic at w=0. The logarithm -ln w is NOT analytic there and is NOT subjected to the Cauchy estimate. It is kept explicitly, including its derivative. The analytic H series is

H(q)=1+sum_{k>=1}(-1)^(k+1)q^k/[k(k+1)].

The complete coefficient jets are

P=1+3Lw+O(w²),
Q=-1+(b-7L)w/2+O(w²).

The following notes replace these qualitative O symbols by uniform bounds on their full tails. Points L=0 are interpreted by the integrable light-parameter boundary limit, not by declaring a divergent individual value finite.
