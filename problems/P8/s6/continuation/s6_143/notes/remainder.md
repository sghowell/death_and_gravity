# Absolute mass-ratio remainder, including endpoints

For 0<r<=1/16 and 0<z<1, put k=r z.
Then

 h(k)=(2k-k^2)/(1-k)^2 <=3k,
 A(k)L(k)=k[-log k]/(1-k)^3
                         <=2k[-log k].

For the first inequality,
3k-h=k(1-5k+3k^2)/(1-k)^2, and its
numerator factor is at least 11/16.
For the second, (1-k)^(-3)<=(16/15)^3<2.
The signs and inequalities are explicit
over the entire parameter interval.

Since 0<log 2<1,
|3-4log 2|<=3. Multiplying by
w=(1-z)^(3/2)/z cancels the factor z
from h and A L. Bound (1-z)^(3/2)<=1.
The first group in the finite integrand
then integrates to at most

 3r[3+2 integral(-log z)dz
           +integral(-log(1-z))dz]=18r.

The second group contributes at most
2r[-log r+1]. Both logarithm moments equal
one. For the second use the explicit real
primitive (1-z)log(1-z)+z and the one-sided
limit z->1 from below; no complex logarithm
branch is part of this interval integral.

Consequently

 |(F_MS-F_MS,0)/(C/Q)|
   <=r[20-2log r]=(20+2log T)/T.

This is a uniform bound on the exact
finite correction, not a numerical fit
at a small mass ratio.

For the actual T=4 times 10^400, the least
integer n with T<=2^n is 1331. Since log 2<1,
the rational correction bound is 2682/T.
The leading finite term satisfies

 47/18+pi^2/12 <71/18<4,

using the elementary pi<4. With Q>=144
and the shared rational Y upper bound,

 |F_MS| <= (12Y_hi/Q_lo^2)[4+2682/T].

The bracket's correction is retained in
the certificate, not rounded to zero.
The actual absolute reference upper bound
is approximately 1.49134255240484e-208.
The finite mass-ratio correction alone
is below 1e-600.
