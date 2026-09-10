# Subtracted radial part and explicit derivative constants

Let a=x(1-x), ell=1-as, L=y+ell, A=y+M-ell,
C=-lambda4+g/(M-s), V0=C+2g/A and c=g/M^2.
One s-channel angular-zero radial integrand, with the fixed local
subtraction already made, is

f(s,y)=V0^2/L^2-C(s)^2/(y+1)^2.

The t-channel angular-zero piece is independent of s. The u-channel
is obtained by s->4-s. Hence its b2 equals the s-channel one:
two channels times the coefficient factor 1/2 times the radial
prefactor 1/(32pi^2). The x interval has length one.

For the actual M>24 and d<5c, throughout the disc,

|C|<5cM, |C'|<=4c, |C''|<=16c/M,
|V0'|<=6c, |V0''|<=18c/M.

The exact derivatives are V0'=g/(M-s)^2-2ga/A^2 and
V0''=2g/(M-s)^3+4ga^2/A^3. Literal symbolic differentiation checks
these and the formula

(V0^2/L^2)'' =
 2[(V0')^2+V0 V0'']/L^2
 +8a V0 V0'/L^3+6a^2 V0^2/L^4.

For 0<=y<=M use |V0|<=20c(y+1), |L|>=y+1/4 and
|(C^2)''|<=192c^2. The second derivative of f is bounded by c^2 times

72/(y+1/4)^2 +720(y+1)/[M(y+1/4)^2]
+240(y+1)/(y+1/4)^3 +150(y+1)^2/(y+1/4)^4
+192/(y+1)^2.

With radial measure y, use (y+1)/(y+1/4)<=4 and y/(y+1/4)<=1.
The integral over this region is at most

c^2 [3624 log(1+4M)+2880].

For y>=M, make the ultraviolet cancellation explicit **before** bounding:

f=4g C/(A L^2)+4g^2/(A^2 L^2)
  +C^2[L^(-2)-(y+1)^(-2)].

Now |A|,|L|>=y. For T=A^(-1)L^(-2),
|T|<=y^-3, |T'|<=3/(4y^4), |T''|<=3/(4y^5).
The first term's second derivative is therefore bounded by

64c^2 M/y^3 +24c^2 M^2/y^4 +15c^2 M^3/y^5.

Its radial integral is (64+12+5)c^2=81c^2.
For U=A^(-2)L^(-2), the second derivative has coefficients
6,-8,6 multiplying a^2 and denominators of total power six.
The second term is bounded by 5c^2 M^4/y^6 and integrates to 5c^2/4.

For E=L^-2-(y+1)^-2, interpolate the denominator between ell and one.
Since |1-ell|=|as|<=3/4 and the entire segment has real part at least y,

|E|<=3/(2y^3), |E'|<=1/(2y^3), |E''|<=3/(8y^4).

Using |(C^2)'|<=40c^2 M and |(C^2)''|<=192c^2,
the third term is bounded by

288c^2/y^3 +40c^2 M/y^3 +(75/8)c^2 M^2/y^4.

Its radial integral is at most (288/M+40+75/16)c^2,
which is bounded by (288+40+75/16)c^2.
The total high-region constant is 6639/16<415.
All power integrals are convergent and checked exactly.

Combining the two regions gives
integral y |f''| dy <=3624 c^2 [1+log(1+4M)].
After both crossed channels and the coefficient/radial factors, this is
less than 13c^2[1+log(1+4M)].

The actual rational parameters obey 1+4M<10^198. The first four positive
Taylor terms of exp(3) sum to 13>10, proving log(10)<3.
Thus log(1+4M)<594<1000 and the leading b2 error is below 14000c^2.

Together with the angular estimate, the complete one-loop error satisfies

|b2_one_loop| <60c^2 M+14000c^2 <10^-6 (4lambda).

The last comparison is exact rational arithmetic, not floating quadrature.
Both contributions are retained. In particular
60c^2 M/(4lambda)=30g(D/M)^3<30g; the actual total relative majorant
is approximately 4.47 times 10^-7. This is an upper bound, not a numerical
evaluation of the actual loop coefficient.
