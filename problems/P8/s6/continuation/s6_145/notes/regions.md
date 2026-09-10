# Symmetric momentum regions with distinct valid soft radii

Let x=q^2,y=l^2, S_q=m^2+x,S_l=m^2+y. Exchange of the two real
loops is a symmetry of the vertex trace (cyclicity and reversal,
or charge-conjugation transposition); p can also change sign.
The full two-point scalar invariant is even. Work on x<=y and
multiply by two. Its MIDDLE region is x<=y<=4S_q; its HIGH region
is y>=4S_q.

On |s-1|<=2, use p=(i sqrt(s),0,0,0), |p|<2. The route prefixes
are below eighteen. For soft scaling on radius sqrt(S)/360 the
Neumann ratio is below 1/2 and each shifted propagator norm is
at most 2/sqrt(S). The degree-four-and-higher Taylor tail is at
most twice the circle majorant divided by that radius to the fourth.

In MIDDLE, the raw four-fermion trace is <=64/(S_q S_l).
Use radius sqrt(S_q)/360. The first subtraction has trace factor
16/S_q times 1/(y S_l), and uses the same radius because K0(l)
is independent of the external soft momentum.
The other subtraction has factor 16/S_l times 1/(x S_q);
use radius sqrt(S_l)/360. Using the first radius for all terms
is unnecessary and loses the required integration control.

The angular average of the unshifted massless chord is 1/max(x,y),
as established in S6.141; the massive scalar chord is smaller.
With common factor B=32*360^4, the raw radial middle term, including
its extra trace factor four, is bounded by

  4 integral_0^infinity x/S_q^3
           log[(m^2+4S_q)/S_q] dx <=4/m^2,

because the logarithm is below log 5<2 and
integral x/S_q^3 dx=1/(2m^2).
The first subtraction has the same integral without factor four.
For the other subtraction, integrate x first:
max(0,y/4-m^2)<=x<=y. The logarithmic ratio is at most five,
both for y<=4m^2 and y>=4m^2. Its bound is also 1/m^2.
MIDDLE therefore contributes at most 6/m^2 on the half-domain,
or 12/m^2 after reflection.

In HIGH the S6.142 resolvent difference estimate applies to
G-K0(l)D2(q,p): its short kernel is bounded by
32sqrt(S_q)/y^(5/2), before couplings/caps.
The two long propagators and trace contribute 16/S_q;
their soft radius is sqrt(S_q)/360. The radial integral is
32 times integral x/S_q^3 dx=16/m^2, in units B.

The complementary K0(q)D2(l,p) keeps radius sqrt(S_l)/360.
Its radial integral, in units B, is

 integral_(4m^2)^infinity y/S_l^3 log[y/(4m^2)] dy
 <= integral_0^infinity y/S_l^3 log(S_l/m^2) dy
 =3/(4m^2).

Thus HIGH contributes 67/(4m^2) on each half, or 67/(2m^2).
Every region, the chord diagonal and both unequal-momentum limits
are controlled. No raw high-region integral is evaluated separately.

The two finite MS anchors each contribute at most
2*integral x/S_q^3 dx=1/m^2 in the same units B.
Combining middle, high and anchors gives 95/(2m^2).
