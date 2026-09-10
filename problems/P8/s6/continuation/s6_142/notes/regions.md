# Real-momentum split and the three legitimate Cauchy bounds

As in S6.141, take real q,l on the long and
short arcs. The first boson endpoint carries
q-l. At the other endpoint the fermions have
momenta q+p and l+p, so the boson difference
is unchanged. Fermion momenta gain external
prefix sums; the boson chord never receives
a complex external shift.

Write x=q^2, y=l^2, S_q=m^2+x, S_l=m^2+y.
Split the real integration domain into

 LOW: y<=4 S_q; HIGH: y>=4 S_q.

The split is independent of the complex
external soft-scaling variable zeta. On the
unit forward disc the prefix norm is below
18. The frozen Dirac/Neumann proof gives
||S(k+delta)||<=2/sqrt(m^2+k^2) when
||delta||_1<=sqrt(m^2+k^2)/20.

In LOW, bound the raw kernel with
R_min=min(sqrt(S_q),sqrt(S_l))/360.
The six fermions and Dirac trace give
256/(S_q^2 S_l). Its massless chord angular
average is exactly 1/max(x,y); a positive
scalar mass only improves it.

The LOW subtraction kernel contains no
external momentum. It can and must use
R_q=sqrt(S_q)/360, independently of l.
Only the four long-arc propagators shift.
The trace bound is 64/(S_q^2 S_l y).
Using R_min here would lose the necessary
outer-momentum decay and is not the proof.

In HIGH, use R_q for both arcs. There
sqrt(S_l)>=2sqrt(S_q), so the short arc is
also controlled. For arbitrary two permitted
short shifts the resolvent identity yields

 ||S(l+d1)S(l+d2)-S(l)^2||
 <=(4||d1||_1+2||d2||_1)/S_l^(3/2)
 <=(3/10)sqrt(S_q)/S_l^(3/2)
 <4sqrt(S_q)/S_l^(3/2).

The actual short arc has one shift zero;
the larger two-shift bound is conservative.

Since |q|<=sqrt(S_q)<=sqrt(y)/2,

 B_b(q-l)<=4/y, B_b(l)<=1/y,
 |B_b(q-l)-B_b(l)|<=12sqrt(S_q)/y^(3/2).

For the last inequality, the numerator
|2q dot l-x| is below 3sqrt(S_q)sqrt(y)
and the denominator is at least y^2/4.
These statements hold for b=0 or 1.

Decompose the paired kernel into the
fermion-product difference times B_b(q-l)
plus S(l)^2 times the boson difference.
Its norm is at most

 16sqrt(S_q)/(y S_l^(3/2))
 +12sqrt(S_q)/(S_l y^(3/2))
 <=28sqrt(S_q)/y^(5/2)
 <32sqrt(S_q)/y^(5/2).

This improved bound is for the paired
difference. The raw HIGH vertex is not
integrated independently. Multiplying four
long-arc propagators and the trace gives
64 times this bound divided by S_q^2.

Gauge caps cost at most four gamma-index
terms; color gives Cf=4/3 and N=6 active
entries. The combined kernel coupling
majorant is N Y^2 c, c=Y+4aCf.
The local MS vertex is bounded separately
by N Y^2(2Y+6aCf)/Q and has only the four
long-arc propagators.

All radii are at least two for m>=720.
Apply the degree-four Cauchy tail bound
2 M/R^4 separately to the three regional
pieces and the MS local term before
integrating. These are bounds on the same
linear soft projection, not different
amplitudes or independent subtraction schemes.
