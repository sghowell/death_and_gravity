# All twelve routings and a decaying logarithm difference

Use the same analytic forward momenta p1,p2,-p1,-p2
as S6.122. Their bilinear squares are -1 and their
squared Hermitian norms are at most 3/2.

For any of the six unordered external pairs at A
let P be their sum. Either complementary leg p_k
can be attached at B; denote the other by p_l.
This gives twelve labelled routings, all checked.

Route the outer light momenta as k and -P-k.
The inner bubble transfer is Q=k-p_k. Combining
the light denominators and setting q=k+xP gives

    (q^2+Delta)^2,  Delta=1-x(1-x)z,
    Q=q+r,  r=-p_k-xP=-(1-x)p_k+x p_l.

Thus r is a convex combination in Hermitian norm,
so |r|^2<=3/2. Direct bilinear calculation gives
r^2=-Delta. On the forward disc, all channels
satisfy |z|<=3, Re Delta>=1/4 and |Delta|<=7/4.

Let y=q^2, a in [0,1/4] and N=2q.r+r^2.
The square estimate gives

    Re[1+a(y+N)] >= 1+a(y/2-3/2)
                  >= (1+ay)/2.

The final margin is at least 1/8 at the worst
parameter endpoint. The unshifted endpoint 1+ay
has an even stronger margin. Their scalar
interpolation therefore stays in the right
half plane. Integrating the logarithmic derivative,

    ln[1+a(y+N)]-ln(1+ay)
      = integral_0^1 aN/[1+a(y+tN)] dt.

This is one anchored analytic logarithm, not
independent uncontrolled branches.

Since |r|<2 and |r^2|<2,

    |N| <= 4sqrt(y)+2,
    a/(1+ay) <= 1/(y+4).

The latter follows from the exact positive
remainder (1-4a)/[(1+ay)(y+4)]. Consequently

    |ln[1+aQ^2]-ln(1+ay)|
       <= [8sqrt(y)+4]/(y+4) = D(y).

D(y) tends to zero at infinity. It is also
below three: with t=sqrt(y), the numerator
of 3-D(y) is
3(t-4/3)^2+8/3, divided by t^2+4.

Therefore

    |I_R(Q)-I_R(y)| <= D(y)/(16pi^2),
    |I_R(Q)| <= [ln(1+y)+3]/(16pi^2).

The decaying first bound, not merely the second
uniform shift estimate, is essential for the
subtracted local core integral.

The same full-diagram first-sheet continuation
and positive parameter bounds as S6.120-S6.123
justify the Feynman shift. Scale external momenta
from zero; the stated margins persist, and the
subtracted terms below have uniform integrable
majorants. Strict margins allow a neighborhood
of the closed disc for Cauchy's estimate.
