# Constructive nonlinear shape and complete contacts

Use Fourier series on the finite2pi L torus. For scalar coefficients
use absolute value and for complex symmetric matrix coefficients use
operator norm. Sum with weight(1+|k|)^8. The triangle inequality gives
submultiplicative weights, so full convolution and matrix products are
continuous Banach-algebra operations. Reality means conjugate coefficients
at opposite momenta. All stated estimates also hold at exponent0.

Define P(k)=I-kk^T/k^2 at k!=0 and P(0)=I. Then Bf(k)=P(k)f(k)
has zero divergence and norm<=||f||. Its trace is Lop f=(2I+Pi0)f,
whose inverse has norm<=1/2. The mean eigenvalue3 is essential.
The orthogonal TT projection is PAP-P trace(PA)/2 at nonzero k and
A-I trace(A)/3 at zero. Its zero sector has five, not two, directions.

For ||tau||<=1/100, set Q=I+tau+Bf. Every symmetric transverse shape
near I decomposes uniquely this way, with tau=P_TT(Q-I) and
f=Lop^-1 trace(Q-I). Write E2(A)=((trace A)^2-trace(A^2))/2.
The exact3x3 determinant identity gives the fixed point

f=-Lop^-1[E2(tau+Bf)+det(tau+Bf)].

Full-convolution bounds give ||E2(A)||<=6a^2 and ||det A||<=6a^3.
Their derivative bounds are12a and18a^2, respectively; these follow
from the trace bound3||A|| and the six full determinant products.
On ||f||<=1/1000, a=11/1000, the map norm is at most
3a^2+3a^3=366993/10^9<1/1000, and its Lipschitz constant is
6a+9a^2=67089/10^6<1/10. Banach contraction proves a unique full
scalar branch in that ball, not just a formal series or projected root.

Reality is preserved by every iteration. Pointwise ||Q-I||<=a<1/8
makes Q positive. It has determinant exactly1 and divergence exactly0.
The weighted algebra embeds into C8 on this finite torus. For smooth
tau, smoothness of f also follows by differentiating the analytic
Banach map along spatial translations. This does not assume that
composition in a same-regularity diffeomorphism group is differentiable.

The complete scalar derivative is K_Q h=cof(Q):Bh.
Since K_Q=Lop plus the derivative of E2+det, the same Neumann
estimate gives ||K_Q^-1||<=1/[2(1-q)]<5/9.
Thus delta f=-K_Q^-1(cof Q:delta tau). Cof Q=Q^-1 on this exact
unit-determinant slice. Every inverse and contraction remains in the
full convolution algebra.

For ||tau0||<=1, define coefficients by
f(epsilon tau0)=epsilon^2 f2+epsilon^3 f3+... .
Trace tau0=0 and the full determinant identity yield

f2=(1/2)Lop^-1 trace(tau0^2),
f3=Lop^-1[trace(tau0 Bf2)-det tau0].

These are coefficients, not derivatives: derivative factors2! and3!
must be restored. The bounds are ||f2||<=3/4 and ||f3||<=33/8.
The second term generates a homogeneous scalar even when tau0 has
zero mean. Dropping it changes det Q.

For |epsilon|<=1/100 the cubic trial lies in the scalar ball.
Put S(t)=t+(3/4)t^2+(33/8)t^3. In the full residual, exact coefficients
through degree3 cancel. The remaining positive majorant is the
degree>=4 part of6S^2+6S^3. Dividing by t^4 and setting t=1/100 gives
173086592988411/2560000000000<100. Applying Lop^-1 and the
contraction inverse yields root error<=50/(1-q)|epsilon|^4
<60|epsilon|^4. This controls the entire nonlinear branch.

Finally gamma=a^2 exp(2v)Q^-1, so det gamma=a^6 exp(6v).
This leaves the independent scalar volume unchanged within this
slice parametrization. It does NOT erase the S258 scalar-volume
correction needed when moving a generic off-slice metric into the slice.
