# Compatible energy faces and the integrable rectangle

Fix a positive angular chart point. Every denominator polynomial has
positive coefficients. For P/Q, its limit as one energy goes to zero
is determined by the minimum power in P and Q. A smaller numerator
power is rejected. A larger one yields zero; equal powers yield the
ratio of leading polynomials. Five elementary rational-limit examples
are checked independently in the tests.

For every component and each held energy, the two sequential limits
onto that energy axis agree by exact cross multiplication. There are
18 such rows per sector/basis case,288 in all. Moreover, writing A
for their common axis value, the verifier proves
 |F-A| <= C(angle)*(sum of the two vanishing energies)*W.
It checks positive denominator coefficients and numerator energy
support, allowing supported midpoint AM-GM monomials if necessary.
At fixed positive angles each collected angular coefficient is finite,
and each denominator coefficient is strictly positive. Thus the
support test gives a finite angle-dependent C. No angular-uniform
Lipschitz constant is claimed here.

Single faces away from axes are ordinary rational limits. Their
nonzero denominators and the joint-axis bounds make these faces
compatible. S318's angular-uniform F=O(W^2) estimate supplies the
zero limit at the origin. The independent mixed derivative bound
is uniform in angles, so the face compatibility need only be pointwise
in angles to pass from a trimmed cube to its boundary by continuity.

The triple fundamental theorem of calculus on [epsilon,a] times
[epsilon,b] times [epsilon,c], followed by dominated convergence,
bounds the alternating eight-corner rectangle by the integral of
721*S/kappa. Define
 I(a,b)=(a+b)log(a+b)-a log(a)-b log(b);
then partial_ab I=1/(a+b), I vanishes on the axes, and
 J=c I(a,b)+b I(a,c)+a I(b,c) has partial_abc J=S.
Therefore ||Delta_a Delta_b Delta_c F||<=721*J/kappa.

For t>=0, the derivative of t-log(1+t^2) is
(t-1)^2/(1+t^2)>=0 and its origin value is zero.
Hence log(1+x)<=sqrt(x), so I(a,b)<=2sqrt(ab).
Consequently J/(abc)<=2(1/sqrt(ab)+1/sqrt(ac)+1/sqrt(bc)).
Also J^2/(abc)<=12(a+b+c), by the three-term Cauchy inequality.
Direct integration over (0,x)^3 gives24*x^2 and18*x^4.
These are bounds for named kernels, not an inclusive event measure.
