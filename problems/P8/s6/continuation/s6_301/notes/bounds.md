# Uniform collinear bound and differentiated dimensional integral

All norms below are Euclidean nuclear norms unless explicitly stated.
They bound both the Frobenius norm and the absolute trace. Let tau=1-t,
s=sqrt(tau), and compare P_t with P_1 in the common R4 embedding.

For rank-one projectors ||P_t-P_1||=sqrt(1-t)=s. A future massive
leg has energy<=2, spatial norm<2 and Doppler denominator>=1/4.
Writing its projected vectors a,b, the rank-one inequality

||a*a^T-b*b^T||_nuclear <= (|a|+|b|)*|a-b|

gives a numerator difference<=8s. The inverse denominator changes by
at most 32tau. Multiplying the boundary numerator bound4 and adding
the first term gives a single massive-current difference<=160s.

For a null leg omega*(1,m), let c=m.n and a=1-c. At every radial value
its nuclear norm is omega*(1+sqrt(t)*c)<=2omega. Thus the difference is
at most4omega. Also 1-sqrt(t)*c >= (1-c)/2 for all c in [-1,1].
The same rank-one inequality and 1-sqrt(t)<=tau give the coarser bound

difference <= 12omega*s/a.

Average the minimum of these bounds on S2. The cap a<=s contributes
at most2omega*s, and the remaining directions contribute
6omega*s*ln(2/s). Consequently

mean difference <= omega*s*[2+6ln(2/s)] <= 14omega*tau^(1/4).

For the final inequality, ln2<1, s<=sqrt(s) and
-sqrt(s)*ln(s)<=2/exp(1)<1. This proof includes arbitrarily collinear
directions and has no minimum separation or multiplicity assumption.

Each massive current has norm<=E+|p|<=4: its numerator/denominator is
E+z-1/(E-z), where z is the projected spatial component.
Each null current has norm<=2omega. Hence the sum has norm and absolute
trace bounded by B=16+2R<=65/4 at all radial points. With four massive
legs the averaged current difference is at most

(640+14R)*tau^(1/4).

The quadratic differences therefore satisfy

|H(t)-H(1)| <= 3B*(640+14R)*tau^(1/4) < 32000*tau^(1/4),
|U(t)-U(1)| <= 2B*(640+14R)*tau^(1/4) < 21000*tau^(1/4).

At the boundary the tensor lives in the two-dimensional transverse
plane, so 0<=H(1)<=B^2<265. Also U(1)<=B^2. At any interior point
|F_e|<=3B^2/2<400 for e>=0; positivity there is not needed.

The exact logarithmic integrals are

integral_0^1 tau^(-3/4)*|ln(tau)|^j dt = j!*4^(j+1).

Thus the radial integral for K1 is absolutely integrable uniformly
over the entire stated family and |K1|<130000. Bounded convergence
alone would not justify this derivative; the quarter-power estimate
is the missing argument.

For 0<=e<=1/8, the digamma series gives
0<psi(3/2+e)-psi(1+e)<=2-2ln2<1.
Therefore A(e)<=exp(e)<2 and |A(e)-1|<=2e.
If I_e denotes the radial H difference integral, then
|I_e|<=128000 and |I_e-I_0|<=512000e, using
1-tau^e<=e*|ln(tau)|. The analogous U difference integral is<=84000.
The exact trace coefficient obeys |c(e)-e/2|<=e^2/2.
Substituting into the exact expression in notes/angular.md gives

|K_e-K0-e*K1|
 <= e^2*[B^2/2+2*128000+2*512000+84000]
 < 1400000*e^2.

S296's unchanged phase p(e) satisfies p<=1, |p'|<4,
|p(e)-1-p'(0)e|<=17e^2/2, p'(0)=EulerGamma-2-lnpi.
It follows that

|p(e)*K_e-K0-e*(K1+p'(0)*K0)| < 2000000*e^2.

These are conservative analytic bounds, not empirical quadrature errors.

## Finite angular-energy measures

For t<1 the null-current kernel is a bounded continuous function of the
radiation direction; the recoil is continuous in its total momentum.
At t=1 the kernel's only discontinuity, as a function of a radiation
direction, is its coincidence with n. A finite positive measure has
only countably many atoms, so for almost every n weak convergence
gives convergence of the current. Its uniform B bound permits angular
dominated convergence. The radial Holder majorants above then permit
dominated convergence of K1 itself. This supplies a continuous finite
positive angular-energy-measure extension of the coefficient and its
bounds, not an interacting quantum state or a functional path measure.
