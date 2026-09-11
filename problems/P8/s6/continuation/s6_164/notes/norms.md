# The same function class, changed support and full heavy inverse

Take the same real Schwartz Psi with Fourier support
in the Euclidean four-momentum unit ball and Fourier
L1 norm at most one. Let U=||Psi||_2. All jets have
sup norm at most one and L2 norm at most U.

The old cubic map obeys ||R_old||_infinity<=C and
||R_old||_2<=C U, with the exact S6.111 coefficient C.
On this class |X|<=4/kappa. A finite coefficient
majorant gives

    Qcap=sum_(j=8)^15 |[X^j]T8| (4/kappa)^j,
    ||Fhat-Fold||_infinity<=Qcap C,
    ||Fhat-Fold||_2<=Qcap C U.

Each finite monomial estimate puts one jet factor
in L2 and the others in L-infinity. No infinite
spacetime volume factor is used.

Let delta bound Qcap C and A=1+C. The old field has
Fourier radius 3, the new difference radius 33, and
the new squared source radius 66. Derivatives do
not enlarge support; products add it. On the joint
source support |p_Minkowski^2|<=4356. The complete
heavy inverse therefore has norm <=1/(M-4356),
provided M>4356. No heavy-propagator series is used.

The quadratic source difference has norm at most
delta(2A+delta)U. Applying Cauchy-Schwarz and the
Fourier derivative bounds gives the following
action-difference coefficients multiplying U^2:

    free:  100 A delta+545 delta^2,
    local: L delta(2A+delta)[A^2+(A+delta)^2]/24,
    heavy: g[2A^2 Jdelta+Jdelta^2]/[8(M-4356)],
           Jdelta=delta(2A+delta), g=G^2.

The free factors are 3*33+1=100 and (33^2+1)/2=545.
The quartic factor is the exact difference of
fourth powers with two factors placed in L2.
The heavy term compares the full bilinear
stationary action with a bounded real multiplier.

For 0<=C,delta<=1/4 the positive higher powers are
bounded, not dropped. Their corner coefficients are

    free/delta <=1045/4<300,
    local/(L delta) <=671/1536<1,
    heavy/[g delta/(M-4356)] <=671/512<2.

Thus the complete difference is at most

    delta[300+L+2g/(M-4356)]U^2.

The old radius 6 source is retained in its original
estimate. Only the new map difference uses the
enlarged radius 66; applying that radius to the
entire old truncation error would lose useful
precision unnecessarily.
