# Closed one-loop flow and uniform actual-hierarchy bounds

Let L0,g0,M0 be the same initial parameters at nu=1, and put

q=(1-3L0 ell)^(1/3)>0,
L=L0/q^3,
g=g0/q^2,
M=M0+(g0/L0)(1-q).

The positive real root is used before its formal endpoint.
Differentiation using dq/dell=-L0/q^2 checks all three ODEs
and every initial value. Composition of two flow steps
multiplies their roots and gives the same heavy mass shift;
all these identities are checked directly.

Define d=L-3g/M, the completed-square quartic margin.
The exact identities

LM-3g=[L0 M0-3g0+4g0(1-q)]/q^3,
dd/dell=3(L-g/M)^2

show that d stays positive and increases for ell>=0 wherever
this positive-root solution exists. The physical potential
has not acquired a new scale dependence: these are tree
parameter coordinates, with the explicit loop-reference
change compensating them at the retained order.

## The declared reference window

The positive exponential Taylor sum through cubic order at
three is 13>10, so ln(10)<3. Together with the ancestor proof
pi>3, this gives for 1<=nu<=10^400

ell=ln(nu)/(16pi^2)<1200/144=25/3<9.

Set epsilon=27L0. The actual rational value is below 1/4.
On ell in [0,9], q^3>=1-epsilon and q>=q^3; hence
L<2L0 and g<2g0. The elementary positive-gap identities in
the code give

0<=L/L0-1<=54L0<10^-202,
0<=g/g0-1<=108L0<10^-202.

For the second bound use q>=1-epsilon and
(1-epsilon)^(-2)-1<=4epsilon on [0,1/4].
The gap numerator is
epsilon(2-7epsilon+4epsilon^2), which is positive there.

Integrating dM/dell=g gives

0<=M-M0<=18g0<3 times 10^-7,
(M-M0)/M0<10^-203.

Since d>0 implies 0<g/M<L/3, its derivative is at most
3L^2<12L0^2. Therefore

0<=d-d0<=108L0^2,
(d-d0)/d0<10^-5.

All last inequalities use the actual exact rational inputs,
including the small completed-square difference d0; that
difference is not rounded away. The pointwise routine
returns rational enclosures for every ell in [0,9], with
the exact running quartic, and exact initial anchors at zero.

The formal endpoint ell=1/(3L0) lies far outside this
window. It is a singularity of the truncated beta ODE,
not a demonstrated all-orders Landau pole, a UV theorem,
or an exclusion of the original ladder row.

No two-loop running, scheme-conversion remainder, high-energy
scattering amplitude, gravitational correction or Regge
contour is bounded by these estimates.
