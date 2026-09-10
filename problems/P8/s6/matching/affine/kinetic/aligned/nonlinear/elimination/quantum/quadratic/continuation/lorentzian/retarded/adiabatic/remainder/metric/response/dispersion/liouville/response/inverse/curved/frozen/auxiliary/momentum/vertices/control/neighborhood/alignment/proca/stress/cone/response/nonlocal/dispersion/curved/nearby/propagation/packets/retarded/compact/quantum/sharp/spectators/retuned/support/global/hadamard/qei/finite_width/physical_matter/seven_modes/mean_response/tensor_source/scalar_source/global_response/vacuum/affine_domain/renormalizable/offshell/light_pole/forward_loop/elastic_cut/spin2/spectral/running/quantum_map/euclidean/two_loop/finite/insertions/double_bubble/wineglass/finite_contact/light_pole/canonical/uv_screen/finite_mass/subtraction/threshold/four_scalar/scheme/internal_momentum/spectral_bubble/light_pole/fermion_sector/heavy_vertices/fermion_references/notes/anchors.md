# Finite scalar and gauge local anchors

Give the fermion denominator Feynman weight x.
With boson mass squared b=1 for scalar exchange
and b=0 for gauge exchange,

    Delta=x m^2+(1-x)b+x(1-x)p_E^2.

After the common regulated momentum shift the
fermion numerator has external weight (1-x)p_E.
At p=0 define J0=integral log(Delta/mu^2),
J1=integral (1-x)log(Delta/mu^2) and
R=integral x m^2/Delta. All x integrals run from
zero to one. The scalar inverse has UV mass
coefficient -Y m/Q and kinetic coefficient
Y/(2Q), with Q=16 pi^2.

At mu=m, put r=b/m^2. Direct primitives give

    J0=-1-r log(r)/(1-r),
    J1=[-3/4-(r-r^2/2)log(r)+r-r^2/4]/(1-r)^2,
    R=1/(1-r)+r log(r)/(1-r)^2.

The r=0 limits are -1,-3/4,1. For 0<=r<=1,
x<=x+(1-x)r<=1 implies J0 in [-1,0] and
J1 in [-3/4,0]. The inequalities
x<=x/[x+(1-x)r]<=1 give R in [1/2,1].
The endpoint logarithms are integrable; a
massless internal gauge line causes no IR
divergence at this off-shell reference.

Keep the full dimensional gauge numerator.
The mass and kinetic factors multiply
(4-2 epsilon)(Ibar-J0) and
(2-2 epsilon)(Ibar/2-J1), respectively.
After the MSbar pole subtraction their finite
parts are -4J0-2 and -2J1-1. At b=0, mu=m
they are two and one half. Early d=4 would
incorrectly lose the constants -2 and -1.

The scalar Yukawa vertex at zero scalar momentum
is the derivative of the proper mass inverse
with respect to its background mass m+y_i Phi.
HOLD mu FIXED during this derivative:

    derivative_m[m J0]=J0+2R.

Only afterwards evaluate mu=m. The scalar finite
vertex ratio is Y(J0+2R)/Q; the gauge ratio is
-6a Cf/Q. Resetting mu=m before differentiation
instead produces +2a Cf/Q for the gauge term,
missing eight units. The independent tests
detect this error. The scalar propagator has
no linear Phi-background change at this vacuum,
so no additional first-background derivative
of its mass is missing here.
