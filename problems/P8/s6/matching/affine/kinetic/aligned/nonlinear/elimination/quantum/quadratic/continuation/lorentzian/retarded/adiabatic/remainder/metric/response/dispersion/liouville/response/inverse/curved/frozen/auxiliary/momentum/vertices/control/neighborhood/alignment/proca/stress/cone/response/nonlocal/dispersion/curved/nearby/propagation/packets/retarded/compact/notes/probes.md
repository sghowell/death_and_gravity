# Explicit smooth compact source and detector

Use the same dimensionless normalized clock, fields and probe units
as S6.91. A comparison in dimensionful parent variables must
restore the corresponding source and canonical normalizations;
none is inferred from the present constants.

Set T=10^-7, s0=-T/4, t0=T/4. The actual S6.91 spacelike
support neighborhoods have radius h=99/(8*10^16). Put

    a=h/10, S_min=T/16, S_max=3T, A_min=10^-8.

Define beta(z)=exp(1-1/(1-z^2)) for |z|<1, zero otherwise.
It is nonnegative, at most one, and greater than one half
on |z|<=1/2: its exponent is at least -1/3 and exp(-x)>=1-x.

## Smoothness and exact normalization

Inside the support,

    beta^(n)(z)=beta(z) P_n(z)/(1-z^2)^(2n),
    P_0=1,
    P_(n+1)=(1-z^2)^2 P_n'
              +[4n z(1-z^2)-2z] P_n.

The recurrence gives a polynomial for every n. At either
boundary put w=1/(1-z^2); any derivative is bounded by a
fixed polynomial in w times exp(-w), which tends to zero.
Thus the zero extension is C-infinity to all orders.
bumps.py checks the recurrence through n=4 as independent
algebraic anchors; the all-orders assertion is the written
induction, not an extrapolation from those checks.
The radial bumps are functions of squared radius, so
they are also smooth at the spatial origin.

Let B1 be integral beta(z) dz over [-1,1], and B3 be
4pi integral_0^1 r^2 beta(r) dr. Neither is set to one.
Use normalized j1(z)=beta(z)/B1 and j3(y)=beta(|y|)/B3.

Since beta>1/2 on the ball of radius 1/2,
B3>pi/12. Also the unnormalized radial bump has squared
L2 norm at most 4pi/3. Consequently

    ||j3||_2^2 < 192/pi <64.

For a positive normal width epsilon, put rho=epsilon/20.
The actual source and detector are

    J(s,y)=rho^-4 j1((s-s0)/rho) j3(y/rho),
    f(t,x)=beta((t-t0)/a) beta(|x_perp|/a)
           beta((x1-S_c(t,s0))/a)
           beta((|x|-S_c(t,s0))/epsilon).

The source is nonnegative, smooth, compact and has spacetime
integral one. Its L1-time/L2-space norm is below 8rho^-3/2.
The detector is also smooth and compact: its spatial support
is separated from the origin, and it is extended by zero
outside a compact interior time interval. The actual
S_c(t,s0), not a frozen straight ray, defines its normal band.

These functions are constructed mathematically from the
validated actual solution. They are not asserted to be a
numerically synthesized pulse or a controlled nonlinear
source/detector apparatus.

## Entire supports stay matter-spacelike

The width chosen later obeys epsilon<a/100. The source
has |s-s0|<rho and |y|<rho, inside the parent's source
neighborhood. On the detector support, |t-t0|<a,
|x_perp|<a and |x1-S_c(t,s0)|<a. Since omega_c<4,

    |x-x0|<a+4a+a=6a<h,
    x0=(S_c(t0,s0),0,0).

Thus every pair in the compact source/detector supports
is matter-spacelike by S6.91, with positive clock-time order.
This is not just separation of the centers. In particular,
the matter delta shell contributes exactly zero to the pairing.

For the source/detector times, 1/4<omega_c<4 and
a+rho<T/4 imply S_min<S_c(t,s)<S_max. All bounds use
the same actual background and physical matter metric.
