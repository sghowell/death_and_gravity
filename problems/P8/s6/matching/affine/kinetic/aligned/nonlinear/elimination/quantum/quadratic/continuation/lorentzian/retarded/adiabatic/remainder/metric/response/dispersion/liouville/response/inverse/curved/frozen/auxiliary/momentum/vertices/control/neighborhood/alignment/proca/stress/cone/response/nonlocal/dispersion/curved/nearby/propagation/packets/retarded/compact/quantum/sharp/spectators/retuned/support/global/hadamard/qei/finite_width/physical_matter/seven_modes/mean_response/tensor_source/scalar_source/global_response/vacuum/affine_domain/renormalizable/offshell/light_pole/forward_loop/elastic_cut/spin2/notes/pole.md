# Fixed-negative-transfer pole subtraction and its finite remainder

Use the propagator and matter-vertex normalization in FORMULATION.md.
For equal-mass external particles let v=(s-u)/2, a=1-t/4 and P1.P2=v/2.
The literal four-dimensional contraction of both full stress tensors is

N=T1.T2-(Tr T1)(Tr T2)/2
 =F(t)^2(v^2-2a^2)+F(t)D(t)a t-3D(t)^2 t^2/8.

Thus the coefficient b2 of v^2 in this t-channel exchange is
-F(t)^2/(M_Pl^2 t). The D form factor cannot change that coefficient.
The code contracts every tensor component in an explicit kinematic
frame and reduces the result to the invariants before differentiation.

Introduce a loop marker hbar and F=1+hbar f, where
f=F_loop(t)-F_loop(0). The one-loop coefficient is -2f/(M_Pl^2 t).
The f^2 contribution is two-loop and is not retained as a resummation.

The tree pole is subtracted at fixed negative t before taking t->0-.
The analytic bound in triangles.md gives, for |t|<=1,

delta b2_t,vertex = -2 F_R'(0)/M_Pl^2 < 0,
|b2_t,vertex(t)-delta b2_t,vertex|
 <=2 F_R'(0)|t|/(3M_Pl^2).

At the actual selected M_Pl=10^400 in light-mass units,

|delta b2_t,vertex| < g/(1296 M^2 M_Pl^2) < 10^-1205.

This is a finite matter vertex correction at order M_Pl^-2.
It is not a forward-limit prescription for a full gravitational
scattering amplitude. Crossed graviton channels, additional matter
graphs at that order, graviton loops and propagator corrections at
higher inverse-Planck order have not been summed or bounded here.
Massless multiparticle cuts require their own infrared treatment.

As a separate normalization diagnostic only, set M=1. For the same
two-distinct-field loop, the exact parameter integral is

J=integral_0^1 [x(1-x)/(1-x+x^2)]^2 dx
 =(45-8pi sqrt(3))/27.

An identical single-field cubic loop instead has the extra symmetry
factor one half, so its slope is g J/(192pi^2). With c2=2b2 this gives
delta c2=-(45-8pi sqrt(3))g/(1296pi^2 M_Pl^2), matching the cited
single-scalar benchmark in light-mass-one units. That diagnostic is
not the actual heavy-mass calibration.

No sign conclusion about a full dispersive gravitational inequality
follows from the negative partial correction. The high-energy
residue, trajectory and contour contributions remain to be derived.
