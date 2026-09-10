# Regional integration and complete quartic primitive bounds

Use Q=16 pi^2 and the radial measures
x dx/Q and y dy/Q. The LOW raw term has
the angular average 1/max(x,y) and the
joint-radius soft denominator min(S_q,S_l)^2.
Its radial integral is bounded by

 I_raw=integral_LOW x y dx dy /
       [max(x,y) S_q^2 S_l min(S_q,S_l)^2].

For y<=x, swap the radial order and use
integral_y^infinity dx/S_q^2=1/S_l.
This gives integral y dy/S_l^4=1/(6m^4).
For x<=y<=4 S_q the inner integral is
log[(m^2+4S_q)/S_q]<=log 5<2.
Its bound is 2 integral x dx/S_q^4=1/(3m^4).
Therefore I_raw<=1/(2m^4).

The LOW subtraction has its own radius R_q.
Its radial integral is

 I_sub=integral x dx/S_q^4
                 integral_0^(4S_q) dy/S_l
 <=integral x dx/S_q^4 [2+log(S_q/m^2)]
 =17/(36m^4).

Here log(1+4S_q/m^2)<=log 5+log(S_q/m^2).
The elementary strict inequality log 5<2
follows from exp(2)>1+2+2+8/6=19/3>5.
The compact beta moments are
I0=integral_0^1 u(1-u)du=1/6 and
Ilog=integral_0^1 -u(1-u)log(u)du=5/36.

For HIGH, the paired estimate gives a
radial integrand proportional to
x/[S_q^(7/2) y^(3/2)]. Its inner integral is

 integral_(4S_q)^infinity dy/y^(3/2)
   =1/sqrt(S_q).

The remaining integral is 1/(6m^4).
The paired-kernel constant 32 is retained.
The MS local anchor has the same outer
moment 1/(6m^4).

With c=Y+4aCf and the common factor
24 times 128 times 360^4 N Y^2/(Q^2 m^4),
the combined difference has coefficient

 c [4(1/2)+17/36+32/6] =c times 281/36.

The local MS term adds (2Y+6aCf)/6.
Since 2Y+6aCf<=2c, the total is bounded by
c times 293/36. Its exact numerical prefactor
is 419948789760000, strictly below 5e14.
Thus

 E_vertex <=5e14 N Y^2 (Y+4aCf)/(Q^2 m^4).

Every positive soft degree n has the same
regional integrals with denominator power
2+n/2, optionally times log(1+t). They
converge for n>0 in both radial limits.
This validates positive-degree projection
and integration before removing the regulator.

The complete S4 subset is Lorentz invariant
after regulated integration. Degree zero is
a local quartic contact; odd degrees vanish.
Degree two has only diagonal/off-diagonal
Gram orbits, related by momentum conservation,
and is constant on the equal-mass shell.
Therefore those degrees have zero b2.
Their cancellation is not asserted pointwise
for an individually routed integrand.

The uniformly integrated projected amplitude
is holomorphic on the unit forward disc.
A second Cauchy estimate bounds b2 by E_vertex.
A finite field/parameter correction on a new
order-two group changes it at order three;
order-two conversions multiplying old
order-one amplitudes remain separate.

At the shared rational bounds E_vertex is
approximately 5.27382729154e-1404.
Adding the disjoint S6.140 and S6.141 groups
gives 1.26571854997016e-1403 for both complete
quartic primitives, relative to the tree
3.16429637492539e-804. No sign is inferred
from these positive error majorants.
