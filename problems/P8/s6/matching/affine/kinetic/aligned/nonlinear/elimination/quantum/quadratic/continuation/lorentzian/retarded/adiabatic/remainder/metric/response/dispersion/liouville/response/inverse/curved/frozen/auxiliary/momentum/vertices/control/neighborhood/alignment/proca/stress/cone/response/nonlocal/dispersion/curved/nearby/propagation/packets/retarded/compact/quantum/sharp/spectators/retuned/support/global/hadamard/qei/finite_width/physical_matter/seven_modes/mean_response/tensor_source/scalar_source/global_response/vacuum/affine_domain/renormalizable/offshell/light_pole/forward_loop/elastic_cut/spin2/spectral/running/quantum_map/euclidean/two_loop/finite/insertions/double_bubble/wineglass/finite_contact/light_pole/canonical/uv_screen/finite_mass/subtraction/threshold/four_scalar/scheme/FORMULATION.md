# P8 S6.133 — complete canonical one-loop reference matching

Original P8 remains OPEN. This checkpoint combines all
one-loop scalar and fermion contributions to the light
scalar pole and forward second coefficient in a single,
explicit reference scheme for GY14-unbroken.

## Reference boundary

Take the prospective MS interaction parameters at
mu=mF=10^200 to be the named L,G,M,a,Y boundary.
M is the heavy scalar mass squared, g=G^2, N=6 and
Q=16pi^2. Physical light mass, H one-point function
and vacuum energy use explicitly fixed local references;
they are not falsely described as purely MS conditions.
The free internal light scalar propagator has mass one.
The large local mass counterterm is paired with the
fermion subgraph, never used as a negative free mass.

Let ell=log(mu^2/m_phi^2)=400log(10). The old complete
scalar amplitude subtracts I0=Ibar_MS+ell and a
previously fixed finite quartic contact sigma.
That canonical scheme is not MS with unchanged numbers.

## Complete field and coupling conversion

Let r=Pi_scalar'(1)>0 in the old inverse convention
s-1+Pi, and fp=f_fermion'(1)>0 in S6.131's continued
Euclidean inverse convention. If zF=-2NY Ibar_MS/Q,
the MS and canonical-star bare field factors through
one formal loop order h are

    Z_MS=1+h zF,
    Z_star=1+h(zF-r+fp).

Matching the bare couplings gives the parameters for
the old-scalar-subtraction plus new-fermion-box
canonical representation:

    delta G_star=-LG ell/(2Q)+(fp-r)G,
    delta M_star=-g ell/(2Q),
    delta L_star=-3L^2 ell/(2Q)-sigma+2(fp-r)L.

The MS bare quartic pole reproduces the inherited
one-loop marginal flow once its field factors are kept.
Both a direct bare-parameter comparison and an independent
total-counterterm comparison verify this conversion.

## Actual one-loop amplitude

With A0=-L+g sum_z (M-z)^-1 and C(z)=-L+g/(M-z),

    A_new = A0+h[A_scalar,old^(1)+A_fermion^(1)
                 +ell sum_z C(z)^2/(2Q)+sigma
                 +2(fp-r) A0] + O(h^2).

All scalar graphs are taken from their complete old
calculation and converted, not transplanted with
unchanged canonical parameters. Every one-loop fermion
box is included from S6.132. There is no direct Phi
gauge vertex or H Yukawa supplying another graph at
this order. The canonical-star representation has unit
light residue and gets no further LSZ copy.

At t=0,u=4-s the known reference tree coefficient is
4lambda. The complete one-loop error, including the
finite scale conversion and both field factors, is
strictly below 10^-6 of it. The resulting formal
one-loop coefficient interval is therefore positive.

## Local vacuum and pole

The H one-point counterterm cancels the stationary
heavy tadpole term. The light scalar mass reference
retains both complete scalar and fermion anchors.
The vacuum-energy reference retains both real scalars,
all fourteen Dirac flavors and the zero scaleless
one-loop gauge/ghost term in this regulator.

The normalized light kernel has mass-one pole and
unit residue, with complex unit-disc remainder
coefficient below 2 times 10^-405. It has no other
pole on that disc and has positive local Phi potential
curvature. No global quantum potential or full
all-momentum/gauge-spectrum result is asserted.

These are complete **formal one-loop** matching results.
They do not inherit the old pure-scalar two-loop bound,
bound all later new-model terms, justify a V high-energy
contour, derive finite-gravity Delta, supply common
bounce-parent B, exclude a full row, or close P8.
