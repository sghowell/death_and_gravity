# Exact failures of two recursive Ward shortcuts

## A bounded conserved outer source does not control an isolated cluster

Put w=1/32 and, for r>0, sn=2r/(1+r^2), cs=(1-r^2)/(1+r^2).
Take q1=w(1,0,0,1), q2=w(1,sn,0,cs), q3=w(1,1,0,0).
Use spatial plus tensors t*t^T-b*b^T with transverse pairs
(x,y), ((cs,0,-sn),y), (y,z), respectively. They have Frobenius
norm sqrt(2). Their total Q has energy3/32 within the original domain.

Let v=(Qx/Q0,1,0,0) and U=v*v^T. Then Q^T eta U=0 and
||U||_F=1+(Qx/Q0)^2<=13/9.
Contract the propagated complete pure-three-soft current with U at
unit canonical coupling. Its four trees (one quartic and three pairs
of cubics) give exactly

A(r)=2*(6r^12-69r^11+16r^10+111r^9-342r^8+796r^7-942r^6
       +1176r^5-904r^4+705r^3-378r^2+161r-48)
     /[9*(r-1)^2*(r^2+1)^4*(2r^2-r+1)].

At r->1, q2.q3=w^2*(r-1)^2/(1+r^2), but Q^2->1/256 and
||U||_F^2->169/81. Nevertheless (r-1)^2*A(r)->2.
Rescaling each physical polarization by1/sqrt(2) and U by9/13
makes every norm at most1 and leaves the nonzero residue9/(13sqrt(2)).
Restoring kappa multiplies this cluster by1/kappa. Thus every fixed
finite kappa still has this angular pole in the named graph class.

The other pair limit r->0 is finite, A->-32/3. Testing only that
direction would have missed the obstruction. The full propagated
root residue matrix at r->1 is

[[-18,12,0,12],[12,-6,0,-12],[0,0,6,0],[12,-12,0,-6]].

The complete pure-three-soft root is itself conserved for symbolic r.
There is no contradiction: after cutting the q2+q3 branch, the remaining
single cubic connecting q1 to the propagated U is NOT a complete
free-leaf hard current. At r=1 its Ward defect is
[1/32,1/24,0,1/32], not zero.

## The actual original hard source also has the isolated-cluster pole

This is not limited to an arbitrary illustrative tensor U. At r=1
take incoming E=5/4 along opposite z directions, ep=77/64, and
dy=sqrt(457)/32. The outgoing massive vectors are
(ep,-Qx/2,+dy,-Qz/2) and (ep,-Qx/2,-dy,-Qz/2).
All four massive shells and total momentum conservation are exact.

The complete47-tree hard root U(Q) is conserved, at both diagnostic
and original parameters. Its contraction with the displayed residue,
divided by kappa, is nonzero in both cases. The exact values, not
floating-point signs, are retained in the report. The corresponding
nested cubic Ward defects are also nonzero.

For nearby r, keep ep=(2E-Q0)/2 and replace dy by
sqrt(ep^2-1-|Qspatial(r)|^2/4). Q has zero y component and dy^2
is positive at r=1. This defines a smooth physical recoil near that
point, with the hard current continuous since its denominators there
are nonzero. Its value at r=1 therefore fixes the nonzero residue
of this actual isolated graph class.

This is NOT a divergence statement for the complete physical
three-real amplitude. Other graphs completing each fixed-pair
remainder supply the Ward cancellations that the isolated class lacks.

## Two off-shell leaves break the free-leaf induction premise

For four soft rays of energy1/32 in the directions +x,-x,+z,-z,
the complete26-tree one-off-shell graviton root is conserved.
Use exact outgoing massive recoil with ep=19/16 and opposite
directions(3/5,0,4/5), preserving incoming E=5/4.

Collapse the third/fourth pair into its propagated off-shell current
and insert it as a leaf in the full434-tree four-Phi/one-h remainder
whose root is the first/second pair momentum. Direct exact evaluation
at the diagnostic parameters gives a nonzero Ward defect.
The report retains that exact vector. This counterexample is sufficient
to disprove an unrestricted multi-off-shell conservation assertion;
it is not labeled an original-parameter physical observable.

The inserted leaf is conserved but is not a free on-shell wave.
Its inverse-kinetic residual is the nonzero pair source. The finite
Noether identity consequently contains additional terms. One cannot
apply the single-root pair bound repeatedly through such an intersection
without estimating those terms and the accompanying complete graph sums.
