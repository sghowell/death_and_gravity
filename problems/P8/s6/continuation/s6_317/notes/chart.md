# Independent finite nonlinear coordinate chart

## Finite labeled ring

Use separate commuting nilpotent variables z_i, z_i^2=0, for the free
waves. For a subset S write z_S=product_(i in S)z_i and
Q_S=sum_(i in S)q_i. At each fixed finite number of leaves this ring
has finite dimension; there is no analytic-series convergence premise.

At unit canonical coupling write g=eta+sum_(S nonempty)g_S z_S exp(iQ_S.x).
The coefficient g_S is2H_S. Consider the formal, possibly complexified
coordinate change

x -> x+i*sum_(S nonempty)Z_S z_S exp(iQ_S.x).

This is a device for exact wave coefficients, not a global spacetime
coordinate or quantum-state construction. Its Jacobian coefficients
are J_0=I and J_S=-Z_S*(eta Q_S)^T.

The pulled-back metric is J^T g(x+iZ) J. For an original target subset T,
translation multiplies g_T by exp(-Q_T^T eta Z). Its coefficient on
a disjoint remaining subset U is evaluated by an anchored partition:
choose the block B containing the least label of U, multiply by
-Q_T^T eta Z_B, and recurse on U minus B. This enumerates every unordered
partition once, exactly reproducing the factorials in the exponential.

Finally sum J_A^T*g_translated,B*J_C over ordered disjoint A,B,C
whose union is S. This keeps both Jacobian factors and all translation
terms. The implementation is independent of the recursive current
projection being checked.

## Unique solution at each degree

Suppose all lower-degree Z are fixed. Let B_S be the coefficient of the
pullback computed using only those shifts. A new Z_S cannot multiply a
nonconstant field or another nonempty shift at degree S: that would
require additional labels or a repeated nilpotent label. It therefore
enters only through the background Jacobian term

g'_S=B_S-(q_S z_S^T+z_S q_S^T),   z_S=eta Z_S.

Here z_S in this equation denotes the shift COVECTOR, not the earlier
nilpotent monomial. The equations g'_(S,0mu)=0 have the unique solution

(z_S)_0=B_(S,00)/(2Q_(S,0)),
(z_S)_i=[B_(S,0i)-(q_S)_i*(z_S)_0]/Q_(S,0).

Every nonempty pure-soft subset has Q_(S,0)>0. Thus the new chart inverse
has only a positive energy denominator, not a new angular-collision
denominator. Existing angular poles in B_S are not bounded by this fact.

For a single field tensor H the associated linear projector Pi_Q is
P^T H P, where P's first column vanishes and its spatial columns are
e_i+(Q_i/Q0)e_0. Consequently

P^T P=diag(0,I_3+v v^T), v=Q_sp/Q0,
||Pi_Q H||F<=||P||op^2 ||H||F<=2||H||F

for future causal Q. The coefficient shift obeys
||Z_S||<=2||B_S||F/Q0 and ||J_S||F<3||B_S||F, by |v|<=1.
These are chart bounds, not a bound on the field B_S.

## Explicit nonlinear degree-three term

Physical spatial TT singleton waves already satisfy the temporal
condition, so their first shifts vanish. For a pair field H_ij,
at kappa=1, let

z0=H00/Q0, zi=(2H0i-qi*z0)/Q0, Z=eta*z.

At triple degree the lower shifts induce

C=sum_(ij,k)[(q_k^T eta Z_ij)*eps_k
            +(eta Q_ij)*(eps_k Z_ij)^T
            +(eps_k Z_ij)*(eta Q_ij)^T].

The new field is Pi_Q(H_123-C), not Pi_Q H_123.
This follows directly by expanding the translation and Jacobians in
J^T g(x+iZ)J. Independent generic-r literal pullback and recursive
calculations reproduce every matrix entry of this formula.

For S316's symbolic cluster, the residue at r=1 of the linear-only
projection is

[[0,0,0,0],[0,2,0,-4],[0,0,6,0],[0,-4,0,0]].

With C retained the residue is zero and the finite limit is

[[0,0,0,0],[0,-5,0,35/6],[0,0,71/2,0],[0,35/6,0,-61/6]].

The independent four-leaf chart has11 nonzero shifts. Its highest
coefficient differs from a linear-only projection, and every one of
its15 nonempty subset coefficients agrees with the new recursion.

## Conserved-source decomposition

Align Q_sp along z, write Q=(W,0,0,W*rho), epsilon=1-rho^2, Q^2=W^2*epsilon.
For a conserved symmetric upper root R, conservation fixes its time
entries by R_0i=rho R_zi and R_00=rho^2 R_zz. Direct trace reversal,
propagation and temporal projection yield

Hxx=-(Rxx-Ryy)/(2Q^2)+Rzz/(2W^2),
Hyy= +(Rxx-Ryy)/(2Q^2)+Rzz/(2W^2),
Hxy=-Rxy/Q^2, Hxz=-Rxz/W^2, Hyz=-Ryz/W^2,
Hzz=(Rxx+Ryy-epsilon*Rzz)/(2W^2).

The entire symbolic matrix identity is checked independently. Only
the transverse trace-free spin-two part keeps the wave denominator.
This suggests a route through angular vanishing of the transverse
Einstein source, but such a general quantitative bound is NOT proved here.
