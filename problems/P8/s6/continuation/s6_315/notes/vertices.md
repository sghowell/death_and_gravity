# General mixed-coefficient construction from the action

Put kappa=1 temporarily and g=eta+2 sum x_i A_i. Work in the formal
power-series ring around zero, where the determinant root is the
branch equal to1. Let x_S=product_(i in S)x_i and r=|S|.
No nilpotency assumption is imposed on the physical field: extracting
square-free coefficients is only a way to label polarizations.

The logarithmic density is
 log sqrt(-g)=1/2 tr log(I+2 eta sum x_i A_i).
Thus its full nonempty mixed coefficient is
 L(S)=(-1)^(r+1)2^(r-1)/r
      *sum_perm tr(eta A1 ... eta Ar).
For D(S)=[x_S]sqrt(-g), D(empty)=1. Differentiating exp(log D) in
one distinguished label i gives the anchored partition identity
 D(S)=sum_(T subset S, i in T) L(T)D(S minus T).
It counts every set partition once, with the block containing i
selected first. In particular, no extra polarization factorial can
be suppressed or inserted.

Let I(S)=[x_S]g^-1. Multiplication by g gives
 I(empty)=eta,
 I(S)=-2 sum_(i in S) I(S-i)A_i eta.
Iterating gives the full permutation expression
 (-2)^r sum_perm eta A1 eta ... Ar eta.
The density inverse is
 K(S)=sum_(T subset S) D(T)I(S-T).

For a field A_i carrying all-outgoing momentum p_i, let C1(i) be the
undoubled linear connection
 C1(i)^r_mn=eta^(ra)(p_(i,m)A_(i,an)+p_(i,n)A_(i,am)
                                  -p_(i,a)A_(i,mn))/2.
The connection of the doubled metric is exactly
 Gamma(S)=sum_(i in S)2 I(S-i)eta C1(i).
This retains the derivative on each separately labeled field.

Write B(M,G,H)=sum_(mnrz) M^mn
 [G^r_mn H^z_rz-G^r_mz H^z_nr].
For r>=3, the stripped Einstein vertex is
 V_r=-1/2 sum_(A,B,C ordered disjoint partition S; B,C nonempty)
             B(K(A),Gamma(B),Gamma(C)).
The minus is the two Fourier derivatives. Restoring canonical fields
multiplies this by kappa^(1-r/2). For r=3 and4 the formula gives the
unchanged S304 cubic and S310 quartic expressions respectively.
Linear de Donder gauge has no higher gauge-fixing vertex.

For a scalar pair of mass-squared m2 the vertex is
 -[p_lower K(S) q_lower+m2 D(S)]/kappa^(r/2).
The heavy/light mixed and contact interactions give -g D(S) and C D(S),
with their matter factorials already differentiated. These formulas
work also for the zero-graviton potential vertices. Kinetic two-point
terms belong in the propagator and are not counted as tree vertices.

## Independent new-order calibration

A separate computation uses the literal4x4 determinant and adjugate
of eta+2tH, expands their ratio and square root, builds Christoffel
symbols by direct index contraction, and extracts the fifth
coefficient. Inclusion-exclusion over all31 nonempty field subsets
polarizes that directional fifth coefficient. It calls neither the
mixed log-density/inverse recursion nor the inherited connection/
bilinear helpers. The new EH5 coefficient is16438 and the fifth
determinant coefficient is-2088 for the recorded tensors.

A further literal-adjugate calculation checks all16 entries of K_4
and the scalar fourth metric vertex; the fourth determinant is620.
These exact nonzero calibrations and recovery of every earlier
implemented metric jet supplement the general coefficient proof.
