# Complete restricted local basis

Two curvature tensors start at two gravitons around flat space. Linear
Ricci, scalar curvature and their derivatives vanish on a real null TT
wave. Therefore only one Weyl tensor, equal here to linear Riemann, can
give a nonzero contact in the stated explicit-curvature class.

Every Weyl index must be saturated by a scalar momentum: all its metric
traces and all contractions with emitted k vanish. With at most six
derivatives total, only four derivatives remain after the curvature.
Thus fewer than six cannot suffice. At exactly six, every remaining
derivative must supply one hard scalar momentum to a curvature slot.
A derivative on Weyl instead supplies k and leaves too few hard slots,
so gives zero. There is no spare scalar dot product or mass-shell trace.
Pairs of Levi-Civita tensors in a parity-even contraction reduce to metric
contractions. Connection or curvature commutators beyond this linear
calculation start at two curvatures. This is an index-counting proof, not
an assumption imported from massless amplitudes.

Every possible contact is consequently a linear combination of
R_abcd p_i^a p_j^b p_k^c p_l^d, including all repeated labels, followed
by the full S4 symmetrization of the four identical scalar fields.
In the retained metric convention this contraction is

Q_ij,kl=epsilon(v_ij,v_kl),
v_ij=(k.p_i)p_j-(k.p_j)p_i.

Set a_i=k.p_i and H_ij=epsilon(p_i,p_j). Conservation and physical TT imply
sum_i a_i=0 and every row sum of H is zero. Without using scalar mass shells,

Q_ij,kl=a_i a_k H_jl-a_i a_l H_jk-a_j a_k H_il+a_j a_l H_ik.

Let T=sum_(i<j)Q_ij,ij. The UNNORMALIZED full 24-label Bose sum is:
0 when either antisymmetric pair repeats a label or the pairs are disjoint;
+4T for identical oriented pairs and -4T for reversed pairs;
-2T when the shared endpoint has the same orientation in both pairs,
+2T for opposite orientation.

Indeed sum_j v_ij=-a_i k, so its physical epsilon contraction vanishes;
summing the shared-edge products gives -2T. Disjoint-edge terms cancel
under reversal of one edge. Identical pairs each occur four times.
All256 words are checked exactly, as are the Bianchi identity, every
curvature trace and every k-slot in a generic null TT frame.
The coefficient matrix has rank1 even before additional physical Gram
constraints. This proves an upper bound, and the independently nonzero
original massive amplitudes prove the lower bound.

The complete operator Phi^2*C_abcd*Hess_ac(Phi)*Hess_bd(Phi)/4 has
four assignments for each of the six unordered pairs, hence contact T.
The coefficient chi therefore has exactly the S336 normalization.
This classifies only this derivative order and observable; it is not a
field redefinition or a full curved-action basis.
