# Full cotangent lift and complete nonzero momentum constraints

All adjoints below use the full real spatial L2 pairing on the fixed
coordinate torus, including symmetric matrix entries and every matter
and vector component. Write c=det(gamma)^(1/3), g=gamma^-1. Then

DQ[h]=c[(trace(g h)/3)g-g h g],
DQ^*[B]=c[(g/3)trace(Bg)-gBg],
A h=div DQ[h].

Let V xi be the ENTIRE Lie action on gamma,W_i,M1,H. On the exact
slice, M_Q=A V is precisely the source-pinned S258 ghost operator.
For ||Q-I||<=1/8 its full nonsymmetric weak form is coercive by3/4
on mean-zero vectors, with inverse L2 norm<=4L^2/3. Constants lie
in its kernel. The range is mean-zero because A is a divergence.
These statements and the derivative-gap local coordinate construction
are retained from S258; no same-regularity group is assumed.

For any ambient tangent delta q, the mean-zero descriptor
xi=M_Q^-1 A delta q makes
H delta q=(I-V M_Q^-1 A)delta q tangent to the slice. Its slice
coordinates are delta v=trace(g delta gamma_h)/6,
delta tau=P_TT DQ[delta gamma_h], and the full horizontal variations
of W,M1,H. The scalar shape derivative in the previous note proves
that these coordinates recover every slice tangent.

For reduced momenta Pi_v,Pi_tau,Pi_W,Pi_M,Pi_H, choose the ambient
extension alpha0 with metric momentum

pi0=(Pi_v/6)g+DQ^*(P_TT Pi_tau),

and unchanged matter/vector momenta. On every slice tangent this
pairs to the complete reduced canonical one-form. Let D0=V^*alpha0.
This is the entire original momentum generator, not just its metric
part; its vector Gauss term is compulsory.

The unique ambient covector extension annihilating the mean-zero
gauge complement is

alpha=alpha0-A^*(M_Q^*)^-1 Pmean0 D0.

Indeed V0^*A^*=M_Q^*, so every mean-zero momentum constraint vanishes.
Any other extension with the same slice pullback differs by A^*lambda;
the invertible adjoint then forces lambda=0. The correction vanishes
on slice tangents. Hence the full canonical one-form pulls back to

integral Pi_v delta v+Pi_tau:delta tau
 +Pi_W dot delta W+Pi_M delta M1+Pi_H delta H.

Taking its exterior derivative gives the canonical cotangent form.
This uses an explicit bounded lift, not an unproved strong Darboux
theorem on an arbitrary infinite-dimensional symplectic space.
For fixed smooth coefficients, A:L2->H^-1,
M_Q^-1:H^-1_mean0->H1_mean0 and V:H1->L2 are bounded; the dual
formula also acts on L2 covectors. Smooth fields give the corresponding
regularity by elliptic bootstrapping.

The equivalent direct shape momentum is obtained by
Ashape=-a^2 exp(2v)Q^-1 pi Q^-1:

Pi_tau=P_TT[Ashape-cof(Q)(K_Q^-1)^*B^*Ashape],
Pi_v=2 pi:gamma.

Both inverse/adjoint contacts are nonlinear and generally nonzero.
A flat TT projection of the original pi is insufficient.

For the ordered constraints (chi,Dmean0), the full second-class block
may be [0,M;-M^*,C] with nonzero C. Its inverse is
[(M^*)^-1 C M^-1,-(M^*)^-1;M^-1,0].
The finite exact and random diagnostics keep M nonsymmetric and C
nonzero. They test signs and adjoints, not a discrete diffeomorphism
algebra. Because M annihilates constants, the cotangent correction
does not alter the three residual translation charges.
