# Full two-point tensors and independent normalization

S6.145 supplies two self-energy words and one vertex word per
scalar/gauge quadratic primitive. Rescale real loop momenta by m,
so S(k)=(1-i slash(k))/(1+k^2), A=1+k^2 and B=1+l^2.
Keep d symbolic and tr 1=4. Euclidean p^2 coefficients are obtained
by the contracted momentum Laplacian divided by 2d, at p=0.

The second derivative of the shifted propagator is

 T(k)=Delta_k S(k)/(2d)
     =-N_k/A^2+4k^2 N_k/(d A^3)+2i slash(k)/(d A^2).

The complete scalar numerator coefficient is the sum of

  2 tr[S_k S_l S_k T_k],

  tr[S_k T_k S_l^2]+tr[S_k^2 T_l S_l]
  -(1/d)sum_nu tr[S_k^2 gamma_nu S_k S_l gamma_nu S_l^2].

The first line includes both marks. The last line's cross term is
essential. Both loop orientations and all word weights are inherited
from the determinant catalog, not added again.

For gauge exchange, insert gamma_mu caps at both ends of the exchanged
boson and sum mu. The self term becomes
2 tr[S_k gamma_mu S_l gamma_mu S_k T_k].
The vertex terms and their cross derivative keep the same order
displayed by tensors.py. The two gauge vertices supply the overall
minus sign relative to the scalar numerator. The group factors
are N Y^2 for scalar and N Y a C_F for gauge.

The trace algorithm expands each numerator into scalar and slash
parts. The even trace is the signed pairing sum of scalar products.
For repeated dummy indices, the contraction graph joins its endpoints;
a closed dummy loop gives d. This retains the d dependence of all
gamma contractions, including the cross term with two summed indices.
No four-dimensional numerator is inserted before pole subtraction.

As independent checks, explicit 4x4 matrix Laplacians in dimensions
three, four and five reproduce each self and vertex tensor. The
momenta in these checks span only the relevant nonchiral subspace.
The raw p=0 result also equals the fixed-mu second mass derivative
of the independently contracted vacuum graphs:

 V_scalar/(NY)=(4m^2-b)I(m^2,m^2,b)-2A_m A_b+A_m^2,
 V_gauge/(a N C_F)=-4m^2 I(m^2,m^2,0)+(d-2)A_m^2.

Here the vacuum expression is differentiated only for active flavors.
Its all-flavor constant is not being evaluated or claimed as a
completed vacuum row. These identities independently test the
relative self/vertex weights and gauge sign.
