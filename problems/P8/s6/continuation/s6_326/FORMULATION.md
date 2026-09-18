# Uniform radiative logarithmic soft coefficient

Keep the original action and parameters, unit scalar mass, kappa=10^800,
and the exact S300 recoil. Incoming energy5/4<=E<=2 and pair-rest-frame
direction u are fixed. Radiation is any finite positive collection
w_j(1,n_j) of total R<=1/8. For an additional soft direction q=(1,n),
use a unit-Frobenius, possibly complex, spatial TT polarization A.
Initially exclude n_j=n; no unique value at those atoms is asserted.

In all-outgoing conventions, set D_i=k_i.q, P_ij=A(k_i,k_j),
B_ij=P_ii D_j/D_i-P_ij, S_ij=B_ij+B_ji and Sbar=sum_i P_ii/D_i.
The explicit coefficient C=F+iG is defined by radiative.components and
notes/conservation.md; every sum over pairs there is unordered.

The attributed one-loop theorem has the form
 H_(N+1),one-loop =
 S0*H_N,one-loop + C_sigma*ln(1/w)*H_N,tree/kappa^(3/2) + O(1)
at fixed other-leg kinematics, with detector-resolution terms retained
separately. It does NOT supply a uniform O(1) remainder over our full
domain. Only the displayed logarithmic coefficient is bounded here.

For all finite radiative states on the specified domain,
 |C_sigma|<320,
 |C_sigma-C_Born|<25000R,
where Born has the same E,u. These are uniform angular approach and
total-energy bounds, not multiplicity-dependent sums of w_j|ln w_j|.
Relative-collinear splitting leaves the coefficient unchanged.
Complex F,G are not separately the real and imaginary parts of C.

For the named Born-coefficient addition
 L_log(w)=C_Born*ln(1/w)/kappa^(3/2)
to the complete47-tree one-real amplitude F1=M5/A0, where A0=Am+AG,
the full signed difference |F1+L_log|^2-|F1|^2 has total variation
below, for0<x<=1/8,
 4e-1187*x*(1+ln(1/x))
 +3000*x^2*(ln(1/x)^2+ln(1/x)+1/2)/kappa^3.
The quadratic term is displayed but is not called one-loop order.
This is not the complete loop-corrected probability or virtual matching.

No unknown hard or evanescent coordinate is selected. The full hard loop,
uniform finite remainder, common-regulator assembly, all-N hard measure,
interacting state, unitarity, absolute Regge and common-parent bounce
remain research obligations. Original V/G/B/P8 remain OPEN.
