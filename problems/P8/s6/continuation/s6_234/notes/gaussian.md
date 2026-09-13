# Exact regulated integration, full Hessian and loop generator

Choose a finite Euclidean lattice with its positive integration weights and symmetric kinetic operators K_phi=-Delta+1 and K_H=-Delta+M_H². Products such as phi² are pointwise lattice products. The same quadratic-form statements apply to continuum test functions with phi² in L2; this does not construct a continuum interacting measure.

Complete the entire heavy square:

S[phi,H]=S_eff[phi]
 +1/2<H-(g/2)G_H phi²,K_H[H-(g/2)G_H phi²]>.

The finite Gaussian change of variables is exact and its determinant has no phi dependence. At fixed regulator and finite volume, G_H is positive and bounded above by1/M_H². Therefore

S_eff>=1/2<phi,K_phi phi>+q integral phi4,
q=-C/24-g²/(8M_H²)=g²(D-1)/[6D²(D+2)]>0.

The regulated integral is integrable. This estimate refers to the given bare classical parameters, not arbitrary regulator-dependent counterterms or a constructed limit. The H determinant would depend on a curved metric, so it cannot be dropped when computing stress.

Vary the nonlocal quartic without commuting multiplication and convolution. Its first derivative is -g² phi(G_H phi²)/2, and its Hessian is

-g² M_(G_H phi²)/2-g² M_phi G_H M_phi.

The contact Hessian adds -C M_(phi²)/2. Independently, the full two-field Hessian at stationary H has heavy block K_H, mixed block -g M_phi, and light block K_phi-C M_(phi²)/2-g² M_(G_H phi²)/2. Its ordered Schur complement is exactly K_phi+W. This proves the full operator identity; finite2/3-site matrices are independent checks, not a substitute for the variation.

The one-light-particle-irreducible effective action follows by doing the remaining light integral. At one loop it is S_eff+Tr Log(K_phi+W)/2, up to the field-independent constant and counterterms. As a FORMAL field Taylor expansion around the positive vacuum, W is quadratic in phi, so its complete degree2 and degree4 terms are

Tr(G_phi W)/2 and -Tr(G_phi W G_phi W)/4.

No large-field convergence or positivity of every background Hessian is assumed. The second term contains contact/contact, contact/exchange and exchange/exchange pieces, including heavy-reducible but light-irreducible diagrams. The light-independent heavy determinant alone does not contain them. The displayed generator still requires counterterm evaluation and momentum-dependent integration before it becomes a renormalized four-point matching calculation.

For the heavy one-point condition, the first light trace in a heavy background gives -g I_phi H/2, where I_phi is the regulated coincident light propagator. Add j1 H with j1=g I_phi/2. The full Gaussian source is J=gphi²/2-j1. Its cross term is +g j1 phi²/(2M_H²). At quadratic order this cancels the -g² I_phi/(2M_H²) self-energy kernel from the local heavy-source term in Tr(G_phi W)/2. Both terms are kept before cancellation. The ordinary quartic tadpole and mixed bubble do not cancel this way.

Analytic continuation uses the Feynman/in-out inverse. Substitution of a retarded kernel in this one-copy functional would produce a different variation and is not the identity proved here. A causal state-dependent stress or expectation value needs a separate in-in construction.
