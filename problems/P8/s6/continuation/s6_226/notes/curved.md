# Actual curvature-coordinate composition and normalization

Use the actual S225 coordinate operatorBc and its formal conformal-density adjoint. Their causal inverse kernels obey

||Y(t,s)||<=(5/2)(t-s)exp(20(t-s)),
||Z(t,s)||<=(5/2)(t-s)exp(108(t-s)).

These distinct kernels and their order are retained. After conjugation by the exponential weight, the scalar majorants are(5/2)tau exp(-(sigma-20)tau) and(5/2)tau exp(-(sigma-108)tau). Since sigma_M>=2m>108, integration of tau exp(-b tau) gives, on the unchanged conformal slabT<=1,

||Y||<= (5/2)/(sigma_M-20)^2,
||Z||<= (5/2)/(sigma_M-108)^2.

These are weighted L2_tH^r bounds uniformly inq, every realr. They use only the actual finite-slab coefficient enclosures; an arbitrary exterior continuation of the metric is unnecessary.

For the specified reference A_V=Bc*(Fdiag+V)Bc, the inverse is

E_V=Y K_V Z.

Its norm is bounded by

125/[4(32+8M)(sigma_M-20)^2(sigma_M-108)^2].

The two products cancel adjacent inverse factors in their actual order. Composition is performed on the same causal distribution graph as S225, with the middle graph now given by notes/inverse.md. The bounded integral operators extend to ordinary source classes by density and the finite-order forward maps. No density or invariance of the FULL S222 quantum-force graph is asserted.

## Time weight and physical density must both be restored honestly

On a finite slab, ||exp(-sigma t)f||<=||f|| and ||u||<=exp(sigma T)||exp(-sigma t)u||. Hence the corresponding unweighted L2 inverse estimate costs an additional exp(sigma_M T). That number may overwhelm the small weighted norm. A small weighted norm is not small backreaction or stability.

The force-normalized curvature reference is

Qbar_V=M_[1/(64pi^2 kappa a^4)] A_V.

Its inverse is E_V M_[64pi^2 kappa a^4], with multiplication on the RIGHT, not throughE_V. The actuala^4<6 gives a weighted physical-reference bound

12000pi^2 kappa /[(32+8M)(sigma_M-20)^2(sigma_M-108)^2].

Kappa remains10^800. Removing the weight still costs exp(sigma_M T), and changing amplitude/dual-source norms to physical metric Frobenius conventions still requires the inherited embedding factors. No physical unweighted smallness is claimed.
