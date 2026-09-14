# Full retained nonlinear spatial gauge: exact scope

## Parent, fields and clock

Pin the S6.257 raw report
ae63ae3c213a5a5cd8b02174ec9eacfd4d19b3d0fb103c05aa221649ead4e05f.
All earlier source functions, their full lapse dependence, fixed profiles,
three vacuum constants, full heavy source, canonical boundaries and
quantum preparations remain unchanged. This is the same retained local
action, not a separately modified candidate. The clock u is held fixed;
no GR constant-mean-curvature time gauge is imposed.

The fields before auxiliary primary reduction are gamma_ij, W_i,
M1, H, N and T, with their canonical momenta. The full spatial generator
includes both auxiliary scalar primaries. On their surface it agrees
with the parent generator. Its equivariance survives the local regular
second-class branch, but no finite anomaly-free regulator is asserted.
The shift multiplier and its primaries for time-dependent descriptors
are separate from the displayed fixed-clock field action.

## Statements

1. The complete canonical Lie pairing gives the full spatial momentum
   generator, including the vector Gauss term, both matter fields,
   auxiliary primaries and an explicit boundary flux. Its continuum
   classical algebra is the spatial Lie algebra.
2. Q = det(gamma)^(1/3) gamma^-1 is a symmetric contravariant density of
   weight 2/3. Under h_phys = C gamma, C > 0, the entire density is
   unchanged, including spatially varying C. For gamma =
   a_hat^2 exp(2v) exp(t), trace t = 0, Q = exp(-t).
3. The full operator M xi = delta_xi(div Q), differentiated before
   imposing chi = 0, equals

       xi.grad chi - chi.grad xi + (2/3)chi div xi
       - Q^jk partial_j partial_k xi^i
       - (1/3)Q^ij partial_j div xi.

   Its complete nonzero-momentum symbol has eigenvalues r,r,4r/3,
   r = k^T Q k > 0, and symmetrizer Q^-1. At Q = I the inverse is
   [I-kk^T/(4k^2)]/k^2, only for k != 0.
4. On a finite flat torus there is a local coordinate slice near the
   homogeneous reference after fixing mean coordinate displacement.
   Translations remain a residual group. On an exact slice with
   ||Q-I||op <= 1/8 the whole weak form is coercive with margin 3/4.
   The mean-zero L2 inverse norm is at most 4L^2/3 for period 2 pi L.
5. Full finite ghost vertices retain both inverse products and second
   operator variations. For an entire gauge-satisfying exponential
   shape family, a 24-dimensional relative determinant has first
   variation zero and second variation 3. Two distinct TT waves have
   a nonzero mixed gauge residual, removed by a displayed coordinate
   correction whose scalar-volume coefficient is 1/32.

The finite determinant and Fourier projections are diagnostics for the
displayed spatial operator, not a spacetime-regulated quantum theory.
No deletion of residual translations, change of fixed state, global
gauge, uniform infrared gap, BRST/covariant quantum measure, physical
counterterm selection, interacting mean or physical cutoff is inferred.

## Evidence

- [Spatial generator, boundary and classical algebra](notes/spatial.md).
- [Whole nonlinear density gauge and principal operator](notes/gauge.md).
- [Local slice, weak coercivity and infrared scope](notes/local.md).
- [Entire finite ghost determinant vertices](notes/ghost.md).
- [Reference TT map and nonlinear gauge-restoration contact](notes/reference.md).
- [Primary-source review and original-problem boundaries](notes/scope.md).
- [Immutable source and read-only validation contract](notes/validation.md).
