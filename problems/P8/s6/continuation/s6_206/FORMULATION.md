# Actual leading endpoint conversion and its boundary

Retain the same massive CD Proca sector, m = 1000, kappa = 10^800, the unit compact time slab, unchanged prepared state and fixed finite covariant prescription. The complete actual-state/unit-W8 correction remains. Work with real symmetric tracefree spatial tensors D and G at a common time t; complex Fourier tensors follow by sesquilinear extension.

Let k = r n, |n| = 1, and l = -k+P. The j = 0 upper endpoint contributes the actual current density

E0 = Re sum_(alpha,beta) A_D#(t,k,l) A_G(t,k,l) / (W_k+W_l),

summed over all nine physical polarization pairs. For fixed P,

E0 = r h(D,G,n,t) + O(1),
h = [4tr(DG)-4 n.(DG+GD)n+3(n.D.n)(n.G.n)]/(8a(t)),

uniformly in angle on the compact CD slab. This is the fixed-mass ultraviolet limit of massive Proca, not the two-polarization Maxwell theory.

Keep the original pair mask chi_pair = 1_(|k|<=K,|-k+P|<=K). Define chi_single = 1_(|k|<=K) only as an explicit comparison mask for calculating a regulator-conversion difference. It does not replace chi_pair in the current or establish a renormalized comparator.

For fixed P != 0, let phat = P/|P|. The first-five-endpoint conversion has leading term

sum_(j=0)^4 (B_j,pair,K - B_j,single,K)
 = -K^3 |P| A(D,G,phat,t) + o(K^3),

A = [18tr(DG)-12(D phat).(G phat)-(phat.D.phat)(phat.G.phat)]
    /(512 pi^2 a(t)).

The shell calculation retains the finite positive-angle grazing strip. The other four endpoints have lower shell order, at most O(K^2), for fixed source time jets. The contact has the identical original one-mode band on both sides of this explicitly defined comparison and contributes zero to this regulator difference.

The angular bracket has eigenvalues 18, 18, 12, 12 and 28/3 in the Frobenius-orthonormal tensor/tensor/vector/vector/scalar basis relative to phat. Thus the conversion has a nonzero |P| cusp, with negative sign for equal nonzero real tensors. At P = 0 the difference is exactly zero.

Both metric factors multiply the coefficient by 4/kappa. This does not turn K^3 growth into a uniform small-response estimate at any finite kappa.

No regulator is changed and no subtraction is performed by this checkpoint. Subleading artifacts, the actual finite/divergent contact and candidate-cell coefficients, and their original fixed covariant matching remain required. No full physical current divergence, model exclusion, matched response, inverse, finite-coupling background, stability, cutoff or original V/G/B closure is established.
