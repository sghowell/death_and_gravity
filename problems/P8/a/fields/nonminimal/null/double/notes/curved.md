# Physical conformal transport with fixed affine normalization

Use the same real massless conformally coupled scalar, signature
(+---), curvature convention and fixed renormalization prescription
as A.20. The full classical improved stress is conformally covariant
off shell, as checked there with arbitrary jets. For differences
between Hadamard states, the state-independent local anomaly cancels:

    Delta T_ab[g]=a^-2 Delta T_ab[eta],
    w_physical=a^-2 w_flat.                              (14)

The reference DOES NOT vanish on curved space. These are difference
identities, not a prescription to discard its curvature terms.
The target's conformal transform need only exist on the Minkowski
slab corresponding to the FLRW patch. The compact proof in plane.md
uses positivity and the Hadamard condition on that slab, not a global
translation-invariant extension of the target state.

In conformal coordinates the physical plane metric is
a^2(deta^2-dz^2), its positive volume is a^2 d eta dz=a dt dz,
and K=a^-2(partial_eta+partial_z). Direct Christoffel contraction
gives g(K,K)=0 and nabla_K K=0. In proper coordinates
K=a^-1 partial_t+a^-2 partial_z. This fixes an affine normalization;
rescaling K by a constant rescales every stress contraction and
bound by its square. It is not a unit timelike vector or the vector
partial_t+a^-1 partial_z, which is in general nonaffine.

Contracting (14) with K twice gives Delta T_KK=a^-6 Delta Tflat_DD.
The physical smeared relative stress therefore becomes

    integral dvol_plane f^2 Delta T_KK
       =integral d eta dz F^2 Delta Tflat_DD, F=a^-2 f.    (15)

Since partial_eta=a partial_t, direct differentiation gives

    F_eta_eta=f_tt-3H f_t+(2H^2-2Hdot)f=L,
    F_eta_z=a^-1(f_tz-2H f_z)=a^-1 M,
    D F=a^-1(f_t+a^-1 f_z-2H f)=a^-1 G.                 (16)

Changing flat measure to physical measure divides by a^2.
Consequently the three derivative weights are a^-2 L^2,
a^-3 L M and a^-4 M^2. In the state term, w_flat=a^2 w_physical
cancels that measure change but (D F)^2 still contributes a^-2.
This is why its physical density is a^-2 w_physical G^2.
Replacing w_physical by its one-sided upper cap is justified by
that nonnegative weight. The derivative density also retains the
pointwise positive-square decomposition transported from (9).

For the ACTUAL scalar reference, A.20 supplies its rho and pressure
with the common hbar/(2880pi^2) factor and independent finite beta_S.
The null contraction equals (rho+pressure)/a^2, and native algebra
gives precisely (4). In particular the rho anomaly alone is not the
null reference and a photon anomaly coefficient is not substituted.
Cosmological and Newton renormalizations follow the fixed A.20
gravitational convention; beta_S is not set by a target state.

The curvature convention gives R_KK=2Hdot/a^2. A cosmological
constant contributes zero to a null contraction. Neither statement
by itself converts the plane average to the null index form. A
single generator uses an affine one-dimensional measure; (15)
uses a physical two-dimensional measure. Additional geometry/state
information and an explicit reduction are required before a
focusing or incompleteness conclusion. No such conclusion is
assigned here, and no proper-width cosmological coefficient is
obtained by simply freezing a(t) in (16).
