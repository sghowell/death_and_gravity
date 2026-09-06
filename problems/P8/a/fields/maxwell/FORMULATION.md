# P8-A.16 — free Maxwell conformal QSEI and conditional focusing dictionary

This is a NEW field-model branch, not a replacement of the scalar field,
state or renormalization prescription in A.1–A.15. It proves a photon-field
inequality on a stated background class and a conditional geometric
dictionary. It does not construct an Einstein–Maxwell semiclassical solution
or assert a cosmological focusing example.

## Geometry, physical field and state domain

Let M=I_eta times R^3 have the smooth metric

    ds^2=a(eta)^2(deta^2-dx^2), a>0,
    dt=a deta, H=a_t/a.

I_eta is any open time interval. There are no conductors, spatial boundaries,
charged matter or interacting-QED dynamics in this theorem. Curves are the
comoving timelike geodesics. The field is one free real Maxwell field with
physical energy density (electric^2+magnetic^2)/2 and the corresponding usual
hbar normalization. Admit every physical Maxwell Hadamard target state on
this strip; no quasifree, homogeneity, zero-mean, SEE or global-state-extension
condition is imposed on the target.

The reference is the conformal transport of the restricted Minkowski vacuum
of this same field. This is an actual positive Hadamard state, not a scalar
reference, instantaneous vacuum or approximate mode solution. The local
positive-type proof does not assume a target state on the strip extends to
all of Minkowski spacetime. A.12 is pinned only for general proper-clock and
Fourier-proof structure; its scalar field coefficient and state are NOT used.

## Explicit new prescription family

Use the FK (+---) convention R=6(Hdot+2H^2), R_00=3(Hdot+H^2),
G_00=-3H^2. The new constant beta_M is defined by

    T_conf = hbar/(2880 pi^2) [62 H3_FK + beta_M I_FK],
    I_FK=(1/sqrt|g|) delta integral sqrt|g| R^2 / delta g^{ab}.

The precise conserved tensor and the source-to-FK sign map are proved in
notes/proof.md. beta_M is symbolic, real and constant. The optional named
zero-type-D specialization is beta_M=0; neither that specialization nor any
other value is inherited from scalar gamma=0. Numerical bound functions
require beta_m as an explicit argument.

Cosmological and Einstein/Newton renormalizations are fixed independently
and displayed in the SEE dictionary, not absorbed unnoticed into beta_M.
Other explicit gravitational curvature terms require separate source
matching. Conformal flatness makes a finite Weyl-squared variation vanish;
the finite curvature-squared ambiguity in this background class is encoded
by the stated I_FK coefficient. Minkowski-vacuum stress is normalized to zero.

For beta_M=0 the physical conformal-vacuum components are

    rho_conf=31 hbar H^4/(480 pi^2),
    E_conf=-31 hbar (H^4+2H^2 Hdot)/(480 pi^2).

The generic finite-beta contribution is retained in stress.reference_jets;
it requires up to Hthird for E_conf. It is not permissible to equate energy
density with EED when the trace anomaly is nonzero.

## Exact all-sampler statement

For real smooth compact proper-time samplers f inside the strip, and every
admitted Hadamard target state omega,

    integral f^2 E_omega dt
      >= -hbar/(8 pi^2) integral |L_H f|^2 dt + integral f^2 E_conf dt,
    L_H f = f''-2H f'+(3H^2/4-3Hdot/2)f.

This is an absolute QSEI in the fixed beta_M prescription, with explicitly
computed local reference term. The coefficient is established, not claimed
optimal. It extends to real H2_0 samplers on any compact proper subinterval
(both zero boundary traces). H2_0 does not mean arbitrary endpoint data.

Given absolute proper-time jet caps H_j for j=0,1,2,3 and an enclosing
proper duration D, the derived nonnegative curvature-envelope form is

    integral f^2 E_omega dt >= -B2 ||f''||^2 - B0 ||f||^2,
    B2=hbar/(8pi^2) A_D^2,
    A_D=1+2H_0 D/3+(3H_0^2/4+3H_1/2)D^2/9,
    B0=hbar/(2880pi^2) V,
    V=186(H_0^4+2H_0^2 H_1)
       +|beta_M|(18H_3+90H_1^2+90H_0 H_2+108H_0^2 H_1).

The rational factor uses pi>3 and is not sharp. beta_M=0 needs only H_0,H_1
in this envelope; other beta values require the displayed higher jets.
Smoothness does not automatically provide any proposed numerical jet cap.

The radiation metric H=1/(2t), t>0, has R=0 and no finite-beta contribution.
An independently checked compact-support IBP control gives

    integral f^2 E_omega dt >= -hbar/(8pi^2)
      integral [f''^2-15 f'^2/(8t^2)+4601 f^2/(1280t^4)] dt.

This is a fixed-background control, NOT a radiation SEE solution. The
negative f'^2 term is part of an exact expanded functional; discarding it
is allowed only in the direction that weakens the lower bound.

## Conditional A.1 match and remaining gates

If the actual state/metric also obey

    G_FK + Lambda g_FK = -kappa(T_Maxwell + T_other), kappa>0,

and E_other>=ell_other uniformly on the required windows, then

    R_UU=Lambda-kappa(E_Maxwell+E_other),
    Q2=kappa B2,
    Q0=kappa B0 + max(Lambda,0) + kappa max(-ell_other,0)

are sufficient nonnegative geometric constants. This coarsening drops only
beneficial source signs. No SEC assumption is inferred for unspecified
additional matter. Extra gravitational tensors may be moved to the source
only with their explicit signs and a proved EED bound.

Use these constants in A.1 only after verifying ALL of its geometric
hypotheses: the extended normal domain [-tau/2,tau], every required local
window, initial R_UU<=rho0<=0 on [0,tau/4], and the sufficient initial K
threshold. The all-comoving photon result is not a bound for arbitrary
boosted normals to another hypersurface. Frozen A.1 numerical optimizer
values are not automatically valid at newly chosen Q2,Q0.

The scaling Q2/tau^2,Q0*tau^2=O(kappa hbar/tau^2) under fixed dimensionless
curvature envelopes is only an algebraic separation of scales. For example,
Hmax*tau<=1 gives |K|tau<=3, whereas A.1's zeta=0 frozen-tail functional has
gradient contribution at least 4. Thus those naive small-Q ratios alone
cannot furnish a focusing example. This negative control is certified.

No interacting-field theorem, exact new SEE existence, global causal
completion, maximal-extension statement, initial-pointwise-premise removal,
"only QEIs" singularity theorem, optimal coefficient or P8 completion is
claimed. The certificate checks algebra and rational envelopes, with the
distributional and geometric implications supplied by written proofs.
