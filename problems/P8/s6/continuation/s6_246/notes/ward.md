# Actual ordered scalar Ward reconstruction and endpoint flux

## Covariant scalar current, not a new state or prescription

Define the contravariant weight-one current density by
\(J_H(D)=\int E^{\mu\nu}h_{D,\mu\nu}\,dt\,d^3x\). In the actual scalar convention,
\[
 E^{00}=-a^3\rho/2=:A,\qquad
 E^{ij}=-aP\delta^{ij}/2=:B\delta^{ij},\qquad E^{0i}=0.
\]
The full reference mean is nonzero and satisfies
\(\rho'+3H(\rho+P)=0\). Both state-dependent and complete finite-counteraction pieces are included.

The standard scalar retarded variation retains the explicit metric variation of the observable as well as its connected stress insertion; see [Hollands–Wald, §4.3, equation (112), and Theorem 5.1](https://arxiv.org/abs/gr-qc/0404074v2). This supplies distributional retarded/contact and conservation context, not the numerical bounds below or an interacting-loop completion. Its existence statement does not select our finite scheme.

Here the scheme is the already fixed full covariant scalar dimensional prescription: vary the entire continued counteraction before the physical limit. Its finite local variation is retained, not readjusted to enforce a desired bound. A compact metric variation has its same-preparation retarded scalar variation; source diffeomorphisms are identity throughout the preparation region. The scalar mode identity in the geometry note independently verifies the complete first Cauchy-data bridge.

## Source covariance and detector conservation are different identities

For \(\xi=(\eta,\chi)\), the full density Lie derivative is
\[
 \mathcal L_\xi E=\xi\cdot\partial E-(\partial\xi)E
                 -E(\partial\xi)^t+(\operatorname{div}\xi)E.
\]
Its components are
\[
\begin{aligned}
 (\mathcal L_\xi E)^{00}&=\eta A'-A\eta'+A\operatorname{div}\chi,\\
 (\mathcal L_\xi E)^{0i}&=-B\partial_i\eta-A\chi_i',\\
 (\mathcal L_\xi E)^{ij}&=\eta B'\delta^{ij}
       -B(\partial_i\chi_j+\partial_j\chi_i)
       +B(\eta'+\operatorname{div}\chi)\delta^{ij}.
\end{aligned}
\]
Every index and the density weight matter.

The ADM chain rule for the actual retarded variation is
\[
 R_H(D,G)=\int h_D:\delta E[h_G]+E:q(D,G).
\]
This does not declare a retarded kernel to be a symmetric single-branch Hessian.

For a source gauge vector identity on the initial germ, local covariance and the same pulled-back state give \(\delta E[\mathcal L_{\xi_G}g]=\mathcal L_{\xi_G}E\). Therefore
\[
 W_s(D,\xi_G;E)=R_H(D,G_{\xi_G})
 =\int h_D:\mathcal L_{\xi_G}E+E:q(D,G_{\xi_G}).
\]
Separately differentiate conservation with a fixed detector vector. The full differentiated flux is
\[
 F^\mu=2\xi_D^\alpha
  \left[\delta E^{\mu\nu}g_{\alpha\nu}
                   +E^{\mu\nu}h_{G,\alpha\nu}\right].
\]
Its final time flux vanishes because \(\xi_D=0\) there, while its initial flux vanishes because the source and its retarded current tangent are zero there. Spatial boundary terms vanish by the compact/Schwartz domain. Thus
\[
 W_d(\xi_D,G;E)=R_H(D_{\xi_D},G)
 =\int -E:\mathcal L_{\xi_D}h_G+E:q(D_{\xi_D},G).
\]
Do not interchange these ordered terms.

## Entire reconstruction

Bilinearity now yields
\[
 R_H(D,G)=R_H(D_{\rm syn},G_{\rm syn})
       +W_s(D,\xi_G;E)+W_d(\xi_D,G_{\rm syn};E).
\]
S6.245 supplies the missing full spatial kernel, with the explicit profile removal in the next note. This is an instantiation by an actual scalar response, not a conditional bound with unknown scalar blocks.

Every upper source endpoint was retained in S6.245. The synchronous source can therefore extend through the observation endpoint. Extend it smoothly outside the slab for compact-variation arguments; causality makes that outside extension irrelevant to the pairing. No sharply varying inside-slab source cutoff is introduced.

The advanced synchronous detector may reach the left endpoint. Since source and output vanish on an initial neighborhood, cutting the detector off inside that neighborhood changes no pairing. Its norm uses no time derivative, so this cutoff adds no derivative cost. Final-neighborhood detector cutoffs and the proved bounded form then give the stated \(L^2\)-time completion.

## Independent full-action controls

The tensor-density Leibniz identity is checked on nonconstant, nondiagonal tensors. For the full covariant constant-density action, both Ward formulas agree with its literal nonlinear ADM Hessian; the detector comparison retains the entire flux divergence. Pure-shift source tests have a nonzero density-Lie contribution and a nonzero chart contribution which cancel. Dropping either fails.

A retarded detector primitive has a nonzero final value in an explicit test and leaves a nonzero differentiated flux. The advanced primitive kills it. These checks distinguish the two endpoint roles rather than treating all compact-support manipulations as interchangeable.
