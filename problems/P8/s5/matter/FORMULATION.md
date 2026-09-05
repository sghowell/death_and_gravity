# S5.5 frozen CD/M1 contract

This is an additive checkpoint. It does not change the original P8(b) question,
the linear classification, the D/M0 interaction results, or the accepted S6
conditional UV framework. The exact sources and report digests used here are
pinned in `src/p8_m1/verify.py`.

## Inputs and conventions

Use the full covariant `CD_matter` witness from
`../../certificates/witness-CD_matter.json`, including its reconstructed F,
Ia completion and backreacting free canonical matter. Signature is (+---),
`X=(partial phi)^2`, and the physical matter metric is the original g.

In dimensionless time `u=t/tau`, `d=1+u^2`, the background is
`a=d^2`, `H=4u/d`, `phi=t`, `chi_dot=1/(10d^6)`. In the unitary clock chart
`X=N^-2`. The covariant tube is `9/10<=X<=11/10`, all real u. Gravitational
normalization is `M^2`; setting `M=tau=1` is only a unit convention.

The new auxiliary spatial metric obeys
`h_ij=T^(-1/2)*hat_h_ij`, `T=1+(N^-2-1)/d^3`. Lapse, shift and chi are kept.
This is a regular field-coordinate change on the tube, not a different physical
matter frame, not a four-dimensional Einstein/disformal parent, and not a
positivity assumption.

## Nonlinear claims admitted here

1. Exact cancellation of quadratic lapse velocities and spatial lapse gradients
   in the full ADM action, with the remaining linear velocity removed by an
   explicitly specified, local canonical boundary primitive.
2. An independent full symmetric-metric-plus-matter Legendre transform, before
   imposing the spatial momentum constraints.
3. Agreement of its nonlinear lapse jets with the earlier independently derived
   background, Theta, Lambda, J and matter-mixing coefficients.
4. Stationary elimination of lapse through cubic order and the resulting
   Hamiltonian through quartic invariant degree. Curvature and momentum
   invariants are retained, not mistaken for fully reduced perturbations.
5. Strict all-time `1/10<d*J<8` and uniform absolute bounds for every recorded
   invariant coefficient in fixed local background units. These are sufficient
   majorants, not sharp coefficients or a physical cutoff estimate.

The lapse implicit-function theorem is local around each finite-time background
point. A uniformly sized nonlinear neighborhood, all-order field analyticity
bounds, a full nonlinear Dirac/PDE analysis and physical amplitudes are not
asserted by the finite-jet calculation.

## Principal-robustness claims admitted here

Require a regular two-derivative principal action with positive-definite kinetic
matrix K and the **same** positive canonical matter diagonal in K and G. The
normal form gives a necessary and sufficient local scalar cone criterion under
exactly those assumptions. Its crossing-safe invariants do not divide by Theta
or Lambda; their use still requires a covering regular chart and `T,J>0`.

For the explicit family
`delta_F=epsilon*J_star(phi)*(X-1)^2/(4*(1+phi^2)^2)` in dimensionless units,
epsilon is an arbitrary fixed finite real parameter. Other covariant functions,
Ia relations and the original free matter action are unchanged. The certified
healthy-causal family is `epsilon>=0`; the interval `-1<epsilon<0` and boundary
`epsilon=-1` are negative controls. All claims use the old background and its
proved geometry, and rederive the changed principal coefficients and tail
corrections. This is not a fresh cosmological reconstruction or loop correction.

Finite-k stability, nonlinear interactions, BKL behavior, all orders, loops,
technical naturalness, UV admissibility and completion remain open. Corrections
that change the matter principal action, introduce additional operators or break
Ia degeneracy need a new reduction. No conclusion about those corrections is
drawn merely from the present matrix theorem.
