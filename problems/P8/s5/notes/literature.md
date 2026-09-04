# S5 primary-source method checks

Checked 2026-09-04. These are external results, not certificates for our
model. Formula adoption requires matching conventions and independent tests.

- [An et al., 2606.03302v1](https://arxiv.org/html/2606.03302v1), June 2026,
  directly follows the A25 benchmark. Sections 3.1--3.3 separate a gamma
  neighbourhood from the usual curvature-variable WKB region and use cubic
  action ratios and a one-vertex three-point function. Its action (3)
  includes Q(phi) and W(phi,chi), so it is not our free canonical M1.
  Its signature is (-+++), unlike this project's (+---). It is useful as
  a crossing-aware diagnostic comparison, not a supplied quartic amplitude,
  all-orders proof, or verification of our D-only model. We do not import
  its numerical cutoff or infer a physical cutoff from a nonlinear gauge
  transformation becoming large.
- [Ageeva and Petrov, 2206.10646](https://arxiv.org/pdf/2206.10646), sections
  3.2--3.3, give an explicit example in which apparently dangerous cubic
  contributions cancel in the tree-level amplitude. This motivates keeping
  off-shell coefficients separate from a physical strong-coupling verdict.
- [Ageeva, Petrov and Rubakov, 2009.05071](https://arxiv.org/pdf/2009.05071),
  section 3 and appendices, extend a Horndeski Genesis power-counting analysis
  to arbitrary order, including the constraint expansion. Their model and
  asymptotic hypotheses are different. Cubic/quartic results here will not
  silently acquire that paper's all-orders scope.
- [Ageeva, Petrov and Rubakov, 2207.04071](https://arxiv.org/pdf/2207.04071),
  sections 2.6.1--2.6.4, distinguish partial-wave normalization, sound speeds
  and scalar/tensor channels. Our luminal witnesses do not permit assuming
  their small-sound-speed channel hierarchy. A future scattering test must
  specify its own normalization and massless-exchange treatment.

UV applicability remains separate. Relevant starting points are
[Grall and Melville, 2102.05683](https://arxiv.org/abs/2102.05683) for
boost-breaking bounds and [Tokuda, Aoki and Hirano, 2007.15009](https://arxiv.org/abs/2007.15009)
for gravitational positivity assumptions. Their titles are not a substitute
for checking dispersion assumptions against a particular time-dependent
background and domain. Neither bound is applied in S5.1.
