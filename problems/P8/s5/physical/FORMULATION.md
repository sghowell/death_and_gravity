# S5.2 acceptance contract — D-only M0 physical reduction

Opened 2026-09-04. This extends S5.1 commit
`2eec52ee3c2ded8136ee445e65ef8a89d52f8eb6` without changing its witness or
the S0--S4 linear classification. The result is a perturbative action
construction through fourth field order, not a nonlinear PDE theorem.

## Claimed domain and deliverable

Use the invertible auxiliary metric and exact lapse-reduced Hamiltonian of
[S5.1](../FORMULATION.md). On a local, flat-FLRW perturbative patch impose
`hat_h=a^2[(1+2 zeta)I+gamma_TT]`. Solve the three spatial momentum
constraints in the full York scalar/vector decomposition, not a scalar
shift truncation. Pull back the canonical one-form including time-dependent
boundaries and compute cubic/quartic kernels for the one physical scalar
and two tensor polarizations.

The implementation accepts arbitrary **rational** 3-momenta and real rational
TT polarizations with total momentum zero and no zero nonempty proper subset
sum. This is an open nonexceptional domain, not a collinear ansatz. The
written rational-operator construction extends to real momenta on that
domain; the code does not claim symbolic multivariate-momentum enumeration.
All integrations by parts assume periodic modes or sufficient spatial decay.
No homogeneous source is silently inverted or set to zero.

Phase kernels have compact time `x=u/sqrt(1+u^2)` in `[-1,1]`; endpoints
mean limits in local background units, not finite cosmic time. Velocity
kernels use either the unitary or gamma canonical chart wherever its exact
kinetic coefficient is finite and positive. Sufficient covering domains are
`x^2>=1/65` for unitary and `x^2<=1/17`, `q>=68` for gamma. In a quartic
velocity kernel that lower bound applies also to every internal pair sum.
It is not a global restriction on which cosmological modes exist.

Canonically normalized velocity kernels include the physical derivative of
the normalization, holding the local unit `ell=tau*sqrt(1+u^2)` constant in
each patch. The quartic Legendre term is retained separately from `-H4`.
Do not identify that term with a cubic-exchange propagator diagram.

## Verification and limits

Acceptance combines written finite-order algebraic lemmas with exact
Gaussian-rational Fourier jets, independent curvature identities, a generic
York inverse identity, the old full quadratic Hamiltonian bridge, a generic
Legendre-series identity, symmetry tests and fail-closed negative controls.
This is adversarial self-review, not independent human review or Lean proof.

The report gives explicit all-time bounds only for its listed phase kernels
at fixed external spatial momenta. It does not bound physical amplitudes
uniformly over momenta, soft channels or energy. S5.3 must assemble exchange
and contact terms, impose a stated frequency/IR/adiabatic domain, handle or
exclude exceptional modes explicitly, and compare a justified perturbative
window with nonzero curvature and coefficient-variation scales.

M1, all-orders/loops, nonlinear stability, boost-breaking positivity and UV
completion are not part of this deliverable. Failure of a chart or a frozen
witness's test is not an exclusion of its entire classification row.
