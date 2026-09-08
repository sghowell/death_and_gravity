# S6.38: specified kinetic promotions of the exact affine lift

Status: scoped exact kinetic-action screens. Original P8 remains OPEN.
This is a new child of the frozen S6.37 action, not an edit to it.

## Fixed data and conventions

The target remains the entire original CD/M1 action with its original
physical metric and free canonical chi, on the same clock tube. Its
exact S6.37 auxiliary lift is the parent input. The parent retains all
64 connection components, with four projective gauge directions.
TRACE_FORM and QUOTIENT_VECTOR preserve that quotient. HOMOTHETIC
instead promotes the trace to a Maxwell field and leaves only its
gradient gauge freedom. Neither torsionlessness nor metric
compatibility is imposed in any of the three actions.

Use the parent source signature (-+++), x=-X_repo, derivative index
last in Gamma^a_bc, and kappa=Gamma-LC(g). Dimensionless M=tau=1
formulas are restored with M^2>0 and tau>0. Physical constant kinetic
couplings zeta_phys are dimensionless. After factoring M^2*tau^2 from
the dimensionless action, the normalized coefficient is
zeta=zeta_phys/(M^2*tau^2), with V_normalized=tau*V_physical.
The original background has phi=u, a=(1+u^2)^2
and chi_dot=1/[10(1+u^2)^6]. Physical momenta are measured in tau^-1.

## Three named actions, not a general kinetic classification

1. HOMOTHETIC adds -zeta_H H_mu_nu H^mu_nu/64, with
   H_mu_nu=R^a_a_mu_nu and constant zeta_H>0.
2. TRACE_FORM adds -zeta_Y tr(Rhat_mu_nu Rhat^mu_nu)/4, with
   Rhat=R-I tr(R)/4 and real constant zeta_Y!=0. The contraction is
   the matrix trace, not a positive norm of matrix entries.
3. QUOTIENT_VECTOR adds -zeta_V F(V)_mu_nu F(V)^mu_nu/4, with
   V_mu=kappa^a_mu_a-kappa^a_a_mu/4 and constant zeta_V>0.
   This is a projectively invariant distortion-trace curl. It contains
   a derivative of the torsion trace and is not a curvature-only term.

The result for the first action is exact spectator reduction.
For the second action the test gives incompatible kinetic signs
of the unconstrained homogeneous spatial STF rank-3 and axial sectors.
For the third, first derive the constrained four-vector Schur action by
eliminating the other 56 quotient components, then retain the complete
metric, clock, free-matter and vector scalar constraints. An isolated
Proca mass or an auxiliary Hessian sign cannot establish health.

The selective vector action preserves the original rolling background,
but its fully reduced configuration-velocity form has two positive
pivots and one negative pivot whenever u!=0 and
q>3*J(u)*(1+u^2)^6/4, for every positive normalized zeta_V.
Here J is exactly the original coupled CD/M1 function, not a refit.
The threshold tends to 3597/3200 as u approaches the crossing; the
punctured proof never divides by Theta at u=0. Thus both TRACE_FORM
and QUOTIENT_VECTOR fail their specified everywhere-healthy parent
tests. The zero kinetic coefficient is a separate auxiliary rank chart.

## Domains and proof obligations

The vector Schur construction uses the entire closed CD tube inside
the parent's open regular domain. Its background and quadratic
linearization use the actual rolling solution, not a flat off-shell
metric probe. Prove exact background preservation and derive the full
linearized source, including coefficient variations and ADM lapse,
shift and spatial-metric variations.

For the scalar health test, q=(tau*k_physical)^2>0 and zeta_V>0.
The ordinary unitary configuration chart is restricted to Theta!=0;
the center Theta=0 requires a separate regular Hamiltonian check. A
negative physical configuration-velocity Schur pivot on a punctured
time interval suffices to reject an everywhere-healthy kinetic parent.
A center momentum-Hessian sign alone is not the claimed proof. Keep
the zeta_V=0 auxiliary limit and every failed Legendre/lapse divisor
separate instead of continuing formulas across rank changes.

## Exclusions from this checkpoint

These are results for three literal kinetic additions only, not a
no-go for all metric-affine theories, tuned curvature contractions,
additional fields, interacting matter or higher-order completions.
No finite-momentum threshold is silently identified with an unhealthy
low-energy frequency, EFT cutoff or heavy-mode gap. A small kinetic
coefficient has a singular auxiliary limit and needs a separate
frequency/domain analysis before an EFT claim.

There is no new vacuum extension, scattering or loop calculation,
Regge/IR estimate, full V/G/B closure, UV completion, original P8
closure, or change to the completed scoped photon and linear results.
All finite algebra is checked symbolically with explicit domains and
written constraint/representation arguments, not Lean formalization.
