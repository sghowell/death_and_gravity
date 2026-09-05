# S5.4 — Applicability audit and the additional UV contract

Date: 2026-09-04. The audit is complete; **positivity is not applied**.
No UV completion is claimed or excluded. Further UV classification needs
a specified admissibility class and matching observable, not a number
inferred from the local tree bound or frozen symbols.

## Exact facts about the selected witness

For a constant clock on a flat metric all clock derivative operators vanish.
The metric equation necessarily requires `F(phi0,0)=0`. For the explicit
polynomial-in-X continuation of the selected D functions,

```
F(u,0)=-(10u^14+62u^12+162u^10+263u^8
         +288u^6+186u^4+34u^2+19)/(1+u^2)^8 < 0.
```

Every numerator coefficient is positive; a separate rational Sturm replay
records the strict sign. Thus **this specific continuation has no finite
constant-clock Minkowski vacuum**. The scalar equation need not be tested
after this necessary condition fails. This does not exclude all UV completions
containing the trajectory. The original contract required a smooth clock tube,
not analytic interpolation to a specified X=0 vacuum. Changing that off-tube
continuation would be a new model/acceptance decision.

Physical geodesic completeness also does not supply an in/out construction.
Here `integral dt/a=pi*tau` is finite, and for fixed nonzero comoving momentum
`k_physical/E_ref=k_comoving/sqrt(1+u^2)` tends to zero at both ends. The
local hard band is not an asymptotic band for one fixed mode. These facts do
not prove nonexistence of alternative cosmological observables or states.

The tensor Planck coefficient is M^2 and tensor exchange is explicitly
present. Gravity has not been decoupled. The hard-transfer criterion excludes
the limits needed for a naive forward-amplitude positivity derivative.

## External theorem hypotheses are not automatic

[Grall and Melville](https://arxiv.org/pdf/2102.05683), section I, use
unbroken rotation/translation kinematics, crossing, unitarity, analyticity
and polynomial boundedness for boost-breaking bounds. Their inflationary
application specifies a subhorizon/decoupling regime. Our finite-time blocks
are not yet matched to that stationary dispersion problem. Principal health
and luminality do not supply the missing analyticity or matching assumptions.

[Tokuda, Aoki and Hirano](https://arxiv.org/pdf/2007.15009), sections 3.2--3.3,
treat the gravitational forward pole together with a high-energy Regge
contribution; strict positivity need not survive an unqualified pole
subtraction. The P8 contract supplies neither that prescription nor an
equivalent gravitational dispersion assumption. Deleting tensor exchange
would change the tested theory.

These are primary-source applicability checks, not imported UV priors. The
code checks the witness-specific algebra, not the external dispersion proofs.

## Required additional specification

A conditional UV classification must first choose:

- the vacuum/stationary matching problem and whether the current off-tube
  continuation may change;
- the symmetry, analyticity, crossing and high-energy boundedness assumptions
  of the intended UV class;
- the graviton-pole/IR prescription and gravitational high-energy assumptions
  of the positivity theorem being tested.

These choices can change which models qualify. Importing a flat-space or
P9 sign now would silently strengthen the research assumptions. Any later
verdict must be conditional on the chosen class; this report keeps
`NOT_APPLIED`. M1 interaction control and P8(a) are separate possible
continuations, not results of this D-only audit.
