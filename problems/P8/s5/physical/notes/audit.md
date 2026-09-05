# S5.2 audit and open S5.3 gate

Date: 2026-09-04. This is local adversarial self-review, not external review.

## Checks and negative controls

- Exact labelled-jet multiplication, factorials, derivative product rule,
  analytic inverse and square root; rational compact-time coefficient field.
- Full coordinate curvature checked against the independent conformal
  identity and a pure coordinate pullback of a flat metric through order 4.
- Generic symbolic three-dimensional York inverse and determinant; exact
  residuals of all three momentum constraints through degree 3 for mixed
  non-collinear inputs. W2/W3 and scalar-generated transverse W2 are not zero.
- Full old scalar quadratic Hamiltonian in both charts, not only its
  principal high-frequency limit; tensor kinetic/gradient terms, absence of
  mixed quadratic terms and vanishing homogeneous scalar tadpoles.
- Removing the time-dependent gamma generator changes the checked
  scalar coordinate-momentum coefficient. Removing normalization derivatives
  changes a checked cubic coefficient; direct uncompactified differentiation
  checks the local-unit chain rule independently.
- Exact unit kinetic normalization for scalar/tensor modes in both charts.
- External-leg permutation and common spatial-rotation covariance. Generic
  rational TT bases retain explicit norms and canonical dual momenta.
- Generic quartic Legendre identity plus independent labelled combinatorial
  factor. Nonzero tensor Legendre channels with four scalar external legs
  are retained, rather than treating the scalar sector as nonlinearly closed.
- Exceptional proper-subset momenta, singular external velocity charts and
  nonpositive internal gamma kinetic denominators fail explicitly. An actual
  nonzero homogeneous constraint source is never discarded.
- Exact report replay rejects altered contact coefficients, source metadata
  or scope limitations. Published S5.1 and linear source hashes are checked.

The report's numerical momenta are exact rationals, not floating-point
samples. Its time-dependent phase kernels are symbolic over the full compact
time interval. Its normalized velocity examples at the bounce and an exterior
point are regression fixtures, not an all-time amplitude certificate.

## Pitfalls resolved during implementation

1. The nonlinear metric determinant, inverse metric and full tensor curvature
   are needed before substituting momentum invariants.
2. Linear momentum constraints are insufficient for quartic interactions;
   the nonlinear vector part is essential even for scalar sources.
3. An exact canonical scalar/tensor gauge avoids silently omitting a nonlinear
   correction to the symplectic form, but the explicit background and gamma
   time-dependent boundaries still have to be kept.
4. Rational TT matrices generally are not unit normalized. Tensor momentum
   source legs and inverse Hessians must use the dual normalization.
5. `L4=-H4(P0)` alone is wrong. The scalar and tensor Legendre terms can be
   sizable and can cancel parts of it; none alone is a cutoff estimator.
6. A varying compact unit is not a physical derivative. The weight correction
   in `D_w` and `D log Z` must be included before calling a mode canonical.
7. A zero pair sum is common in centre-of-mass kinematics, not an exotic
   numerical corner. The present API rejects it; no assertion about its
   limiting amplitude follows.

## Next: S5.3, rather than declaring S5 complete

The reduced interactions now supply the input to a physical control test.
Construct a local frequency-space kernel including time-dependent canonical
normalization and the background measure, with an explicit adiabatic error
budget. Assemble scalar and both tensor cubic-exchange channels together
with quartic contacts. Check chart-equivalent observables, not equality of
off-shell individual vertices.

Start on a nonexceptional frequency/angular domain with all internal chart
conditions stated. Treat genuine massless propagator poles separately from
constraint inverses. Then either give a controlled extension/IR prescription
for exceptional channels or state their exclusion in the final result.
Do not call the gamma-chart finite-q denominator a physical strong-coupling
scale. Fixed-momentum boundedness and the free parameter M tau alone do not
give a quantitative perturbative energy window.

Use the nonvanishing `E_ref=1/(tau sqrt(1+u^2))` and the earlier curvature/
coefficient-variation bounds, not |H|, to compare scales. All-time coverage
must include both infinite tails and chart overlaps. The appropriate initial
target is explicitly tree-level and through quartic order. All-orders,
loops, nonlinear stability, M1 and UV assumptions require later work.
