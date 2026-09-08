# S6.52: dimensional local subtraction and matching controls

Original P8 remains OPEN. Continue the S6.50/51
retained-vector calculation without changing its action or its frozen
physical subtraction prescription. Derive the local order 0/2/4
adiabatic integrals in D=3-2epsilon_DR spatial dimensions, retaining
the D-1 transverse multiplicity and dimension-dependent canonical
friction. Their epsilon dependence contributes finite local terms.

First test the zero-order finite lapse term against S6.47 and the
ordinary-Proca derivative poles against lapse variations of the
S6.49 local curvature counterterms. Determine the actual extra clock
mass-insertion terms. A successful local comparison is not by itself
full covariant finite matching, a Hadamard-state theorem or V/G/B.

## Precisely certified local statements

Direct D-dimensional action variation fixes the lapse/pressure forms
and canonical rates. The generic WKB coefficients reduce exactly to
the frozen actual-clock ones at D=3. Dimensional radial Gamma
integration retains both the D-1 transverse multiplicity and the
dimension-dependent canonical friction through the finite term.
Taking that count to four dimensions too early loses 2+4b_N in
the finite zeroth-order lapse term, in the common normalization.

The actual clock-sensitive order 0/2/4 lapse poles have an independent
covariant reconstruction, retaining the scalar-gradient part of the
constrained vector propagator. The finite zeroth-order lapse term
agrees with the already frozen S6.47 potential and its N derivative.

For the separate ordinary-Proca control a_N=b_N=0, varying the
dimensional counterterm action before taking the four-dimensional
limit supplies finite terms that restore the Ward identity. The
matched local energy terms, with respective prefactors m^4,m^2,1
and common normalization 1/(64pi^2), are

    3ell-5/2,
    H^2(6ell-10),
    (ell+2)(6H^2H'+2HH''-H'^2), ell=log(m^2/mu^2).

The independent finite local heat action gives the same values.
Component-wise radial finite parts alone fail the conservation
diagnostic. This ordinary-Proca diagnostic is not substituted for
the actual clock-dependent finite energy.

## Boundary

All statements concern local first variations and the specified
dimensional prescription. A component-wise Laurent finite part is
explicitly distinguished from a covariantly matched stress component.
The finite clock-mass derivative counterterms are not fully matched;
their continuation must preserve S6.47's finite potential and jets.
The ordinary vector component's Ward identity is not imposed on the
actual component that exchanges with the clock equation.

No all-order Hadamard state, all-time quantum control, full corrected
bounce/constraints/cones, interacting cutoff, higher-loop error,
finite-gravity Regge bound or V/G/B closure follows. See the
[written derivation](notes/proof.md). No frozen action or finite
matching coefficient is changed. The report is exact symbolic
verification with a written local proof, not proof-assistant formalization.
