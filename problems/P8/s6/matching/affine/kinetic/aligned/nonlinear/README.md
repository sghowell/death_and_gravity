# Local nonlinear auxiliary-constraint theorem

[S6.45](FORMULATION.md) derives the complete scalar-clock ADM and
Hamiltonian constraints of the unchanged source-aligned S6.42 action.
Both lapse derivative cancellations retain their required boundary
terms. The joint lapse/normal-vector Jacobian is nonsingular through
the compact bounce interval, including the center.

The [written proof](notes/proof.md) gives seven physical modes including
free chi on a nonempty local nonlinear canonical neighborhood. It retains
the full vector, shear and spatial-diffeomorphism sectors and does not
infer an extra lapse boundary datum from a velocity-form equation.

This does not establish nonlinear stability, a quantitative neighborhood
size, other-background causality, higher-order or quantum errors, a
cutoff, vacuum/finite-gravity V/G/B conditions, or original P8 closure.

Read-only replay with the local P8 source roots on PYTHONPATH:

    python -m p8_affine_nonlinear.verify --check

Tests include independent dense Legendre solves, omitted-boundary and
Gauss-term controls, singular-pivot controls, exact compact bounds,
source-manifest and report-mutation checks. All ancestors stay frozen.
