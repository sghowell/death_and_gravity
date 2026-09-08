# Exact vector elimination and nonlinear source remainder

[S6.46](FORMULATION.md) derives the ordered, exact classical vector
elimination identity for unchanged S6.42. Its [operator proof](notes/operators.md)
keeps boundary conditions, retarded state and the formal quantum
determinant distinct. A formal resolvent is not an inverse bound.

The [source proof](notes/source.md) obtains explicit cubic remainders
for the full actual nonlinear source and its first spatial derivative,
plus a fifth-order substitution error for the first local curl density.
Nontrivial source preparation survives the full lapse map. Higher
Fourier harmonics and variable vector coefficients are not suppressed
by assertion or removed from the full source.

These are exact classical identities and specified-field source bounds,
not full nonlinear evolution, nonlocal-action or quantum error bounds,
a cutoff, V/G/B conditions, or original P8 closure.

Read-only replay with local P8 source roots on PYTHONPATH:

    python -m p8_aligned_elimination.verify --check
