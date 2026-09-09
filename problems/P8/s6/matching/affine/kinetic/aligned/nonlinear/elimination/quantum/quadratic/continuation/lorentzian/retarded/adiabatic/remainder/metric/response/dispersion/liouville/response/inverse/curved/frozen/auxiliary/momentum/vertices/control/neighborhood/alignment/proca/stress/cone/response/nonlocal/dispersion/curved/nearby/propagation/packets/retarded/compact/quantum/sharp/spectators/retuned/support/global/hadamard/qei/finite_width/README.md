# Actual finite-width relational quadratic bound

S6.101 constructs a fresh globally transported generalized Hadamard
reference state with explicit initial graph bounds. It gives a fully
specified finite-width upper bound on the original relational
test-energy reference functional for real compact samplers in
|u|<1/100. The final formula is in [notes/bound.md](notes/bound.md).

The low-frequency proof uses a new canonical chart regular at k=0.
The high-frequency proof retains the sixth-order graph's exact
remainder and the all-order state's initial tail, using whole-interval
ball jets rather than time samples. Its matrix Fourier estimate uses
positive frequency sums and an opposite-sign resolvent, not a small
clock-minus-matter denominator.

[The formulation](FORMULATION.md) freezes the scope. The read-only
certificate module is p8_proca_finite_width.verify, with --check
comparing a full ancestor/source replay to the saved report.
Scientific SymPy is unchanged. The ball computations use the existing
python-flint installation at an explicitly scoped 160-bit precision.

The result is deliberately nonoptimal and concerns the actual
quadratic relational test energy, not the full physical stress.
It does not close original P8.
