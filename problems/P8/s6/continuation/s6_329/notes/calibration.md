# Exact and independent tests

Checks derive the frozen CDF, Levy exponent, finite-cutoff count,
omitted-energy mean and variance, independent-band energy sum,
conditional moments and original numerical margins. Domain tests
reject floating, symbolic, negative and out-of-range scalar inputs.

Tests preserve exact positive values at index1/10^800 and very small
resolution; Monte Carlo cannot resolve those rare scales and is not
used as evidence. Independent inverse-CDF quadrature at moderate
calibration indices checks the conditioned moment formula, and direct
integrals check the omitted low-energy Laplace exponent bound.
Poisson second moments retain their variance in addition to the
square of the mean.

Actual original S328 coefficient values in real/complex TT frames,
including collinear directions and split atoms, are used as marks.
A deterministic total-energy example rejects replacing a calorimetric
event by individual-particle cuts. The proof of the infinite cloud
and marked Lp limit is the written positive-band and dominated
convergence argument, not any finite sample or numerical quadrature.
