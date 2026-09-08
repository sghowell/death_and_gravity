# Prepared on-shell linear light family

[S6.44](FORMULATION.md) constructs a nonzero solution of the actual
regular CD/M1 linear light equations with the preparation required
by S6.43. Exact initial Cauchy jets and continuous bounds hold through
the entire declared bounce window, including the center.

The actual quadratic source is nonzero. A real cosine's output
momentum and source-jet envelope meet the leading retarded-response
theorem, with an explicit conservative rational initial amplitude.
See the [written proof](notes/proof.md).

This is linear light evolution and its leading heavy response, not
a full nonlinear parent solution, higher-order error bound, nonlinear
constraint/stability theorem, loop/cutoff or V/G/B UV verdict. Original
P8 remains open and no frozen physical action changes.

Read-only replay with the local P8 source roots on PYTHONPATH:

    python -m p8_aligned_onshell.verify --check

Tests independently reconstruct the moving canonical map and Cauchy
recurrence, check physical ADM readouts, exact compact bounds,
nontrivial preparation, Fourier output and response transfer, and
reject report mutations and inexact or invalid inputs.
