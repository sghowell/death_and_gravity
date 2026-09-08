# P8: exact moving-clock covariance and physical stress readout

Original P8 remains open; the completed scoped P8(a) objective
is unchanged. There is no user-intervention blocker. This follows
the [acoustic-reduction audit](assessment-2026-09-08-p8-acoustic-vector-reduction.md).

## New result

[S6.71](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/dispersion/liouville/response/FORMULATION.md)
gives the exact covariance and tangent transformation from the
original physical modes to the two prepared acoustic clocks.
It preserves the covariance determinant and includes the
initial-zero time shifts in both the generator and the readout.

The two source-dependent acoustic shifts are different. They
need not vanish at the final endpoint even when the physical
source has compact support. Zero varied initial data do not
erase the original unvaried covariance or its squeezing.

The retarded mode Green function has the correct normalization
and sign and is independent of the normalized mode basis.
The variance-response kernel still depends on the selected
state; an independent squeezing fixture changes it by -3/2.

All four physical energy/pressure readouts and contact rows
are recovered, including the moving output volume, source
vertices, second mass contacts and covariance-history terms.
The history term differentiates the covariance, not an
arbitrarily enlarged full output observable.

The dimensional transverse and longitudinal pump jets are
retained before the finite limit. Mapping the same physical
subtraction through the acoustic coordinates leaves S6.68's
finite C10-to-C0 response unchanged. No independent scalar
renormalization or new C2-to-C0 stress bound is introduced.

## Next work and boundary

S6.72's exact isolated massive causal-inverse certificate has
been generated and hash-checked; its fresh ordinary, CLI and
full regression runs are underway. The actual curved no-loss
normal form and a prepared homogeneous tree-plus-retained-
Gaussian inverse have a separate draft proof and focused
tests underway. Neither is part of the frozen S6.71 claim.

Quantum stability/cones, nonlinear neighborhood estimates,
spatial/independent-state response, other loops, interactions,
cutoff, finite Wilson matching and common-parent V/G/B remain
separate obligations. A representation of the response is
not itself an inverse or a physical stability theorem.

## Verification

The report pins 15 sources and fully rebuilds S6.70. It checks
46 named identities containing 82 scalar entries, 13 audit
gates and 40 rejected inputs.

The focused science suite passes 18 tests in 41.17 seconds.
The ordinary suite passes 46 tests in 1329.67 seconds without
the broad GCD adapter. The independent read-only CLI passes.

The full P8 regression passes 5705 tests in 2135.77 seconds,
without excluding a frozen checkpoint. The adapter passes
128 original tuple comparisons and records 13731 domain
fallbacks and 6509 exact descents. It collected before S6.72
tests were added.

Report SHA-256:
8c1f2be8713fee678e977a07941510692cfd2f98a984b36db0e6cb180dbeb4c3.

The evidence is exact symbolic/rational computation with
written prepared transport, physical readout and unchanged-
regulator arguments, not proof-assistant formalization.
