# P8: exact acoustic vector reduction and prepared potential bounds

Original P8 remains open; the completed scoped P8(a) objective
is unchanged. There is no user-intervention blocker. This follows
the [metric dispersion and chart audit](assessment-2026-09-08-p8-metric-dispersion-chart.md).

## New result

[S6.70](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/dispersion/liouville/FORMULATION.md)
constructs positive sector-dependent acoustic times directly
from the physical transverse and longitudinal Hamiltonians.
The canonical transformation preserves the Wronskian and
retains the factor required for the physical stress readout.

The longitudinal mode potential is a scalar principal potential
plus an exact momentum-dependent remainder. The potential has
a continuous zero-momentum limit, but its canonical mode map
is used only at positive momentum. On I=[-1/2,1/2] the remainder
is nonnegative and bounded by
88*m^2/(k^2+m^2).

Prepared lapse/scale variations retain both the time-rate
variation and the fixed-acoustic-time history shift. The
physical source chart is invertible, and a rational continuous
bound gives

||delta R_k||_C0 <= 18000*m^2/(k^2+m^2)*||source||_C2.

This is a mode-potential bound, not a C2-to-C0 bound on the
renormalized stress. The original selected state, nonzero
unvaried initial mixing and physical finite prescription remain
unchanged. No off-clock ordinary-Proca Hadamard or cone result
is transferred through the acoustic coordinate change.

## Next work and boundary

The exact moving-clock covariance and physical-readout bridge
has a generated, hash-checked S6.71 certificate and is undergoing
fresh replay and full regression. Separately, a first-sheet
and causal-kernel inverse for the isolated exact massive flat
dispersion block is in development. Neither is included in
the frozen S6.70 claim.

The actual curved no-loss stress remainder and full coupled
causal inverse remain open, as do spatial/independent-state
response, quantum stability/cones, interactions, other loops,
cutoff, finite Wilson matching and common-parent V/G/B.

## Verification

The report pins 14 source files and fully rebuilds S6.69.
It verifies 30 named identities containing 36 scalar entries,
12 audit gates and 93 rejected inputs.

The focused science suite passes 17 tests in 86.57 seconds.
The ordinary suite passes 45 tests in 1366.21 seconds without
the broad GCD adapter. The independent read-only CLI passes.

The full P8 regression passes 5659 tests in 2031.46 seconds,
without excluding a frozen checkpoint. The adapter passes
128 original tuple comparisons and records 13567 domain
fallbacks and 6509 exact descents. This run collected before
S6.71 tests were added.

Report SHA-256:
daab6d10f0dfa4190153fa63f5bb6eed9829c1ae73a0704eac34efb5c5b0da10.

The evidence is exact symbolic/rational computation with written
canonical, moving-clock and continuous interval arguments,
not proof-assistant formalization.
