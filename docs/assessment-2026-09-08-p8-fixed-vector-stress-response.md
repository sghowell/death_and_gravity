# P8: controlled fixed-vector-stress homogeneous response

Original P8 remains open. This follows the
[clock-cone margin audit](assessment-2026-09-08-p8-clock-cone-margin.md).
The action, vector state, finite matching prescription and
adopted closure contract are unchanged.

## New result

[S6.57](../problems/P8/s6/matching/affine/kinetic/aligned/margin/response/FORMULATION.md)
turns the certified vector stress bounds into an actual
first-order homogeneous response estimate on [-1/2,1/2].

The literal local action is restricted before eliminating a
nonzero spatial mode. Physical metric variation gives
F_v=3p and F_n=3delta*p-rho, retaining the lapse-dependent
spatial metric. The clock source is fixed by the same Ward
identity. The regular Hamiltonian has no inverse Theta or
spatial momentum, and the lapse is solved rather than imposed
as an additional zero initial datum.

A canonical normalization gives a two-component system with
matrix norm below 3/2 and forcing below 49eta0. Exact
continuous estimates control the entire interval:

    |v|,|4a³(p_v+3ell*s)|<=131eta0,
    |n|<=94eta0, |n'|<=48000eta0+4eta1,
    |s|<=44eta0.

The preserved free-matter momentum is the joint quantity
a³(p_m+3ell*v), not p_m separately.

The exact original physical point-chart lift stays in the
timelike tube. Its fractional scale-factor change is
<=450eta0 and its physical Hubble correction is
<=(98000eta0+8eta1)/tau. It retains opposite Hubble signs
at u=+-1/4 and hence an interior scale-factor minimum.

At L=10^12,R=1000, the actual vector value and derivative
bounds are below 10^-14. Thus this specified approximate
metric changes its scale by less than 4.5*10^-12 and its
Hubble rate by less than 10^-9/tau.

## Remaining work

This is the exact solution of the stated fixed-source linear
problem. The nonlinear chart lift is not an exact quantum
solution, and no full semiclassical residual or feedback
contraction is claimed.

The actual causal second functional variation, its contact
terms, renormalization and compatible initial covariance
remain to be controlled. An initial-state audit now shows why
simply copying old Cauchy data is not sufficient on the
response metric; its per-mode variation and state-preparation
calculation is undergoing separate checks. Corrected cones,
interacting cutoff and V/G/B remain open. There is no
user-intervention blocker.

## Verification

The report pins 12 sources and fully rebuilds S6.56.
It checks 22 named identities comprising 23 scalar entries,
36 continuous proof checks and 34 rejected inputs.
Focused science passes 8 tests in 1.35 seconds. The ordinary
suite passes 43 tests in 464.21 seconds without the broad
GCD adapter; the independent read-only CLI passes.

The full regression through S6.57 passes 5051 tests in
966.01 seconds without excluding a frozen checkpoint.
The adapter passes 128 original tuple comparisons and records
5879 domain fallbacks and 6509 exact descents. This run
collected before the S6.58 working tests were added.

Report SHA-256:
a200c10d71ce6771c8ce921e120e0a00da13e53658ba75860d6e31d5737ae3b7.
The evidence is exact symbolic/rational verification with
written continuous proofs, not proof-assistant formalization.
