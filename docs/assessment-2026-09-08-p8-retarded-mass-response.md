# P8: selected-state retarded mass response

Original P8 remains open. This follows the
[physical-signature pole audit](assessment-2026-09-08-p8-physical-signature-quadratic-pole.md).
The original clock, selected state and fixed scalar profiles are unchanged.

## New result

[S6.64](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/FORMULATION.md)
derives the bare selected-state mass-response kernel on the unchanged
physical background. It includes the local temporal-constraint contact,
which is present even for a covariant mass source linear in its amplitude.

The actual two transverse and one longitudinal modes give explicit
TT, TL, LT and LL integrands for arbitrary external spatial momentum.
An independent four-index Wick contraction gives the same result.
The temporal-spatial longitudinal phase and the regular isotropic
zero-momentum limit are retained.

The homogeneous Kubo commutator agrees with the exact canonical
covariance response for general symbolic transfer and covariance
matrices. In flat space, the bubble plus contact reproduces the
exact vacuum-energy susceptibility and all three frozen timelike
frequency UV residues. Omitting the contact fails both checks;
its pole error is nonzero.

Sources vanishing on an initial neighborhood leave the selected
initial covariance unchanged for this restricted perturbative response.
Ordinary-Proca CCR imply physical-background causal support away
from coincidence. No finite momentum cutoff, arbitrary varied
Cauchy state or finite off-clock cone is thereby certified.

## Remaining work

The homogeneous finite local adiabatic subtraction component is
now in its separate replay and regression. Its result is not imported
into this frozen bare-kernel claim. Work on the subtracted exact
nonlocal response and its momentum bounds is underway.

Common-dimensional finite matching, the remaining integral, full
metric and second mass-source blocks, coupled feedback/stability,
cones, interactions, cutoff and V/G/B remain open. There is no
user-intervention blocker.

## Verification

The report pins 16 sources and fully rebuilds S6.63. It verifies
65 named identities containing 147 scalar entries, 13 algebraic audit
checks and 140 rejected inputs. The focused science suite passes
14 tests in 2.09 seconds.

The ordinary suite passes 45 tests in 535.66 seconds without the
broad GCD adapter. The independent read-only CLI passes. The full
regression through S6.64 passes 5383 tests in 1181.24 seconds, with
no frozen checkpoint excluded. The adapter passes 128 original
tuple comparisons and records 6823 domain fallbacks and 6509 exact
descents. This run collected before the S6.65 tests were added.

Report SHA-256:
17a31b89b503e42852dc6b60c986e615c02c582cd4b72599b28774f8e1d235a5.
The evidence is exact canonical, polarization, covariance and
Gamma-residue computation with written causal arguments, not
proof-assistant formalization.
