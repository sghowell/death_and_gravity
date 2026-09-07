# Sources, sign conventions, and claim boundary

Read on 2026-09-07. Paths below are relative to this gate unless an absolute
repository path is stated. No source below asserts the present parent is UV
complete, or supplies its missing original C/D matching data.

## Immutable action and preparation lineage

The actual physical action, conventions and background are the pinned
S6.20 VARIABLE ones in `../../FORMULATION.md` and
`../../src/p8_variable_beta/canonical.py`. The convention is `+---`,
`R_B=-6(D_T H+2H^2)`, Einstein action `-M^2 R_B/2`, with `u=T/tau`.
S6.21 (`../`) fixes the actual canonical/physical Cauchy and source maps.
Its uniform matched-transfer momentum restriction is **not** imported into
the present `0<=K<=4` coefficient proof. The mass pole isolation is derived
again independently of K.

S6.23 (`../prepared/`) establishes its analytic source-free two-dimensional
sector and symplectic/physical Wronskian. The new symplectic source projection
uses that theorem, but is not identification of a general forced solution
with that sector. The read-only replay pins S6.23 report SHA
`a8d89126980dd8ec458b3161bc8b397f74211736f02e64d4e9e29b95e43a0ade`
and replays its S6.21/S6.20 lineage. The adopted original contract remains
`problems/P8/s6/FORMULATION.md`, not a renamed local-EFT criterion.

## Exact special function, not an imported approximation

[NIST DLMF 15.10.1–2](https://dlmf.nist.gov/15.10) gives the hypergeometric
equation and its local Gauss-function solution. We substitute the stated
parameters into that equation and derive the real-potential Jost Wronskian
from its normalized series at zero argument. The single solution at that
argument is valid despite the integer difference between a and b; no pair
of solutions at infinity is used. The retarded jump and the finite-delta
Liouville factors are derived here rather than read off a scattering formula.

## Decoupling language audit — a different physical model

[Achucarro et al., arXiv:1205.0710v2](https://arxiv.org/pdf/1205.0710),
equations (9)–(14), distinguishes propagating frequency separation from a
field-space mass and requires an adiabatic low-energy regime when replacing
the heavy equation algebraically. Its two-scalar inflation action is not
our two-tensor parent. We use this source only to police terminology:
neither its frequency/cutoff formulas nor its EFT existence conclusion is
transferred here. Our proven elimination retains a retarded memory kernel;
the source RMS and B-based quantities are explicitly only the defined
diagnostic scales, not substitutes for that paper's normal-mode frequencies.

## What is new and what stays open

The finite-source norm estimate, fresh continuous K-inclusive coefficients,
causal feedback bound, actual physical map, and projected-source control are
the calculations in this gate. They do not require a new background or
source prescription. The result is local in its specified time domain but
nonlocal as a time-response functional. A local derivative expansion with
a controlled class of approximately low-temporal-band sources, general
initial states and omitted-operator errors remains unproved. No numerical
cutoff is manufactured, and no original S6/P8 or UV closure is asserted.
