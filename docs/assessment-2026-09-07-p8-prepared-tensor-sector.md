# P8: a convergent prepared tensor sector, not yet a sourced EFT

Continuation of the [physical-response and constraint checkpoint](assessment-2026-09-06-p8-physical-response-and-regular-constraints.md).
The unchanged variable-coupling parent now has a proved analytic two-mode
tensor sector with its actual physical normalization and source-loading
bounds. Original P8 remains open: selecting a source-free sector does not
establish causal low-energy reduction for general matter sources, the
original C/D operator match, or the adopted vacuum/finite-gravity UV gates.

## What is now proved

The [S6.23 prepared-sector formulation](../problems/P8/s6/matching/variable/response/prepared/FORMULATION.md)
and [complete proof](../problems/P8/s6/matching/variable/response/prepared/notes/proof.md)
retain the full coupled tensor equations of the pinned S6.20 action.
Their physical domain is `0<delta=c-2<=10^-9`, `|u=T/tau|<=1/100`,
and `0<=K=(tau*k_com)^2<=4`. The apparent heavy pole is cleared by the
exact substitution `Q=Dq`, `D=(2+delta)(1+u^2)^4-2`.

A coefficient-space inverse and an explicit contraction below one half
construct two exact joint-analytic solutions, with light initial data
`(l,l')=(1,0)` and `(0,1)` at the bounce. The complex proof radii are
`|u|<=1/20`, `|delta|<=1/400`; they do not describe physical negative-delta
parents. This is a convergent solution-space theorem, not only a formal
inner expansion or a sampled propagation calculation.

The restricted conserved symplectic normalization obeys
`|Omega-1|<=3delta^3/5`. The actual g-metric Wronskian is between `3/5`
and `2` throughout the physical slab. Consequently the exact closed
source-free g equation on this selected sector has a positive kinetic
normalization between `9/20` and `11/6`.

At the center, its zero-derivative coefficient G satisfies

```
|G(0)/K - 1 - (14/33)delta| <= 48001delta^2.
```

The ratio at K=0 is defined by proved analytic continuation: a constant
common physical tensor is an exact solution there. The coefficient
differs from the locked-field approximation's `1+(4/5)delta+O(delta^2)`.
This comparison concerns a prepared source-free response at finite K;
it is not a characteristic-cone or causality theorem.

## Fixed preparation versus a retuned source

The endpoint comparison retains S6.21's actual imported band `1<=K<=4`
and fixed physical slices `u=+-1/100`. Its analytic inclusion differs
from the limiting regular-light inclusion by at most `200delta`.
The full fixed-slice output error for fixed limiting light data is
therefore bounded by `8600delta` times the input norm. This improves the
upper rate from order `delta^(1/3)` to order delta; it does not establish
a sharp or nonzero leading leakage coefficient.

For the same fixed g-only source used in S6.21, the bound is

```
delta [8600 ||target|| + 126000000 integral |sigma_0| du].
```

The fixed-source estimate does not give that protocol the specially
prepared center expansion. A separately specified delta-dependent source
loads the analytic sector exactly. Its support remains strictly inside
`(-1/50,-1/100)`. A pole-cancelled loading formula retains four derivatives
of an explicit smooth cutoff and gives finite source and retuning bounds.
The deliberately coarse supremum constant is `321137485366608000` per
unit l1 target amplitude. This is a linear-probe estimate, not a claim of
small finite-amplitude stress or of a low temporal-frequency band.

## Remaining completion gates

The positive kinetic factor belongs to the selected source-free sector,
not to an unrestricted sourced, spatially local EFT. A general incoming
state or source is not assumed to select that sector. A controlled causal
source/state prescription, a justified frequency window and omission
errors are still needed before interpreting an effective light operator.

The scalar/vector constraint result from S6.22 remains distinct from a
uniform physical scalar propagation theorem; its center-frozen-map warning
must be resolved by the ongoing unreduced-action and exact-symbol work.
No provisional scalar coefficient is promoted by this tensor checkpoint.

A successful parent would still need a common-action/domain and actual
matter-frame match to the original C/D witness, appropriate gap/cutoff,
omitted-operator and loop estimates, and the adopted vacuum plus
finite-gravity positivity analysis. Named-parent results cannot decide an
entire DHOST row's UV status. The completed scoped photon objective is
unchanged. None of these remaining mathematical tasks currently requires
a new user choice or external permission.

## Verification record

The frozen [report](../problems/P8/s6/matching/variable/response/prepared/certificates/analytic-prepared-sector.json)
has SHA-256
`a8d89126980dd8ec458b3161bc8b397f74211736f02e64d4e9e29b95e43a0ade`.
It pins 18 sources and records 33 exact residuals, 38 strict continuous
margins, 27 rejected-domain controls, 52 independent inverse-polynomial
checks and independent Fraction replays of all norm/source calibrations.
The separately authored physical-action and norm audit contains 26 tests.

The full child suite passes **81 ordinary tests in 163.91 seconds**.
Its fresh standalone CLI also passes with `PYTHONHASHSEED=0`, without an
arithmetic adapter. The full-tree exact-arithmetic-adapted regression
passes **2,533 tests in 531.57 seconds**, after 128 startup equivalence
checks; its counters are 6,509 exact descents and 3,157 domain fallbacks.
All frozen source hashes, Ruff and default whitespace/EOF checks pass.
The unfinished microlocal
child is excluded from that checkpoint regression and commit manifest;
unrelated P4/P9 work remains untouched.
