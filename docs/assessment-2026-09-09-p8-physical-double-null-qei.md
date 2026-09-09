# P8: an actual physical double-null inequality

Original P8 remains open; no user intervention is needed.
This follows the [null-line state-cap obstruction audit](assessment-2026-09-09-p8-null-state-cap-obstruction.md).

## New certified result

[A.22](../problems/P8/a/fields/nonminimal/null/double/FORMULATION.md)
derives a finite nonoptimal inequality for the actual physical
null stress of a real massless scalar, smeared over a timelike
two-plane, for 0<=xi<=1/2. It retains the full nonminimal improvement
and a one-sided relative Wick-square upper cap on that plane.

For xi=1/6 and the stated normalized compact profile, a rational
boost gives the explicit flat-space lower bound

    integral dt dz f^2 <:T_DD:>
      >= -14117 hbar/(1536pi^2 delta_plus^3 delta_minus)
         -4 Phi_*^2/delta_plus^2.

Both widths are positive. The quantum cost diverges when the
second width shrinks to zero, consistently with A.21. This is not
an optimal bound or the square of a transversely smeared field.

The conformal FLRW transport retains physical plane volume,
the specified affine null normalization, the actual scalar
reference stress and its finite beta_S term. A finite-plane
inequality alone is not a single-ray focusing theorem.

## Verification

The losslessly saved report pins 17 sources and fully rebuilds A.21
and its ancestry. There are 39 named identities, 51 scalar entries,
18 proof gates, 39 rejected inputs and mutation controls for all
20 report fields.

Native science passes 95 tests in .29 seconds.
Independent ordinary native replay passes 139 tests in 8.46
seconds; the separate read-only CLI also passes.
Scientific SymPy is unchanged.

The complete regression passes 8494 tests in 2513.93 seconds,
with all 453 captured files present and unchanged. Path-list SHA-256:

    4ee3078a75ea20abcb27d86f99fc3c72165e3ef890530d126c07bc7c7a89faef

The unchanged exact GCD adapter passes 128 original tuple comparisons,
with 32031 domain fallbacks, 7340 exact descents and 4 mixed fallbacks.
This snapshot predates A.23's added tests.

Report SHA-256:

    7f37845da53396f5bbc171c43572a07eb10c493c135247b58fe255784c50b15b

The positive-type plane restriction, Sobolev limits and conformal
QFT arguments are written and source-pinned, not proof-assistant
formalized or independently peer reviewed. Exact 21-file staging
verifies every staged source and report hash plus audit updates.

## Continuing work

A.23 now has a frozen 19-source report, 75 native science tests,
123 ordinary tests and a separate CLI passing. Its complete
455-file regression is running. It derives a conditional global
FLRW null endpoint within affine duration 2tau using the actual
SEE, finite-plane QEI and an outgoing index identity. Its exact
worst-case margin exceeds one half. It does not assume a
homogeneous quantum state or a single-null QEI.

A.23 also provides distinct controls: a future-complete geometry
satisfying all geometric caps, and an actual thermal scalar SEE
state satisfying the past conditions. The former is not a small-
source SEE witness; the latter does not prove the conditional
future caps. Root closure is not assigned before its full replay.

Private S6.102 work reconstructs the physical normal matter density
and pressure with the second-order lapse response retained. At the
bounce, 420 native identities across all 105 independent zero-output
phase pairs pass. The tensor and Proca sectors have nonzero
constraint-induced matter-energy contributions. These are classical
observable checks, not a conserved renormalized source or quantum
backreaction solution.

Original P8(b)'s full source/Ward completion, interactions,
omitted operators and actual V/G/B UV matching remain open.
