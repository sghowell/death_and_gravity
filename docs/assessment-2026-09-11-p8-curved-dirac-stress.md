# P8 continuation: absolute specified curved free Dirac tensor

S6.171 controls the absolute one-loop tensor of the specified free Dirac
sector on the actual CD geometry. Original P8(b), V/G/B and P8 remain OPEN.

## Result and reference choices

The separately named GY14-SAT8-MR/EC-N0 prescription retains the saturated
mass reference and specifies the dimensionally continued Euler/Weyl
curvature pole basis, with no added finite R^2 term in that basis.
A fixed zero-field Newton reference subtracts this free sector's
soft-metric contribution at the chosen mass scale. These are recorded
finite choices, not a fit to the bounce energy or a complete interacting
Newton dictionary.

All42 species and both actual S6.170 Hadamard states are retained.
Exact twenty-frame diagonal and off-diagonal remainder estimates give

    energy remainder < 1e410,
    pressure remainder < 1e416.

The full D-dimensional metric variation is taken before the finite limit.
This retains the curvature improvement and evanescent contributions.
A different dimensional pole basis produces an explicit finite R^2
difference; its nonzero bounce energy is a negative control.
The retained Euler contribution is conserved through H=0.

The referenced local two-derivative pieces obey absolute bounds below
1e594 for energy and 1e595 for pressure. Combining the entire paired
potential with these local, geometric and exact-state terms gives

    |rho| < 1e789, |P| < 1e789,
    both below 1e-11 times the specified kappa=1e800 scale.

This is a prescribed free tensor in the comoving orthonormal frame,
not an interacting parent stress, a relative error against the zero
bounce density, or a controlled corrected background/response.

The curvature calculation checks the primary spinor heat-kernel
coefficients in [arXiv1703.00908v2](https://arxiv.org/pdf/1703.00908v2)
with the repository's own pole normalization, curvature-sign dictionary
and full lapse/scale variation. External continuum results and written
analytic estimates are not claimed as formalized theorems.

## Verification

Final private science passed350 tests in25.52 seconds.
Fresh repository science passed350 in25.68 seconds.
Ordinary replay passed375 in2028.60 seconds; independent native CLI
replay passed. The captured complete P8 snapshot passed29183 tests
in3335.55 seconds, with final exit code0.

All595 captured test files remained present and unchanged.
Path-list SHA-256:

    fb2a2059711fbf2c0e9a0ae9d9974f3d0f22e16596527babf6bdee796ebeeb13

The full-only exact GCD adapter passed128 original tuple comparisons.
Final counters:34090 domain fallbacks,7340 exact descents,88 mixed
fallbacks. Native, direct science, ordinary and CLI used original SymPy.

The native93415-character report was transferred losslessly in eight
chunks and independently checked against19 frozen scientific hashes.
Its20 fields record58 named identities,58 scalar entries,32 proof gates,
9 controls and185 rejected inputs. Report SHA-256:

    38ce293e50e209addb4bd4aa19f739be8a4dfc54ece5ac3f1890fc6bba449d60

Independent checks include complete twenty-frame projectors, finite-angle
and radial estimates, nonzero-regulator Gamma integrals, fixed-comoving
pressure differentiation, covariant local variations, and full/reduced
Euler lapse and scale variations. Actual tiny coherences remain nonzero
at4800 and5200 digits and agree under refinement. Full mass-profile and
double-zero Newton-reference tests use850 digits. These supplement the
written bounds; they are not validated numerical integrations.

An initial unsimplified exact zero was simplified privately before freeze.
Subsequent lint/import/format corrections did not change a bound.
No prior frozen source or report byte was changed.
Exact23-file staging excludes later continuations and unrelated P4/P9.

## Next frontier

S6.172's native and repository-science checks establish a separately
named minimal-Einstein canonical-clock null-equation obstruction.
Its ordinary and independent CLI replays have passed; its full replay
is running. It is not a whole-row or full interacting-parent exclusion.

S6.173 has native and repository-science checks for the actual target's
leading clock sources, including its nonminimal curvature and derivative
terms and their constraint relation. Its fresh replays are running.
These calculations show why vacuum-only matching and a small prescribed
fermion tensor do not establish a common parent bounce.

The next construction must address the leading clock action, physical
matter metric, constraint structure and common vacuum/clock domain.
Full interacting state/stress, cutoff and higher-loop bounds, quantum
matching and controlled background/response remain open, as do the
vacuum contour/truncation and fixed finite-gravity IR/Regge estimates.
Scoped P8(a) and A.20-A.23 remain unchanged. No user-choice blocker
has been identified.
