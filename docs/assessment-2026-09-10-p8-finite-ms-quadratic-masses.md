# P8: finite MS masses of both quadratic primitives

Original P8 remains OPEN. This follows the
[finite-slope checkpoint](assessment-2026-09-10-p8-finite-ms-quadratic-slopes.md).
No user intervention is required for the current research.

## Complete regulated mass references

[S6.147](../problems/P8/s6/continuation/s6_147/FORMULATION.md)
computes the zero-momentum MS masses of the same complete scalar
and gauge quadratic primitive rows. Dimension-symbolic vacuum
derivatives independently reproduce the direct tensor reduction.
Both self-energy placements, the vertex word and all assigned
proper fermion mass, kinetic and Yukawa counterterms are retained.
Fixed-scale differentiation precedes setting mu=m.

The finite massless-exchange anchors are -56 in NY^2 m^2/Q^2
units and +40 in NYaC_F m^2/Q^2 units. Their regulated double
and simple poles cancel only after the complete forest is paired.
Finite epsilon-times-pole contributions are kept; the pi-squared
terms cancel exactly. The gauge exchange is actually massless.
The scalar massless exchange is only an auxiliary reference.

Replacing that auxiliary scalar by the physical mass-one scalar
uses a regulated integral over its squared mass, including the
nonzero whole-fermion-cycle quartic counterterm. The leading
difference is NY^2[78+32 log(m^2)]/Q^2; its explicit remainder
is bounded by NY^2[175+50 log(4m^2)]/(4Q^2 m^2).
The finite difference is about 5.11704625009938e-410.

The combined absolute zero-momentum mass bound is approximately
2.18964072398307e-12. The previously bounded slope and nonlocal
remainder are retained in continuing to s=1; the resulting
on-shell mass bound is below 1e-10. These are bounds on two
assigned primitive rows, not the complete canonical pole.

## Independent verification

The report pins 18 source/proof/test files and 20 fields:
49 named identities, 49 scalar entries, 23 proof gates,
9 controls and 51 rejected inputs. Final private science passed
186 tests in 24.05 seconds; fresh repository science passed
186 in 23.97 seconds. Ordinary replay passed 211 tests in
2342.45 seconds. Independent native command-line replay passed.
The complete P8 snapshot passed 22666 tests in 3261.07 seconds.

All 547 captured test files were present and unchanged. Path-list SHA-256:

    1c4b70c0d91a29033a22771c0a1aae575bfdb204a721c9a610e9e2fee126bf35

The full-run adapter passed 128 original tuple comparisons.
Counters were 34031 domain fallbacks, 7340 exact descents and
88 mixed fallbacks. This snapshot predates S6.148. Only full
regression uses the separately audited exact GCD adapter;
native, ordinary, CLI and direct science use unmodified SymPy
with interpreter-only allowances.

The 33124-character native report was transferred losslessly,
and all eighteen source hashes were independently verified.
Report SHA-256:

    f7dd5bec638832ee8df808daf51ac8efe0cab88f0aafd618ec30b3226dafc0e0

Independent tests include explicit Dirac matrices in dimensions
three, four and five, complex-regulator finite-part extraction,
fixed-scale field/coupling variation, massive triangle integrals,
the quartic-counterterm tadpole, regulator-first scalar mass
integration, finite remainders and mutation controls.
All private cleanup preceded the final private run and freeze.
No frozen scientific source or report changed.
Written analytic proofs and exact tests are not formalization or
independent peer review. Exact 22-file staging excludes later
continuations and unrelated P4/P9 changes.

## Remaining work

The older scalar insertion's finite local references have passed
private, repository, native, ordinary and command-line checks;
their full regression is running. Both fermion vacuum rows have
passed private, repository and native checks; their independent
replays are running. Neither is counted as published here.

Other parameter/counterterm and canonical conversions, the
complete matched two-loop pole/error, remaining vacuum/reference
assembly, truncation control, V contours, finite-gravity G and
common-parent B remain open. No all-orders UV construction is
added to the adopted finite-EFT/necessary-positivity scope.
Scoped P8(a) and A.20-A.23 are unchanged. Original P8(b) and
original P8 remain OPEN.
