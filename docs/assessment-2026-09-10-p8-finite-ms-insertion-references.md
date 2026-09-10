# P8: finite MS local references of the older fermion insertion

Original P8 remains OPEN. This follows the
[quadratic mass checkpoint](assessment-2026-09-10-p8-finite-ms-quadratic-masses.md).
No user intervention is required for the current research.

## The paired insertion's actual local terms

[S6.148](../problems/P8/s6/continuation/s6_148/FORMULATION.md)
completes the finite outer MS mass, slope and assigned heavy-source
references of scalar_Phi2_W1_F0. The inner fermion subgraph retains
its full physical scalar mass/residue subtraction. Its nonlocal
part was already bounded; that subtraction did not determine or
erase the finite outer local references.

The stationary heavy-field source is included once. It cancels
the apparent -g/M piece in the local tadpole coefficient, leaving

    f_MS(s)=L T_insert,MS/2 - g F_MS(M;s).

The inserted tadpole has two exact beta anchors and a finite
remainder bounded by 24/T, where T=4m^2. In units C/Q, C=2NY/Q,
the finite part is

    T(13/4+pi^2/8) - 47/9 - pi^2/6 + delta_T.

The heavy bubble uses both independent ratios 1/T and M/T.
Its leading finite anchor is 47/18+pi^2/12 in the same C/Q
units. Its remainder is bounded by
18/T+2(M/T)[log(T/M)+1]. Nonzero simple-pole residues in the
remainder integrals are paired with the finite regulator
prefactors before MS subtraction. The physical scalar mass
and actual heavy mass are not replaced by an auxiliary limit.

The resulting local mass bound is approximately
8.96020706191107e-13; the separate heavy-bubble mass term is
bounded by 2.22227357686286e-216. The on-shell mass is below
1e-10, and the slope is below 5.55568394215716e-617.
The assigned heavy source can be large in absolute units:
its bound is 4.77877709968590e188. The corresponding source/M
bound is 2.44673387503918e-9, below 1e-6. This is not the
complete vacuum/source or canonical pole calculation.

## Independent verification

The report pins 18 source/proof/test files and 20 fields:
47 named identities, 47 scalar entries, 23 proof gates,
9 controls and 105 rejected inputs. Final private science passed
227 tests in 26.75 seconds; fresh repository science passed
227 in 26.93 seconds. Ordinary replay passed 252 tests in
2156.51 seconds. Independent native command-line replay passed.
The complete P8 snapshot passed 22918 tests in 3234.73 seconds.

All 549 captured test files were present and unchanged. Path-list SHA-256:

    92ea4373d894a96da7698fb2b3ea4e456df2a8e33897546f48c920ceb9e851ed

The full-run adapter passed 128 original tuple comparisons.
Counters were 34045 domain fallbacks, 7340 exact descents and
88 mixed fallbacks. This snapshot predates S6.149. Only full
regression uses the separately audited exact GCD adapter;
native, ordinary, CLI and direct science use unmodified SymPy
with interpreter-only allowances.

The 41315-character native report was transferred losslessly,
and all eighteen source hashes were independently verified.
Report SHA-256:

    b7740228c719dede347afa67d3b5853c677c984ba4a15bfa4209836e38c8408d

Independent tests cover nonzero-regulator beta/spectral
integrals, complex finite-part extraction, both mass ratios,
joint on-shell derivatives, exact source elimination and source
cross-order ownership. Explicit Gamma recurrences, rather than
an unreduced symbolic simplifier result, verify the anchors.
All private fixes preceded the final private run and freeze.
No frozen scientific source or report changed.
Written analytic proofs and exact tests are not formalization or
independent peer review. Exact 22-file staging excludes later
continuations and unrelated P4/P9 changes.

## Remaining work

The two fermion vacuum rows have passed private, repository and
native checks; their independent replays are still running.
The isolated fixed-contact conversion is private research.
Neither is counted as published here.

Other parameter/counterterm and canonical conversions, the
complete matched two-loop pole/error, remaining vacuum/source
assembly, truncation control, V contours, finite-gravity G and
common-parent B remain open. No all-orders UV construction is
added to the adopted finite-EFT/necessary-positivity scope.
Scoped P8(a) and A.20-A.23 are unchanged. Original P8(b) and
original P8 remain OPEN.
