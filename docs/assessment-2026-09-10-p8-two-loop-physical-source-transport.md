# P8: full two-loop physical-source derivative-map transport

Original P8 remains OPEN. This follows the
[fixed-order reference assembly](assessment-2026-09-10-p8-complete-second-heavy-source.md).
No user intervention is required for the current research.

## Full action and the specified physical source

[S6.160](../problems/P8/s6/continuation/s6_160/FORMULATION.md)
transports the named canonical GY14 Phi pole, b2, vacuum and
stationary-H results through two loops using the literal
S6.111 cubic derivative map F(Psi)=Psi+R[Psi].

The full parent action, all fixed parent counterterms and
the Jacobian are transformed together. The source is
J F(Psi)+J_H H, not J Psi. The map coefficients are fixed
numerically at the named canonical interaction boundary.
The same D-dimensional lift, with the required dimensional
coupling factors, is used before every finite-part operation.

This defines the pullback of the canonical renormalized
generating functional. It does not introduce an independent
ordinary-Psi MS composite-operator prescription. In particular,
no numerical bound on that operator's second residue is claimed.

The full generated scalar action contains degrees six, eight,
ten and twelve as well as the quadratic and quartic pieces.
The cubic-Phi Yukawa term, H-dependent terms and generated
parent-counterterm interactions are retained. The H Gaussian
elimination commutes with this H-independent change of Phi.

For E physical scalar source endpoints, parent action valence
sum N, vertex count V, counterterm loop weight j and total
map degree p, the connected effective loop count is

    L_eff = 1-E/2+(N-2V)/2+j+p.

At loop two this permits map degree three for four physical
Phi sources, degree two for two sources and degree one for
the vacuum or a single H source. The 42 enumerated support
tuples are necessary possibilities, not a count of actual
connected mixed-field graphs. All higher generated terms
remain in the action even when they enter beyond this order.

## Regulated transport and diagnostics

Finite-dimensional Jacobians are retained in independent
Gaussian and coupled-map tests. General Gaussian moment
recurrences cancel through map degree three. Deliberately
omitting the Jacobian, transformed source, sextic/octic
scalar terms or cubic-Phi Yukawa term produces nonzero
diagnostic defects at the predicted orders.

In the continuum, the local derivative-map Jacobian yields
scaleless polynomial ghost momentum integrals in dimensional
regularization. This statement is made within the joint
regulated action/source prescription, not by separately
discarding poles and then multiplying finite parts.

A symbolic ordinary-field overlap cancels from the on-shell
physical amplitude when all source and amputation factors
are included. That algebra does not provide a numerical
second ordinary-Psi overlap or a global inverse map.

The transported physical-reference bounds remain unchanged:

    total through-two-loop b2 relative allowance <1e-6,
    second b2 relative allowance                 <1e-7,
    unit-disc Phi pole coefficient              <1e-18,
    complete second vacuum reference            <1e595,
    complete second stationary-H source         <1e189.

These are the existing fixed-order canonical bounds, not
new physical truncation estimates.

## Independent verification

The native report pins 18 source/proof/test files and 20
fields: 81 named identities, 81 scalar entries, 23 proof
gates, 9 controls and 101 rejected inputs. Final private
science passed 221 tests in 21.88 seconds; fresh repository
science passed 221 in 22.13 seconds. Ordinary replay passed
246 tests in 2086.02 seconds. Independent native CLI replay
passed. The complete captured P8 snapshot passed 25552
tests in 3203.88 seconds with final exit code 0.

All 573 captured test files were present and unchanged.
Path-list SHA-256:

    a96e3d1a4d866107be66d06860866c5bd6c2eddb5a5739f5069e562430562aea

The full-run adapter passed 128 original tuple comparisons.
Counters were 34104 domain fallbacks, 7340 exact descents
and 88 mixed fallbacks. This snapshot predates S6.161.
Only full regression uses the separately audited exact
GCD adapter. Native, direct science, ordinary and CLI retain
unmodified SymPy with interpreter-only allowances.

The 68249-character native report was transferred losslessly
in six chunks; all eighteen source hashes were independently
verified. Report SHA-256:

    e54116c62e54a0587bcb1e9dbb312284a905434a6db8276fdabef2842cd4dea5

Cleanup preceded final private verification and byte freeze.
No frozen scientific source or report changed. The exact
checks and written regulated arguments are not formalization
or independent peer review. Exact 22-file staging excludes
later continuations and unrelated P4/P9 changes.

## Remaining work

The complete second low-energy Phi elastic cut has passed
native and direct science; independent regressions are still
running. A bound matching the full finite-kappa analytic
target on the earlier restricted classical function class
is in private development. Neither is counted as published here.

Ordinary-Psi composite normalization, physical finite-EFT
truncation, global V contours and cuts, finite-gravity G and
common-parent B remain open. No global inverse, rolling-state
transport or all-orders UV construction is inferred or added
as a new requirement. Scoped P8(a) and A.20-A.23 are unchanged;
original P8(b) and original P8 remain OPEN.
