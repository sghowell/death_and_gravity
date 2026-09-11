# P8: a free in/out state with complete transition remainder and finite energy

Original P8 remains OPEN. This follows the
[clock-transparent polynomial map](assessment-2026-09-10-p8-clock-transparent-polynomial-map.md).
No user intervention is required for the current research.

## A specified quadratic state and observable

[S6.165](../problems/P8/s6/continuation/s6_165/FORMULATION.md)
fixes the free Dirac operator in Minkowski space with

    M_+/-(t)=mF +/- Delta s(t/tau),
    s(x)=x/(1+x^8)^(1/8),
    mF=10^200, Delta<=3e197, tau=10^-100.

The clock map supplies the same kinematic scalar argument.
It does not make that clock a parent background solution.
Gauge interactions and higher-loop fermion effects are
not part of this explicitly quadratic calculation.

Positive asymptotic masses and integrable mass tails give
modewise in/out evolution. The free in-vacuum defines a
quasifree state. The observable is its free out-particle
Hamiltonian energy per spatial volume, normal ordered
against the out-vacuum. No finite total energy or global
infinite-volume Fock implementer is asserted.

## The complete transition remainder is retained

For each helicity,

    H_p=p sigma1+M(t)sigma3,
    omega=sqrt(p^2+M(t)^2),
    v=p Mdot/(2omega^2).

The instantaneous-basis interaction equations preserve
the mode norm. Only odd transition-coupling Dyson terms
contribute to beta. With E=sqrt(p^2+m0^2), m0=.99mF,

    eta=integral |v|dt<=pDelta/E^2<1/100,
    |beta-I1|<=sinh(eta)-eta<=eta^3/5.

This is a convergent remainder for all higher quadratic
transition terms, not an all-Feynman-loop claim.

Two boundary-controlled integrations by parts in the
first transition integral use the full varying frequency.
The exact profile estimates are

    integral s'=2,
    integral |s's''|=1,
    integral (s')^3<=2,
    integral |s'''|<=36.

They give the complete bound

    |beta_p|<=A+B,
    A=5pDelta/(tau^2 E^4),
    B=p^3Delta^3/(5E^6).

Both terms fall as p^-3. Neither the higher terms nor
their interference allowance is discarded.

## Complete momentum integration

There are N=6 active color/flavor components and four
Dirac spin/particle-antiparticle factors. The twelve
inert flavors have constant masses and zero flat-space
production. Using omega_out<=2E and
(A+B)^2<=2A^2+2B^2 gives

    rho_out<=8N/pi^2 [
        5Delta^2/(tau^4 m0^2)+Delta^6/(225m0^2)].

The two radial integrals are exactly 1/(5m0^2) and
1/(9m0^2). With pi^2>9 the actual allowances are

    first-transition part: 2.4487297214569941843e396,
    higher-transition part: 1.7630853994490358127e783,
    complete out-particle energy: below 1e784.

The uniform transition amplitude is below 1e-9.
The energy allowance divided by the named
kappa*m_Phi^4 reference scale is below 1e-16.
This is not a relative error against the bounce
density, which vanishes at H=0.

Finite out-particle energy alone does not establish
Hadamard regularity, transient local stress, curved
production or backreaction.

## Independent verification

The native report pins eighteen source/proof/test files
and twenty fields: 48 named identities, 48 scalar entries,
22 proof gates, 9 controls and 154 rejected inputs.
Final private science passed 262 tests in 22.47 seconds;
fresh repository science passed 262 in 22.77 seconds.
Ordinary replay passed 287 tests in 2229.71 seconds.
Independent native CLI replay passed. The full captured
P8 snapshot passed 27239 tests in 3251.43 seconds with
final exit code 0.

All 583 captured test files were present and unchanged.
Path-list SHA-256:

    e179fbdef16e88b6fea5e02fce8542a627d2e6c3e3b5dd56b3242648ea976e4e

The full-run exact GCD adapter passed 128 original tuple
comparisons. Final counters: 34099 domain fallbacks,
7340 exact descents and 88 mixed fallbacks. This snapshot
predates S6.166. Only full regression uses the adapter;
native, direct science, ordinary and CLI use unmodified
SymPy with interpreter-only allowances.

The 24150-character native report was transferred
losslessly in three chunks. All eighteen source hashes
were independently verified. Report SHA-256:

    b9b4c42ca2a054b86797a78e4e6b3e53658288cee87e9cf2b65579015cccd6d0

Independent tests integrate full profile and radial
functions, solve finite-interval mode equations with
step refinement and norm checks, and resolve nonzero
higher transition terms. Constant phase, zero momentum,
constant mass and opposite-profile controls are retained.
The actual extreme-parameter bound is analytic and
rational, not inferred from numerical mode sampling.
These diagnostics are not validated integration or
formalization.

Cleanup preceded final verification and freeze. Exact
22-file staging excludes later continuations and
unrelated P4/P9 work. No published scientific source
or report changed.

## Remaining work

S6.166's explicit Hadamard-theorem application has passed
native and direct science; its independent replays are
running. A sharper exact-frame state-energy estimate has
passed private tests and is frozen, but awaits its
native replay rebuild. It is not counted as published.

A state of the free flat quadratic operator is not the
interacting curved common-parent state. Absolute local
stress, physical cutoff and loop errors, quantum target
matching, global V contours/cuts and finite-gravity G
remain open. Scoped P8(a) and A.20-A.23 are unchanged;
original P8(b) and original P8 remain OPEN.
