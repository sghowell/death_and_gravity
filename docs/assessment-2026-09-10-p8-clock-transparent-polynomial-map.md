# P8: a polynomial clock-transparent map preserving the full-target match

Original P8 remains OPEN. This follows the
[protected Yukawa profile](assessment-2026-09-10-p8-protected-yukawa-mass-profile.md).
No user intervention is required for the current research.

## One map with two explicit local limits

[S6.164](../problems/P8/s6/continuation/s6_164/FORMULATION.md)
replaces the earlier source coordinate by the separately named
polynomial map

    Fhat(Psi)=Psi+q8(X) R_old(Psi), X=(partial Psi)^2/kappa,
    q8=1-T8, T8'=51480 X^7(1-X)^7.

The full fifteenth-degree polynomial is

    T8=6435X^8-40040X^9+108108X^10-163800X^11
       +150150X^12-83160X^13+25740X^14-3432X^15.

At the vacuum, q8-1 is divisible by X^8. The new field-map
difference starts at field degree nineteen, so the named
vacuum/H/Phi2/Phi4 data through two loops are unchanged,
including the full dimensional action, counterterms,
Jacobian and physical source. This is not a result for
an independently prescribed ordinary-Psi composite source.

At X=1, q8 and its first seven derivatives vanish.
For the covariant clock Psi=sqrt(kappa)t, Fhat=Psi.
Its first seven field/metric variations agree with the
identity map there. The eighth variation is nonzero.

This aligns the protected mass profile's clock argument
and resolves the specific naive source-map derivative
screen exposed in S6.163. It is a kinematic statement:
the map does not make the canonical polynomial parent
solve the target bounce background equations.

## Full classical action comparison on the same restricted class

The map is polynomial, so its Fourier support remains
finite. On the original real Schwartz/Fourier unit class,
the old map support is three, the new map support is
thirty-three, and the new quadratic heavy source support
is sixty-six. The full heavy resolvent therefore has
the strict denominator margin M-4356.

This finite support is important. An infinite analytic
gate cannot simply inherit this same momentum-domain
bound.

The original cubic norm coefficient satisfies C_R<1e-403.
The actual gate and map-difference estimates are bounded
by the deliberately looser exact allowances

    |q8-1| allowance 1e-6390,
    map-difference allowance delta=1e-6793.

The unexpanded free/local/heavy action comparison is bounded
by delta times

    300+L+2g/(M-4356).

Its coefficient is below 1e-6790. Adding it to the previous
full finite-kappa analytic-target estimate preserves

    |S_polynomial[Fhat(Psi),H_stationary]-S_full_target[Psi]|
       <1e-800 ||Psi||_2^2

on that SAME restricted class. The proof retains the full
heavy inverse and integrates powers of the field norm,
not a constant over infinite spacetime. It does not bound
the action or stress on the rolling clock domain.

The dimensional gate argument uses the required mu^(2epsilon)
factor in X_D before finite parts. No independent finite
reference shift is hidden in the new coordinate.

## Verification and negative controls

The native report pins eighteen source/proof/test files
and twenty fields: 61 named identities, 61 scalar entries,
25 proof gates, 9 controls and 180 rejected inputs.
Final private science passed 301 tests in 23.16 seconds;
fresh repository science passed 301 in 23.20 seconds.
Ordinary replay passed 326 tests in 2282.51 seconds.
Independent native CLI replay passed. The full captured
P8 snapshot passed 26952 tests in 3259.75 seconds with
final exit code 0.

All 581 captured test files were present and unchanged.
Path-list SHA-256:

    dbdc8fa8f90880f60b56b4a668123723eb2455e90b44bd3e005452e402d74716

The full-run adapter passed 128 original tuple comparisons.
Final counters: 34087 domain fallbacks, 7340 exact descents,
88 mixed fallbacks. This snapshot predates S6.165.
Only full regression uses the audited exact GCD adapter;
native, direct science, ordinary and CLI retain unmodified
SymPy with interpreter-only allowances.

The 155298-character native report was transferred
losslessly in thirteen chunks. All source hashes were
independently verified. Report SHA-256:

    2684caa2dd2764c3fcc525c31b571188b323ca8c2cc81da0ce8eb9b2e69bd1b0

Independent finite-Fourier calculations evaluate the full
stationary action with exact rational arithmetic. A control
at M=4356 detects the new source's genuine frequency-66
pole even though the old source remains below threshold.
Other controls retain the nonzero eighth clock variation,
mixed field/metric variations, reflection identity and
higher positive powers of the norm allowance. Periodic
diagnostics are not called Schwartz functions.

Cleanup and formatting preceded final private verification
and immutable source freeze. Written proofs and exact
replay are not formalization or independent peer review.
The exact 22-file publication excludes later continuations
and unrelated P4/P9 changes.

## Remaining work

S6.165's specified free flat in/out state and S6.166's
Hadamard-theorem application have passed native and fresh
direct science; independent replays are running. The next
private derivation seeks a sharper exact transition bound
and a local in/out energy DIFFERENCE, not absolute stress.

The new coordinate and restricted action comparison do
not establish a common rolling parent, an interacting
curved state, a cutoff, quantum target matching or physical
finite-EFT errors. Global V contours/cuts, finite-gravity G
and common-parent B remain open. Scoped P8(a) and A.20-A.23
are unchanged; original P8(b) and original P8 are OPEN.
