# P8: protected Yukawa mass profile with unchanged named two-loop data

Original P8 remains OPEN. This follows the
[full analytic-target classical match](assessment-2026-09-10-p8-full-analytic-target-classical-match.md).
No user intervention is required for the current research.

## A separately named finite-EFT candidate

[S6.163](../problems/P8/s6/continuation/s6_163/FORMULATION.md)
introduces GY14-SAT8. The two active opposite Yukawa masses use

    M_+/- = mF +/- y f_R(Phi),
    f_R(z)=z/[1+(z/R)^8]^(1/8), R=10^300.

The other twelve Dirac flavors are unchanged. This is a separately
named nonpolynomial EFT, not the previous globally renormalizable
model. Its positive real root is smooth and real analytic for all
real arguments; it is not entire in the complex field plane.

For real z, |f_R(z)|<R, 0<f_R'(z)<=1 and
|f_R''(z)|<=9/R. The actual fixed parameters give

    mF=10^200, |y|<3e-103,
    0.99mF < M_+/-(z) < 1.01mF

uniformly over real z. This is a pointwise mass enclosure, not a
propagator, state or interacting spectral-gap theorem.

## Named two-loop references remain unchanged

The first new interaction is proportional to Phi^9 times a fermion
bilinear. Its weighted field excess, together with the connected
loop identity, excludes it from the named vacuum, H-source,
two-Phi and four-Phi coefficients through two loops. The proof
retains induced counterterm forests and map-field degree; it does
not discard self-contractions. The dimensional continuation uses
R_D=mu^(-epsilon)R and y_D=mu^epsilon y consistently.

This statement is deliberately not about every external field
count at two loops. Independent determinant and finite Gaussian
controls show a nonzero ten-Phi change already at one loop.
The lower-field finite references are fixed, not silently retuned.

At zero classical fermion field the changed interaction vanishes
identically, so the preceding restricted classical full-target
action comparison is unchanged. No global marginal UV-flow claim
is transferred to the nonpolynomial EFT.

## A clock-argument problem is exposed, not hidden

For an assumed scalar argument z(t), the derivative bounds imply

    |Mdot|<=|y| |zdot|,
    |Mddot|<=|y| (|zddot|+9|zdot|^2/R).

The diagnostic z=sqrt(kappa)t has small corresponding mass
variation ratios. But the original literal derivative source map
does not give that argument on the flat clock. Its slope at zero
is approximately 2e200. Using the actual Yukawa lower bound, the
naive mapped-clock first-derivative ratio exceeds 1e97.

That is a failure of this small-derivative screen for this map.
It is not a particle-production proof, a whole-row exclusion or
a failure of every possible common parent. The separately frozen
S6.164 polynomial clock-transparent map addresses this argument
mismatch kinematically; it does not supply parent background
equations or a state.

## Verification

The immutable native report pins 18 source/proof/test files and
20 fields: 52 named identities, 52 scalar entries, 29 proof gates,
9 controls and 224 rejected inputs. Final private science passed
372 tests in 22.97 seconds; fresh repository science passed 372
in 22.95 seconds. Ordinary replay passed 397 tests in 2331.15
seconds. Independent native CLI replay passed. The complete
captured P8 snapshot passed 26626 tests in 3361.75 seconds, with
final exit code 0.

All 579 captured test files were present and unchanged.
Path-list SHA-256:

    03532d58ef7844242487efd7d0d0beda700675009ca9392fac6e4681e1385035

The full-run exact GCD adapter passed 128 original tuple
comparisons. Its final counters were 34093 domain fallbacks,
7340 exact descents and 88 mixed fallbacks. This snapshot
predates S6.164. Only full regression uses that audited adapter;
native, direct science, ordinary and CLI use unmodified SymPy.

The 72880-character report was transferred losslessly in seven
chunks and all eighteen source hashes were independently checked.
Report SHA-256:

    56a114490c2334adcdfbfb14e1a2add42e7f7e769e51b024ff91141ec63074cb

Independent tests cover the complete real profile, complex branch
limits, loop/forest counting, nonzero high-field determinant
changes, the dimensional dictionary, and both direct and mapped
clock derivative controls. Cleanup preceded final private tests
and byte freeze. Exact replay and written proofs are not
formalization or independent peer review.

Publication stages exactly the 22 files belonging to this
checkpoint, excluding later continuations and unrelated P4/P9 work.

## Remaining work

S6.164 and the free flat-clock Dirac production estimate S6.165
have passed native and direct science; their fresh independent
replays are running. The next private calculation checks
all-order switching tails and short-distance state regularity.

The pointwise protected mass is not the interacting B-state
dictionary. A common rolling parent, curved/local stress,
cutoff and physical finite-EFT errors, global V contours/cuts
and finite-gravity G remain open. Scoped P8(a) and A.20-A.23
are unchanged; original P8(b) and original P8 are not closed.
