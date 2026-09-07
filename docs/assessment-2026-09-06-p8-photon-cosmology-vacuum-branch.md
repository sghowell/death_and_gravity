# P8: cosmological-strength photon theorem and a direct matching obstruction

Recorded 2026-09-06, after the
[history-focusing/finite-response checkpoint](assessment-2026-09-06-p8-photon-focusing-finite-response.md).
The original problem statement and all ancestor certificates remain
unchanged. P8 as a whole is not closed.

## New results and their scope

| Gate | Established | Not inferred |
|---|---|---|
| [A.18](../problems/P8/a/fields/maxwell/focusing/cosmology/FORMULATION.md) | A photon-QSEI timelike incompleteness theorem with robust radiation/matter-era history and cosmological-strength constants | Observed-universe applicability, unrestricted fields/spacetimes, or a classical comparison being a quantum solution |
| [S6.12.COMPOSITE](../problems/P8/s6/matching/composite/vacuum/FORMULATION.md) | Exact tree-level symmetric vacuum branch and direct CD/M1 operator/background matching obstruction | All-branch, all-parent, loop, cutoff or UV-positivity classification |

## A: a nontrivial cosmological calibration

The [A.16–A.18 chain](../problems/P8/a/fields/maxwell/focusing/cosmology/notes/closure.md)
now supplies a physical photon inequality, a timelike incompleteness
implication and a cosmological-strength specialization. The source is an
actual free Maxwell Hadamard state with its explicit finite prescription;
the geometry is a stipulated smooth global spatially flat FLRW product.

The actual short history lies in a positive-radius C3 neighborhood of
time-reversed power laws with exponents between radiation and dust. It
has a quantified contraction rate, with no initial pointwise energy
sign inserted as a separate premise. This calibrated neighborhood does
happen to satisfy initial timelike convergence; it is not advertised
as an SEC-violating example. No future pointwise SEC is imposed.

An important independent audit rejected the first attempted future
envelope. Reusing the tight past bounds would already force a Taylor
contradiction by time tau/20, without the QSEI. That was not a
nontrivial cosmological quantum-energy calibration.

The final proof instead separates past and future caps:

    past d=(21/10,9,70,800),
    conditional future c=(4,128,16384,1048576).

A compact seventh-degree smoothstep followed by positive mollification
constructs, for every reference exponent, a smooth future-complete
metric inside the past C3 tube and obeying the final future caps.
It becomes static at positive scale factor. Thus these geometric
assumptions alone do not imply incompleteness. This comparison is not
an allowed small-source Maxwell SEE solution; the quantum/source
inequality is essential to the contradiction.

The exact affine quantum cost is `C0+Cbeta*abs(beta_M)`, with each
coefficient below 13 million. Under

    delta*(1+abs(beta_M))<=10^-8,
    sigma=tau^2*(Lambda_++kappa*max(-ell,0))<=5,

the strict focusing margin is at least `47079/350000>1/8`.
No finite prescription is silently set to zero. The proper-time
requirement is `tau/t_P>=10^4*sqrt((1+abs(beta_M))/pi)`, in seconds.
The cosmological-constant bound is stated independently of this
Planck ratio and retains the actual versus reference Hubble-rate
difference in the unanchored tube.

The future caps are conditional on extension and remain hypotheses,
not conclusions inferred from the short history. A hypothetical smooth
normal reaching tau contradicts the QSEI/index identity. The result
is timelike incompleteness of the specified spacetime, not curvature
blow-up or inextendibility in every larger one.

This reaches a precisely stated photon/global-flat-FLRW version of
P8(a)'s realistic-field and cosmological-strength objective. Its root
status is **SCOPED OBJECTIVE COMPLETE (CERTIFIED)**: free-Maxwell
timelike incompleteness in global spatially flat FLRW, with an
all-Hadamard QSEI, explicit geometric history/conditional-extension
and finite-prescription/source hypotheses, and robust
cosmological-strength constants. This qualified status supersedes
the earlier checkpoint-time open assessments; no frozen report is
rewritten and no unrestricted theorem is claimed. A new
exact SEE solution, interacting QED or observed-universe reconstruction
is not silently added as a prerequisite to this conditional theorem.
Conversely, the theorem is not promoted to arbitrary spacetimes or
fields. Any additional ordinary matter needs its own explicit energy
assumptions; the pure-Maxwell statement has no such additional source.
It is this pure-Maxwell specialization that has only a QSEI as its
energy input. Neither a generic initial SEC assumption nor an
unproved energy inequality for a different field is inserted.

## B: genuine relative-field elimination in the physical matter frame

The new audit keeps `mathcalG=g_eff` physical and defines the relative
matrix by the Cayley transform of the positive square root. The exact
inverse is

    g=mathcalG*(I+Delta)^2/4,
    f=mathcalG*(I-Delta)^2/4.

Equal Einstein coefficients and exchange-symmetric potential coefficients
make the full action even in Delta. Every prescribed matter/probe
coupling depending on mathcalG is independent of Delta. Consequently
Delta=0 is an exact off-shell stationary branch, including dynamical
light sources and all derivatives of the metric map.

For the specified beta2 potential the branch action is exactly

    integral sqrt|mathcalG|[-M^2*R/4-3*M^2*m^2/8+L_light].

With canonical first-derivative input matter it induces no C/D DHOST
operators at tree level. This is a branch statement, not a loop
cancellation or proof that higher operators initially inserted into
L_light disappear. Relative-field loops and other stationary states
remain possible.

A proportional Minkowski candidate requires the same light potential
to have value `-3*M^2*m^2/8`, a stationary point and a nonnegative
canonical mass Hessian. The constant vacuum source shifts all five
Hassan–Rosen coefficients. Including those shifts and the physical
clock gives `M_eff^2=M^2/2` and `m_FP^2=m^2/4`. Omitting either
dictionary gives the wrong mass. The physical source has no massive
residue on this branch, although the massive five-polarization sector
has positive kinetic normalization and nonzero mass.

The pinned composite family has one internal scalar; CD/M1 has a clock
and an additional free rolling scalar. A separately named canonical-light
extension supplies the correct light content without relabeling the old
one-scalar background. On its zero-relative branch the old CD/M1
coefficient gaps at the bounce are `abs(F2_X)=1/2` and `abs(A3)=1` in
the pinned normalization. A claimed correcting remainder must account
for those gaps in the specified physical metric and clock dictionary.

There is also a physical-background obstruction. Canonical NEC matter
makes the branch Hubble rate nonincreasing, whereas the CD target has
opposite endpoint values `+-8/(5*tau)`. The maximum endpoint error
therefore cannot be below `8/(5*tau)`. Repairing the exact CD null
equation while retaining its exact free M1 charge needs residual
magnitude at least `(801/100)*M_eff^2/tau^2` at the bounce.

The connected stationary-branch theorem explicitly assumes a common
boundary/state prescription and an invertible full constrained
functional derivative. A positive mass or algebraic potential Hessian
is not that inverse. No rolling inverse norm is claimed. S6 requires
a common matching construction, not a physical-time trajectory from
vacuum to bounce; the proof does not add the latter requirement.

## What remains before original P8 closure

The 32-row linear-principal DHOST classification, its exact matter
coupling and its subluminality contract remain complete. The new
vacuum result rules out one source-preserving matching route directly;
it does not turn that failure into an exclusion of the surviving rows.
The adopted UV branch still needs applicable, controlled vacuum/bounce
matching and finite-gravity positivity information for a common model,
or a theorem excluding the relevant remaining class under its stated
hypotheses. Healthy vacuum data alone are insufficient.

No full UV completion, generic nonlinear stability theorem or interacting
Standard Model result is substituted for the original obligations.
These open items remain research work, not a missing credential,
approval or unanswered user preference. Original P8 is not finished.

## Frozen replay inventory

| Report | SHA256 | New hashed files | Exact symbolic residuals | Targeted tests |
|---|---|---:|---:|---:|
| [A.18](../problems/P8/a/fields/maxwell/focusing/cosmology/certificates/cosmological-calibration.json) | `b03767fb43c49a0f6c5ff34307355c4cd37a75930f61e07f4826c06101e4d23d` | 20 | 35 | 61 |
| [S6.12.COMPOSITE](../problems/P8/s6/matching/composite/vacuum/certificates/composite-vacuum.json) | `b06bfa4c45aa4e7ad5b0bd4a2c5c928d9616534b69e4c765a2d752966c78462e` | 15 | 33 | 33 |

Both standalone read-only `--check` commands pass, including all pinned
ancestor replays. A.18 reconstructs the full core data with stdlib
Fraction arithmetic; S6.12 additionally checks four coefficientwise
Fraction polynomial identities. The separately authored root audit
suites contain 13 and 11 tests respectively. All 35 new source hashes
were independently checked, as were local links in the 12 Markdown
files in this checkpoint. No source-hashed ancestor was edited.

The checkpoint-wide regression passes **1,753 tests in 419.59 seconds**
using the [opt-in exact GCD runner](p8-exact-regression-runner.md):

```sh
PYTHONHASHSEED=0 .venv/bin/python scripts/p8_exact_regression.py problems/P8 -q -o faulthandler_timeout=60 --ignore=problems/P8/a/fields/maxwell/thermal --ignore=problems/P8/s6/matching/trimetric
```

The two ignores exclude subsequent children, not tests belonging to this
checkpoint. The adapter made 6,509 exact univariate real-domain descents
with both cofactor identities checked on every call, after 128 comparisons
against the untouched Gaussian GCD implementation. Its ten independent
boundary/restoration tests pass. This is an adapted exact-arithmetic run;
the ordinary combined runs were stopped after repeated legacy SymPy
slowdowns and are not reported as passes. The ordinary isolated legacy
quadratic bridge, both new targeted suites and both standalone certificate
CLIs pass separately. P8-wide Ruff and `git diff --check` also pass.
