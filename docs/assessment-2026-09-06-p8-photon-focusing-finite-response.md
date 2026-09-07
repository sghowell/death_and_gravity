# P8: history-based photon focusing and finite physical tensor response

Recorded 2026-09-06, following the
[photon-QSEI/asymptotic-response checkpoint](assessment-2026-09-06-p8-maxwell-asymptotic-response.md).
These are two new conditional theorems. Original P8 remains open; the
frozen problem statement and existing certificate ancestors are unchanged.

## What is now established

| Gate | New result | Retained boundary |
|---|---|---|
| [A.17](../problems/P8/a/fields/maxwell/focusing/FORMULATION.md) | Free-photon QSEI implies global-flat-FLRW timelike incompleteness without an initial pointwise SEC premise | Explicit contraction history, future geometry bounds conditional on extension, actual SEE and signed source budget |
| [S6.11.COMPOSITE](../problems/P8/s6/matching/composite/perturbations/cones/matching/asymptotics/response/FORMULATION.md) | A finite-positive-parameter physical TT probe disproves the specified locked-action response within a quantified relative-error contract | Selected composite parent, fixed interval and spatial band; no general matching, cutoff or original-row verdict |

Certificate SHA256 values:

    A.17: 22cfba547368b420767bbe87068c5b45fb556992c9e37a389091525f3e2ef62a
    S6.11.COMPOSITE: 34fa2ee35c2b33557e84b8fc181a202bed4793ed00d6aa8b76964ebbdc5dbe12

A.17 has 18 source hashes and 25 exact identities. S6.11 has 16 source
hashes, 28 exact identities and three separate Fraction polynomial
identities, plus exact interval-derivative records and comparison gates.
The independent physical/covariant audits contain ten and twelve tests
respectively. Their continuous-time, distributional and geometric
implications are written proofs, not proof-assistant formalizations.

Validation: **1,659 P8 tests passed in 482.64 seconds**, with
`PYTHONHASHSEED=0` and `faulthandler_timeout=60`. The two slow-test
diagnostic stack dumps came from unchanged S5 symbolic algebra, not
test failures. The two targeted suites pass 43 and 33 tests. Both
standalone read-only certificate replays, full P8 Ruff, all 34 new
source hashes, eleven Markdown files' local artifact links and diff
checks pass. Only the 39 verified checkpoint paths are included;
unfinished follow-ons and unrelated P4/P9 changes are excluded.

## A: the initial energy-sign premise is removed, not disguised

The actual field is the free Maxwell field of A.16, in an arbitrary
Hadamard state on a smooth global product `I_t times R^3`. Its complete
anomaly, finite renormalization parameter, cosmological constant and
additional-source dictionary are retained. No state homogeneity or
quasifree condition is added.

For matched rising/falling cubic samplers u and g on a short past and
proposed future segment, the exact FLRW square identity gives

    J[g]+K <= 3(||g'||^2+||u'||^2)+Q_Maxwell[u+g]
              -3 integral_past (u'-Hu)^2.

Here `u+g` denotes the joined sampler on their disjoint interiors, not
two simultaneously supported functions. The past is genuinely present
and has `H<=-h/tau`; it supplies a lower bound on the last square without
a pointwise Ricci/SEC assumption. The future bounds permit inverse-power
growth in `tau-t` and are imposed only on a normal which actually reaches
tau. The quadratic endpoint zeros make their sampler costs integrable.
Nothing is integrated through a singular or absent endpoint.

With `delta=kappa*hbar/(8*pi^2*tau^2)` and the explicit dimensionless
additional-source allowance sigma, the exact sufficient margin is

    3h+39h^2*r/35-18/5-delta*(C_P+C_F)-13*(1+r)*sigma/35.

If this is nonnegative, a smooth normal reaching tau would have both
`J+K<=0` and `J+K=3 integral_future(g'-Hg)^2>0`. The global FLRW time
coordinate then proves a finite upper bound on every future timelike
curve from the initial slice, without importing a compact-Cauchy theorem
onto noncompact spatial R^3. This is incompleteness of the stipulated
spacetime, not curvature blow-up or inextendibility in every larger one.

At `r=1/100`, `h=3/2`, past/future caps `(2,4,16,96)`, `|beta_M|<=1`
and `delta<=10^-8`, the quantum coefficient is below 13 million. The
zero-source margin exceeds 3/4 and sigma<=1 leaves more than 1/3. The
Planck dictionary makes the quantum restriction compatible with
macroscopic or cosmological durations. This dimensional calibration
does not verify an observed universe or classical radiation history.

A smooth mollified-quintic Hubble history meets the geometric bounds,
violates the initial pointwise SEC and becomes static in the future,
so its future timelike and null geodesics are complete. It cannot satisfy
the permitted small-quantum/small-source Maxwell SEE. This independently
checked comparison shows why the quantum/source input is load-bearing;
it is not presented as an actual allowed state/solution witness.

## B: a physical source with a finite error bound

Use the unchanged composite beta2 parent and physical proper time
`u=mT`. The prescribed spatial band is `kcom/m in [1/16,1/8]`; the pulse
`p=u^2(1-u)^2` is C1 when extended by zero. It has no temporal bandlimit.
The probe is a separately conserved external TT stress, not a TT mode
silently assigned to the internal homogeneous canonical scalar.

The literal u-action source is `J_u=m*M^2*epsilon*sigma*p/4`, and its
physical stress is `Pi=-2m*J_u/Ae^3`. Both time/volume factors and both
physical composite projection weights are checked independently. The
initial background root is exactly `sqrt(1+3*(1+epsilon)^3)`; replacing
it by its limit 2 would change the pinned finite-member density.

In the limit the normalized response kernel is a sum of two causal
kernels, `G_L+sqrt(v(u))*G_H*sqrt(v(s))`. The relative contribution is
positive and at least the light contribution for this nonnegative
pulse. Its positivity is proved for the evolving equation, not from a
frozen mass sign. All canonical-normalization derivatives are retained.

Exact outward-rational interval differentiation bounds the full finite
background and flux equations. Continuation and Duhamel estimates give
uniform normalized output errors below 1/600 for
`0<epsilon<=2^(-187116)`. At u=1 the full-minus-locked gap is at least
1/75 and the relative omission is at least one third. A separate
`sigma<=2^(-67162)` gate keeps each individual metric amplitude at most
1/128. These are deliberately unoptimized sufficient constants, stored
as exact binary exponents. The absolute probe and composite response
vanish with epsilon*sigma; the normalized mismatch does not.

Prepared source-free equal-metric value/velocity data provide a positive
control: full and locked predictions agree within 1/600 under the same
parameter gate. Thus the failure concerns the stated physical source,
initial data and response contract, not every use of a light-only model.

The tiny parameter range is not a cutoff or weak-curvature certificate.
The singular two-metric endpoint, fixed time window, internal reconstructed
potential and untested scalar/vector sectors remain explicit. A different
parent or a reduced theory with a different operator, memory or source
map needs its own test.

## Original closure obligations

The original P8(a) asks for realistic-field incompleteness with constants
usable at cosmological scales. A.17 is a real advance toward that target:
it supplies the field-specific implication and removes the old initial
pointwise energy sign. The remaining assessment must test the geometric
and source hypotheses in a meaningful cosmological regime, without
mistaking Planck-scale arithmetic for an application. A new exact SEE
solution, interacting QED, or a theorem of curvature inextendibility is
not automatically an extra prerequisite to every singularity theorem.

For P8(b), the frozen 32-row linear classification and its matter/speed
contract remain complete. The selected UV branch still needs controlled
vacuum-to-bounce matching before vacuum positivity can decide surviving
rows. Neither a failed selected parent nor a healthy vacuum alone
resolves that question. Full UV completion is not substituted for the
weaker, original positivity-consistency obligation.

These unresolved items are research obligations, not a missing user
choice, credential or approval. The project is not finished or closed.
