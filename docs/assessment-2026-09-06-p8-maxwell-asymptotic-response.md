# P8 after the photon QSEI and large-asymmetry response gates

Recorded 2026-09-06, following the
[actual-extension/light-reduction checkpoint](assessment-2026-09-06-p8-extended-qsei-light-reduction.md).
Original P8 remains open. The new results address a realistic free field
and an actual time-dependent parent response, respectively; neither
completes the corresponding original research question.

## New results

| Gate | Established | Boundary |
|---|---|---|
| [A.16](../problems/P8/a/fields/maxwell/FORMULATION.md) | Exact all-Hadamard photon QSEI on smooth flat-FLRW conformal strips, with complete anomaly/prescription and source dictionaries | No removed initial pointwise premise, new SEE solution or cosmological incompleteness theorem |
| [S6.10.COMPOSITE](../problems/P8/s6/matching/composite/perturbations/cones/matching/asymptotics/FORMULATION.md) | Actual limiting homogeneous tensor growth, physical projection and fixed-forward-interval convergence | Singular normalization limit; no uniform finite-amplitude, finite-band, long-window or general matching verdict |

Certificate SHA256 values:

    A.16: f59d63524bf1d98998a928fbe5139eef45f3703ebd1e61cbdf470e738ce94b94
    S6.10.COMPOSITE: bfb0b74a2d50a508a5e5fc346b59a31958c9f9a8ccfad621d4caa34883e90eef

A.16 hashes 19 source files and checks 27 exact identities against a
separate Fraction-only curvature and envelope engine. S6.10 hashes 13
sources and checks 43 symbolic identities plus six Fraction identities.
Targeted suites pass 57 and 28 tests respectively, including 10 and 9
independent covariant/physical audits. Both standalone report replays
pass. The distributional, positivity and continuous-time implications
remain written proofs, not proof-assistant formalizations.

Validation: **1,583 P8 tests passed in 415.09 seconds**, using
`PYTHONHASHSEED=0` and `faulthandler_timeout=60`. Both diagnostic
slow-test stack dumps came from unchanged S5 symbolic algebra; neither
was a test failure. Full P8 Ruff, both standalone certificate replays,
32 new source hashes, local Markdown artifact links and diff checks
pass. Neither an unfinished follow-on nor an unrelated P4/P9 change is
part of this checkpoint.

## A: an exact photon inequality with its physical choices visible

This is a separate free-Maxwell model, not a replacement of the previous
minimal-scalar field or its finite prescription. Four-dimensional
conformal invariance transports the restricted flat vacuum to an actual
Hadamard reference. The positive-type proof works locally on the strip
for arbitrary target Hadamard states; no target global extension,
homogeneity, quasifree or semiclassical-equation hypothesis is needed.

The exact proper-time operator is

    L_H f=f''-2H f'+(3H^2/4-3Hdot/2)f.

The lower bound is the conformal-reference EED minus
`hbar/(8pi^2)*norm(L_H f)^2`. Its two photon polarizations, conformal
measure and H2 endpoint conditions are checked independently. There is
no scalar scattering-history duration restriction and no optimality
claim. The reference stress is
`hbar/(2880pi^2)*(62H3_FK+beta_M I_FK)`; beta_M remains a required
explicit input to numerical curvature envelopes. Zero type D is only a
separately named specialization, not an inherited scalar choice.

Energy density is not EED when the trace anomaly is nonzero. At de Sitter
the reference density is positive and its EED negative. On the radiation
metric, R=0 removes the finite-beta term but not the anomaly: the exact
absolute IBP coefficient is `4601/(1280t^4)`, with the negative first-
derivative term retained. The radiation metric is a control, not a new
solution of the photon SEE.

The source audit uses the flat Maxwell QEI of
[Fewster–Pfenning](https://arxiv.org/pdf/gr-qc/0303106) and
[Pfenning](https://arxiv.org/pdf/gr-qc/0107075), with the conformal anomaly
law of [Herzog–Huang](https://arxiv.org/pdf/1301.5002). A literal derivative-
tensor formula in [Markowicz et al.](https://arxiv.org/pdf/2006.06884)
does not pass the stated FLRW conservation conversion and is explicitly
not imported; the finite-beta tensor is independently defined by the
conserved FK inverse-metric variation. The precise source-use boundary,
not a broad judgement about that paper, is part of the source-hashed proof.

The conditional Einstein dictionary retains Lambda and every additional
source with its own sign and bound. Under dimensionless curvature caps,
its quantum constants scale as `kappa*hbar/tau^2`. This arithmetic does
not supply a focusing example: the naive `Hmax*tau<=1` cap permits
`|K|tau<=3`, below A1's zero-initial-credit gradient threshold of at least 4.
The initial pointwise Ricci premise, extended normal domain and sufficient
contraction still require a separate theorem or verified application.

## B: a full limiting equation, not just a center coefficient

In the specified beta2 family, set epsilon=1/y0 and use physical time
u=mT. The exact scaled ODE has the forward limiting solution

    z=12 cosh(u)+2 sinh(u)-11,
    v=z'/z, j=11u/z, R=v-j, eta=3(R^2-1).

For u>=0, z>=1, R>0 and `1<=v^2<=140/19<9`. The density may become
negative; the finite-parameter canonical null source remains positive.
The required potential is reconstructed separately for each finite
member and is not the old free M1 action.

The exact zero-charge homogeneous relative equation retains its full
normalization curvature. Its coefficient `D_H=N_R-m_alg^2` is at least
1/4 along the whole limiting orbit. The positive Green comparison gives
`H(u)>=cosh(u/2)` for the normalized data H(0)=1,H'(0)=0. This is an
actual time-dependent solution bound, not an inference from a frozen
mass sign. Zero heavy data still give zero homogeneous response.

The physical composite projection tends to `-2sqrt(v)H/M` and has a
nonzero tidal jet. On every fixed finite forward interval, a compact
positive-root tube and Gronwall estimates give regular sufficiently
small positive-parameter backgrounds and uniform normalized-response
convergence. No numerical parameter threshold or interchange with a
growing time interval is asserted.

The relative canonical normalization is O(epsilon). Unit canonical data
therefore do not represent a uniformly small two-metric perturbation.
Multiplying the column by epsilon times a small linear amplitude keeps
the individual metrics small on a fixed interval, but sends its absolute
composite amplitude to zero. Normalized growth ratios survive. This is
not a finite-amplitude nonlinear instability or a general EFT exclusion.

## Next obligations

On A, the useful next question is whether the photon inequality yields
a focusing criterion without A1's initial pointwise SEC premise, while
retaining any needed geometric-history or continuation assumptions
explicitly. A known anomaly-driven exact cosmology, or a small Planck
ratio by itself, is not a replacement for that implication. Interacting
fields remain beyond the new free-photon result.

On B, a prepared light initial state and a physical composite-metric
source need not excite the relative tensor in the same proportion.
A nonzero-momentum, retarded source-response calculation must retain
both normalizations and control finite-parameter errors before it can
support a finite-band matching verdict. Such a result would still be
about a stated parent family, not all possible UV parents of the original
DHOST witnesses.

The original task requires realistic-field cosmological incompleteness
and a stated-basis bounce classification with controlled matter and UV
positivity applicability. Full UV completion is not being substituted
as a stronger closure requirement, but a healthy vacuum or a failed
selected parent alone does not decide the surviving bounce rows.
No missing credential, approval or unresolved user choice currently
blocks the authorized calculations. P8 is not finished or closed.
