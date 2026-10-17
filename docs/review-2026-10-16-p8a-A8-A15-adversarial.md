# Independent adversarial review: P8(a) rows A.8–A.15

Date: 2026-10-16. Reviewer: a separate agent with no part in authoring
P8(a). It was given the code, the ledger rows and the level definitions,
but not the earlier triage in
[review-2026-10-16-p8a-adversarial.md](review-2026-10-16-p8a-adversarial.md).
The review was read-only; an mtime scan shows nothing was written to the
repository. Base path below:
R = problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual.

**Overall verdict:** A.8–A.15 are internally correct but physically
confined. They are off-shell perturbative estimates plus exact SEE
solutions lasting at most about 0.74 Planck time, bounded by the classic
Planck-rate runaway, with after-the-fact flat-space-strength QSEIs on one
constructed spacetime and every reported focusing test failing. They do
not move P8(a) toward a cosmological incompleteness theorem.

## 1. Replays

All runs used `PYTHONDONTWRITEBYTECODE=1` and `-p no:cacheprovider`.

- All eight `--check` commands pass (exit 0).
- pytest passes for six suites: A.8 63, A.9 44, A.10 74, A.12 54, A.13 40,
  A.14 45.
- **A.11 and A.15 fail to load in a worktree** with a pytest collection
  error. Their conftest files put every `src` folder under problems/P8
  (382 of them) on the import path:
  `R/existence/preparation/tests/conftest.py:5-8` and
  `R/existence/preparation/continuation/extension/qsei/tests/conftest.py:4-7`.
  In the worktree that reaches a 1045-character path and raises
  `OSError [Errno 63] File name too long`. From the main checkout, whose
  path is 32 characters shorter, they pass: 75 and 37. Whether they replay
  depends on the checkout path length.

## 2. Cross-cutting findings

1. **Durations are Planck-scale, and no row or note says so.**
   - With κħ = 8π t_P² and δ = 10⁻¹⁴: t_star = 1.49×10⁵ t_P and
     T0 = 2.97×10⁵ t_P.
   - The SEE domains are about 3.7×10⁻⁵ t_P for A.11/A.12 and about
     0.74 t_P for A.14/A.15.
   - A.13's comparator durations 10⁻⁵ and 10⁻⁴ are about 7.4 and 74 t_P.
   - A.13's pole p0 = 1.231×10⁶ is an e-folding time of 0.60 t_P: the
     classic Planck-rate runaway of the SEE, generally treated as outside
     semiclassical validity. A.14 spans only about 1.2 of its e-folds
     (e^{p0L} = 3.4).
   - References: Horowitz PRD 21, 1445 (1980); Suen PRL 62, 2217 (1989);
     Simon PRD 43, 3308 (1991); Parker–Simon PRD 47, 1339 (1993);
     Flanagan–Wald PRD 54, 6233 (1996); Anderson–Molina-París–Mottola PRD
     67, 024026 (2003).
2. **The QSEIs point the wrong way for a singularity theorem.** A.12 and
   A.15 fix one constructed metric and quantify over target states. An
   incompleteness theorem needs an inequality valid for unknown SEE
   solutions, and the SEE state's own effective energy density is already
   fixed by the geometry. Both constants are the Minkowski value
   ħ/(16π²) times root², with roots 1+2.7×10⁻⁹ and 1+1.1×10⁻⁶, rounded up
   to 2.
3. **The certificates cannot fail on their decisive conclusions.**
   `--check` only compares regenerated JSON with the stored JSON, so a
   literal value always passes.
   - A.14 `verify.py:77-80,90` hard-codes `True` for same data and
     constraint, Hadamard state, full SEE, smoothness, and agreement with
     A.11.
   - A.13 `verify.py:84` hard-codes `full_kernel_nonnegative: True`.
   - Negative controls that are literals or prose strings: A.9
     `verify.py:106-111`, A.12 `scattering.py:103-105`, A.13
     `resolvent.py:95-98`, A.14 `verify.py:94-99`.
   - Tautologies counted as identities: A.14, 6 of 11
     (`regularity.py:23-25`: `exp(σL)exp(−σL)−1`, `n−n**1`, distributivity;
     `weighted.py:57`; `map.py:39-42`); A.11, 8 of 18
     (`regularity.py:62-65`, including `(u′v+uv′)−(u′v+uv′)`;
     `actual_map.py:91-93,105-107`). A.9 counts the Leibniz rule four times
     (`mode_bounds.py:58-61`). A.10's `phase_derivative_transfer` is the
     product rule, not the Dyson transfer it is named after.
   - By a rough classification of test names and bodies, about 218 of 420
     test items only mutate reports or hashes, or reject bad input types
     and domains.
4. **The mathematics checks out.** The reviewer re-derived the steps listed
   per row and confirmed the A.11 and A.13 numbers independently with
   mpmath.
5. **The hypothesis classes are non-empty.** A.11 and A.14 construct their
   solutions. Nothing relies on content A.1–A.7 do not state, but A.8–A.11
   inherit A.7's written finite-amplitude constants (b_j, C2).

## 3. Verdict and level per row

| Row | Verdict | Recommended level |
|---|---|---|
| A.8 | GAP (presentation) | CERTIFIED, conditional on A.7's constants |
| A.9 | PROVEN AS STATED | WRITTEN_PROOF |
| A.10 | PROVEN AS STATED | WRITTEN_PROOF |
| A.11 | GAP (presentation) | WRITTEN_PROOF |
| A.12 | GAP (presentation) | WRITTEN_PROOF |
| A.13 | PROVEN AS STATED | CERTIFIED |
| A.14 | GAP (presentation; the certificate overstates) | WRITTEN_PROOF |
| A.15 | GAP (presentation) | WRITTEN_PROOF |

- **A.8.** The reconstruction identities are symbolic, the bounds rational,
  and the kernel lemmas routine. Verified: the anomaly coefficients
  R²/(1152π²) and (Riem²−Ric²)/(2880π²), the Euler term
  ρ = ħH⁴/(960π²) (equal to A.3), and the 1922 bound. "<10⁻¹⁷ of classical
  scales" hides that the residual is up to about 1% of the quantum source
  itself. The renormalization length λ = 4√2 t_star moves with the
  background, and on FLRW a λ-shift equals a γ-shift, so the ε = 1 clause
  describes a family of schemes, not one theory.
- **A.9.** Re-derived the I0 formula, the 3/2 constant, f0, f1, f2,
  C0root < 2.109, the UV moment 1/3 and the full-Parseval IR bound. The
  metric is off shell.
- **A.10.** Checked the derivative transfer (which needs U(0)=0), the R₂
  kernels, the bound 18z²e^{2z}, the shared-history term, the keyhole
  residue and cut sign, and the inverse norm.
- **A.11.** Verified the u/(60δ) trace normalization,
  d(a) = −19/60 − log(a/2)/2, the constraint propagation C′ = a⁴hD_trace,
  and the center and pair constants. The exact inverse norm is 0.061
  against the bound 0.21; the regularity induction is sound. But the
  unforced SEE holds on only about 3.7×10⁻⁵ t_P; the Einstein term
  contributes only about 7.5×10⁻⁸ of the contraction constant (< 0.567),
  which comes from the log/local operator; and ‖X*‖ < 7.2×10⁻⁸ still allows
  u′ to be about 10³ times the prepared metric's |ū′|. It is an
  explicit-constant, massless instance of the Meda–Pinamonti–Siemssen
  fixed-point scheme (Ann. Henri Poincaré 22, 3965 (2021), proved there for
  massive fields).
- **A.12.** Verified the forward/backward split, the Parseval factors and
  the reference numerator ≥ δ/54. The index failure is trivial:
  J·T0 ≳ 2×10¹⁰ against the ≤ 3/4 needed.
- **A.13.** The pole count and positivity are routine and the brackets
  rational. Independent checks: the damped pair has Re p = −1.21×10⁶ with
  Im w = 6.07 < 2π; exact unweighted norms are 0.24, 1.56×10⁴ and
  2.1×10⁵² at durations 10⁻⁶, 10⁻⁵ and 10⁻⁴; the weighted resolvent norm
  R_f(σ) = 0.222 < 240/769.
- **A.14.** Reproduced the weighted gains (α, α², 18α³), KG ≈ 8.34×10⁻⁶
  and η < 9×10⁻⁸. The e^{σL} cost of returning to unweighted norms caps
  the method, and its a priori bounds cannot fix the sign of the quantum
  effective energy density.
- **A.15.** Verified the N−N̄ decomposition, the loss of −0.0237, the
  (2/5)L⁴ penalty and root < 1+2×10⁻⁶. The reference effective energy
  density may be negative at about 1.2×10⁴ times the A.3 magnitude.

## 4. Wording fixes recommended

- A.9–A.12, A.14 and A.15: change CERTIFIED to WRITTEN_PROOF.
- Every row: state its durations in Planck units.
- A.8: add "residual ≲1% of the quantum source (<10⁻¹⁷ of classical
  terms); λ∝t_star; threshold t_star≈1.5×10⁵ t_P".
- A.11: add "L=10⁻¹⁰ η_star; unforced SEE on ≈3.7×10⁻⁵ t_P; the Einstein
  term contributes 7.5×10⁻⁸ to the contraction".
- A.12: add "the coefficient 2 is the Minkowski value rounded up; no
  focusing step uses this QSEI".
- A.13: add "comparator durations ≈7.4 and 74 t_P; the positive pole is
  the known Planck-rate runaway (e-folding 0.60 t_P)".
- A.14: write "extend to L=10⁻⁶ η_star (source-free ≈0.74 t_P, ≈1.2
  runaway e-folds)", and drop the literal-True entries or label them as
  prose conclusions.
- A.15: replace "improving to at least 2/5" with "Q2/τ²≥2/5, versus 10⁻⁶
  in A.1's calibration (≈560 t_P would be needed); the reference effective
  energy density may be negative".
