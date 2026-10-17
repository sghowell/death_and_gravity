# Independent adversarial review: P8(b) classification (rows P8-0 to P8-3)

Date: 2026-10-16. Reviewer: a separate agent with no part in authoring
P8(b). It read the formulation, notes S0–S3, the certificates and the
source before reading the authors' self-review. The review was read-only;
the repository was unchanged before and after.

**Overall verdict: the chain survives.** Every load-bearing coefficient
agrees with an independent covariant derivation done by a different
method. The no-go is sound, and holds under weaker hypotheses than
stated. All three witnesses replay exactly under that derivation, with
all-time sign proofs. What remains are presentation-level gaps.

## 1. Method

The reviewer wrote an independent SymPy derivation. It expands the full
covariant action √−g[F + K□φ + F2R + ΣAiLi], plus the free χ field, in
signature (+,−,−,−), using arbitrary fields α, β, ζ(t,x) rather than a
single Fourier mode. It compares the result with the claimed reduced action
using Euler operators, which ignores total-derivative differences. Its Σ
and Θ match the published KYY11 formulas for the cubic and quartic
Horndeski special cases.

## 2. Verdicts

**P8-0: PROVEN AS STATED.**
- The G_S and F_S constraint reduction re-derives, and
  det(K−G) = −Yq²(f−g)² holds.
- Ḣ(0) = (ε+3)/(6ετ²) = 1/375 and a1(0) = 13/(12cosh²(1/10)) re-derive.
- The raw A25 v2 HTML equation source shows (10g) D = 2φ̇(F2X−A1) and
  (39) Λ = 1−3g1. The printed D would give Λ = 1−g1, so A25 contradicts
  itself.

**P8-1: GAP (presentation-level).**
- Every stated identity is true. The expansion gives GT = −2(F2−XA1) and
  FT = −2F2, the fully covariant tensor check. The background curvature is
  R = −6(Ḣ+2H²).
- The unreduced scalar action matches, with Δ = (2A1−4F2X−A3)/(2T) and
  Λ = −2F2+4F2X+FTΔ. The same Θ appears in both the 6Θnv̇ term and the
  shift terms, and the matter couplings match. Euler residuals are 0 for
  both M0 and M1.
- The gap: the certificate never machine-derives the gradient matrix of
  the regular crossing chart (the "b chart"). `gamma.py:61-68` builds
  G_b11 from hand-typed `before_ibp` and `derivative` expressions, types
  in G_b12 = l/2 and G_b22 = 1/2 directly, and its residual
  `regular_gradient_IBP` only compares two hand-written formulas.
- Independent re-derivation of the b-chart Lagrangian from H0 + R²/J
  (q = k²/a²) gives these order-q terms: ḃb: −2qTΘ/Λ;
  b²: q(FTΘ²/Λ² − 2HTΘ/Λ − 2ṪΘ/Λ); bs: −ql; s²: −q/2. ḃs and ṡb are only
  order 1, so there are no order-q gyroscopic terms. The unitary-chart M1
  reduction also has none at order k².
- So G_b = [[Θ²FS/Λ², l/2], [l/2, 1/2]] and G_b11 = −TΘ̇/Λ at Θ = 0 are
  correct. The tilted-slicing degeneracy check, which supplements the
  unitary-gauge one, was not re-derived.

**P8-2.C, P8-2.D, P8-2.CD: PROVEN AS STATED.**
- Only the covariant functions recorded in the witness JSON files were fed
  into the independent derivation. All background residuals are 0,
  including χ backreaction for CD. Dropping the matter leaves exactly
  χ̇²/2 = 1/(200d¹²), so the backreaction is real.
- Θ, Λ and J are identical to the certificates, and GS = FS. CD satisfies
  the exceptional relation exactly, and a³χ̇ = 1/10.
- Independent SymPy real-root isolation (not the repo's FLINT code): J,
  Θ/t and T−1/4 have no real roots; FS's only pole is at t = 0; Λ's zeros
  sit at |t| ≈ 0.474, 0.510 and 0.380, outside the |t| ≤ 1/4 crossing
  chart; G_b11(0) = 15/4, 2 and 6.
- The certificates prove genuinely all-time statements: signs by Sturm on
  t² ∈ [0,∞), and tail bounds by Bernstein on a compactified variable.
  Only trivial written lemmas remain (a ≥ 1 implies completeness; the
  field-redefinition bookkeeping), and both check.
- Caveat: the "conventional tails" are exponential-potential scaling
  tails, not free-scalar tails. The potential term is the same order as
  H²: lim d·(F0−FX+FXX/2) = −1/4, −10, −44, against lim d·H² = 1/4, 4, 16.
  So potential energy ≈ (3p−1) × kinetic energy, i.e. ½, 5 and 11 times.
  The power laws are p = ½, 2 and 4, all above 1/3. In the contracting
  tail, homogeneous anisotropy therefore grows relative to H² even at
  linear order. The existing BKL/infrared exclusion covers this, but it is
  not stated next to "stable".

**P8-2.reg: PROVEN AS STATED.**
- The 457 tiles are contiguous on [−500, 500]; the 9 crossing-chart tiles
  lie on [2.20, 4.39]. The smallest certified lower bounds are positive:
  P ≥ 2.95e−8 and J−P ≥ 9.96e−9.
- Spot-checked tail constants hold: |e| ≤ xE, |e_x| and |e_xx| ≤ 2xE, the
  bounds on g, max x(1+x²)e^{−x} ≈ 1.51 < 3, P0 ≥ 1/(4000x²), and the
  gaps (100+1/3) and (100+7/25).
- The derivation gives Σ_derived − Σ_(47) = −2g1[(Θ+ġ1)˙ + 3H(Θ+ġ1)],
  exactly as claimed against the raw A25 (47)/(3.24) text. It reproduces
  A25 (55). Equation (28)'s sigmoid placement in the raw HTML matches the
  repo.
- Strong-gravity model re-derived: Θ = 1/(2d³), Σ = 0, GS = FS = 3T,
  ξ = 2I², ∫aFT dt = π²/2.

**P8-3: PROVEN AS STATED.**
- The case split is exhaustive: every component of {Θ≠0}, plus Θ ≡ 0
  handled separately.
- The proof is stronger than written. With Λ > 0 at a finite left
  endpoint of a Θ > 0 component, ξ → +∞ there, which contradicts ξ
  increasing. So the endpoint step needs neither strict health at Θ = 0
  nor the b chart; the b chart is needed only for the Θ ≡ 0 case.
- M0/M1 dictionary steps re-derived. C disabled: Δ = (T−FT)/(FT+3T) ≥ 0,
  so Λ = FT(1+Δ) > 0. D disabled: the mismatch factor is
  2X(FT+XA1)(A1−2F2X)/T, forcing A1 = 2F2X and Λ = T > 0.
  Reparametrizing the clock preserves F2X ≡ 0, A3 ≡ 0 and A1 ≡ 0; the
  induced terms feed only K and F.
- The minimal supports follow from the set of E rows being closed upward.
- However, `classification.json`'s rows, minimal supports and counts are
  generated from a hard-coded rule (`exclusions.py:46-58`,
  `verify_all.py:108-109`). That file restates the theorem rather than
  checking it.

**Conventions** are consistent throughout: X = φ̇² with no ½; the Einstein
term −R/2 gives GT = FT = 1; the Riemann convention is MTW's. A search for
P9 conventions finds only disclaimers.

## 3. Replay

- **pytest:** 83 tests collected, 83 passed in 74.74 s, exit 0, run with
  `PYTHONDONTWRITEBYTECODE=1` and `-p no:cacheprovider`.
- **`p8.verify_all --check`:** "P8(b): S0 plus 6 certificate replays
  passed; 32 scoped sector-row verdicts", exit 0, 47 s.
- The environment matches the certificate's pins (Python 3.11.14, SymPy
  1.14.0, flint 0.9.0, uv.lock hash). No in-scope file has changed since
  commit 417724fc.

## 4. Wording fixes recommended

1. **s3 §1:** "independently derived identity ξ′ = a(FS+FT)" is just the
   definition of FS. The real derived content is that FS is the
   unitary-chart gradient coefficient.
2. **s3 §1:** add the monotonicity endpoint argument, and state that the b
   chart is needed only for Θ ≡ 0.
3. **s1 §3:** "derived action with its boundary included" is true only for
   K_b. Either have `gamma.py` derive G_b or reword.
4. **s2-rational-witnesses §3:** replace "potential contribution → 0" with
   its size relative to H² (scaling tails, potential ≈ (3p−1) × kinetic).
5. **README / P8-3:** next to "stable", say the contracting tails have
   p > 1/3, so anisotropy grows at k → 0.
6. **Witness JSON "coverage" fields:** "core [−2,2]" and "unitary
   |t| ≥ 1/8" misdescribe the evidence. The Sturm proofs are global, and
   the unitary chart is valid for all t ≠ 0.
7. **Ledger P8-3 "finite lattice enumeration":** say the rows are
   generated from the proved rule.
8. **literature-digest:** D1 ("remains open") and D2 ("not replayed yet")
   are stale.
9. **s2-benchmarks:** "(100+1/3)/d²" means 100·(J0−P0); say so.

Editing any of these notes changes the source hashes stored in
`classification.json`, so that certificate must be regenerated afterwards.
