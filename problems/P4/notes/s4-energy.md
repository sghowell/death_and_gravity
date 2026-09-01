# P4 S4-2 — Theorem C route (b1): the direct energy identity and its certificate

Status: S4-2 done (2026-08-31). Code (new; nothing pre-existing modified):
`src/p4/validated/{tc_energy,tc_monotone}.py`; certificates
`results/theorem_c/{tc_energy,tc_monotone}.json`; tests
`uv run pytest problems/P4/tests/test_tc_reduction.py problems/P4/tests/test_tc_energy.py -q`
(15 tests). Inputs: `s4-reduction.md` (the exact (χ, η) reduction), `s2-theorem-b.md`
§§2–3 (class, σ(κ), E ⟺ eigenvalue), theorem-c-plan.md §4.2, GHJS arXiv:2509.12435
Lemmas 3.7–3.9 (the mechanism), the A1/A2/A3 certified background.

## 1. The normal form and the identity

On the constraint surface the fluid pair closes (S4-1, exact): χ′ = aχ + bη,
η′ = cχ + dη, and where b ≠ 0 the scalar equation χ″ − Tχ′ − Uχ = 0 holds with
T = a + d + b′/b, U = a′ + bc − a(T − a). The Liouville substitution v = χ/f,
2f′/f = T (GHJS Lemma 3.7's move) gives the **normal form**

    v″ = 𝒬(x; κ) v,     𝒬 = U + T²/4 − T′/2,

with 𝒬 an *exact rational function of κ with real x-dependent coefficients*.
Its κ²-density is the characteristic discriminant of the fluid block,

    q₂(x) = −det C + (tr C)²/4 = (3 − V²)/𝔴 + (tr C)²/4 = (μ₊ − μ₋)²/4 > 0,

a **sum** of positives (hyperbolicity: μ± real distinct). The "square-completion
margin" Hq = (3 − V²) + D(trC)²/4 of `s4-reduction.md` §3, negative on the sonic
shoulder (−0.354, 0), is the *difference* — it does not obstruct this route; the
shoulder is only S4-3's (strip) problem.

**Identity.** Let κ = ε + iτ with ε ≥ c₀ and let p be a nonzero solution in the
class of `s2-theorem-b.md` §3.8 ((i) analytic at the sonic point, (ii) on Σ,
(iii) regular at the centre, (iv) sonic gauge N_p(0) = 0). With χ = ℓ(p),
v = χ/f: multiplying v″ = 𝒬v by v̄ and integrating over (−∞, 0),

    ∫ |v′|² + ∫ 𝒬 |v|² = [v̄v′]₋∞⁰ = 0                    (§2 for the boundary terms).

Writing 𝒬 = q₂κ² + q₁κ + q₀ + R(x; κ) (κ-polynomial part + proper rational
remainder; poles only at κ = A(x), −F_N(x) and the three roots of the Bc-cubic,
all with Re ≤ 2 — condition (d) below):

    Im-part:  0 = ∫ Im𝒬 |v|²,   Im𝒬 = τ(2εq₂ + q₁) + Im R,  |Im R| ≤ |τ| B₁,
    Re-part (τ = 0):  0 = ∫ |v′|² + ∫ (q₂ε² + q₁ε + q₀ + R(ε)) |v|²,  |R(ε)| ≤ B₀,

where B₀ ≥ sup_{Re κ ≥ c₀} |R|, B₁ ≥ sup |∂R/∂κ| (maximum modulus: the sups sit
on the line Re κ = c₀; |Im R(ε+iτ)| = |∫₀^τ Re R′(ε+is) ds| ≤ |τ|B₁). Hence the
pointwise conditions, for all x ∈ (−∞, 0),

    (a) q₂ > 0,  (b) 2c₀q₂ + q₁ − B₁ > 0,  (c) c₀²q₂ + c₀q₁ + q₀ − B₀ > 0,
    (d) all κ-poles of 𝒬(x; ·) have Re ≤ r̄ = 2 < c₀,

give: (b) ⟹ Im𝒬/τ > 0 for every ε ≥ c₀ (q₂ > 0 makes 2εq₂ + q₁ increasing in ε)
⟹ τ = 0; then (b)+(c) ⟹ 𝒬(x; ε) > 0 for all ε ≥ c₀ ⟹ v ≡ 0. So **χ ≡ 0**;
b ≠ 0 (condition (d): the Bc-cubic has no roots with Re ≥ 2) gives η ≡ 0, so
ℓ(p) = m(p) = c(p) = 0, i.e. p(x) ∈ span g(x) pointwise (S4-1: ker(ℓ, m)|_Σ =
span g). Both p and g solve the linear system and P is invertible on the open
cone (SΔ̃ ≠ 0 ⟸ the certified 𝔴 = −D > 0 row of §6), so p = const·g; the gauge condition (iv) then forces const·N₀(κ − κ̄) = 0,
and κ ≠ κ̄ (Re κ ≥ c₀ > κ̄) gives p = 0 — contradiction. **No eigenvalues with
Re κ ≥ c₀ in the class; hence E(κ) ≠ 0 there** by the E = 0 ⟺ eigenvalue
equivalence of `s2-theorem-b.md` §3.8, whose three ingredients are κ-independent
or hold for all κ with Re σ(κ) < 0 — automatic on Re κ ≥ 0 (plan §1).

## 2. Function space and boundary terms

*Sonic end.* An admissible p is analytic at x = 0; χ = Δ̃[(N′+κN)V_p − V′N_p]
has the factor Δ̃ = 4SWD with D(u(0)) = 0 **exactly** (A1 closed forms), so
χ = x·(analytic), χ, χ′ bounded. The transform weight obeys
xT(x; κ) → 1 + σ(κ) = (1−σ₀)(1 − κ): certified enclosure (closing box, §3)
xT = τ₀ + τ₁κ + O(2.3e−5) with τ₀ = 0.90989405…, τ₁ = −0.90989403…
(= ±(1−σ₀), σ₀ = 0.090105970507920(5e−16)). Hence |f| ~ |x|^{(1+Reσ)/2} and
|v| = O(|x|^{(1−Reσ)/2}), |v̄v′| = O(|x|^{−Reσ}) → 0 since Re σ(κ) =
−σ₀ − (1−σ₀)ε < 0; ∫|v′|² and ∫|𝒬||v|² converge at 0 for the same reason
(exponents −1−Reσ > −1). *Centre end.* Class (iii) bounds the scaled
perturbations, so χ, η = O(1) as x → −∞ (the covector weights Δ̃(N′+κN) ~ Ce⁻ˣ
compensate V_p ~ ṽ_p eˣ etc.); T → T∞ with **T∞ = −3 exactly** (the two
χ-branches at the centre have t-exponents {0, −3} and the normal form has no
first-order term, forcing exponent-sum 0 = −T∞ − 3; float: 𝒬 → T∞²/4 = 9/4 ✓),
so |f| ~ e^{−3x/2} → ∞ and v, v′ → 0 exponentially; all integrals converge.
The interchange limits defining [v̄v′] therefore vanish; no uniformity in κ is
needed (the identity is per-κ).

## 3. The certificate (`tc_energy.py`), c₀ = 6

The conditions (a)–(d) are certified in outward-rounded ball arithmetic on boxes
covering **x ∈ [−4.5, 0]**; the far tail x < −4.5 is Lemma D (§4). Structure:

* **Exact algebra** (fmpq_mpoly, verified at import): SD = Dp·SDr, D2 = Dp³·D2r,
  Bc = Dp²·Bcr, Ac = Dp²·Acr, Ae = Dp·Aer, Be = Dp·Ber, TN1 = Dp⁵·TN1r,
  QNm = Dp¹⁰·QN_red (Dp = the sonic polynomial D), and the reduced identity
  U1 = Dp⁶(Bcr·UN3r − Acr·EBc·D2r). 𝒬 = QN/(P³·SDr³D2r²Bcr²) with QN assembled
  from the small reduced atoms only — no hidden Dp-power cancellations.
* **The field identity.** QN_red(u₀) = 0 **exactly in ℚ(V₀, √3)** at the KHA
  closed-form sonic point, identically in V₀ and κ (mode ``field``, ≈ 2 min of
  exact fmpq_poly arithmetic; re-verified this session). This is the order-11
  vanishing of QN at x = 0 and the engine of the closing box.
* **Jets.** Every box is evaluated in order-2 jets anchored at the box midpoint
  (value and d/dx tight via the exact D1-mpolys; only the curvature slot sees
  the width); the κ-division QN/QDr runs in jet arithmetic, so the quotient
  q₂, q₁, q₀ and the remainder coefficients stay tight. On the tube the atoms
  are evaluated in the scaled variables (n, w, v, t) with explicit t-powers.
* **B₀, B₁** by τ-boxed line-suprema on Re κ = c₀ (48 boxes to τ = 12,
  geometric to TTAIL = 64, second-order recentring, factor-wise |QD| lower
  bounds) + a coefficient tail for |τ| ≥ TTAIL; B₁ = min(direct |R′|-numerator
  bound, Cauchy shift sup_{Re ≥ c₀−1.2}|R|/1.2). Pole condition (d) by
  Routh–Hurwitz on the shifted Bcr-cubic plus the A(x), −F_N(x) balls.
* **Regions.** (i) tube steps with x − h ≥ −4.5 (certified A3+A2 tube,
  sub-boxes of half-width from ~|x|/1500 down to ~|x|/45000, adaptive); (ii) the sonic ladder
  [−0.05, −10⁻⁸] (certified A1 series, adaptive bisection); (iii) the **closing
  box [−10⁻⁸, 0]**: by the field identity and the mean value theorem
  QN_red = x·QQ with QQ = ∇QN_red·u′ hull and D = x·E_D (D(u(0)) = 0 exact),
  E_D ∈ 7.5605725(1) = D₁, so x²𝒬 = QQ/(E_D³·SDr³D2r²Bcr²) has nonvanishing
  factors up to and including x = 0. Known answer (all three contained):
  x²𝒬 → (σ(κ)²−1)/4, i.e. q̃₂ ∋ (1−σ₀)²/4, q̃₁ ∋ σ₀(1−σ₀)/2, q̃₀ ∋ (σ₀²−1)/4.
  Margins at the closing box (c₀ = 6): q̃₂ ≥ 0.207, (b) ≥ 2.17, (c) ≥ 7.24.

At the sonic end the conditions are tight where the σ-threshold lives: (c)
scales as ((σ(ε)² − 1)/4)/x², positive iff ε > 1 — the GHJS "Re λ > 1"
threshold transferred. The certified inequality is tightest at the deep end of
the tube (x ≈ −4.4: mB-margin 2.9e−4 against q₂ ≈ 2.4e−4, i.e. ≈ 1.2 q₂) and in
the mid-tube x ≈ −1…−2; the working c₀ = 6 is set by those interval margins,
not by the analysis (float margins stay positive down to c₀ ≈ 4.5).

## 4. Lemma D (deep tail x ≤ −4.5) — analytic, float-validated, NOT yet certified

**Statement.** For x ≤ −4.5 and Re κ ≥ 6: the pole condition (d) holds (all
κ-poles of 𝒬(x;·), including the Bc-cubic roots, have Re ≤ 2 — float: they
cluster at κ = 1 with sup Re = 1.00 down to x = −16), Im𝒬(x;κ)/τ > 0 (τ ≠ 0)
and 𝒬(x; ε) > 0. *Evidence and structure.* With t = eˣ: q₂/t² → 1.9618,
q₁/q₂ → −1, q₀ → 9/4 (exactly, from T∞ = −3), and the rational remainder has
all its poles clustering at κ = 1 with residues O(t²); float line-sups give
B₀ ≤ 1.2e−6, B₁ ≤ 1.6e−6 at x = −8 against margins (2c₀q₂+q₁) ≥ 2.6e−6 and
(c) ≥ 2.25 − o(1); on the 0.1-grid over [−16, −4.5] at c₀ = 6 the relative
margin of (b) stays ≥ 0.39 and of (c) ≥ 0.99 (improving like t² beyond −16). *Why not certified here:*
the t²-smallness of the remainder is a coefficient-relational cancellation
(it vanishes modulo the fixed-point ideal (2NV+1), not term-by-term), so every
generic interval evaluation of the division remainder carries an absolute junk
floor that crosses the t²-margins near x ≈ −5.5 (certified boxes still pass at
x = −4.3 with sub-boxes of half-width |x|/44000, and fail at −5.8 at any width;
the lemma-region has t² ≤ 1.2e−4). The identified closing
computation — exact reduction of the t⁰/t¹-parts of QN_red and of the
denominator modulo (2NV + 1) in the centre-scaled variables (t-power labels
then carry the margins for all t) — is inherited by S4-3.

## 5. Certified results (c₀ = 6, RBAR = 2; JSONs under results/theorem_c/)

| region | boxes | min (a) q₂ | min (b) | min (c) | status |
|---|---|---|---|---|---|
| tube x ∈ [−4.5, −0.05] | 14520 | 2.42e−4 | 2.91e−4 | 2.249 | PASS |
| sonic ladder [−0.05, −1e−8] | 16666 | 0.207 | 2.33 | 7.33 | PASS |
| closing box [−1e−8, 0] | 1 | 0.207 | 2.17 | 7.24 | PASS |

**Certified statement (P4-8 draft).** In the class of `s2-theorem-b.md` §3.8,
there is no eigenvalue with Re κ ≥ 6 — hence E(κ) ≠ 0 on {Re κ ≥ 6} —
**conditional on Lemma D** (whose conditions are the same pointwise
inequalities on x ≤ −4.5, float-validated with relative margin ≥ 0.39). All other
ingredients (identity, class, boundary terms, conditions on [−4.5, 0]) are
certified or exact.

## 6. P4-3 promotion: the tube sign certificate (`tc_monotone.py`) — PASSES

All in ball arithmetic on (u, u′)-enclosures (no singular denominators):

| claim (closed cone) | tube (344 steps ×4) | sonic [−0.05, 0] | centre t ∈ (0, e⁻⁸] |
|---|---|---|---|
| v_rel′ > 0 | ≥ 1.36e−4 | ≥ 0.478 | ≥ 0.404 (of v_rel′/t) |
| −ρ̂′ > 0 | ≥ 7.69e−6 | ≥ 3.25 (of −ρ̂′/ρ̂) | ≥ 82.7 (of −ρ̂′/t²) |
| v_rel > 0 | ≥ 1.36e−4 | ≥ 0.552 | ≥ 0.404 (of v_rel/t) |
| 1/3 − v_rel² > 0 | ≥ 0.0273 | ≥ 3.03e−3 (interior) | ≥ 0.333 |
| −D > 0 (𝔴 > 0) | ≥ 0.338 (of −t²D) | D′ ≥ 7.44 & D(0) = 0 exact | ≥ 1.53 (of −t²D) |
| D′ > 0 (so 𝔴′/𝔴 < 0) | ≥ 6.33 | ≥ 7.44 | ≥ 3.06 (of t²D′) |

The last row certifies GHJS's commutator sign 𝔴′/𝔴 = D′/D < 0 **globally** on
the open cone (S4-1 §3 observed it in floats; S4-3 needs it only near x = 0).
Two exact inputs close the endpoints: D(u₀) = 0 and v_rel(0) = 1/√3 (closed
forms), so v_rel < 1/√3 and w = 1/3 − v_rel² > 0 on the *open* cone follow from
v_rel′ > 0; and w′ = −2v_rel v_rel′ < 0. At the centre the k = 1 coefficient of
the ŵ-series vanishes exactly (the order-1 recursion is homogeneous; checked),
so −ρ̂′/t² → 2|ŵ₂| + ŵ∞(A−1)/t²|₀ ≈ 82.7 > 0 certifies ρ̂′ < 0 on the whole
open tail. V_R (= KHA's V) remains non-monotone (one interior zero,
x = −0.2528) — float remark only, `results/ec_monotonicity.json`.

## 7. Draft CLAIMS rows (the main session owns the ledger)

* **P4-3** (promote CONJECTURED → CERTIFIED): "Monotonicity on the closed sound
  cone: v_rel ↑ (to 1/√3), ρ̂ ↓, w = 1/3 − v_rel² ↓ and > 0 on the open cone,
  and the reduction weight 𝔴 = −D = 3(N+V)²w > 0 — certified as an Arb tube
  sign certificate, together with D′ > 0 (𝔴′/𝔴 < 0, the GHJS commutator sign)
  (tc_monotone.json: A3 tube + A1 sonic series incl. x = 0 + A2 centre tail,
  t-factored); V_R non-monotonicity stays a float remark."
* **P4-8** (new, CERTIFIED*): "Direct-energy exclusion: no eigenvalues (class
  of P4-6) with Re κ ≥ 6, hence E ≠ 0 there; certified interval conditions on
  x ∈ [−4.5, 0] (tc_energy.json; exact reduction identities + the ℚ(V₀,√3)
  field identity at the sonic point; GHJS-style normal form) — *conditional on
  Lemma D (deep tail x ≤ −4.5, same pointwise conditions, float-validated,
  closing computation identified; see s4-energy.md §4)."

## 8. S4-3 handoff

* (b2) strip: Re κ ∈ [0, 6], |Im κ| ≥ b₁ with target **b₁ ≤ 14** (then Theorem
  B covers the rest of Ω; the plan's R has |Im| ≤ 14). Inherits: the Hq < 0
  sonic shoulder (−0.354, 0) ⟹ GHJS commutator (m·𝔴′/𝔴 with 𝔴′/𝔴 < 0 now
  CERTIFIED via P4-3, D′ > 0 rows); the certified normal-form machinery of tc_energy
  (jets, line-sups, closing box) reusable verbatim.
* Lemma D closing computation (§4): exact mod-(2NV+1) reduction of the t⁰/t¹
  parts at the centre (the certified region already reaches −4.5; the floor is
  ≈ −5.5 at any sub-box width, so only the exact route closes the tail).
* Numbers: σ₀ = 0.090105970507920(5e−16); D₁ = E_D = 7.5605725(1); T∞ = −3;
  xT-affine (0.9098940, −0.9098940) ± 2.3e−5; q₂ → (1−σ₀)²/4 at the sonic end,
  → 1.9618t² at the centre; float c₀-floor ≈ 4.5, certified c₀ = 6.

## 9. Reproduction

```
PYTHONPATH=problems/P4/src uv run python -m p4.validated.tc_energy all <tube.json> 6.0   # ~1-2 h
PYTHONPATH=problems/P4/src uv run python -m p4.validated.tc_energy field                 # ~2 min, exact
PYTHONPATH=problems/P4/src uv run python -m p4.validated.tc_monotone <tube.json>         # ~4 min
uv run pytest problems/P4/tests/test_tc_reduction.py problems/P4/tests/test_tc_energy.py -q
```
