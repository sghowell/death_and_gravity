# P2 S4 report: part (i) at P_B, uniformly in l ≥ 20

Status: S4 exit, 2026-10-18; revised 2026-10-19 after the independent review. Branch `p2/s4-part-i`,
from main at `dc5bfd27`.
Pre-registration: `notes/s4-preregistration.md`, written at `d00f9a3c` before any production run.
- Amendments 1–3 (requested by the coordinator) were committed at `c8bd7436`, Amendment 4 at
  `8f189b18` and Amendment 5 at `e639d96e`, each before the production runs it affects.
- Amendment 6, at `a16b52cb`, is post hoc and disclosed as such. It defines when K4's float secant
  search has found a zero (§3). A dated correction (2026-10-19) adds that the launch-time code had a
  rule, under which K4 fails, and gives the actual secant mechanism (review S2).

Proof: `notes/s4-proof.md`.

**Review.** `docs/review-2026-10-18-p2-s4-adversarial.md` (2026-10-18/19) found the theorem proven as
stated, with nothing blocking. The user approved its fixes, applied on this branch on 2026-10-19:
- **S1.** The S closing box's truncation term now bounds P_N over the complex ρ the contour visits.
  - All S leaves were regenerated at `68de8184`: the two closing leaves for S1, and the rest for W4's
    outward printing.
  - The closing margins agree with the first run to 8 digits (§1).
- **S2.** A dated correction to Amendment 6 (§3, §4).
- **S3.** The simply connected union of the strips and the checked windows, in W1/W2 (proof §§2–3).
- **S4.** The verdict and ledger wording (here, §7, FORMULATION v0.6).
- **S5.** The pilot discrepancy (§4).
- **W1–W8.**
  - Errata appended to the pre-registration (W1, W5).
  - The proof's W̃ minimum and L tiling (W2, W3).
  - Bounds printed outward, and L's summary taken from ball endpoints (W4).
  - N2's h-range (W6) and T's target-limited margin (W7), both in this report.
  - The first K4R launch's logs restored to `results/part1/` (W8).

**Verdict: PASS, l₀ = 20.**
- **The theorem.** At P_B, for every l ≥ 20 and every ω with Re ω ≥ 0 and −d_δ ≤ Im ω ≤ 0
  (d_δ = 0.55κ₋), W_l(ω) ≠ 0.
- **P2-1 at P_B.** With the conjugation symmetry and P2-0, W_l has no zero in S_δ for every l ≥ 20.
  This contains P2-1w's statement (its windows lie in that half-strip at l = 20, 50, 100). The level
  is WRITTEN_PROOF:
  - the inequalities are certified in Arb;
  - W1 (vertical-strip analyticity) and W2 (Jost continuation and decay on the curved throat
    contours) are non-routine written arguments, as in P2-1w;
  - the Im ω > 0 part is at P2-0's level.
- **Every l < l₀ = 20 is counted at P_B**, on P2-3B's boxes R′_l (S3 counted on R′_l, not on D7's
  R_l). So the S3 fallback was not needed.
- **Status.** No kill criterion fired. Compute was ≈ 3.9 core-h of the 48 core-h cap (projection 9).
  The review fixes added ≈ 0.5 core-h.

## 1. What was certified

| regime | domain (every h = 1/(l+½) ≤ h₀ = 2/41) | certificates | certified margin |
|---|---|---|---|
| L | 0 ≤ Re ω ≤ 8κ₋ | `L_PB_plus.json`, `L_PB_c.json` (1760 ω-boxes of side κ₋/20, none refined) | Re y₊ ≥ 0.1245, Re y_c ≤ −0.0804 at the exits |
| S | 6κ₋ ≤ Re ω ≤ 0.9Ω_P/h | `S_PB_plus.json` (47 ζ-leaves), `S_PB_c.json` (45) | min q̂ = 0.0130 (+), 0.0091 (c) |
| T | 0.8Ω_P/h ≤ Re ω ≤ 2.0Ω_P/h | `T_PB.json` (24 ζ-boxes, 37 146 cells) | E ≥ 0.1003·a*/h on the contour, for every h ≤ h₀ |
| A | Re ω ≥ 1.5Ω_P/h | `A_PB.json` (ζ_A plus W7) | E ≥ 6.955·a*/h at Re ω = ζ_A/h, monotone above |

Here a* = √(|U″|/2) = 2Ω_Pγ₀ = 0.029969, and ν ∈ [−d_δ, 0] throughout. Also certified:
- `bulk_PB.json`: U is monotone on each side of r_P; W̃ ≥ −99/12500; C_W = 0.03168.
- `cover_PB.json`: the cover inequalities and the tilings.

### L
- The extremes over all boxes, taken from the stored balls' endpoints and rounded outward (W4), are
  Re y₊ ≥ 0.12456385 and Re y_c ≤ −0.080453488.
- The exits lie at h²τ_a ≤ h₀²·8 = 0.019 from r₊, and at h²τ_a^c ≤ h₀²·5529/64 = 0.206 from r_c.
- Every box certified at the registered side κ₋/20. The quadrisection (down to κ₋/80) was never
  used.
- **Arithmetic.** The scaled Frobenius recurrence ran in affine arithmetic, in the deviations of ω
  and of h², on 16 sub-intervals of h² ∈ [0, h₀²] (Amendment 5).
  - Plain Arb rectangles wrap by up to √2 per complex multiplication. Over a box at Re ω ≈ 3κ₋ they
    lost all accuracy, and at Re ω = 8κ₋ the series did not converge in interval arithmetic.

### S
- **Leaves.** Per side: the closing box ζ ∈ (0, 0.02Ω_P] plus the 44 registered ζ-boxes of width
  0.02Ω_P.
  - On the + side, [0.02, 0.04] and [0.06, 0.08] (Ω_P units) were halved once.
  - No ε-box needed halving.
- **The smallest margins** are near the top: 0.0142 (+) and 0.0098 (c) on [0.88, 0.9], and 0.0130 (+)
  and 0.0091 (c) on [0.86, 0.88] and [0.82, 0.84].
- **The closing boxes.** Design w = 2, y_a = 1/20, φ_∞ = −π/3. The near-horizon model has
  ρ_max/R = 0.48 (+) and 0.58 (c), with order 14.
  - **Regenerated after S1** at `68de8184`, with the truncation term bounding P_N over complex ρ.
    | | mu_cells (outward) | mu_tail | review's prediction |
    |---|---|---|---|
    | + | [0.41468628, 0.41468629] | [1.5341632, 1.7320506] | 0.41468629; 1.5341632 |
    | c | [0.35342094, 0.35342095] | [1.5341632, 1.7320506] | 0.35342094; 1.5341632 |
  - The first run printed 0.41468629 and 0.35342094 (rounded to nearest). The values are unchanged
    to 8 digits.
  - All 92 S leaves were regenerated at that commit. Every design, X_a, k, target and cell partition
    is as before. The printed bounds moved outward by at most one unit in the 8th significant
    digit: at most 1.0e−7 relative, and up to 1e−5 absolute for the exit margins of order 10²–10³
    (review §9, FU1).
- **The other leaves.**
  - Designs from the registered grid plus the added half-widths w = 2, 3 (Amendment 5). w = 4 was
    never selected.
  - The validated contour integration found every a-priori box.
  - The tail discs are ≤ 6.2·10⁻⁵.
  - Every W1 check passed. The + contours lie in Re x < Re r*(0); the c contours cross the two
    special lines only within S2's checked segments.

### T and A
- **T.** The worst cells sit at s ≈ ±2.4 for ζ ≤ 0.95Ω_P, and at s = 0 above that. The kinetic
  coefficient c ≥ 0.1235, and E₀ ≥ 0 on every cell.
  - **The margin 0.1003a* is target-limited** (review W7). It is set by the bisection target a*/10,
    not by the design: cells are refined only until they reach it. The float minimum of the same box
    bound ([0.8, 0.85]Ω_P) is 0.139a*, and the certified values stay well below the float ones up to
    ζ ≈ 1.45Ω_P.
  - The registered √ ground state was replaced by the linear one (Amendment 4). The float h-uniform
    margin dropped from 0.175a* to 0.156a* at ζ = 0.8Ω_P.
- **A.** E₀ ≥ 0.0109 everywhere. The monotonicity condition holds:
  φ ≤ −0.2 < −arctan(d_δh₀/ζ_A) = −0.00364.

## 2. Method

As pre-registered, with the amendments listed above. In outline:
- **T.** S2's contour, with the exact splitting h²E = E₀ + hE₁ + h²E₂. The per-cell bound
  E₀/h₀ + E₁ + h₀·min(0, E₂) ≥ μ (with E₀ ≥ 0) holds for every h ≤ h₀, and endpoint evaluation in
  (ζ, ν) is exact.
- **A.** The same contour, with φ₀ = −0.2 and monotonicity in Re ω.
- **S.** Half-line coercivity on C³ contours below each turning point, in the normalised form
  (h enters only through ζ and ε). Validated curved-contour integration and mean-value cells. A
  near-horizon Taylor model closes ζ → 0.
- **L.** The τ-scaled Frobenius series at both horizons, uniform in h² ∈ [0, h₀²]; the exit rule;
  then Riccati invariance through the forbidden bulk.
- The four are joined by the cover (§1 of the proof).

## 3. Acceptance criteria and controls (`controls_PB.json`, `replay_PB.json`)

All pass. K-unsound did not fire.

| item | result |
|---|---|
| **K1** | S2's geometry certificate replays (A0). |
| **K2** | L's scaled-Frobenius Jost values, continued to r_P, give W_l agreeing with zsol and arbw at all 12 points (l ∈ {20, 30, 50} × ω ∈ {0, 2 − 0.2i, 5 − 0.55i, 8 − 0.3i}κ₋). All balls overlap; relative radii ≤ 1.5·10⁻³⁸ (L), < 10⁻⁴⁹ (zsol, arbw). |
| **K3** | At 3 points in each overlap (L∩S, S∩T, T∩A), both adjacent certificates hold. |
| **K4** | The float acceptance check, not the kill criterion of that name (W5). At 200 deterministic points (64 L, 52 S, 37 T, 47 A; l ∈ [20, 400]), the float W_l is nonzero everywhere. Under the launch-time rule K4 fails (61 of 200); under Amendment 6's criterion it passes. See below. |
| **N1** | The T checker on boxes holding the PS fundamental mode does not certify, at l = 20, 50, 100: IN (±γ₀/50) gives μ/a* = −0.39, −0.070, −0.070, and STRAD gives −0.42, −0.100, −0.100. The STRAD value −0.100 is exactly E₁(0)/a* at ν = −1.05γ₀: the checker's threshold sits at ν ≈ −0.95γ₀, 5% short of the mode. |
| **N2** | The h-box [1/50.5, 1/20.5] with the ζ-box [Ω_P − 0.001κ₋, Ω_P] and the ν-box [−1.69, −1.68]κ₋ does not certify (μ/a* = −0.38). The T checker has no lower bound on h, so it ran on (0, h₀], a superset of the registered h-box (W6). A certification there would have been the stronger claim, so the control is valid. |
| **N3** | The L continuation path at l = 1, 2 does not certify the NE-mode boxes: Re y₊ ∈ [−0.232, −0.222] and [−0.316, −0.293], a clean sign failure. The same checker certifies mode-free boxes at half the depth (Re y₊ ≈ 0.146 and 0.164). |
| **N4** | The S checker with the exit in the allowed region does not certify: Ũ(x_a) − 1 ≈ −0.62 (+) and −0.66 (c), and min q̂ ≈ −0.81 and −0.92. |
| **N5** | L's production path with the exit in the allowed region does not certify, in 4 of 4 boxes. In the two boxes away from the real axis, both the y-sign test and the exit rule fail. In the two with ν ∈ [−0.05, 0]κ₋, only the exit rule fails: there Re y ≈ Im ω + O(U) keeps its sign in the allowed region. So the exit rule (the forbidden bulk) is essential to L's argument, as the proof says. |
| **Replay** | After the review fixes, a **full** replay at `68de8184` (`verify full`, 9 min) agrees exactly: every one of the 37 005 T leaf cells and 1536 A cells, both tails, all 92 S leaves recomputed from their stored designs and targets, all 1760 L boxes, the bulk and the cover. The first, 5% replay at `0e3c40f7` had also agreed. |
| **Tests** | `tests/test_part1.py`, with the S1 test added: fast 9 pass and 2 skip; `P2_SLOW=1` 11 pass. The whole P2 suite: fast 86 pass and 17 skip (35 s); with `P2_SLOW=1`, 103 pass, with 1 warning (the pre-existing overflow in S2's m3 Riccati design; 662 s). |

**K4 in detail.**
- **The raw outcome.** The secant search (≤ 12 steps, seeded at each point) stopped inside H⁺ at
  61 points, all at large |ω| (46 A, 15 T; Re ω up to 4214κ₋).
- **None is a zero.** At every one, |W(root)|/|W(seed)| ∈ [1.0000008, 1.0000013]. The rigorous zsol
  ball excludes 0 at all 61, with zsol's Frobenius term cap raised for this check only.
- **The rule changed after launch.** The code at launch (`35201758`; K4 was launched from
  `ca66270b`) counted every secant stop in H⁺ as a zero: `assemble_K4` set
  `pass = … and not any(root_in_H_plus)`. **Under that rule K4 fails, 61 of 200.** Amendment 6
  replaced the rule after three such stops among the first 14 points, and its first text omitted
  this (correction of 2026-10-19, review S2).
- **Why the search stopped: overshoot and return** (the review's traces at three stops).
  - The first secant step overshoots to Im ω ≈ −7.5, −13.0 and −6.9 κ₋, where |W| ≈ 10¹⁵ to 10⁴⁸.
  - The next step returns to the previous iterate to 16 digits.
  - The step after that is ≤ 3·10⁻¹⁶ relative, so the step-only test fires at the perturbed seed
    x₀(1 + 10⁻⁶) + 10⁻⁹i.
  - The phase of W turns only slowly there (d arg W/d Re ω = −3.35, −1.87, −3.68). The first text's
    "fast-turning phase" explanation was wrong.
- **K4 carries almost no weight.** It is a float sanity check, and the proof does not use it.
- **A side finding.** The first residual run stopped at zsol's 6000-term cap. Its launch record is
  kept, and the cap was raised as above.

## 4. Deviations from the pre-registration (all dated amendments)

| # | what | when | why |
|---|---|---|---|
| A1–A3 | W1 per contour (`w1_check`), N5, the ζ → 0 wording | before production | coordinator's go message |
| A4 | T's ground state is linear in s, not √(q₀/c) | before production | √ of an interval at s = 0 |
| A5 | S's compact C² family, w grid extended by {2, 3, 4}; S refinement and cell details; the closing box's near-horizon model; L's affine arithmetic and 16 h²-sub-intervals; N3 and N4 set-up; K4 sampling; replay sampling | before the L, S and control runs | development showed each was needed (§1) |
| A6 | K4: a secant stop counts as a zero only if W is small there | after the first K4 results (post hoc) | a step-test stop is not a zero; the launch-time code counted every stop in H⁺ as one, so K4 fails 61/200 under it |
| A6 correction; errata W1, W5 | the launch-time rule and the secant's overshoot-and-return mechanism (S2); Amendment 4's a* − a_d is 3.0·10⁻⁹ (W1); the two meanings of "K4" (W5) | 2026-10-19, after the review | review findings S2, W1, W5 |

Two pilot numbers in §9 of the pre-registration did not survive:
- **The S margin near the top was 0.087 in the pilot; S4 certifies 0.009–0.014 there.**
  - The pilot's 0.087 is a leading-order margin: min_t q₀/|Y′|² at ε = 0, normalised by |Y′|². The
    certificates bound the full q̂, which includes the O(ε) terms, with ε up to h₀/ζ ≈ 0.26–0.29 near
    the top. The first text attributed the whole drop to the contour family (review S5).
  - The drop has three parts (the review's `pilot_compare.py` on the stored top-leaf designs):
    1. The compact contour family (and the ζ-box): 0.087 → 0.025–0.032, at leading order.
    2. The O(ε) terms at the l = 20 corner of each leaf → 0.015–0.018 (full q̂, float).
    3. The certification loss of the mean-value cells → 0.009–0.014 (certified).
- The registered grid alone gives float margins ≤ 0.06 (+) and < 0 (c) at small ζ. Amendment 5 added
  half-widths 2, 3 and 4.

## 5. Compute

| item | CPU |
|---|---|
| pilot float studies (before the pre-registration) | ≈ 0.05 core-h |
| development: float design studies (S family, bend widths, closing-box exits), L arithmetic prototypes, S/T/L sample runs, control trials (a stopped N4 run of about 5 min and an N3 run of 3.5 min), timing | ≈ 0.8 core-h (estimated; not logged per run) |
| production: T 2.4 s, A 0.1 s, S 213 s (both sides), L 337 s, bulk 1 s | 0.15 core-h |
| K4: 200 float points | 2.00 core-h (31 min wall, 4 workers) |
| K4 residual and zsol checks (Amendment 6), incl. the stopped first launch and a timing test | 0.61 core-h |
| controls (N1–N5, K1–K3), replay, tests (fast and slow), float and zsol cross-checks after production | ≈ 0.25 core-h |
| **S4 total** | **≈ 3.9 core-h** (cap 48; pre-registered projection about 9) |
| review fixes (2026-10-19): S regenerated (211 s), full replay (≈ 540 s), controls (≈ 190 s), the whole P2 test suite fast and slow (≈ 690 s), closing-box and S1-test checks | ≈ 0.5 core-h (cap 2) |

At most 4 workers ran at once, all at nice 10. Launch records are in `results/part1/launch_*.json`.
The regime certificates cost far less than projected: S about 0.06 core-h against 3, and T seconds
against 0.3. The float K4 sample was the largest item.

## 6. What S4 does not establish

- **Independent review.** Done: `docs/review-2026-10-18-p2-s4-adversarial.md` (2026-10-18/19)
  found the theorem proven as stated. Its fixes are applied (header).
  - The review's scripts are to be committed by the reviewer under
    `problems/P2/certificates/review_s4/`.
- **P_A.** Not attempted (§10 of the pre-registration).
- **The assembly** (P2-4, P2-5, the ledger rows): S5.
- **CERTIFIED for P2-1.** That would need W1 and W2 machine-checked (2D validated continuation of
  r(x)).

## 7. Ledger text (paste-ready; for the coordinator, not applied here)

This is the review's row (`docs/review-2026-10-18-p2-s4-adversarial.md` §7), as the user approved,
with one change: the evidence column also cites `problems/P2/certificates/review_s4/`, where the
reviewer will commit the review's scripts. The l = 20 count clause is kept. The row quotes no
closing-box number, so S1's regeneration does not change it.

**P2-1 (new row):**

| P2-1 | P2 | At P_B (ΛM² = 3/50, Q = (999/1000)Q_max, δ = 1/20), uniformly in l ≥ l₀ = 20: W_l(ω) ≠ 0 for every ω with Re ω ≥ 0 and −(½+δ)κ₋ ≤ Im ω ≤ 0. With the conjugation symmetry W_l(−ω̄) = conj W_l(ω) and P2-0 (Im ω ≥ 0), W_l has no zero in S_δ for any l ≥ 20; this contains P2-1w. With P2-3B, every l < l₀ is counted at P_B, on P2-3B's boxes R′_l (the assembly P2-4 is S5). Not P_A | WRITTEN_PROOF. Four regimes cover {l ≥ 20} × the closed half-strip; with h = 1/(l+½) ≤ h₀ = 2/41 and ζ = h Re ω, the cover's inequalities and tilings are certified in certificates/part1/cover_PB.json. **L** (Re ω ≤ 8κ₋): τ-scaled Frobenius series at both horizons (r = r_h ± h²τ), in Arb affine arithmetic over 1760 ω-boxes and 16 sub-intervals of h² ∈ [0, h₀²]; Re y₊ ≥ 0.1245 and Re y_c ≤ −0.0804 at the exits (y = Φ′/Φ); then Riccati invariance through the forbidden bulk. **S** (6κ₋ ≤ Re ω ≤ 0.9Ω_P/h): half-line coercivity on C³ contours below each turning point, uniform in h through ζ and 1/Re ω; validated curved-contour integration; a near-horizon Taylor model with a Rouché disc and Cauchy remainder for ζ ≤ 0.02Ω_P; min q̂ = 0.0130 (+) and 0.0091 (c) over 47 and 45 ζ-leaves. **T** (0.8–2.0 Ω_P/h): S2's rotated contour with the exact splitting h²E = E₀ + hE₁ + h²E₂ and a linear ground state; E ≥ 0.1003·2Ω_Pγ₀/h for every h ≤ h₀. **A** (≥ 1.5Ω_P/h): E ≥ 6.955·2Ω_Pγ₀/h at 1.5Ω_P/h, and monotone in Re ω above. Machine-checked (Arb, 128-bit): every cell, tail, exit and order condition, the bulk facts (U monotone on each side of r_P; W̃ ≥ −0.00792) and the W1 crossing checks. Written (notes/s4-proof.md): W1, the vertical-strip analyticity of r(x) and the simply-connected union of the strips and S2's checked windows (extending S2 Lemma 1); W2, the Jost continuation and decay on the curved throat contours and the near-horizon map (extending S2 Lemma 3); and routine W3–W9. Controls: N1–N5 do not certify (N1, N2: PS modes at l = 20, 50, 100 and an h-box holding the l = 20–50 modes; N3: NE modes at l = 1, 2 by the continuation path, while mode-free boxes certify; N4, N5: exits misplaced into the allowed region); K2: W_l from the L code agrees with zsol and arbw at 12 points (relative radii ≤ 1.5e−38); K4 (float only, not used by the proof): W_l ≠ 0 at 200 points; the secant search stopped in H⁺ at 61 of them, none a zero (residual and zsol); the launch-time rule counted such stops as zeros, and Amendment 6 (post hoc) replaced it. Pre-registered (notes/s4-preregistration.md) with six dated amendments. Independent review 2026-10-18: proven as stated. Every certificate was replayed in full; reviewer code sharing nothing with p2 confirms every T, A, S and L bound (each certified bound lies below the reviewer's value; L at actual l from 20 to 10⁵; both closing boxes by the exact near-horizon map); S3's counting pipeline with the P2-2 threshold independently finds no zero of W₂₀ in R′₂₀; the closing-box truncation term was corrected after the review. The Im ω > 0 part is at P2-0's level | problems/P2/src/p2/part1.py, part1_T.py, throat.py, lowfreq.py, affine.py; problems/P2/certificates/part1/; problems/P2/notes/s4-preregistration.md, s4-proof.md, s4-report.md; problems/P2/tests/test_part1.py; docs/review-2026-10-18-p2-s4-adversarial.md; problems/P2/certificates/review_s4/ | 2026-10-18 |

**P2-1w (optional note, from the review).** Append to its statement column: "Contained in P2-1
(2026-10-18)." Its level and evidence are unchanged.

## 8. Commits

| commit | subject |
|---|---|
| `d00f9a3c` | the pre-registration (before any production run) |
| `c8bd7436` | Amendments 1–3 (coordinator) |
| `8f189b18` | Amendment 4: linear T ground state |
| `e639d96e` | Amendment 5: L and S methods fixed in development, before their production runs |
| `35201758` | the checkers (part1, part1_T, throat, lowfreq, affine) and tests; the production commit |
| `ca66270b` | S assembly checks the tail regions lie in the strips (W2) |
| `776aa7ef` | certificates T, A, S±, L±, bulk, cover |
| `a16b52cb` | Amendment 6 (post hoc): when K4's secant search has found a zero |
| `53fa7080`, `0e3c40f7` | K4's residual and zsol checks as a worker; zsol's term cap raised for that check |
| `279878d4`, `b24e3087` | the cover checks L's tiling box by box |
| `56f8d43c` | point lookups at contour nodes (no certificate changes) |
| `5b2b3a41` | controls and replay |
| `c8bf85bf` | proof, report and FORMULATION v0.6 |
| `8761ae07`, `930587a6`, `0e329217` | the report's commit table and wording |
| `2a72dd6b` | the independent review (reviewer's commit; fast-forwarded onto this branch) |
| `68de8184` | review fixes S1, W4: the closing box's truncation term over complex ρ; outward printing; `verify full`; the S1 test |
| `0b883e50` | S leaves regenerated at `68de8184`; L re-assembled; cover, controls, full replay |
| (this commit) | review fixes S2–S5, W1–W3, W5–W8: the Amendment 6 correction and errata, proof, report, FORMULATION v0.6 |

- **Which commit each certificate was computed at:**
  - **bulk, T, A and L's 1760 boxes:** at `35201758`, by workers launched there (see the launch
    records).
  - **L_PB_plus.json and L_PB_c.json** were re-assembled at `68de8184` from those rows, with the
    summary extremes taken from ball endpoints (W4).
  - **S's 92 leaves (both sides, including the two closing boxes):** regenerated at `68de8184`, by
    workers launched there (`launch_S+_20261019T130338.json`, `launch_Sc_20261019T130338.json`).
    The first S rows (at `35201758`) are kept, uncommitted, in `results/part1/archive-35201758/`.
  - **cover, controls and replay:** at `68de8184`.
- **The "+dirty" labels.** Some S+ and L rows record `35201758+dirty`. Those rows were written while
  an edit to `throat.assemble` was uncommitted (it became `ca66270b`).
  - The workers had loaded the code at launch, and never call `assemble`.
  - The replay at `0e3c40f7` reproduced the sampled S and L leaves exactly, and the full replay at
    `68de8184` reproduces all L boxes, including the two dirty ones. (The S rows are now all clean
    `68de8184` rows.)
- The first cover, controls and replay ran at `279878d4` and `0e3c40f7`.
- **What the later code commits change:** only assembly checks, the controls' evaluation, and a
  lookup that production never reaches.

## 9. Reproduction (cwd `problems/P2/src`)

```
python -m p2.part1 bulk                      # bulk_PB.json (<1 s)
python -m p2.part1 launch L 2                # detached workers (nice 10), launch record in results/part1/
python -m p2.part1 launch S+ 1; python -m p2.part1 launch Sc 1
python -m p2.part1 launch T 1; python -m p2.part1 launch A 1
python -m p2.part1 assemble L; python -m p2.part1 assemble S+; python -m p2.part1 assemble Sc
python -m p2.part1 assemble T; python -m p2.part1 assemble A
python -m p2.part1 cover                     # cover_PB.json
python -m p2.part1 launch K4 4               # float sample (about 30 min wall)
python -m p2.part1 controls                  # controls_PB.json (N1-N5, K1-K4)
python -m p2.part1 verify full               # replay_PB.json (every cell, leaf and box; about 9 min)
python -m pytest problems/P2/tests/test_part1.py -q        # from the repo root; P2_SLOW=1 for the slow set
```
