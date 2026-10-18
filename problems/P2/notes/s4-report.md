# P2 S4 report: part (i) at P_B, uniformly in l ≥ 20

Status: S4 exit, 2026-10-18. Branch `p2/s4-part-i`, from main at `dc5bfd27`.
Pre-registration: `notes/s4-preregistration.md`, written at `d00f9a3c` before any production run.
- Amendments 1–3 (requested by the coordinator) were committed at `c8bd7436`, Amendment 4 at
  `8f189b18` and Amendment 5 at `e639d96e`, each before the production runs it affects.
- Amendment 6, at `a16b52cb`, is post hoc and disclosed as such. It defines when K4's float secant
  search has found a zero (§3).

Proof: `notes/s4-proof.md`.

**Verdict: PASS, l₀ = 20.**
- **The theorem.** At P_B, for every l ≥ 20 and every ω with Re ω ≥ 0 and −d_δ ≤ Im ω ≤ 0
  (d_δ = 0.55κ₋), W_l(ω) ≠ 0.
- **P2-1 at P_B.** With the conjugation symmetry and P2-0, W_l has no zero in S_δ for every l ≥ 20.
  The level is WRITTEN_PROOF:
  - the inequalities are certified in Arb;
  - W1 (vertical-strip analyticity) and W2 (Jost continuation and decay on the curved throat
    contours) are non-routine written arguments, as in P2-1w;
  - the Im ω > 0 part is at P2-0's level.
- **P2-3 at P_B is complete.** S3 certified l ≤ 19 (P2-3B), and l₀ = 20, so the S3 fallback was not
  needed.
- **Status.** No kill criterion fired. Compute was ≈ 3.9 core-h of the 48 core-h cap (projection 9).

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
- **The closing boxes.** Margins are 0.414 (+) and 0.353 (c), with design w = 2, y_a = 1/20,
  φ_∞ = −π/3. The near-horizon model has ρ_max/R = 0.48 (+) and 0.58 (c), with order 14.
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
- **S.** Half-line coercivity on C² contours below each turning point, in the normalised form
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
| **K4** | At 200 deterministic points (64 L, 52 S, 37 T, 47 A; l ∈ [20, 400]), the float W_l is nonzero everywhere; see below. |
| **N1** | The T checker on boxes holding the PS fundamental mode does not certify, at l = 20, 50, 100: IN (±γ₀/50) gives μ/a* = −0.39, −0.070, −0.070, and STRAD gives −0.42, −0.100, −0.100. The STRAD value −0.100 is exactly E₁(0)/a* at ν = −1.05γ₀: the checker's threshold sits at ν ≈ −0.95γ₀, 5% short of the mode. |
| **N2** | The h-box [1/50.5, 1/20.5] with the ζ-box [Ω_P − 0.001κ₋, Ω_P] and the ν-box [−1.69, −1.68]κ₋ does not certify (μ/a* = −0.38). |
| **N3** | The L continuation path at l = 1, 2 does not certify the NE-mode boxes: Re y₊ ∈ [−0.232, −0.222] and [−0.316, −0.293], a clean sign failure. The same checker certifies mode-free boxes at half the depth (Re y₊ ≈ 0.146 and 0.164). |
| **N4** | The S checker with the exit in the allowed region does not certify: Ũ(x_a) − 1 ≈ −0.62 (+) and −0.66 (c), and min q̂ ≈ −0.81 and −0.92. |
| **N5** | L's production path with the exit in the allowed region does not certify, in 4 of 4 boxes. In the two boxes away from the real axis, both the y-sign test and the exit rule fail. In the two with ν ∈ [−0.05, 0]κ₋, only the exit rule fails: there Re y ≈ Im ω + O(U) keeps its sign in the allowed region. So the exit rule (the forbidden bulk) is essential to L's argument, as the proof says. |
| **Replay** | Agrees exactly: 1878 T cells and 78 A cells (5% plus each box's worst cell, plus the tails), 3 S leaves per side recomputed from their stored designs, 88 L boxes, the bulk and the cover. |
| **Tests** | `tests/test_part1.py`: fast 8 pass and 2 skip (8 s); `P2_SLOW=1` 10 pass (195 s). |

**K4 in detail.**
- **The raw outcome.** The secant search (≤ 12 steps, seeded at each point) stopped inside H⁺ at
  61 points, all at large |ω| (46 A, 15 T; Re ω up to 4214κ₋).
- **None is a zero.** At every one, |W(root)|/|W(seed)| ∈ [1.0000008, 1.0000013]. The rigorous zsol
  ball excludes 0 at all 61, with zsol's Frobenius term cap raised for this check only.
- **Why the search stopped.** Its relative-step test fires where W's phase turns fast. The criterion
  for a found zero was defined post hoc (Amendment 6), after three such stops among the first 14
  points.
- **A side finding.** The first residual run stopped at zsol's 6000-term cap. Its launch record is
  kept, and the cap was raised as above.

## 4. Deviations from the pre-registration (all dated amendments)

| # | what | when | why |
|---|---|---|---|
| A1–A3 | W1 per contour (`w1_check`), N5, the ζ → 0 wording | before production | coordinator's go message |
| A4 | T's ground state is linear in s, not √(q₀/c) | before production | √ of an interval at s = 0 |
| A5 | S's compact C² family, w grid extended by {2, 3, 4}; S refinement and cell details; the closing box's near-horizon model; L's affine arithmetic and 16 h²-sub-intervals; N3 and N4 set-up; K4 sampling; replay sampling | before the L, S and control runs | development showed each was needed (§1) |
| A6 | K4: a secant stop counts as a zero only if W is small there | after the first K4 results | the secant's step test stops where W's phase turns fast |

Two pilot numbers in §9 of the pre-registration did not survive:
- The S margin near the top was 0.087 in the pilot. The pilot used a softplus contour that is not
  exactly real at the exit. The compact family certifies 0.0098–0.0142 there.
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

At most 4 workers ran at once, all at nice 10. Launch records are in `results/part1/launch_*.json`.
The regime certificates cost far less than projected: S about 0.06 core-h against 3, and T seconds
against 0.3. The float K4 sample was the largest item.

## 6. What S4 does not establish

- **Independent review.** The S4 certificates have not been reviewed. The reviewer should attack
  first:
  - the S closing box (the near-horizon Taylor model, Rouché disc and Cauchy remainder);
  - W2 on the curved throat contours;
  - the affine arithmetic of L.
- **P_A.** Not attempted (§10 of the pre-registration).
- **The assembly** (P2-4, P2-5, the ledger rows): S5.
- **CERTIFIED for P2-1.** That would need W1 and W2 machine-checked (2D validated continuation of
  r(x)).

## 7. Ledger text (paste-ready; for the coordinator, not applied here)

**P2-1 (new row; P_B only).** P2-1w stays: it is a separate, reviewed result at three l.

| P2-1 | P2 | At P_B (ΛM² = 3/50, Q = (999/1000)Q_max, δ = 1/20), uniformly in l ≥ l₀ = 20: W_l(ω) ≠ 0 for every ω with Re ω ≥ 0 and −(½+δ)κ₋ ≤ Im ω ≤ 0. With the conjugation symmetry W_l(−ω̄) = conj W_l(ω) and P2-0 (Im ω ≥ 0), W_l has no zero in S_δ for any l ≥ 20. With P2-3B (l ≤ 19), P2-3 is complete at P_B (the assembly P2-4 is S5). Not P_A | WRITTEN_PROOF. Four regimes cover {l ≥ 20} × the closed half-strip; with h = 1/(l+½) ≤ h₀ = 2/41 and ζ = h Re ω, the cover's inequalities and tilings are certified in certificates/part1/cover_PB.json. **L** (Re ω ≤ 8κ₋): τ-scaled Frobenius series at both horizons (r = r_h ± h²τ), in Arb affine arithmetic over 1760 ω-boxes and 16 sub-intervals of h² ∈ [0, h₀²]; Re y₊ ≥ 0.1245 and Re y_c ≤ −0.0804 at the exits (y = Φ′/Φ); then Riccati invariance through the forbidden bulk. **S** (6κ₋ ≤ Re ω ≤ 0.9Ω_P/h): half-line coercivity on C² contours below each turning point, uniform in h through ζ and 1/Re ω; validated curved-contour integration; a near-horizon Taylor model for ζ ≤ 0.02Ω_P; min q̂ = 0.0130 (+) and 0.0091 (c) over 47 and 45 ζ-leaves. **T** (0.8–2.0 Ω_P/h): S2's rotated contour with the exact splitting h²E = E₀ + hE₁ + h²E₂ and a linear ground state; E ≥ 0.1003·2Ω_Pγ₀/h for every h ≤ h₀. **A** (≥ 1.5Ω_P/h): E ≥ 6.955·2Ω_Pγ₀/h at 1.5Ω_P/h, and monotone in Re ω above. Machine-checked (Arb, 128-bit): every cell, tail, exit and order condition, the bulk facts (U monotone on each side of r_P; W̃ ≥ −0.00792) and the W1 crossing checks. Written (notes/s4-proof.md): W1, the vertical-strip analyticity of r(x) (extending S2 Lemma 1); W2, the Jost continuation and decay on the curved throat contours and the near-horizon map (extending S2 Lemma 3); and routine W3–W9 (half-line coercivity, the C¹ ground state, Riccati invariance, the scaled-Frobenius identification, monotonicity, the splitting, the cover). Controls: N1–N5 do not certify (N1, N2: PS modes at l = 20, 50, 100 and an h-box holding the l = 20–50 modes; N3: NE modes at l = 1, 2 by the continuation path, with a clean sign failure, while mode-free boxes certify; N4, N5: exits misplaced into the allowed region); K2: W_l from the L code agrees with zsol and arbw at 12 points (relative radii ≤ 1.5e−38); K4: 200 float points, W_l nonzero, and no float zero in H⁺ (61 secant stops at large ω were non-zeros by residual and by zsol; criterion defined post hoc, Amendment 6). Replay exact on a 5% sample. Pre-registered (notes/s4-preregistration.md) with six dated amendments; the sixth, a post-hoc definition of K4's float zero, is disclosed. Not yet independently reviewed. The Im ω > 0 part is at P2-0's level | problems/P2/src/p2/part1.py, part1_T.py, throat.py, lowfreq.py, affine.py; problems/P2/certificates/part1/; problems/P2/notes/s4-preregistration.md, s4-proof.md, s4-report.md; problems/P2/tests/test_part1.py | 2026-10-18 |

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
| (this commit) | the report's commit table |

- **The regime certificates** (bulk, T, A, S±, L±) were computed by workers launched at
  `35201758`; see the launch records.
- **The "+dirty" labels.** Some S+ and L rows record `35201758+dirty`. Those rows were written while
  an edit to `throat.assemble` was uncommitted (it became `ca66270b`).
  - The workers had loaded the code at launch, and never call `assemble`.
  - The replay at `0e3c40f7` reproduces the sampled S and L leaves exactly.
- The cover was recomputed at `279878d4`, and the controls and replay ran at `0e3c40f7`.
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
python -m p2.part1 verify                    # replay_PB.json
python -m pytest problems/P2/tests/test_part1.py -q        # from the repo root; P2_SLOW=1 for the slow set
```
