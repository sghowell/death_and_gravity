# Independent adversarial review: P2 stage S2 (part-(i) barrier-top probe, route M3)

Date: 2026-10-17. Reviewer: a separate agent with no part in authoring P2 S2.

Scope: branch `p2/s2-part-i` at `6b933896`. The four S2 commits are:
- `0d115788`: pre-registration;
- `c4159e7f`: route M3;
- `f950b34a`: the KAT-6 test fix;
- `6b933896`: certificates, proof note and report.

The review read the author's files without changing any of them. The reviewer's scripts, outputs and
README are in `problems/P2/certificates/review_s2/`; the script names below refer to that directory.

**Overall verdict:** no flaw makes the certificate unsound.
- **The proposed row is proven as stated.** The window statement holds at l = 20, 50 and 100:
  - every step of the analytic chain is correct;
  - the Arb implementation is sound;
  - an independent re-implementation sharing no code with `m3.py` reproduces every margin.
- **The level should be WRITTEN_PROOF.** The proposed "CERTIFIED (Arb) … with the written
  implication at WRITTEN_PROOF pending review" is not a README level. Two of the written steps are
  non-routine prose: the continuation into the complex sectors and the decay of the Jost solutions on
  the rotated line.
- **GO for S3 follows from the pre-registered rule.** The projected l₀ = 2, however, is the formal
  output of a mis-specified fit. It has no predictive content and is described too strongly.
- **Nothing is blocking.** The should-fix items concern disclosure, wording, one immaterial float
  leak and test coverage.

## 1. Verdicts

| item | verdict |
|---|---|
| 1a. QNM ⇒ decaying solution on the rotated line; sector analyticity (A0) | PROVEN. The direction used (W_l = 0 ⇒ an L² solution) holds at both ends for every ω in each window. Lemma 3 has two exposition gaps (S5) |
| 1b. Tails, Jost decay, boundary terms | PROVEN. Independent tail discs: 0.00904 ≤ 0.00956 and 0.345 ≤ 0.351 |
| 1c. Coercivity identity | PROVEN. Every sign and constant was checked in sympy |
| 1d. Ground-state step | PROVEN (sympy). Its hypotheses hold: c ≥ 0.1199; G′ is continuous, piecewise linear and compactly supported |
| 2. Arb implementation | SOUND. It covers every ω in each window through the ω² hull, not sample points; the cells cover [−64, 32], with separate tail checks; all enclosures are rigorous. There is one float leak (m3.py:398–399) of relative size 7e-17, which is immaterial |
| 3. Negative control A2 | VALID BUT WEAK. Branch (a) runs the same code path. Branch (b) never reaches the checker. A sharper control designed near the mode (§5) fails on every box containing the mode, by −0.10 to −0.17 μ_HO. It certifies a mode-free box 0.1γ₀ above the mode at l = 20 and 50 |
| 4. P2-2 overlap | JUSTIFIED at all three l, with certified margins 4.30, 2.99 and 1.008 κ₋. a_100 replays (1.881 < 2). At l = 100 the actual deviation is 0.012, against the certified bound 0.881 |
| 5. Kill criterion, l₀ fit, M2 → M3 | GO SUPPORTED. The projection is overstated in the wording (S2). M3 was recorded before any production run, as git confirms. Two pre-registered criteria were implemented differently without disclosure (S3) |
| 6. Replays | ALL REPRODUCE. `p2.m3 verify` OK; pytest 47 passed / 5 skipped; P2_SLOW=1: 52 passed |
| 7. KAT-6 test fix | CORRECT |
| 8. Level | WRITTEN_PROOF. The P2-2 amendment (a_100) is acceptable with its evidence path stated |

## 2. Findings by severity

### Blocking

None.

### Should-fix before the merge

- **S1. The ledger level** (s2-report.md:137).
  - "CERTIFIED … the written implication … at WRITTEN_PROOF level, pending independent review" puts
    two levels in one row. The README gives a row one level: CERTIFIED only when "any written lemmas
    are routine (otherwise the row is WRITTEN_PROOF)".
  - Lemmas 1 and 3 of s2-m3.md are the M3-specific analysis that s0-gate §4 named as M3's main risk.
    They are correct, but they are not routine.
  - Use the WRITTEN_PROOF text of §10.

- **S2. The l₀ projection and the loss attribution** (s2-report.md:89–97; FORMULATION.md v0.4, the
  "S2 result" bullet).
  - The certified m(l) = 0.528, 0.568, 0.560 is the sum of design losses with opposite l-dependence
    (§4, table). Fitting m_∞ − C/(l+½) through three such points predicts nothing:
    - with the design fixed, m(l) eventually falls like −√l rather than tending to 0.577;
    - nothing was computed below l = 20.
  - "The barrier-top regime poses no obstruction at any l ≥ 20, and on this fit at any l ≥ 2"
    overstates three certified points.
  - "The certificate loses about 0.11–0.13 μ_HO to its discretization (spline ground state,
    mean-value cells)" is half right. The mean-value cells cost < 1e-4 μ_HO. The losses are:
    - the 5% design slack;
    - one ground state for the whole ω² hull;
    - the fixed node spacing of the spline.
  - Recommended wording: "projected l₀ = 2 is the formal output of the pre-registered fit".

- **S3. Two pre-registered criteria were implemented differently, without disclosure.**
  - **A0** (s2-preregistration.md:68–71, 90–94) was to be "a fan of validated rays", "ray spacing below
    half the certified analyticity radius" and "the funnels are certified".
    - The implementation is a different argument: the Re r\* identity with Σ A_i = 0, and two
      validated vertical segments.
    - That argument is valid and simpler, but the report presents it as the A0 that was
      pre-registered.
    - `Ray.analyticity` (m3_contour.py:158) is the unused remnant of the planned fan.
  - **A3** (s2-preregistration.md:100–101) compares μ_cert with "the float best margin λ_min(l)".
    - The code computes the float λ only at the low corner (m3.py:535, 657).
    - The worst corner is the low right corner, (a₂, −β): float 0.648, 0.661 and 0.667 μ_HO against
      0.661, 0.681 and 0.679 at the low corner.
    - A3 still passes: μ_cert lies below the float λ at all four corners and along the bottom edge
      (`indep_float_spectrum.py` B).
  - The pre-registration's design-study figure of 0.58 at l = 20 (s2-preregistration.md:38–42) is not
    reproduced by the committed design, which gives 0.648 by the reviewer's computation.
  - All three points should be recorded in s2-report §2.

- **S4. A2's strength is overstated** (s2-report.md:66–69; s2-m3.md:148–153).
  - "The method therefore cannot certify a true QNM" overstates one negative test.
  - Branch (b) always fails in the design step (Riccati overflow, m3.py:212), so it never runs the
    checker.
  - Branch (a) fails by −0.11, −0.064 and −0.095 μ_HO. A checker optimistic by less than that would
    pass it too.
  - Reword the claim, and add the near-mode control of §5 (or the reviewer's script) to the tests.

- **S5. Gaps in the written proof** (s2-m3.md:80–91 and 138–139). All the missing facts are true and
  were checked here; they should be written in.
  - **Lemma 3** should state the convergence radii of the Frobenius factors: g₊ converges for
    |r − r₊| < r₊ − r₋ = 0.0970 > 0.00956, and g_c for |r − r_c| < 4.77 > 0.351.
  - **Lemma 3** should also say that on {Re x ≤ −47.69} the map r(x) = r₊ + ψ(e^{2κ₊x}) is
    single-valued.
    - The image of the contour winds round r₊: arg(r − r₊) = −3.17 at s = −64.
    - So the representation Φ₊ = c·r·g₊·e^{−iω(x − c₊ − R̃₊)} must be read as a function of x.
  - **The tail check (§5)** uses two facts:
    - |φ′(s)| ≤ |φ′(S)|, because |φ′| decreases for |s| > (3/5)^{1/4}s₁ = 4.40;
    - φ(s) lies between φ(S) and φ_∞.
  - **The citation of the S1 LG identification and β < κ₊ (s2-m3.md:80–81) is unnecessary:**
    W_l = 0 already gives R₊ ∝ R_c.

- **S6. A float round trip in `tail_bound`** (m3.py:398–399).
  - `float(t.lower())` rounds to the nearest double. The added −1e-30 is far below one ulp, since the
    terms are of size 5–40.
  - The left disc radius therefore comes out 6.8e-17 (relative) below its rigorous value:
    0.0095571080701487445 against 0.0095571080701487452 (`float_leak_tail_bound.py`).
  - The effect is immaterial: with the disc enlarged by 1e-12, every tail bound is unchanged (13.195,
    96.10, 403.42 on the left; 1.068, 22.58, 112.30 on the right).
  - Fix it with exact Arb min/max, as in `indep_cells_arb.tail_disc`.

- **S7. Test coverage** (test_m3.py:59–70).
  - `test_stored_controls_and_overlap` and `test_summary_fit_and_verdict` only read stored flags.
  - The l = 100 P2-2 criterion is replayed by `p2.m3 verify` but by no pytest, and not by
    `p2.lg_bounds verify`, which reads only `certificates/lg/`.
  - A2 is re-run only at l = 100, under P2_SLOW.

- **S8. The P2-2 amendment must name its evidence path.** It should cite
  `certificates/m3/lg_PB_l100.json` and its replay command (§10).

### Wording

- **W1.** s2-m3.md:136: "the hull adds no overestimation in ω²".
  - Interval evaluation over the rectangle is exact, but the rectangle is larger than the image of
    the window.
  - In floats this costs about 0.02–0.03 μ_HO (min of E over the hull corners against the true window
    corners).
- **W2.** s2-m3.md:158: "the first PS overtone layer" should be "the PS fundamental (n = 0) layer at
  Im ω ≈ −γ₀". The first overtone sits near −3γ₀; the reviewer's complex-scaled spectrum finds it at
  −5.05κ₋.
- **W3.** FORMULATION.md:3 still says v0.3, while the log ends at v0.4 (line 399). The v0.3 entry
  updated the header, so this one should too.
- **W4.** The commit table at s2-report.md:146–149 omits `f950b34a`. It also lists "KAT-6 replay
  hygiene" under `6b933896`, whose subject does not say that.
- **W5.** The docstring at winding.py:20 still gives the pre-fix bound 2Mr/(R − r)²; the code uses
  the corner distance √2·r. This is left over from S1.

## 3. Item 1: the analytic chain

### 1a. QNM ⇒ decaying solution; the sectors

**Lemma 1 (sector analyticity).** Each step checks.
- **The Re r\* identity.** Along any continuation, d/dt Re r\*(r(t)) = Re x′(t), because r\*′ = 1/F.
  The roots r_i are equilibria of dr/dx = F and are never reached. Hence Re r\*(r(x)) = Re x.
- **Σ A_i = 0.** This holds because r²/Δ = O(r⁻²). So Re r\* → 0 as r → ∞, and Re r\* → Re r\*(0)
  as r → 0.
- **Roots and logarithmic branch points are not singularities of r(x).** The four roots of Δ
  (including r_n < 0) are the logarithmic branch points of r\*(r). They are equilibria of dr/dx = F,
  reachable only as Re x → ±∞, since A_i log|r − r_i| diverges. The r ↔ r\* map therefore has only
  two finite-x obstructions:
  - r = 0: a cube-root branch point of r(x), since r³ ≈ 3Q²(x − x₀); V_l has a pole there;
  - r = ∞: a simple pole of r(x); V_l has a double pole there.
- **Where the continuation can fail.** A vertical continuation that fails at a finite height must
  leave every compact subset of ℂ∖{0}. It then tends to 0 or to ∞, which forces
  Re x ∈ {Re r\*(0), 0}.
  - On the Riemann sphere, r = ∞ is a regular point of the ODE.
  - V_l, however, has a pole there (V ~ (2Λ²/9)r²), so excluding it is necessary and correct.
- **Σ₋.** It lies left of x_P = −2.43797 < Re r\*(0) = −2.29664, with margin 0.14133 (reproduced).
- **Σ₊.** It meets the two special lines at heights 0.14133 and 2.43797.
  - The author validates both segments. Re r\*(0) is a ball, but `real_point`'s interval Newton
    encloses r(X′) for every X′ in it, so the exact segment is covered.
  - Reviewer's float values: |r| ∈ [1.98262, 1.98276] and [2.45138, 2.45197], inside the author's
    validated [1.961, 2.004] and [2.412, 2.491].
- **The nearest singularities are far away.** Continuing the same two lines from the real axis, the
  singularities lie at:
  - Im x = ±7.947 (r → 0; this equals π(A₋ + A₊));
  - Im x = ±16.82 (r → ∞; this equals −πA_c = 16.821). The float event at |r| = 10³ fires at 16.771,
    which is (3/Λ)/|r| = 0.05 earlier.
  - Both are far outside the sectors. A fan of 180 vertical lines covering both sectors meets no
    singularity (`indep_geometry.py`).
- **The Cauchy–Riemann step is correct.** G = ∂_X R − F(R) solves ∂_Y G = iF′(R)G with G(X, 0) = 0.

**Lemma 2 (the contour image).** The independent check agrees closely.
- The reviewer enclosed r(s) by a different method: Krawczyk on the branch-tracked equation
  r\*(r) = x_P + e^{iπ/4}s.
- At all 3,074 points of spacing 1/32 on [−64, 32], these balls overlap the author's Taylor-model balls.
  The author's radii are ≤ 1.5e-23, and the midpoints agree to 2e-24.
- The Taylor model itself is correct:
  - the a-priori box is a valid Picard enclosure;
  - the componentwise Lagrange remainder uses r_{N+1}(B), which is valid because the ODE is
    autonomous;
  - the flow derivative is ∂r/∂r₀ = F(r(t))/F(r₀);
  - the step lookups by float are safe, because cells are split at multiples of 1/16 in exact
    rationals and never cross a step.

**Lemma 3 (QNM ⇒ decaying u).**
- If W_l(ω) = 0 then R₊ ∝ R_c, and Φ₊ continues through Σ₋ ∪ Σ₊.
- On {Re x ≤ −47.69} the Frobenius representation holds, with the radii given in S5. It gives
  |u| ≲ e^{s(Re ω + Im ω)/√2} (sympy item 10), and the same bound for u′.
- The cosmological end is identical.
- Re ω + Im ω ≥ a₁ − β > 0 is certified at every ω of each window, not only at the centre. The decay
  rate is at least 61κ₋ (from 86.82κ₋/√2).

**Match with D1/D1′.**
- P2-1w is a statement about zeros of W_l (D1), and uses only the implication "W_l = 0 ⇒ an L²
  solution on the contour".
- D1′ (WRITTEN_PROOF since S1) is not used. It is needed only to read "QNM-free" as "no localized
  resonance", at D1′'s level.

**Float cross-check of the sign and direction conventions** (`indep_float_spectrum.py` A). The
reviewer computed the eigenvalues of −d²/ds² + iV_l on the rotated line (finite differences,
h = 0.004).
- **The PS fundamentals are reproduced.** The S1 values are recovered to ≤ 3.1e-6 κ₋ at
  l = 20, 50 and 100.
- **The first overtones appear** near Im ω = −5.05κ₋ ≈ −3γ₀.
- **No eigenvalue lies in any window.** The nearest is the PS fundamental, 1.13κ₋ below.
- So complex scaling along this contour reveals the QNMs, with the conventions of Lemma 3.

### 1b. Tails and boundary terms

- **Tail discs.** The reviewer bounded them with the triangle inequality in Arb (`indep_geometry.py` T):
  - the discs: 0.00904 (left) and 0.345 (right), against the author's 0.00956 and 0.351;
  - the start points, by Krawczyk: |r(−64) − r₊| = 0.00575 and |r(32) − r_c| = 0.187;
  - a float continuation to s = −400 and s = 200 stays inside the discs, with maxima 0.00575 and
    0.187.
- **Boundary terms.**
  - In §3, [−u′e^{−iφ}ū] → 0 because u and u′ decay exponentially.
  - In Lemma 4, cgg′|v|² vanishes identically outside supp G′. The nodes span [−2.99, 2.96],
    [−1.85, 1.84] and [−1.31, 1.30], and G′ is exactly 0 at both end nodes.

### 1c. Coercivity identity

Checked in sympy (`indep_identities_sympy.py` items 1–2):
- the pointwise identity e^{−iφ}(−u″)ū = −(e^{−iφ}u′ū)′ + e^{−iφ}|u′|² − iφ′e^{−iφ}u′ū;
- the bound Re(−iφ′e^{−iφ}u′ū) ≥ −|φ′||u′||u|;
- Young's step with ε₁ = η|φ′|, which leaves the remainder |φ′|(√η A − B/(2√η))² ≥ 0.

So J(u) ≥ ∫ c|u′|² + q|u|², with c = cos φ − η|φ′| and q = Re(e^{−iφ}Q) − |φ′|/(4η). The signs and
constants are as written.

### 1d. Ground-state step

- **The identity is correct.** c|u′|² + q|u|² = cg²|v′|² + E|u|² + (cgg′|v|²)′, with
  E = q + c′G′ − cG′² + cG″, checks in sympy (item 3).
- **Its hypotheses hold:**
  - c ≥ 0.1199 on every cell (reviewer's run) and c ≥ 0.1328 on the tails;
  - G′ is continuous and piecewise linear, with exact rational nodes and compact support;
  - cg′ is absolutely continuous;
  - E = q is bounded on the tails, and u, u′ ∈ L².

## 4. Item 2: the implementation

- **ω enclosure.**
  - E and E′ depend on ω only through ω², and affinely.
  - Both are evaluated on the rectangle W2 ⊇ {ω² : ω ∈ 𝒲_l}: Re ω² ∈ [a₁² − β², a₂²] and
    Im ω² ∈ [−2a₂β, 0]. The endpoints are taken from outward-rounded ball endpoints.
  - The mean-value bound combines inf over W2 of E(s_m) with sup over W2 × cell of |E′|.
  - So μ_cert holds for every ω in the window, not at sample points.
- **Coverage.**
  - The cells partition [−64, 32] exactly. Each popped cell is either accepted or has both halves
    pushed.
  - The two tails are checked separately.
- **Enclosures.**
  - V_l and V_l′ are evaluated on acb balls, with Δ in factored form on the certified roots.
  - r(s) comes from the validated Taylor model, confirmed independently above.
  - G′ is exact, and G″ is exact on each cell.
  - The φ-derivatives use `arb_series` on intervals.
  - sympy confirms every coded formula: E′, c′, c″, q′, Q′, V, V′ and F (items 4–8).
- **Parameters.**
  - They are rebuilt by `params.Point`, which is the S1 code that writes and verifies `params.json`.
  - The reviewer's 60-digit mpmath values lie inside the `params.json` balls.
- **Float leaks.**
  - (i) `tail_bound`: see S6.
  - (ii) m3.py:545 picks the minimum of three bounds by a float key. A tie within one ulp does not
    occur (0.32, 1.07 and 13.2).
  - (iii) The float step lookup is safe, as explained under Lemma 2.
  - (iv) The construction `(lo+hi)/2 + arb(0,1)·(hi−lo)/2` (for example `_iv`, m3.py:128) is ball
    arithmetic throughout.
    - arb(0, 1) has radius 1 + 2⁻²⁹, so the construction widens outward and never misses an endpoint.
    - The sibling problem's pattern, arb(mid, rad) built from rounded floats, does not occur in S2.
- **Independent re-check** (`indep_cells_arb.py`). It shares no code with `m3.py`:
  - parameters parsed from `params.json`;
  - r(s) from point and interval Krawczyk enclosures;
  - V_l from the expanded polynomial;
  - derivatives by forward dual numbers;
  - the reviewer's own cover, hull and tail discs.

  Its results:

  | l | E ≥ (reviewer) | μ_cert (author) | tails left / right (reviewer) | c ≥ |
  |---|---|---|---|---|
  | 20 | 0.3245483 | 0.3245472 | 13.0 / 12.3 | 0.11999 |
  | 50 | 0.8599116 | 0.8599105 | 94.9 / 90.7 | 0.11999 |
  | 100 | 1.6853896 | 1.6853881 | 398.7 / 382.2 | 0.12000 |

  - Bisection stopped once a cell reached the author's μ_cert, so the reviewer's minima sit just
    above it. No cell stayed below.
  - At 148 sample points the two implementations agree on E and E′ to ≤ 1.5e-7. That is the outward
    widening of the hull radius.
- **Where the margin goes** (float; `indep_float_spectrum.py` B and C, plus a scan of the node
  endpoints).
  - Each certified μ equals, to five digits, the float minimum of E at the left end of a spline
    interval, where cG″ jumps.
  - In units of μ_HO:

  | l | float λ at the low corner | × 0.95 (design slack) | min of E at interval midpoints (one ground state for the hull) | min of E at node ends (spline sawtooth) | μ_cert |
  |---|---|---|---|---|---|
  | 20 | 0.661 | 0.628 | 0.551 (−0.077) | 0.528 (−0.023) | 0.528 |
  | 50 | 0.681 | 0.647 | 0.616 (−0.031) | 0.568 (−0.048) | 0.568 |
  | 100 | 0.679 | 0.645 | 0.630 (−0.015) | 0.560 (−0.070) | 0.560 |

  - The hull loss falls like 1/l. The sawtooth grows like √l at the fixed node spacing of 0.01.
  - This explains the non-monotone m(l), and tells S4 what must scale:
    - node spacing ∝ (l+½)^{−1/2};
    - window width ∝ the oscillator quantum.

## 5. Item 3: the negative control A2

- **Branch (a)** runs the A1 code path (`certify_l` → `Checker.run` → tail checks) on an explicit box,
  using the window's ground state. It gives μ = −0.068, −0.097 and −0.285.
- **Branch (b)** designs a ground state for the box. Its Riccati sweep overflows (the RuntimeWarning
  in the slow tests, m3.py:212). The design fails, and the failure is counted as "not certified"
  (m3.py:597–598), so the checker never runs.
- **A2 catches only gross errors.** The window's ground state was designed 1.1κ₋ above the mode, so
  its failure there is expected for any checker optimistic by less than about 0.06–0.11 μ_HO.
- **A sharper control** (`sharp_negative_control.py`) runs the author's checker with a ground state
  designed at Im ω = −0.88γ₀:

  | l | float margin at the design point | box around the mode | box straddling Im ∈ [−1.05, −0.95]γ₀ | mode-free box Im ∈ [−0.90, −0.86]γ₀ |
  |---|---|---|---|---|
  | 20 | +0.101 | −0.102 (fails) | −0.131 (fails) | +0.017 (certifies) |
  | 50 | +0.116 | −0.113 (fails) | −0.143 (fails) | +0.007 (certifies) |
  | 100 | +0.119 | −0.142 (fails) | −0.172 (fails) | −0.022 (fails) |

  (μ in units of μ_HO; boxes of real half-width γ₀/50 about the mode's real part; 416 s of CPU in all.)

  - Every box containing the mode fails, with a near-optimal design, by a clear margin.
  - At l = 20 and 50, a mode-free box 0.1γ₀ above the mode certifies. So the checker's verdict changes
    within about 0.1γ₀ of a true QNM. The A2 failure is therefore not just an artefact of a poor
    design.
  - At l = 100 the mode-free box does not certify, because the 0.07 μ_HO spline loss of §4 exceeds the
    remaining float margin.
- **The float margin vanishes exactly at the modes.** It is −0.019, −0.004 and −0.001 μ_HO at the PS
  mode, which reproduces the author's figures. The complex-scaled eigenvalues of §3 sit at the same
  points.

## 6. Item 4: the P2-2 overlap

`indep_p22_overlap.py` checks the overlap with the reviewer's own window, built from the
`params.json` balls:

| l | a_l source | a_l/κ₋ | window top/κ₋ | certified a_l < top | margin |
|---|---|---|---|---|---|
| 20 | certificates/lg/lg_PB_l20.json (S1) | 106.640 | 110.935 | yes | 4.296κ₋ |
| 50 | certificates/lg/lg_PB_l50.json (S1) | 255.507 | 258.500 | yes | 2.994κ₋ |
| 100 | certificates/m3/lg_PB_l100.json (S2) | 503.434 | 504.442 | yes | 1.008κ₋ |

- **The thresholds match.** The m3 certificates carry exactly the a_l of these files.
- **The l = 100 threshold is sound.**
  - `lg_bounds.py` is unchanged since S1 (`a78a3057`).
  - Its criterion replays at the stored (a_100, r_X) as 1.88101480199 < 2.
  - P2-2's bound is uniform in a ≥ a₀ (s1-p2-2-lg.md, Lemma 3), so a single comparison per l closes
    the gap above the window.
- **The thin margin is not a soundness issue.** It is a certified inequality between an exact dyadic
  and a ball.
- **The actual Wronskian is far inside the bound.** An independent float integration along the real
  axis gives |W e^{iΘ}/(2iω) − 1| ≤ 0.0124 at l = 100 for Re ω ∈ {a_100, 1.01a_100, 1.05a_100,
  1.2a_100} and Im ω ∈ {0, −β/2, −β}, against the certified bound 0.881. This matches the S1
  review's 0.004–0.015 for l ≤ 50.

## 7. Item 5: kill criterion, fit and route switch

- **The rule.** GO if A0–A4 hold and the projected l₀ ≤ 100 (s2-preregistration.md:114–118).
  - All criteria hold, and l₀ = 2 (the C/(l+½) fit) or 5 (the C/(l+½)² fit).
  - The verdict would be GO under any reasonable reading: the certified margin at l = 100, where the
    kill rule bites, is 0.56 μ_HO.
- **The projection.** It is the formal output of the rule, and is computed as the rule says
  (⌈C/m_∞ − ½⌉).
  - §4 shows that the fit is mis-specified.
  - The report should say so, and drop "no obstruction at any l ≥ 20" (S2).
- **The M2 → M3 switch.** It happened before production, as far as git can show:
  - the pre-registration (`0d115788`, 10:11:04) records M3 first, and no M3 code had been committed
    at that point;
  - the M3 code followed at 10:30:22 (`c4159e7f`);
  - every certificate records commit `c4159e7f` with no `+dirty` flag;
  - the pre-registration has not changed since;
  - no M2 code exists;
  - the float design study came before the pre-registration and is disclosed there.
  - The A0 and A3 implementation changes (S3) are the only undisclosed deviations.

## 8. Item 6: replays

All replays were run at nice 10; the outputs are in `review_s2/replay_*.txt`.

| command | result |
|---|---|
| `python -m p2.m3 verify` (cwd `problems/P2/src`) | `m3 replay: OK`, 61.6 s. Every window's A1 verdict and μ_cert reproduce identically; the l = 100 LG criterion replays |
| `pytest problems/P2/tests -q` | 47 passed, 5 skipped, 11.6 s |
| `P2_SLOW=1 pytest problems/P2/tests -q` | 52 passed, 1 warning (the A2 design overflow), 199 s |

After all runs, `git status` showed no change to any tracked file. The KAT-6 and m3 tests no longer
overwrite certificates.

## 9. Item 7: the KAT-6 test fix (`f950b34a`)

- `run_kat6` gains `write=True`, and `test_kat6_rerun` now calls `run_kat6(write=False)`.
- The fix is correct and minimal. The command-line path still writes the certificate.
- The slow run confirmed that the certificate is no longer rewritten.
- The stale docstring at winding.py:20 is W5.

## 10. Item 8: the level, and paste-ready ledger text

**The proposed level is not coherent under the README.**
- One row carries one level.
- CERTIFIED requires the decisive steps to be machine-checked and any written lemma to be routine,
  "otherwise the row is WRITTEN_PROOF".
- A row whose written implication still awaits review is therefore not CERTIFIED, and "pending
  review" is not a level.

**WRITTEN_PROOF fits.**
- The written proof (s2-m3.md) is complete up to the S5 additions, and this review has checked it.
- Every constant and bound it uses is machine-checked:
  - the margins;
  - c ≥ 0;
  - x_P < Re r\*(0);
  - the two special segments;
  - the tail discs;
  - the decay condition;
  - the separation of the tail discs from the other roots.

**CERTIFIED does not fit.** Lemma 1 (continuation via the Re r\* identity and the two obstruction
lines) and Lemma 3 (Jost solutions on the spiralling tails) are non-routine prose. The precedent is
P2-0, which was rated WRITTEN_PROOF in S1 because its regularity and integration steps were prose.
The combined QNM-free statement also inherits P2-0's level for Im ω > 0.

**P2-1w (new row; replaces the text of s2-report §6):**

| P2-1w | P2 | At P_B (ΛM² = 3/50, Q = (999/1000)Q_max, δ = 1/20) and for l = 20, 50, 100: W_l(ω) ≠ 0 for every ω with (l+½)Ω_P − 8γ₀ ≤ Re ω ≤ (l+½)Ω_P + 6γ₀ and −(½+δ)κ₋ ≤ Im ω ≤ 0 (the barrier-top window; Re ω/κ₋ ∈ [87.37, 110.94], [234.94, 258.50], [480.88, 504.44]). Certified coercivity margins μ/μ_HO = 0.528, 0.568, 0.560, with μ_HO = 2(l+½)Ω_Pγ₀. With P2-2 (a_l below the window top: a_20, a_50 from the S1 certificates; a_100 = 503.43κ₋ from certificates/m3/lg_PB_l100.json) and P2-0, no zero of W_l lies in {Re ω ≥ (l+½)Ω_P − 8γ₀} ∩ S_δ at these l. Three discrete l only, so not part (i): there is no uniformity in l, and the sub-barrier and low-frequency regimes are open | WRITTEN_PROOF. Machine-checked (Arb, 128-bit): the ground-state inequality E ≥ μ on s-cells covering [−64, 32] and on both tails, for the whole ω² hull of each window; c ≥ 0.119; the contour enclosure; x_P < Re r*(0); the two special vertical segments; the tail discs; the decay condition Re ω + Im ω > 0. Written (notes/s2-m3.md): the sector analyticity of r(x) (Lemma 1); the continuation and decay of the Jost solutions on the rotated line (Lemma 3); the coercivity and ground-state integrations by parts (§3, Lemma 4). Independent review 2026-10-17: proven as stated. The margins were reproduced by an independent Arb implementation, the identities were checked in sympy, and the complex-scaled float spectrum reproduces the PS modes and has no eigenvalue in the windows. The Im ω > 0 part is at P2-0's level | problems/P2/src/p2/m3.py, m3_contour.py; problems/P2/certificates/m3/; problems/P2/notes/s2-m3.md, s2-preregistration.md, s2-report.md; problems/P2/tests/test_m3.py; docs/review-2026-10-17-p2-s2-adversarial.md; problems/P2/certificates/review_s2/ | 2026-10-17 |

**P2-2 amendment (accept).** Two changes to the existing row:
- **Statement:** after "for every 0 ≤ l ≤ 50", insert "and, at P_B, also l = 100
  (a_100 = 503.43κ₋)".
- **Evidence:** append "problems/P2/certificates/m3/lg_PB_l100.json (S1 code, unchanged since
  a78a3057; run at c4159e7f; criterion 1.881 < 2; replayed by `python -m p2.m3 verify`, not by
  `p2.lg_bounds verify`); independent review 2026-10-17: replayed; actual deviation 0.012 at l = 100".

The level is unchanged: CERTIFIED on the lower strip.

## 11. Before the S2 branch is merged into main

1. **Ledger text.** Use the §10 text for P2-1w and the P2-2 amendment. Correct s2-report §6 to match,
   or mark it superseded by this review.
2. **The projection (S2).** In s2-report §3 and FORMULATION v0.4:
   - describe l₀ = 2 as the formal pre-registered output;
   - replace "any l ≥ 20" with "at l = 20, 50 and 100";
   - replace the loss attribution with the decomposition of §4.
3. **Pre-registration deviations (S3).** Record the A0 and A3 implementation changes and the
   design-study discrepancy in s2-report §2.
4. **Proof note and report wording (S4, S5, W1, W2).**
   - Add the Frobenius-radius and single-valuedness sentences to s2-m3.md Lemma 3.
   - Add the |φ′| monotonicity and the range of φ to §5.
   - Drop the S1 identification citation.
   - Fix "first PS overtone layer" and "the hull adds no overestimation".
   - Soften the A2 claims in s2-m3.md §6 and s2-report §2.
5. **The float leak (S6).** Either fix m3.py:398–399 with exact Arb min/max (the window certificates
   replay unchanged), or record the leak and its size in the report.
6. **Housekeeping (W3, W4).** Bump the FORMULATION header to v0.4, and correct the s2-report §7 commit
   table.
7. **Optional (S7).** Add a pytest replay of `lg_PB_l100.json` and a near-mode negative control.

Under the repository's rules the merge into main must be run by the user.
