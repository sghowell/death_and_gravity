# Independent adversarial review: P2 stage S4 (part (i), the large-l exclusion at P_B, uniformly in l)

Date: 2026-10-18/19. Reviewer: a separate agent with no part in authoring P2 S4.

Scope: branch `p2/s4-part-i` at `0e329217`, from main `dc5bfd27`. The S4 commits are:
- `d00f9a3c` (10:31:36): the pre-registration;
- `c8bd7436` (10:39:09), `8f189b18` (10:52:04), `e639d96e` (11:52:57): Amendments 1–3, 4 and 5;
- `35201758` (11:53:14): the checkers and tests (the production commit);
- `ca66270b` (11:56:52): an assembly check (W2 tail strip);
- `776aa7ef` (11:57:38): the regime certificates;
- `a16b52cb` (12:01:46): Amendment 6 (post hoc), with the K4 residual check;
- `53fa7080`, `0e3c40f7`, `279878d4`, `b24e3087`, `56f8d43c`: K4 worker, zsol cap, cover tiling, point lookups;
- `5b2b3a41`, `c8bf85bf`, `8761ae07`, `930587a6`, `0e329217`: controls and replay, proof, report, FORMULATION v0.6.

The review read the author's files without changing any of them, and read the author's gitignored
working files (`results/part1/*.jsonl`) in place, read only.

**Where the review's scripts are.** The coordinator asked for them under
`problems/P2/certificates/review_s4/`. Creating that directory was **denied by the tool permission
system** (auto-mode classifier, reason "Interfere With Workloads"), and under the review's rules the
denial was not worked around. So this document is the only file the review adds to the repository.
The scripts and their JSON outputs are in the reviewer's session scratch directory,
`/private/tmp/claude-501/-Users-seanhowell-dev-research-death-and-gravity/5e36ca60-7404-4050-9ab7-cadfa3093ba6/scratchpad/rv/`,
which is not durable. §5 lists every script, what it does and its result, precisely enough to redo
each check. A second denied call is listed in §6.

**Overall verdict:** the S4 theorem is proven as stated, and P2-1 at P_B with l₀ = 20 holds at
WRITTEN_PROOF. Nothing is blocking. One machine-checked error term in the S closing box is not
rigorous as committed (S1); with it corrected, the review's re-run gives the same margins, so the
proof stands on the corrected computation, and the certificate should be regenerated before the
merge.
- **The theorem holds.** At P_B, for every l ≥ 20 and every ω with Re ω ≥ 0 and
  −(½+δ)κ₋ ≤ Im ω ≤ 0, W_l(ω) ≠ 0. With the conjugation symmetry and P2-0 (Im ω > 0), W_l has no
  zero in S_δ for any l ≥ 20.
- **How this was checked.**
  - Every analytic step was re-derived: the cover, W1–W9, the near-horizon model, the affine
    arithmetic.
  - Every certificate was replayed **in full** with the author's code at the clean head (not a 5%
    sample): 24 T boxes and A, 92 S leaves, 1760 L boxes. All agree field for field, including the
    three rows stamped `+dirty`.
  - The certified inequalities were re-evaluated by reviewer code that shares nothing with `p2`
    (float/mpmath), and every certified bound lies below the reviewer's value: in all 24 T boxes, A,
    all 92 S leaves (both closing boxes by the exact near-horizon map, no Taylor model) and all 1760
    L boxes (35 200 values at actual l = 20, 21, 100, 10⁵, plus 1872 more at l = 33 and 1000 on 52
    boxes).
  - **An independent rigorous count at l = 20**, S4's hardest l, with S3's reviewed pipeline and the
    certified P2-2 threshold a_20, shares no code or argument with S4: the winding number of W₂₀
    around ∂R′₂₀ is 0 (ball 0 ± 7.4e−4, 223 segments, \|W₂₀\| ≥ 4.07 on the right half). With P2-2 and
    P2-0, W₂₀ has no zero in S_δ, independently of every S4 regime.
- **One machine-checked bound is not rigorous as written (S1).** In the S closing box, the truncation
  error of the near-horizon Taylor model bounds |P_N(ρ)| on the *real* segment, while the contour
  visits *complex* ρ. Re-run with a rigorous bound, every cell bound moves by at most 8.2e−5 and
  both closing-box margins are unchanged to 8 digits (0.41468629, 0.35342094). The result stands; the
  certificate must be regenerated.
- **Level: WRITTEN_PROOF**, as the author proposes. W1 and W2 are non-routine prose (as in P2-1w), and
  the Im ω > 0 part is at P2-0's level.
- **Should-fix items:**
  - S1: the closing-box truncation term (above);
  - S2: Amendment 6's disclosure is inaccurate on two points (the launch-time code did define the
    criterion, under which K4 fails 61/200; and the stated mechanism is wrong). The new criterion itself
    is legitimate, and K4 carries almost no weight;
  - S3: W1/W2 as written do not state the fact that makes the c-side continuation consistent (the
    union of the strips and the checked windows is simply connected);
  - S4: ledger wording ("P2-3 complete" needs P2-3B's R′_l; the review sentence; P2-1w);
  - S5: the pilot discrepancy is mis-attributed (the pilot number was a leading-order, ε = 0 quantity).

## 1. Verdicts

| item | verdict |
|---|---|
| A1. Pre-registration against execution | FOLLOWED, with six dated amendments. l₀ = 20, h₀, every regime boundary, the controls N1–N4 and the kill criteria are as registered. Every design constant in the code at `35201758` matches the pre-registration as amended. The original text is unchanged (only "(none)" under Amendments was replaced) |
| A2. Amendments 4 and 5 before the runs they affect | YES. A4 10:52:04; A5 11:52:57; code 11:53:14; launches L 11:53:20, S± 11:53:21, T 11:54:40, A 11:54:50 (launch records, all clean `35201758`); the first certificate rows finish 11:53:21. A5 precedes production by 23 s; the development runs before it (≈ 0.5 core-h, not logged per run) are disclosed |
| A3. Amendment 6 | LEGITIMATE in substance, DISCLOSED INACCURATELY (S2). None of the 61 stops is a zero: the reviewer's independent code gives D = 1 at 13 of them, and traces of three secant runs show the real mechanism. K4 is a float sanity check and carries almost no evidential weight for the theorem |
| A4. `+dirty` provenance | CONFIRMED. Three rows are dirty: S+ [22/25, 9/10] (11:56:40) and L [159/20, 8] × [−1/2, −9/20], [−11/20, −1/2] (11:56:38). The edit became `ca66270b`, which touches only `throat.assemble`; the workers load their modules at launch and never call it. All three rows reproduce field for field at the clean review head, whose production code differs from `35201758` only in `assemble` and a point lookup that gives identical results for cell midpoints |
| B5. The cover (W9) | COMPLETE for every integer l ≥ 20 and every ω ∈ H⁺, including Re ω = 0, Im ω = 0 and −d_δ, the regime boundaries, ζ → 0 and Re ω → ∞. Re-checked from the stored files with exact rationals and independent constants, plus 2000 random (l, ω) lookups with l up to 10⁶: every tiling, flag and overlap holds, and no point is uncovered |
| B6. W1 (vertical-strip analyticity) | PROVEN. The continuation exists at every height on each open strip, and Cauchy–Riemann holds on unbounded strips. The c contours cross Re x = −2.29664 only on their real part (height 0) and Re x = 0 at height ≤ 0.0537 (reviewer's recomputation, matching `w1_check`). The written argument omits the simply-connected union (S3) |
| B7. W2 (Jost continuation and decay) | PROVEN, including the closing box: the near-horizon map, the Rouché disc and the Cauchy remainder are correct, and the constants check independently. Decay needs only Re ω + Im ω > 0, which holds on every S, T and A box |
| B8. W3 (half-line coercivity) | PROVEN. The boundary term, the signs on both throats (Re y₊ > 0, Re y_c < 0), the first-order term and the zero-kinetic Young weight all check; q̂ is as coded |
| B9. T, W4, W8 | PROVEN. The splitting re-derived (and checked numerically to 1e−14); the per-h argument is correct; the linear ground state is C¹ with compact support, which is all W4 needs |
| B10. A and W7 | PROVEN. ∂q/∂Re ω = 2\|ω\| sin(arg ω − φ) ≥ 0 on all of A, at every h ≤ h₀ and ν ∈ [−d_δ, 0]; checked numerically on 9 (l, ν) pairs |
| B11. L, W5, W6, `affine.py` | PROVEN. The τ-scaled recurrence, the majorant tail and the affine forms are rigorous on a line-by-line reading; 160 020 randomized enclosure tests of `affine.py` pass; the 16 sub-intervals cover [0, h₀²] |
| B12. Level and claims | WRITTEN_PROOF. s4-report §6 does not overclaim, apart from "P2-3 complete" (S4) |
| C13. Recomputation | ALL REPRODUCE (§5) |
| C14. Point checks | 150 random (l, ω), l ∈ [21, 1743], all regimes and both strip edges: W_l ≠ 0 everywhere by an independent float method that finds the l = 20 and l = 50 PS modes to 2.5e−8κ₋ and 5e−7κ₋ |
| C15. Controls | N1, N2, N4, N5 exercise the production checkers; N3 the L continuation path (as disclosed). N5's two near-axis boxes fail only through the exit rule: acceptable (§5) |
| C16. Tests | Fast: 85 passed, 17 skipped (31.6 s). P2_SLOW=1: 102 passed, 1 pre-existing warning (S2's m3 overflow), 639.6 s |

## 2. Findings by severity

### Blocking

None. Checked without finding a flaw that defeats the theorem:
- the cover and its tilings, with exact rationals and independent constants;
- W1 (existence of the vertical continuation, Cauchy–Riemann on unbounded strips, the crossings);
- W2 (the Frobenius representation on the tail half-planes, the decay, the near-horizon map);
- W3 (boundary term, signs, first-order term, Young weight, q̂);
- T (the splitting, the per-h bound, the box monotonicity, E₀ = s²Ẽ₀ at 0, the tails, c ≥ 0);
- A (W7, the lower edge);
- L (the τ-recurrence, the majorant, the affine arithmetic, the exit rule, the bulk, W5);
- the closing box (Rouché, Cauchy remainder, ρ₁, the model Ũ; S1 is a defect in one error term, not in
  the result).

### Should-fix before the merge

- **S1. The closing-box truncation term evaluates P_N on the real axis** (throat.py:964–975 and 1034).
  - `certify_closing` writes Ũ = e^{2Y}P(ρ)/P(ρ_T), with ρ = ρ_T e^{2Y(t)} **complex** along the
    contour. The model term is correct: `model()` uses `nh_P_seg(nh, acb(rho1) * e2Y)`.
  - The error term is not. `err_cell` bounds |P_N(ρ)| from below (in `e1`) and from above (in `Uabs`)
    through `PN = nh_P_seg(nh, acb(rho1 * e2t))`, the real segment {λρ₁e^{2t} : λ ∈ [0, 1]}. The tail
    term `UL` does the same at e^{2t_L}.
  - For complex ρ these are not bounds. The Taylor coefficients have mixed signs (side +: g₁, …, g₆ > 0,
    g₇, …, g₁₂ < 0; side c alternates in pairs), so neither end of the real-segment enclosure bounds
    |P_N| on the arc the contour visits. For example, at arg ρ = π (reached on the bend, where
    2η = −π), P_N ≈ g₁ − g₂|ρ| lies below the segment's lower end g₁ (g₂/g₁ = 14.9 on side +).
  - **Why it is harmless here.** Where |ρ| = ρ_T e^{2t} is not tiny (near the exit), η and so arg ρ are
    near 0; where arg ρ is large (the bend's far end and the π/4 ray), |ρ| ≤ ρ₁e^{2t} is tiny and
    P_N ≈ g₁ whatever the argument.
  - **Effect on the committed numbers: none.** The reviewer re-ran the author's closing-box check with
    only these two terms replaced by rigorous ones: the larger of the disc lower bound
    |g₁| − Σ_{n≥2}|g_n||ρ|^{n−1} and the complex-segment enclosure `nh_P_seg(nh, rho1 * e^{2Y(cell)})`,
    and the smaller of the two upper bounds. On the production partitions (161 cells per side):
    - the error term grows by at most a factor 1.0040 (+) and 1.0000046 (c);
    - no cell bound drops by more than 8.2e−5 (+) and 2.0e−6 (c);
    - mu_cells stays 0.41468629 (+) and 0.35342094 (c), and mu_tail stays 1.5341632.
  - The exact near-horizon model, with no Taylor truncation at all, gives 0.4596 (+) and 0.3970 (c)
    (§5), so the result is not in doubt.
  - Fix: use the complex enclosure (or the disc bounds) in `err_cell` and in `UL`; regenerate both
    closing-box leaves. The review's corrected run predicts identical stored values.
  - Failure scenario: a later design with a larger bend angle at the exit side, or a horizon where
    g₂/g₁ is larger, puts P_N's complex values far from its real-segment values. The error term is then
    underestimated and a cell can be certified on a false bound, with no control or test noticing:
    no control targets the closing box.
- **S2. Amendment 6 is disclosed inaccurately** (s4-preregistration.md, Amendment 6; s4-report.md §3
  and §4).
  - **"No criterion for a found zero."** The amendment says the registered *wording* had none. That is
    true of the text, but the code at launch had one, and the disclosure omits it. At `35201758` (and
    at `ca66270b`, where K4 was launched), `assemble_K4` sets
    `pass = ... and not any(r["root_in_H_plus"])`: every secant stop in H⁺ fails K4. Under the
    launch-time rule, K4 fails with 61 of 200 points. Amendment 6 replaced that rule after the first
    results.
  - **The stated mechanism is wrong.** The amendment and the report say the stop happens because "the
    phase of W turns fast (d log W/dω is large), so the step f/W′ is tiny". The reviewer traced the
    author's secant (`qnm_float.Solver.find`, dps 20) at three of the stops (`k4_secant_trace.py`):
    - d arg W/d(Re ω) at the seeds is only −3.35, −1.87 and −3.68;
    - the first secant step overshoots to Im ω = −7.49κ₋, −13.0κ₋ and −6.91κ₋, where
      \|W\| = 1.8e15, 1.1e48 and 1.8e17;
    - the next step returns to the previous iterate to 16 digits (f₂ ≫ f₁ makes x₃ ≈ x₁);
    - the step after that is at most 3e−16 relative, and the step-only test fires at the perturbed
      seed x₁ = x₀(1 + 10⁻⁶) + 10⁻⁹i, which is why every stop lies within 10⁻⁶ of its seed.
  - **The new criterion is nonetheless legitimate.** A secant stop is a zero only if W is small there;
    the 61 stops are not zeros: the residual ratio is 1.0000008–1.0000013, the rigorous zsol balls
    exclude 0, and the reviewer's independent code gives D = |y₊ − y_c|/(|y₊| + |y_c|) = 1 at 13 of
    them, including all three named in the amendment.
  - **K4's weight.** Little. W_l ≠ 0 at 200 random points is generic for an analytic function, and a
    12-step secant from random seeds would rarely find a zero even if one existed. The proof does not
    use K4.
  - Fix: append to Amendment 6 (dated) and to s4-report §3: the launch-time code's rule and K4's
    verdict under it; the actual mechanism. Keep the new criterion.
  - Failure scenario: a reader of the ledger's "criterion defined post hoc" assumes no rule had been
    fixed in advance. In fact a fixed rule existed and the check failed under it; that fact, and the
    reason the change is still harmless, belong in the record.
- **S3. W1/W2 omit the fact that makes the c-side continuation consistent** (s4-proof.md §§2–3).
  - The c contours start at a real exit in Π₋ or Π₀, cross into Π₊ through the window at Re x = 0, and
    W2's tail representation lives on {Re x ≥ X_L} ⊂ Π₊. The proof shows analyticity on each strip,
    with "each strip simply connected", and continuation "through" the windows.
  - What is needed is that the **union** D = Π₋ ∪ Π₀ ∪ Π₊ ∪ N(window at X₀) ∪ N(window at 0) is simply
    connected (its complement in the Riemann sphere is four rays joined at ∞), so that r and Φ_c are
    single-valued on D. Then Φ_c continued along the contour equals the D-function, which on the tail
    half-plane equals the Frobenius representation (∗).
  - The fact is true and two lines long, and the needed neighbourhoods of the windows exist by
    continuous dependence of the vertical continuation on X (with Cauchy–Riemann as in W1). Equivalently:
    the region between the real axis and the contour is a curvilinear trapezoid inside D.
  - Add it to W1 or W2. Failure scenario: without it, the written argument proves analyticity strip by
    strip but not that the branch reached along the contour is the one represented on the tail, which
    is exactly the non-routine step that keeps P2-1 at WRITTEN_PROOF.
- **S4. Ledger wording** (s4-report.md §7 and the verdict; FORMULATION v0.6).
  - "With P2-3B (l ≤ 19), P2-3 is complete at P_B." P2-3B counts on R′_l, not D7's R_l, and its own
    row says so. Write "every l < l₀ = 20 is counted at P_B, on R′_l (P2-3B)".
  - P2-1 implies P2-1w's statement (all three windows lie in H⁺, l ≥ 20). Say so; P2-1w can stay as a
    separate result with independent certificates.
  - Replace "Not yet independently reviewed". Use the row in §7.
  - Failure scenario: a reader takes "P2-3 complete" as D7's R_l with η = κ₋/10 and W_l the P2-2
    threshold, which is not what was counted.
- **S5. The pilot discrepancy is mis-attributed** (s4-report.md §4: "The pilot used a softplus contour
  that is not exactly real at the exit").
  - The pilot's 0.087 (s4-preregistration §9) is the **leading-order** margin min_t q₀/|Y′|², at ε = 0
    and normalised by |Y′|². The certificates bound q̂, which includes the O(ε) terms with ε up to
    h₀/ζ ≈ 0.26–0.29 near the top.
  - For the stored compact-family designs on the top leaves (`pilot_compare.py`): leading order
    0.025–0.032; full q̂ (float vertex minimum) 0.015–0.018; certified 0.009–0.014.
  - So the drop from 0.087 has three parts: the contour family (and the ζ-box), 0.087 → 0.025–0.032 at
    leading order; the O(ε) terms at the l = 20 corner of each leaf, → 0.015–0.018; and the
    certification loss of the mean-value cells, → 0.009–0.014. The report names only the first.
    Correct the sentence.
  - Failure scenario: a planner for P_A, where the margins are thinner, budgets from the leading-order
    pilot number and is surprised by the ε terms.

### Wording and nits

- **W1.** Amendment 4: "a* = 0.02996861191… (a* − a_d ≈ 3.7·10⁻¹¹)". The values are
  a* = 2Ω_Pγ₀ = 0.0299686149113 and a_d = 0.0299686118960, so a* − a_d = 3.0e−9. The best 30-bit
  dyadic below a* is 32178555/2³⁰. Harmless: any a_d < a* is admissible.
- **W2.** s4-proof §0: "The true minimum is −0.0079101, near r = 3.94". That is the certified lower
  bound over 6000 cells. The float minimum is −0.0078910 at r = 3.9461. W_bound = 0.00792 holds either
  way (margin 2.9e−5 against the true minimum).
- **W3.** s4-proof §1, "L's leaves have total area 8 × 0.55": since `279878d4` the cover checks the
  tiling box by box. Update the sentence.
- **W4.** Printed bounds. `throat.lohi` rounds to nearest at 8 digits, so a printed lower bound can sit
  above the arb endpoint (e.g. the c leaf [43/50, 22/25] prints mu_tail 1.3962179; the true tail
  value is 1.39621789). `L_PB_plus.json` "min_Re_y_plus_lower" is the minimum of *midpoints*. Every
  pass/fail decision uses the arbs, so this is presentation only; round outward.
- **W5.** "K4" names both the float acceptance check (§4) and a kill criterion (§7) of the
  pre-registration.
- **W6.** N2 is registered on h ∈ [1/50.5, 1/20.5]; the T checker has no lower h bound, so it runs on
  (0, h₀]. That is a superset (a certification would be stronger), so the control is valid; say so.
- **W7.** s4-report §1: T's "certified margin 0.1003·a*" (box [0.8, 0.85]) is set by the bisection
  target a*/10, not by the design: cells are refined only until they reach it. The float minimum of
  the same box bound is 0.139a*, and the certified values stay well below the float ones up to
  ζ ≈ 1.45Ω_P. Correct as a certified statement; say that it is target-limited.
- **W8.** The first K4R launch's logs (`results/part1/K4R_w*_20261018T124005.log`) are absent from the
  author's working directory although its launch record names them. Gitignored, so nothing in the
  repository refers to them.

## 3. Pre-registration integrity (scope A)

### A1. The pre-registration against execution

| registered (`d00f9a3c` + Amendments) | executed (`35201758`, certificates) | verdict |
|---|---|---|
| Target P2-1 at P_B, l₀ = 20, h₀ = 1/20.5; Im ω > 0 by P2-0; Re ω < 0 by symmetry | as registered | ✓ |
| Cover: Ω_L = 8κ₋, Ω_S = 6κ₋; ζ_SU = 0.9, ζ_TL = 0.8, ζ_TU = 2.0, ζ_A = 1.5 (Ω_P) | `part1.py:39–42`, `throat.py:44–46`; `cover_PB.json` | ✓ |
| T: φ_∞ = −1.4375, s₁ = 7, p = 4, η = ½, ε_g = 0.05, taper [2, 4], depth ≤ 12; linear Γ (A4); target a*/10 (A5) | `DESIGN_T` | ✓ |
| A: φ₀ = −0.2, Γ = 0, monotone above ζ_A | `DESIGN_A`, `monotonicity_A` | ✓ |
| S: π/4 tail; grid w ∈ {0.15, 0.3, 0.6, 1} ∪ {2, 3, 4} (A5), y_a, φ_∞, p_c, p_w as registered; ζ-boxes 0.02Ω_P halved to 0.0025Ω_P; ε halved ≤ 4 levels; compact C³ family, Young weight, cells, closing box (A5) | `throat.py:43–61` and the code | ✓ (w = 4 never selected) |
| L: boxes κ₋/20 halved to κ₋/80; h² ∈ [0, h₀²] in 16 pieces; affine forms; τ grid 1/64 (A5) | `lowfreq.py:41–54` | ✓ (no box needed halving) |
| Controls N1–N5, K1–K4, replay, tests | `controls_PB.json`, `replay_PB.json`, `test_part1.py` | ✓; K4's criterion changed post hoc (S2) |
| Kill criteria | none fired; ≈ 3.9 core-h of 48 | ✓ |

**The original text is intact.** `git diff d00f9a3c 0e329217` on the pre-registration removes one line
only, "(none)" under Amendments. Amendments 1–6 are appended, each dated.

**Development before Amendment 5.** The checkers were developed between 10:52 and 11:53 on sample
boxes (≈ 0.5 core-h, disclosed but not logged per run), and the methods were fixed in A5 23 s before
production. That is acceptable here: the certificates are proofs, any design that certifies is valid,
and nothing whose change could favour a verdict (l₀, the regime boundaries, the controls, the kill
criteria) moved.

### A2. Timing

| event | time (2026-10-18) |
|---|---|
| A4 (T ground state) | 10:52:04 |
| A5 (S, L methods) | 11:52:57 |
| production code `35201758` | 11:53:14 |
| launches: L (2 workers), S+, Sc / T / A | 11:53:20 / 11:53:21 / 11:54:40 / 11:54:50 |
| rows finished (author's JSONL): Sc / S+ / L / T / A | 11:53:21–11:53:38 / 11:53:21–11:56:40 / 11:53:21–11:56:38 / 11:54:41–44 / 11:54:51 |
| `ca66270b`; K4 launch | 11:56:52; 11:57:02 |
| certificates `776aa7ef`; A6 | 11:57:38; 12:01:46 |

Every launch record has commit `35201758` without `+dirty` (the flag is computed from
`git status --porcelain -- .` in `src/p2`), so the tree was clean when the workers loaded their code.

### A3. Amendment 6

See S2 for the disclosure. On the substance:
- **Definition or weakening?** A definition. A secant iteration that stops on its step test has not
  found a zero unless the function is small there. The amendment requires both a small residual and a
  rigorous zsol ball containing 0 to call a stop a zero; that is the natural reading of "finds a QNM".
- **Was the run changed?** No: the points, the solver and its settings are as launched; only the
  evaluation of the stored stops changed.
- **Independent check.** At 13 of the 61 stops (every sixth, plus the three named in the amendment),
  the reviewer's float code (`indep_points.py`, below) gives D = 1.000 to 15 digits. A zero of W_l has
  D = 0, and the same code finds the l = 20 PS mode with D = 2.6e−6.
- **Why the secant stops.** S2 above: overshoot to a point with enormous \|W\|, return, zero step.
- **Weight.** K4 is a float consistency check that is not part of the proof; its pass or fail says very
  little about the theorem.

### A4. The `+dirty` rows

- Which rows: the author's working files hold 45 S+ rows (44 clean, 1 dirty), 45 Sc (clean), 1760 L
  (1758 clean, 2 dirty), 24 T and 1 A (clean). No duplicates, no failures, no re-runs.
- The edit: `ca66270b` is the only record of it. It changes `throat.assemble` (a W2 tail-strip check
  and a tail-disc maximum), which the workers never call; nothing else in `src/p2` changed between the
  launches and that commit.
- From `35201758` to the review head, `lowfreq.py`, `part1_T.py` and `affine.py` are unchanged.
  `throat.py` changes in `assemble` and in `CurveContour.find` for a == b (a point at a step node),
  which returns the same step as before for every point strictly inside a step, as all cell midpoints
  are.
- **Recomputed at the clean head:** S+ [22/25, 9/10] and both L boxes reproduce every stored field
  (design, X_a, k, exit, w1_check, eps-box bounds, tail disc, cell counts; τ_a, Re y bounds, terms, M).
  So do all other 89 S leaves and 1758 L boxes (§5).

## 4. Mathematics (scope B)

### B5. The cover (W9)

Fix an integer l ≥ 20, so h = 1/(l+½) ∈ (0, 2/41], and ω = Re ω + iν with Re ω ≥ 0, ν ∈ [−d_δ, 0],
d_δ = 0.55κ₋ = 11/20 κ₋.
- **Re ω ≤ 8κ₋**: L. Its 1760 boxes tile [0, 8] × [−11/20, 0] (κ₋ units) exactly, and each holds
  for all h² ∈ [0, h₀²] ∋ h². This includes Re ω = 0 and both strip edges.
- **6κ₋ ≤ Re ω ≤ 0.9Ω_P/h**: S. Then ζ = h Re ω ∈ (0, 0.9Ω_P] (ζ ≥ 6κ₋h > 0) and
  ε = 1/Re ω ≤ min(1/(6κ₋), h₀/ζ). Each side's leaves tile [0, 0.9] (Ω_P units) and, in each leaf,
  the ε-boxes tile [0, ε_max] with ε_max ≥ min(1/(6κ₋), h₀/ζ_lo). The closing leaf (0, 0.02] takes
  ζ → 0, at any Re ω ≥ 6κ₋ (Amendment 3). The + and c leaves are different partitions; each side is
  used with the leaf containing ζ, whose ε_max bounds ε.
- **0.8Ω_P/h ≤ Re ω ≤ 2Ω_P/h**: T (24 boxes tiling [0.8, 2.0]).
- **Re ω ≥ 1.5Ω_P/h**: A, which reaches Re ω → ∞ by W7.
- **Overlaps.** 6 < 8; 0.8 < 0.9; 1.5 < 2; and 0.9Ω_P/h₀ = 90.75κ₋ > 8κ₋, so S is non-empty and
  reaches L for every h ≤ h₀.
- **Decay.** Re ω + Im ω > 0 on S (6κ₋ > d_δ), T (0.8Ω_P/h₀ = 80.7κ₋) and A.
- **Reductions.** E₊ ∪ E_c lies at Im ω ≤ −κ₊ = −0.8229κ₋ < −d_δ. Re ω < 0 follows from
  W_l(−ω̄) = conj W_l(ω); Im ω > 0 from P2-0.

`cover_recheck.py` re-derives all of this from the stored files with exact rational endpoints and
the reviewer's own constants:
- T's 24 boxes tile [4/5, 2] exactly, all with h₀ = 2/41, ν_low = −d_δ and the decay flag; A sits at
  3/2 with φ₀ = −1/5 < −arctan(d_δh₀/ζ_A) = −0.00364.
- On each side the S leaves (47 and 45) tile [0, 9/10]; in every leaf the ε-boxes tile [0, ε_max] and
  ε_max ≥ min(1/(6κ₋), h₀/ζ_lo); every exit, order, η/φ-at-exit and W1 flag is set; the + contours
  and tails have Re x < X₀, the c tails X_L > 0, and the c crossings sit at heights ≤ 0 (X₀) and
  ≤ 0.0537 (Re x = 0), at most 0.211 from the exit in t.
- The 1760 L boxes are distinct, on the κ₋/20 grid, inside [0, 8] × [−11/20, 0] with total area
  8 × 11/20, the same on both sides, and all pass.
- All eight overlap and pole inequalities hold.
- 2000 random (l, ω), l log-uniform up to 10⁶ and Re ω up to 50Ω_P/h (20% of them below 9κ₋), each
  found a covering certificate piece by exact lookup (L 405, S 520, T 136, A 1102 hits; none
  uncovered).

### B6. W1: vertical-strip analyticity

- **Existence at every height.** Along the vertical continuation from a real X,
  d/dY Re r*(r) = Re(i) = 0, so Re r*(r(X + iY)) = X. Re r* is continuous at r = 0 with value
  X₀ = Re r*(0) = −2.29664 and tends to Σ A_i log\|r\| → 0 as r → ∞ (Σ A_i = 0; reviewer: 5e−49). So
  for X ∉ {X₀, 0} the solution of dr/dY = iF(r) cannot enter a neighbourhood of 0 or of ∞. It cannot
  reach a root of Δ either (equilibria, and Re r* → ±∞ there). It therefore stays in a compact subset
  of ℂ∖{0}, where F is analytic, and exists for all Y.
- **Joint analyticity on unbounded strips.** The S2 argument is local: for each (X, Y), continuous
  dependence gives a C¹ solution near (X, Y), G = ∂_X R − F(R) solves ∂_Y G = iF′(R)G with G(X, 0) = 0,
  and Cauchy–Riemann follows. Unboundedness plays no role.
- **The windows.** S2's validated segments are at the exact lines: `real_point` receives X₀ as a ball
  (radius 4.9e−12) and encloses r for every X in it. The segments run to heights 0.14133 and 2.43797,
  with \|r\| ∈ [1.961, 2.004] and [2.412, 2.491].
- **The crossings actually used** (`indep_S.py`, recomputed from the stored designs): on the c side
  the contours meet Re x = X₀ (when X_a < X₀) on their real segment (height 0), and Re x = 0 at height
  at most 0.0537 (leaf [22/25, 9/10]), at t-distance at most 0.22 from the exit. These match the
  stored `w1_check` fields. All + contours have Re x ≤ X_a < X₀; all tail regions satisfy X_L < X₀
  (+) or X_L > 0 (c).
- **Single-valuedness.** Along one contour the continuation is trivially single-valued. What must be
  shown is that it agrees with the tail representation; that needs the simply-connected union (S3).
  Missing from the text, true, and short.

### B7. W2: the Jost solutions on the contours

- **Tail representation.** On {Re x ≤ X_L} ⊂ Π₋ (side +), the first-exit argument with
  Re r*(r) = Re x keeps r in the disc \|r − r₊\| ≤ d (stored d ≤ 6.2e−5; c side ≤ 3.6e−7), starting
  from the real points, which lie in the disc by the same argument. There
  Φ₊ = c·r·g₊(r)·exp(−iω(x − c₊ − R̃₊(r))) is analytic in x (radius of g₊: 0.097 ≫ d) and agrees with
  Φ₊ on a real half-line, hence on the half-plane.
- **Decay.** On the + ray x = X_a + ((1+i)t − t_a)/k, \|e^{−iωx}\| = const·e^{(Re ω + Im ω)t/k} → 0 as
  t → −∞, and u′ = x′Φ_x has the same exponential. On the c ray, Φ_c ~ e^{iωx} and
  x = X_a − ((1+i)t − t_a)/k give the same rate. Re ω + Im ω ≥ 5.45κ₋ on S, so the rate is uniform on
  every box, including ε → 0 (faster) and ζ → 0 (k = κ_h, same rate).
- **The closing box.** With s = ±(r − r_h), ρ(s) = s·f(s) = e^{±2κ_h(x − c_h)} and f analytic and
  non-zero on \|s\| ≤ d_r (d_r = 0.06 < r₊ − r₋ = 0.097 and d_r = 1 < r_c − r₊ = 4.77):
  - Rouché on \|s\| = d_r with \|sf\| ≥ m_f·d_r = R gives a unique, simple, analytic inverse ψ on
    \|ρ\| < R with values in the disc; ψ is real on [0, R).
  - On the half-plane H = {\|e^{±2κ_h(x − c_h)}\| < R}, r̃(x) = r_h ± ψ(e^{±2κ_h(x − c_h)}) is analytic,
    equals r on the real part of H, hence solves dr/dx = F(r) on H and is the continuation of r there.
    H lies in Π₋ (side +: H = {Re x < −21.45}) or Π₊ (side c: H = {Re x > 16.97}). The stored check,
    r*(r_h ± d_r) < X₀ (+) or > 0 (c), implies this, since ρ(r_h ± d_r) ≥ m_f d_r = R on the real
    axis.
  - G(ρ) = U(r_h ± ψ(ρ)) is analytic on \|ρ\| < R with \|G\| ≤ M_G, so \|g_n\| ≤ M_G/R^n and the stated
    Cauchy remainder holds. W̃(r_h) = 0 gives \|W̃\| ≤ M_W\|ρ\|/R (the code's bound is weaker, valid).
  - G is increasing on [0, ρ₁] because U is monotone on (r₊, r_P) (resp. (r_P, r_c)) and ρ(r) is
    increasing; with G(ρ₁) ≥ ζ₁² this gives ρ_T ≤ ρ₁.
  - The contour stays in H: \|ρ\| = ρ_T e^{2t} ≤ ρ₁e^{2t_a} = 0.477R (+) and 0.579R (c).
  - Independent check (`nh_check.py`): m_f, R, M_G, M_W are valid (sampled extremes 0.6399 ≥ 0.6323,
    \|sf\| ≥ 0.0406 ≥ 0.0379, max\|U\| 0.00687 ≤ 0.00737, max\|W̃\| 0.00123 ≤ 0.00148 on the + side;
    0.866 ≥ 0.337, 0.874 ≥ 0.337, 0.00931 ≤ 0.0121, 0.0123 ≤ 0.0199 on the c side); g₁ and g₂ agree to
    all printed digits with an FFT Cauchy integral; G(ρ₁)/ζ₁² = 1.0046 (+) and 1.00026 (c).
  - The defect is the error term only (S1).

### B8. W3: half-line coercivity

- u = Φ∘x satisfies −u″ + (x″/x′)u′ + Qu = 0 with Q = x′²(V − ω²). Multiplying by e^{−iφ}ū,
  integrating over (−∞, t_a] (the boundary term at −∞ vanishes by W2) and taking real parts gives the
  identity in s4-proof §4 exactly.
- Young with weight cos φ/κ_tot: \|κ\|\|u′\|\|u\| ≤ cos φ\|u′\|² + \|κ\|²\|u\|²/(4 cos φ). The kinetic term
  ∫cos φ\|u′\|² is cancelled exactly, which is allowed (no ground state is used), and
  q = Re(e^{−iφ}Q) − κ_tot²/(4 cos φ). cos φ ≥ ½ (φ ∈ [−π/3, 0]); the code also certifies cos φ > 0
  per cell.
- At t_a: φ = 0, x′ = s/k, so Re(u′ū) = (s/k)\|Φ\|² Re y. Positivity gives Φ(x_a) ≠ 0 and
  Re y₊ > 0 (s = +1), Re y_c < 0 (s = −1).
- Normalised form: x′² = Y′²/k², x″/x′ = Y″/Y′, V − ω² = (Re ω)²[U/ζ² − 1 + 2iε\|ν\| + ε²(ν² + W̃)],
  so q = (Re ω/k)² q̂ with q̂ as coded (`throat.parts`, `Params.combine`). q̂ is affine in each of
  1/ζ², ε\|ν\|, ε²ν², ε² separately, so the interval combination is a valid bound over the box.
- The contour: m₀ (coefficients M0–M3) is the stated C³ smoothed max(0, u): m₀(±1), m₀′(±1),
  m₀″(±1) = 0, m₀‴(±1) = 0 all check, so η is C³ and the exit is real for t ≥ w.
- Bulk: each side's exit gives Ũ − 1 > ε_max²W_bound ≥ ε²W_bound for its own leaf; by unimodality of
  U, Re(V − ω²) = (Re ω)²[Ũ − 1 + ε²(ν² + W̃)] > 0 between the exits; X_a(+) < x_P < X_a(c).

### B9. T, W4 and W8

- **Splitting.** With V_l = h⁻²U + W̃ (L = h⁻² − ¼; sympy) and ω = ζ/h + iν,
  Re(ie^{−iφ}(V − ω²)) = h⁻²Re(ie^{−iφ}(U − ζ²)) + h⁻¹·2ζν cos φ + Re(ie^{−iφ}(W̃ + ν²)), and with
  G′ = Γ/h the ground-state terms give −cΓ²h⁻² + (cΓ)′h⁻¹. This is the stated E₀, E₁, E₂. The
  reviewer's float code evaluates S2's E from V_l directly and matches E₀ + hE₁ + h²E₂ to 1.0e−14.
- **Per-h bound.** E = h⁻¹[E₀/h + E₁ + hE₂] ≥ h⁻¹[E₀/h₀ + E₁ + h₀ min(0, E₂)] when E₀ ≥ 0 and h ≤ h₀.
- **Box monotonicity.** Re(ie^{−iφ}) = sin φ ≤ 0 (φ ∈ [−1.4375, 0]) and cos φ > 0, so E₀ is
  non-decreasing in ζ², E₁ non-increasing in ζ and non-decreasing in ν, and E₂ non-increasing in ν².
- **At s = 0.** E₀(0) = 0 exactly; (U − U_P)/s² ∈ ½ conv U_ss([0, s]) since U_s(0) = 0, and
  sin φ/s² = sinc(φ)·φ_∞s²/(s₁⁴ + s⁴) (φ₀ = 0 in T). The code uses only E₀ ≥ 0 on those cells.
- **Ground state.** Γ = 0.95·a_d·s·θ(s), θ a C¹ cosine taper: Γ is C¹ with compact support, so G′ is
  absolutely continuous and the boundary term in W4 vanishes outside \|s\| ≤ 4. Cells never straddle
  \|s\| = 2, 4 or 0, so E₁′ (which needs Γ″) exists on each cell. c ≥ 0.1235 on every cell.
- **Tails.** \|s\| ≥ 64, 32, Γ = 0: E₀ ≥ ζ_lo²\|sin φ(S)\| − sup\|U\|, E₁ ≥ 2ζ_hiν_lo cos φ(S),
  E₂ ≥ −sup\|W̃\| − ν² − \|φ′(S)\|/(4η), all valid with S2's discs.

### B10. A and W7

- With Γ = 0, E = q, and ∂q/∂Re ω = −Re(ie^{−iφ}·2ω) = 2\|ω\| sin(arg ω − φ).
- On A, arg ω ∈ [−arctan(d_δh/ζ), 0] ⊂ [−0.00364, 0] (since Re ω ≥ ζ_A/h ≥ ζ_A/h₀), and
  φ ∈ [φ_∞, φ₀] = [−1.4375, −0.2]. So arg ω − φ ∈ (0, 1.44) ⊂ (0, π) and the derivative is ≥ 0 at
  every s (tails included), every h ≤ h₀ and every ν ∈ [−d_δ, 0].
- Hence E(Re ω) ≥ E(ζ_A/h) ≥ μ_A/h for every Re ω ≥ ζ_A/h; c is ω-independent.
- At the lower edge the per-h bound needs E₀ ≥ 0, certified (0.01087).
- The pre-registered form "E ≥ h⁻²(m_A − h·2ζ_Ad_δ − h²max\|E₂\|)" was implemented as T's per-cell
  bound, which is the same inequality with cell-wise constants.

### B11. L, W5 and W6

- **The τ-equation.** With r = r_h + sh²τ, Δ = h²D(τ) and Lh² = 1 − h²/4, the radial equation becomes
  D²R_ττ + sDE R_τ + (ω²r⁴ − (1 − h²/4)D)R = 0, polynomial in τ and h². Re-derived; `tau_polys` builds
  exactly these coefficients.
- **The recurrence.** A₂ = B₁ = d₁², C₀ = ω²r_h⁴, so the indicial factor is d₁²n(n + 2ρ) with
  ρ = −iω/(2κ_h); n + 2ρ ≠ 0 since Im ω > −κ_h. The lags (≤ 6) and the majorant g_j are those of
  arbw, each factor decreasing in n for n > 2\|ρ\|.
- **Identification (W6).** For each actual h the τ-series is the Frobenius series of the Jost solution
  in a rescaled variable; R_τ/R is independent of the factor (h²)^ρτ^ρ.
- **Radius.** (r₊ − r₋)/h₀² = 40.8 against τ_a ≤ 8, and (r_c − r₊)/h₀² = 2006 against
  τ_a^c ≤ 86.4; the verified q·τ < 1 enforces convergence anyway.
- **y.** y = h²D/r³ + s(D/r²)R_τ/R, as coded, with s = −1 at r_c.
- **Exit rule and bulk.** h⁻²U = D/r⁴ ≥ M ≥ max(4\|ω\|², C_W) at both exits and, by unimodality, in
  between; W̃ ≥ −W_bound ≥ −C_W/4. So Re(V − ω²) ≥ M − M/4 − M/4 + ν² > 0.
- **W5.** (Re y)′ = Re(V − ω²) − (Re y)² + (Im y)² > 0 wherever Re y = 0, so Re y₊ > 0 propagates
  forward and Re y_c < 0 backward; a zero of Φ₊ would force Re y₊ → −∞ first. Then y₊ ≠ y_c between
  the exits and W_l ≠ 0.
- **Bulk facts** (`bulk_and_affine.py`): Δ′r − 4Δ = −2(r² − 3r + 2Q²) and the splitting in sympy;
  3 − r_P = 1.04127 < r₊; U′ > 0 / < 0 on 2000-point grids on each side; min W̃ = −0.0078910.
- **`affine.py`.** The semantics is pointwise: for the true (ω, h²), with dw = ω − wc and
  dh = h² − m (m the exact midpoint, inside the ball hc), each quantity equals c₀* + c₁*dw + c₂*dh + E*
  with c_i* in the balls and \|E*\| ≤ e. Sums, products (remainder \|c₀\|e′ + \|c₀′\|e + dev·dev′) and
  the inverse (1/c₀ − L/c₀² + E′ with \|E′\| ≤ e/\|c₀\|² + s²/(\|c₀\|²(\|c₀\| − s)), \|c₀\| ≥ a > s) preserve
  it. The randomized test: 300 random programs (sums, products, quotients, inverses, squares, scalar
  multiples) over random affine inputs, at 60 random admissible points each, in mpmath at 50 digits:
  160 020 checks, no violation, largest distance/dev 0.9992.
- **h² pieces.** [0, h₀²] = ∪_{j<16}[jh₀²/16, (j+1)h₀²/16]; each certificate is the conjunction.

### B12. The level and the claims

- **WRITTEN_PROOF** (README "Ledger levels"). The inequalities are machine-checked in Arb, but W1 and
  W2 are non-routine prose (the M3 continuation analysis, exactly as in P2-1w), and the Im ω > 0 part
  rests on P2-0's written argument. CERTIFIED would need W1 and W2 machine-checked.
- **P2-1 and P2-1w.** P2-1 contains P2-1w's statement (its windows lie in H⁺ at l = 20, 50, 100).
  P2-1w's certificates are an independent method on the same contour; keeping the row costs nothing.
- **"P2-3 complete at P_B with P2-3B."** True in P2-3B's form: every l < l₀ = 20 is counted, on R′_l.
  The statement needs that qualifier (S4).
- **s4-report §6** disclaims the review, P_A, the assembly and CERTIFIED. Apart from S4's wording,
  nothing is overclaimed. With P2-0, P2-1 (l ≥ 20), P2-2 and P2-3B (l ≤ 19), every ingredient of P2-4
  at P_B now exists; the assembly is S5's and this review does not assert it.

## 5. Recomputation and tests (scope C)

### C13. Independent recomputation

All scripts are in the reviewer's scratch directory `…/scratchpad/rv/` (see the note at the top),
each with its JSON output of the same name; run from that directory with the repository's
`.venv/bin/python <script> <out.json>` (`count_l20.py` takes an output directory and a worker count;
`replay_all_SL.py` takes `S` or `L` first). "Independent of" names the code a check does not use.

| check | independent of | result |
|---|---|---|
| `indep_params.py`: P_B constants from the exact definitions (mpmath, 50 digits) | all `p2` | κ₋, κ₊, κ_c, γ₀/κ₋ agree with P2-par; X₀ = −2.2966368932, x_P = −2.4379660346, a* = 0.0299686149 |
| `indep_TA.py`: r(s) on S2's contour by DOP853 (rtol 1e−13), s-spacing 1/256 on [−64, 32], tails to −160/80; F per T box from the corners; S2's E from V_l directly at l ∈ {20, 21, 25, 37, 50, 100, 400, 2000, 10⁵} | all `p2` | All 24 T boxes: certified ≤ float minimum (e.g. 0.1003 vs 0.1394 a* on [0.8, 0.85]; 0.2964 vs 0.2965 on [1.95, 2]); E₀ ≥ 0 on the grid; c ≥ 0.1329 on the tails. 54 direct samples: hE ≥ μ_cert(box) everywhere, smallest ratio 1.0136. Splitting matches to 1.0e−14. A: certified 6.9553 a* vs float 6.9568; E₀ ≥ 0.010874; q non-decreasing along 41 Re ω values for 9 (l, ν) pairs |
| `indep_S.py`: every one of the 92 S leaves; contour rebuilt from the stored (X_a, k, design), s = ±(r − r_h) integrated as log s (relative accuracy), t-spacing 1/512, tail to t = −40; q̂ at the 16 vertices of each leaf's parameter box | all `p2` | Every certified mu_cells and mu_tail lies below the float minimum (tolerance: the 8-digit printing, W4). Thin leaves: + [0.82, 0.84] 0.01333 vs 0.01650; + [0.86, 0.88] 0.01303 vs 0.01591; + [0.88, 0.9] 0.01424 vs 0.01718; c [0.82, 0.84] 0.00917 vs 0.01477; c [0.86, 0.88] 0.01177 vs 0.01782; c [0.88, 0.9] 0.00984 vs 0.01670. Exit margins and order as stored |
| same, closing boxes: the actual contour x_T(ζ) ± Y/κ_h for ρ_T = λρ₁, λ ∈ {1, 0.7, …, 10⁻⁵}, with U~ = U(r)/U(r_T) exact; and λ → 0 (U~ = e^{2Y}, W̃ = 0) | all `p2`, and the Taylor model | float minima 0.4596 (+) and 0.3970 (c) against certified 0.4147 and 0.3534; tails ≥ 1.534163 (certified 1.5341632); exit margins ≥ 59.2 |
| `nh_check.py`: m_f, the Rouché bound, M_G, M_W by dense sampling; g_n by an FFT Cauchy integral with ψ solved by Newton; G(ρ₁) | all `p2` | §4 B7: all stored constants valid, g₁, g₂ equal |
| `closing_err_check.py`: the author's closing-box check with a rigorous error term (S1) | — | margins unchanged; largest cell drop 8.2e−5 |
| `indep_L.py`: Jost y at the stored exits by a Frobenius series in r itself (mpmath, 30 digits), and the exit rule on 201 h² values per box. (a) 52 boxes (near Re ω = 0, 8κ₋, Im ω = −d_δ, and random), 6 points per box, l ∈ {20, 21, 33, 100, 1000, 10⁵}; (b) `indep_L.py out.json all`: all 1760 boxes, 5 points per box (corners and centre), l ∈ {20, 21, 100, 10⁵} | all `p2` | (a) 1872 values, (b) 35 200 values: none violates the stored bound; smallest margins 9.8e−4 (Re y₊) and 5.6e−5 (Re y_c); the exit rule holds on every box |
| `bulk_and_affine.py` | all `p2` (bulk); `affine.AF` tested | §4 B11 |
| `replay_author.py`, `replay_all_SL.py`: the author's code at the clean head | — (reproducibility) | All 24 T boxes and A recomputed in full from their stored partitions; all 92 S leaves (205 CPU-s) and 1760 L boxes (335 CPU-s) recomputed; every produced field equals the certificate |
| `count_l20.py`: S3's reviewed counting pipeline (`p2.part3`, zsol box bounds, Cauchy discs) at l = 20 on R′_20 = [−a_20, a_20] × [−b, e], a_20 from the CERTIFIED P2-2 certificate `lg_PB_l20.json`; working files redirected to the scratch directory | all S4 code and arguments | Winding 0 (sum over H: 0 ± 2.3e−3; /π: 0 ± 7.4e−4), certified; 223 segments, K ≤ 104, \|W₂₀\| ≥ 4.07 on H, 2008 CPU-s. a₂₀ = 76112707/2²⁴ (106.64κ₋), b = 6280879/2²⁸, e = 4567911/2³⁰. part3, zsol, params and the lg certificates are unchanged on the S4 branch. With P2-2 and P2-0: no zero of W₂₀ in S_δ |
| `cover_recheck.py` | all `p2` | §4 B5: every check holds; no random point uncovered |
| `k4_secant_trace.py`: the author's float secant (`qnm_float`, dps 20, 12 steps) re-run step by step at K4 stops 0, 11, 23, with the phase rate of W at each seed | — (diagnosis) | S2: overshoot to Im ω ≈ −7…−13κ₋, return, zero step; phase rates −3.35, −1.87, −3.68 |
| `k4_stops_indep.py`: D at 13 of the 61 stops with `indep_points.py`'s code | all `p2` | D = 1.000 at all 13 |
| `pilot_compare.py`: leading-order and full q̂ for the stored top-leaf designs | all `p2` | S5 |
| `indep_points.py` | all `p2` | §5 C14 |

### C14. Point checks

`indep_points.py` (no `p2` code): Riccati in the tortoise coordinate, dr/dx = F(r),
dy/dx = V_l − ω² − y², integrated (DOP853, rtol 1e−11) from r_h ± h²/4, with start values from the
mpmath Frobenius series of `indep_L.py`, to x_P from both sides. W_l = 0 iff y₊ = y_c there; the
statistic is D = \|y₊ − y_c\|/(\|y₊\| + \|y_c\|).
- **Positive control.** At the l = 20 and l = 50 PS fundamental modes (below the strip), D = 2.6e−6 and
  1.1e−4, and a secant on y₊ − y_c seeded 0.36κ₋ away converges to the modes within 2.5e−8κ₋ and
  5.0e−7κ₋.
- **150 random points**, l log-uniform in [20, 2000] (drawn: 21–1743), regimes L 36, S 31, T 21,
  A 32 and 30 far above the barrier (3–10 Ω_P/h), 15% on Im ω = −d_δ and 10% on the real axis:
  D ≥ 0.99999999981 at every point.
- This is a gross check only: D saturates at 1 away from a zero, so random points cannot find an
  isolated zero. The certified inequalities themselves were checked pointwise in C13 (T: hE at actual
  l; S: q̂ at the box vertices, which bound every actual (ζ, ε, ν); L: Re y at actual l).

### C15. Controls

| control | path | assessment |
|---|---|---|
| N1 (IN ±γ₀/50, STRAD) at l = 20, 50, 100 | `certify_box_abs` → `run_box`, the production T cell checker at h₀ = 1/(l+½) | good; STRAD's −0.100a* is exactly E₁(0) at ν = −1.05γ₀ (2Ω_Pν + 0.95a*), so the checker's border sits 5% above the mode layer |
| N2 | same checker on (0, h₀] (W6) | good |
| N3 (NE modes, l = 1, 2) | the L continuation path, not production (disclosed) | good as a positive check of the y-sign logic: Re y₊ ∈ [−0.232, −0.222], [−0.316, −0.293] at the modes; mode-free boxes at half depth certify |
| N4 | `certify_zbox`, the production S checker, exit moved into the allowed region; ε-levels 0, depth 2 | good; q̂ < 0 on an interval near the exit (−0.81, −0.92) and Ũ(x_a) − 1 < 0 |
| N5 | `LChecker.box`, the production L checker with the exit at D/r⁴ ≥ (Re ω)²/4 | acceptable. In the two boxes with ν ∈ [−0.05, 0], Re y keeps its sign at an allowed-region exit (Re y ≈ ν + O(U) > 0), so only the exit rule rejects. That is correct physics: the sign test alone means nothing in the allowed region, and the argument needs the forbidden bulk, which N5 shows is enforced |

Not exercised by any control: the S closing box (S1 lives there) and A's monotonicity. This review's
checks cover both (exact closing model; W7 numerically).

### C16. Tests

All runs at nice 10 in the review worktree at `0e329217`; `git status` stayed clean.

| command | result |
|---|---|
| `pytest problems/P2/tests -q` | 85 passed, 17 skipped, 31.6 s |
| `P2_SLOW=1 pytest problems/P2/tests -q` | 102 passed, 1 warning (the pre-existing S2 m3 Riccati overflow), 639.6 s |

`test_part1.py` ran 8 tests fast (2 slow ones skipped) and all 10 with P2_SLOW=1, as the report
states. It covers the splitting against S2's E, one box per regime, both closing boxes, the affine
forms, K2 at one point, and N1 (l = 20), N3, N5 when slow. It does not test the closing-box error term
against complex ρ (S1).

## 6. Compute and tool notes

| item | CPU |
|---|---|
| pytest fast and P2_SLOW=1 | 0.19 core-h (669 s) |
| reviewer-only float/mpmath checks (T/A, all S leaves, closing model, near-horizon constants, L on all boxes, 150 points, K4 stops, bulk, affine, cover, pilot) | 0.26 core-h (≈ 920 s; the 150 points take 592 s, all L boxes 266 s) |
| replays with the author's code (all T/A, all S, all L, dirty rows, closing boxes, corrected closing error term) | 0.16 core-h (581 s) |
| K4 secant traces (author's float solver) | 0.05 core-h (194 s) |
| rigorous count at l = 20 (S3 pipeline, 3 workers) | 0.56 core-h (2008 s) |
| **total** | **≈ 1.22 core-h** (cap 4) |

All runs at nice 10. At most 3 worker processes ran at once, except for one 2-second diagnostic
(printing the near-horizon coefficients) that overlapped the 3-worker count.

**Denied tool calls.**
1. `mkdir -p …/p2-s4-review/problems/P2/certificates/review_s4 && cd … && git status --short | head`:
   denied (auto-mode classifier, "Interfere With Workloads"). Not retried, not worked around; the
   review's scripts and outputs therefore live outside the repository (top of this document).
2. Stopping a background shell of the reviewer's own (an `ls` that hung; `ls` is `eza` here):
   denied for the same reason. Two such `eza` processes of the reviewer's remain, sleeping at 0% CPU.

## 7. The ledger

**Level of P2-1: WRITTEN_PROOF.** The decisive inequalities are machine-checked in Arb, every
certificate replays exactly, and reviewer code sharing nothing with `p2` confirms every certified
bound. W1 and W2 are
non-routine written arguments (with S3's addition), and the Im ω > 0 part is at P2-0's level. The row
below assumes S1 is fixed and the two closing-box leaves regenerated; it quotes no closing-box number,
so it needs no change if their values move in the last digit.

**P2-1 (new row; replaces the text of s4-report §7):**

| P2-1 | P2 | At P_B (ΛM² = 3/50, Q = (999/1000)Q_max, δ = 1/20), uniformly in l ≥ l₀ = 20: W_l(ω) ≠ 0 for every ω with Re ω ≥ 0 and −(½+δ)κ₋ ≤ Im ω ≤ 0. With the conjugation symmetry W_l(−ω̄) = conj W_l(ω) and P2-0 (Im ω ≥ 0), W_l has no zero in S_δ for any l ≥ 20; this contains P2-1w. With P2-3B, every l < l₀ is counted at P_B, on P2-3B's boxes R′_l (the assembly P2-4 is S5). Not P_A | WRITTEN_PROOF. Four regimes cover {l ≥ 20} × the closed half-strip; with h = 1/(l+½) ≤ h₀ = 2/41 and ζ = h Re ω, the cover's inequalities and tilings are certified in certificates/part1/cover_PB.json. **L** (Re ω ≤ 8κ₋): τ-scaled Frobenius series at both horizons (r = r_h ± h²τ), in Arb affine arithmetic over 1760 ω-boxes and 16 sub-intervals of h² ∈ [0, h₀²]; Re y₊ ≥ 0.1245 and Re y_c ≤ −0.0804 at the exits (y = Φ′/Φ); then Riccati invariance through the forbidden bulk. **S** (6κ₋ ≤ Re ω ≤ 0.9Ω_P/h): half-line coercivity on C³ contours below each turning point, uniform in h through ζ and 1/Re ω; validated curved-contour integration; a near-horizon Taylor model with a Rouché disc and Cauchy remainder for ζ ≤ 0.02Ω_P; min q̂ = 0.0130 (+) and 0.0091 (c) over 47 and 45 ζ-leaves. **T** (0.8–2.0 Ω_P/h): S2's rotated contour with the exact splitting h²E = E₀ + hE₁ + h²E₂ and a linear ground state; E ≥ 0.1003·2Ω_Pγ₀/h for every h ≤ h₀. **A** (≥ 1.5Ω_P/h): E ≥ 6.955·2Ω_Pγ₀/h at 1.5Ω_P/h, and monotone in Re ω above. Machine-checked (Arb, 128-bit): every cell, tail, exit and order condition, the bulk facts (U monotone on each side of r_P; W̃ ≥ −0.00792) and the W1 crossing checks. Written (notes/s4-proof.md): W1, the vertical-strip analyticity of r(x) and the simply-connected union of the strips and S2's checked windows (extending S2 Lemma 1); W2, the Jost continuation and decay on the curved throat contours and the near-horizon map (extending S2 Lemma 3); and routine W3–W9. Controls: N1–N5 do not certify (N1, N2: PS modes at l = 20, 50, 100 and an h-box holding the l = 20–50 modes; N3: NE modes at l = 1, 2 by the continuation path, while mode-free boxes certify; N4, N5: exits misplaced into the allowed region); K2: W_l from the L code agrees with zsol and arbw at 12 points (relative radii ≤ 1.5e−38); K4 (float only, not used by the proof): W_l ≠ 0 at 200 points; the secant search stopped in H⁺ at 61 of them, none a zero (residual and zsol); the launch-time rule counted such stops as zeros, and Amendment 6 (post hoc) replaced it. Pre-registered (notes/s4-preregistration.md) with six dated amendments. Independent review 2026-10-18: proven as stated. Every certificate was replayed in full; reviewer code sharing nothing with p2 confirms every T, A, S and L bound (each certified bound lies below the reviewer's value; L at actual l from 20 to 10⁵; both closing boxes by the exact near-horizon map); S3's counting pipeline with the P2-2 threshold independently finds no zero of W₂₀ in R′₂₀; the closing-box truncation term was corrected after the review. The Im ω > 0 part is at P2-0's level | problems/P2/src/p2/part1.py, part1_T.py, throat.py, lowfreq.py, affine.py; problems/P2/certificates/part1/; problems/P2/notes/s4-preregistration.md, s4-proof.md, s4-report.md; problems/P2/tests/test_part1.py; docs/review-2026-10-18-p2-s4-adversarial.md | 2026-10-18 |

If the l = 20 count is not wanted in the ledger, drop that clause; nothing else depends on it.

**P2-1w (optional note).** Append to its statement column: "Contained in P2-1 (2026-10-18)." Its level
and evidence are unchanged.

**FORMULATION v0.6: corrections.**
- "With P2-3B (l ≤ 19), P2-3 is complete at P_B" → "With P2-3B, every l < l₀ = 20 is counted at P_B,
  on P2-3B's boxes R′_l".
- "C² contours" → "C³ contours" (as built; C² is all W3 needs, so optional).
- After the merge: the header ("not yet reviewed") and "Open. The independent review of S4" should
  record this review and S1's correction.

**P2-3B.** No change. P2-1 removes P2-3B's caveat "l₀ comes from S4" in substance; the row can stay.

## 8. Before the S4 branch is merged into main

1. **S1.** Fix `err_cell` and the tail `UL` in `throat.certify_closing` (complex enclosure or disc
   bounds); regenerate the two closing-box leaves and re-assemble `S_PB_plus.json`, `S_PB_c.json`;
   re-run `cover` and `verify`. Add a fast test that compares the error term with complex ρ, or a
   closing-box negative control.
2. **S2.** Append a dated correction to Amendment 6 and fix s4-report §3–§4: the launch-time rule
   (`assemble_K4` at `35201758`) and K4's verdict under it; the overshoot-and-return mechanism.
3. **S3.** Add the simply-connected union (or the trapezoid argument) to W1/W2 in s4-proof.md.
4. **S4.** Use the §7 row; correct the s4-report verdict, §7 and FORMULATION v0.6 wording.
5. **S5.** Correct s4-report §4 on the pilot discrepancy.
6. **W1–W8** as convenient.
7. The review's scripts: the owner decides whether they should be committed under
   `problems/P2/certificates/review_s4/` (their creation was denied in this session).

Under the repository's rules the merge into main must be run by the user.
