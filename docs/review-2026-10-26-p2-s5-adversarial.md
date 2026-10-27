# Independent adversarial review: P2 stage S5 (the assembly P2-4 and the conversion P2-5 at P_B)

Date: 2026-10-26. Reviewer: a separate agent with no part in authoring P2 S5.

Scope: branch `p2/s5-assembly` at `aaeb62dc`, from main `3dd0b2bc`. The S5 commits are:
- `e2c75aa6` (18:24:16): the plan, `notes/s5-plan.md`;
- `7aa9fbd0` (18:30:09): `src/p2/assembly.py` and `tests/test_assembly.py`;
- `942f89c1` (18:30:32): `certificates/assembly/assembly_PB.json`;
- `a04a4323` (18:39:56): `notes/s5-assembly.md` and FORMULATION v1.0;
- `adac7b31`, `aaeb62dc` (18:42:10, 18:42:28): `notes/s5-report.md`.

The review worked in its own worktree (`.claude/worktrees/p2-s5-review`, branch `p2/s5-review`) and read
the author's files without changing any. It adds this document and `problems/P2/certificates/review_s5/`
(scripts, outputs and a README with the replay commands). No tool call was denied.

**Sources.** Fetched as TeX from arXiv on 2026-10-26 and read by the reviewer, not taken from the author's
quotes:
- **[HV17]** arXiv:1512.08004**v2** (posted 2020-05-27, "the published version"), `cauchyhor-r.tex`, TeX
  date line "December 25, 2015. Final revision: July 12, 2017.", sha256 `a0d3302f…33cf2be`;
  J. Math. Phys. 58, 081509 (2017). It is the file the author cites: every line number in
  `s5-assembly.md` matches it.
- **[CF17]** arXiv:1607.01018**v2** (posted 2017-10-23), `RN_deSitterRevised6.tex`, sha256 `dcaab0e0…a37464d453`;
  Ann. Henri Poincaré 18, 3371–3398 (2017), DOI 10.1007/s00023-017-0592-z (Crossref).
- Read: HV17 §1 (including §1.1–1.2) and §2.1–2.7 in full; CF17 §1–3, Assumption 2.2 and the
  propagation to the Cauchy horizon (§4, "Propagating the energy estimates…"). Theorem numbering was
  re-derived by counting the environments: HV17 Def 2.1, Prop 2.4, Prop 2.12, Thm 2.13, Lemma 2.14,
  Prop 2.17, Cor 2.19, Thm 2.24; CF17 Assumption 2.2, Thm 3.1. All as the author has them.

**Overall verdict:** P2-4 and P2-5 at P_B are proven as stated, P2-4 at **WRITTEN_PROOF** and P2-5 at
**WRITTEN_PROOF, conversion cited**. Nothing is blocking.
[Update 2026-10-26: the author's fixes (`3f2a17aa`) were re-checked and all hold; see §9.]
- **The decisive question is settled as the author reads it.** HV17's α-hypothesis (§2.5, ℓ. 1313) is about
  "the forward problem in Ω₂₃", and the next sentence says in so many words that the global problem may have
  other resonances in that half-space. That is P2-D1′'s notion, so P2-4 with P2-D1′ supplies the hypothesis
  exactly and P2-5 is not conditional. The reliance on HV17 Lemma 2.14(1)'s "more general" clause, which
  HV17 proves in outline, is legitimate: the outline is routine, and §4.C6 below writes it out.
- **P2-4's cover of S_δ is complete for every l**, and the assembly check is a faithful machine version of
  its bookkeeping. Reviewer code sharing nothing with `p2` re-checks every seam (280 checks, P_B recomputed
  from its definitions): b − (½+δ)κ₋ = 4.546·10⁻¹⁰ against an enclosure width of 2.7·10⁻⁵⁸ for (½+δ)κ₋;
  the a_l identities; the 2351 discs; the re-sums; l₀ = 20; and a 20 295-point sampled cover of S_δ for l
  up to 10⁶. 22 further mutations of the inputs each make the author's check fail on the right item.
- **Three should-fix items, all wording or precision in the ledger texts:**
  - S1: the decay rate e^{−αt₀} at α = d_δ exactly needs the contour argument of H7 finished;
  - S2: the two-ended / bifurcation-sphere scope is stated as "not uniform up to B", which says both too
    little and too much;
  - S3: the P2-5 row should name the HV17 steps that are outlines and the source of the order-1 claim at 0.

## 1. Verdicts

| item | verdict |
|---|---|
| A1. Plan before code; plan followed | YES. The plan (`e2c75aa6`, 18:24:16) precedes the code (`7aa9fbd0`, 18:30:09) and the certificate (18:30:32), and is unchanged at `aaeb62dc`. Deliverables, checks C1–C9, controls N1–N8 (N8 as N8a/N8b), acceptance and budget are as registered. The stop rule was correctly not triggered |
| A2. Certificate at a clean commit; replay | YES. `commit` = `7aa9fbd0` with no `+dirty` (the flag covers `src/p2`); all 54 inputs' sha256 equal the blobs on main `3dd0b2bc` and on disk; `python -m p2.assembly verify` → identical |
| B3. Cover of S_δ for every l | COMPLETE, including Im ω = 0, Im ω = −d_δ, Re ω = ±a_l, the l = 19/20 seam at h₀ = 2/41 exactly, the regime boundaries and 0 ∈ int R′₀. E₊ ∪ E_c lies below S_δ |
| B3. The assembly check | FAITHFUL for the bookkeeping; by design it trusts each ingredient certificate's own replay. N1–N8 exercise their items; 22 more mutations (M1–M22) all caught; a wrong a_l, a missing l and a wrong winding are each caught |
| B4. Zeros → resonances; ω = 0; β | CORRECT. S_δ ⊂ {Im ω > −min(κ₊, κ_c)} and misses E₊ ∪ E_c; at 0 only l = 0 contributes, m₀(0) = 1, so pole order 1, one chain of length 1, constant state (no Jordan block). β ≥ 0.55 (not strict), 0.55 ≤ β ≤ 0.90938098 |
| B5. Level of P2-4 | WRITTEN_PROOF, the minimum; no ingredient level misreported |
| C6. Resonance notion | LOCALIZED (Ω₂₃), confirmed from the TeX; identical to P2-D1′'s. Lemma 2.14(1)'s general clause: outline, routine, citable (S3: say so in the row). Global poles on the contour line: not treated by HV17, harmless, but H7's second route weakens the decay (S1) |
| C7. Other HV17 hypotheses | ALL MATCHED: non-degenerate RNdS, Q > 0, dimension 4, unique photon sphere in (r₊, r_c); massless scalar; α/γ₀ = 0.3268; β₁ = 1/κ₋, exponent 21/20 > 1. No hidden smallness assumption for RNdS in HV17 §2 |
| C8. CF17 | MATCHED. Price law with p = α on CF17's H⁺ = H⁺_A (the component in the v chart, adjacent to i⁺, crossed by HV17's H_I); t₀ − v constant on H⁺; T = ∂_v there. 2 min{p, κ₊} = 1.1κ₋ > κ₋, and κ₋ < 2κ₊ is implied by d_δ < κ₊ (W1) |
| C9. Overclaiming | Rows, FORMULATION v1.0 and the .tex text are accurate apart from S2 (scope at B) and W6–W7. The CF17 bibitem is correct |
| D10. Tests | Fast: 90 passed, 17 skipped (37.5 s). `P2_SLOW=1`: 107 passed, 1 pre-existing warning (685 s). `p2.assembly verify`, `p2.params verify`, `p2.lg_bounds verify PB 19`, `p2.part3 verify` (5% sample) all OK |

## 2. Findings by severity

### Blocking

None. Checked without finding a flaw that defeats P2-4 or P2-5:
- the decomposition S_δ = U ∪ H₊ ∪ H₋ and every seam (§4.B3), recomputed independently;
- the zero-to-resonance translation, multiplicity at 0 and the D3 infimum (§4.B4);
- the resonance notion in HV17 §2.5, Lemma 2.14 and the proof of Prop 2.17 (§4.C6);
- every other hypothesis of HV17 Thm 1.1/§§2.5–2.6 and CF17 Thm 3.1, including the Price-law conversion
  (§4.C7–C8).

### Should-fix before the merge

**S1. The decay rate at α = d_δ exactly rests on a contour argument that H7 leaves half-done.**
(`s5-assembly.md` §2.2 H7; the P2-5 row's "abs(u′) ≤ Ce^{−αt₀}"; FORMULATION v1.0's P2-5 bullet.)
- HV17's proof of Prop 2.17 shifts the inverse Mellin transform to Im σ = −α and expands over "the
  resonances of 𝒫 in −α < Im σ < −α̃" (ℓ. 1338). It silently assumes that the global P̂(σ)⁻¹ has no pole on
  the line. HV17 can afford this, since its α is existential; P2-5 fixes α = d_δ, so it must not.
- H7 offers two routes.
  - Route 2 ("run at α̂ ↑ α") keeps the Sobolev conclusion: H^{1/2+α̂/κ₋−0} for every α̂ < α is the same as
    H^{21/20−0}. But it gives decay e^{−α̂t₀} for each α̂ < α only, not the row's e^{−αt₀}. H7 does not say so.
  - Route 1 (a global pole on the line is not an Ω₂₃ resonance, so its singular part vanishes in r > r₁) is
    correct, but integrating through a pole needs one more step: subtract the singular parts at the finitely
    many poles on the line before integrating. The remainder is holomorphic near the line with the same
    high-energy bounds, and the subtracted terms vanish in r > r₁.
- **Simplest complete fix** (one sentence).
  - Ω₂₃ resonances are discrete, with finitely many in Im σ ≥ −γ for every γ < γ₀. They are among HV17's
    global resonances (restriction identity, §4.C6), of which there are finitely many in such strips
    (ℓ. 1530: "only finitely many resonances in any strip Im σ > −α > −γ₀").
  - None lies in Im σ ≥ −d_δ except 0 (P2-4). So some α′ ∈ (d_δ, γ₀) leaves 0 as the only Ω₂₃ resonance in
    Im σ ≥ −α′, and some α″ ∈ (d_δ, α′] has no global pole on Im σ = −α″.
  - HV17 at α″ gives decay e^{−α″t₀} ≤ e^{−d_δ t₀} in t₀ > 0 and regularity
    H^{1/2+α″/κ₋−0} ⊂ H^{21/20−0}: the row's statement at α = d_δ.
  - The same argument gives β > 0.55 strictly at "conversion cited" level. Nothing needs it, and the
    author's choice not to claim it can stand.
- Failure scenario: a reader checks the row's e^{−αt₀} with α = 0.55κ₋ against the note, meets route 2,
  and concludes that only e^{−(0.55κ₋−ε)t₀} is established; or reads route 1 as integrating through a pole.

**S2. The two-ended setting: "not uniform up to the bifurcation sphere" says both less and more than is
proved.** (FORMULATION v1.0, "Reading of P2-5"; `s5-report.md` §6 P2-5 row and §7 .tex items (c)–(d);
`s5-assembly.md` §2.1 (ii) and §2.5.)
- What is proved:
  - (a) for smooth data on HV17's one-ended H_I, u ∈ H^{1.05−0} in D(H_I), near the part of the Cauchy
    horizon to its future, towards i⁺;
  - (b) for smooth data on a complete Cauchy surface, the same near every point of the Cauchy horizon other
    than the bifurcation sphere B, by Killing translates of H_I and their mirror images (s5-assembly H2 —
    correct), with constants depending on the point;
  - (c) nothing at B.
- CF17 adds nothing at B either:
  - its estimate (mainEst) holds for 0 ≤ u₁ ≤ u₂ (ℓ. 810), the region next to i⁺ adjacent to the
    event-horizon component H⁺_A on which the Price law is assumed (the v chart, ℓ. 771). In CF17's
    orientation, u → +∞ at H⁺_A (ℓ. 526);
  - by translation, it gives bounded energy and C⁰ ∩ H¹_loc near every point of CH⁺_A other than B.
- The texts' "not uniform up to the bifurcation sphere" suggests a non-uniform statement at B, where there is
  none.
- The frozen §2 template reads "Every solution of □_gφ = 0 with smooth Cauchy data on RNdS … φ ∈ H¹_loc up to
  the Cauchy horizon". The v1.0 "reading" narrows it. That is a scope change, and the revision log should
  record it as one.
- Fix:
  - in FORMULATION v1.0, call it a narrowing of the §2 P2-5 template and state (a)–(c);
  - in the P2-5 row, replace "for two-ended data, not uniformly up to the bifurcation sphere" by "for data on a
    complete Cauchy surface: near every point of the Cauchy horizon other than its bifurcation sphere B
    (constants depend on the point); nothing at B";
  - in (ii), say that the extension is to the Cauchy-horizon component adjacent to the event-horizon component
    carrying the Price law, away from B;
  - in the .tex, change (c) and (d) as in §7.
- Failure scenario: a reader of the ledger or the .tex takes "H¹_loc up to the Cauchy horizon" for two-ended
  data to include B. That is the point an argument about extendibility across the whole Cauchy horizon would
  need.

**S3. The P2-5 row should say which HV17 steps are outlines, and where the order-1 claim at σ = 0 comes
from.** (`s5-report.md` §6, P2-5 row.)
- The conversion leans on two steps HV17 does not prove in detail:
  - Lemma 2.14(1)'s general clause (ℓ. 1259; "follows along the same lines", ℓ. 1266). This is what removes
    the global non-Ω₂₃ resonances in r > r₁ in the proof of Prop 2.17 (ℓ. 1342);
  - Lemma 2.14(2)'s order 1 at σ = 0, taken from Schwarzschild–de Sitter and Kerr–de Sitter references
    (ℓ. 1268), together with the closed upper half-plane, where HV17 recalls only the Im σ > 0 argument
    (ℓ. 1264).
- At P_B, P2 supplies the second itself: D4 (m₀(0) = 1), P2-0 (Im ω ≥ 0, including the real axis) and
  P2-D1′.
- The first is routine. §4.C6 gives the argument: a restriction identity, hyperbolic uniqueness in
  r₁ < r < r₂ − 2δ, and an induction over the Laurent coefficients. Citing it at "conversion cited" is
  acceptable.
- Fix: add to the row's level column: "HV17's removal of the global (non-Ω₂₃) resonances uses its
  Lemma 2.14(1), whose general clause HV17 proves in outline (checked by the S5 review); HV17's order 1 at
  σ = 0, taken there from Schwarzschild–de Sitter references, is supplied at P_B by D4, P2-0 and P2-D1′".
  The row in §7 does this.
- Failure scenario: an auditor of the "conversion cited" level assumes every cited step is proved in full in
  HV17.

### Wording and nits

- **W1. CF17's "substantive part"** (`s5-assembly.md` §2.4 F3; `s5-report.md` §4).
  - 2κ₊ > κ₋ is implied by d_δ < κ₊: 2κ₊ > 2d_δ = 1.1κ₋. P2-4 uses d_δ < κ₊ already (C4).
  - With p = d_δ, CF17's whole condition therefore carries nothing beyond d_δ < κ₊ and δ > 0. For the same
    reason N8b cannot isolate C8.km_lt_2kp: any κ₊ < κ₋/2 also breaks C4 (N8b fails 45 checks).
  - Say so. As the note says, the CF17 route's content is the Price law with p > κ₋/2.
- **W2. The controls table** (`s5-assembly.md` §1.4).
  - The "fails on" column lists a subset. N3 also fails `C2.l5.lower_strip_in_box_and_P2_2`; N7 also fails
    `C3.cover_replay_identical`; N8a also fails `C1.overlap.gamma0`; N8b fails 45 checks.
  - Add "among others". The JSON is complete.
- **W3. Name of H_I.** "HV17's Cauchy surface H_I" (P2-5 row, `s5-assembly.md`): H_I is a one-ended partial
  Cauchy surface (HV17 writes "'Cauchy' surface" in §1.2). Say "initial surface".
- **W4. Source versions.**
  - `s5-plan.md` calls the HV17 file the "final revision of 2017-07-12", and `s5-assembly.md` quotes its
    date line. It is arXiv v2, posted 2020-05-27 ("the published version"), whose TeX date line reads
    "Final revision: July 12, 2017".
  - Record the arXiv version for both sources (CF17 is v2, 2017-10-23).
- **W5. An HV17 typo, harmless here.**
  - §2.3 (ℓ. 997) places the complex absorption in {r₁ − δ < r < r₂ − δ}. That is inconsistent with "a
    small neighborhood of πL₀ ∪ πΓ_*", which lies in r₀ ≤ r < r₁, and with "modifications … in r ≤ r₁"
    (ℓ. 992). The intended set is {r₀ − δ < r < r₁ − δ}.
  - The Lemma 2.14(1) argument needs 𝒬 = 0 on Ω₂₃, which the intended support gives. S5 depends on nothing
    else here, but a reader checking Lemma 2.14 should know about it.
- **W6. FORMULATION v1.0 status wording.** "The ladder at P_B is complete up to P2-5" should read
  "assembled up to P2-5 (P2-4 and P2-5 pending review and approval)" until the rows are entered.
- **W7. The .tex v1.3 entry.**
  - Item (c), "lies in H^{1.05−0}, in particular in H¹_loc, at the Cauchy horizon", needs S2's localization.
  - The optional Status sentence ("this region has since been established") should say at what level: a
    written proof with machine-checked certificates.
- **W8. Redundant or weak checks.**
  - `C2.l*.lower_strip_in_box_and_P2_2` repeats `a_identical` and `b_ge_d`.
  - `C1.Lambda_M` is a string match on the conventions field; the real check is the recomputed `Point`.
  - Both are harmless.
- **W9. Authoring order.** The plan says it was "written … before any S5 code". Commit order confirms the
  plan came first (6 minutes before 668 lines of code and tests), but authoring order cannot be checked.
  Nothing in S5 could favour a verdict post hoc: the check is deterministic bookkeeping, and the plan fixes
  the controls.

## 3. Integrity (scope A)

### A1. The plan against execution

| registered (`e2c75aa6`) | executed | verdict |
|---|---|---|
| Deliverables 1–5: `assembly.py` (run/verify), the fast test, `s5-assembly.md`, `s5-report.md` (rows, .tex text), FORMULATION v1.0 | all present | ✓ |
| Inputs: committed certificates only, pinned | 54 files, sha256 in the certificate, equal to main `3dd0b2bc` | ✓ |
| C1–C9 as listed in §3 of the plan | 366 named checks, each in the listed group | ✓ (C1 adds a string match for Λ, M; W8) |
| N1–N8, each failing on its named item | N1–N7, N8a, N8b, as specified (N1 at ⌊d_δ·2²⁸⌋/2²⁸ = 6280878/2²⁸, path moved) | ✓ |
| Acceptance: check passes; controls fail; suite green; proof from listed rows only; every §4 item quoted and matched; no frozen convention changed; rows and .tex proposed only; only `problems/P2/`; no AI attribution | met. The v1.0 "Reading of P2-5" narrows the §2 template without saying so (S2) | ✓ with S2 |
| Stop rule (global notion → stop) | not triggered, correctly | ✓ |
| Budget | ≈ 0.05 core-h | ✓ |

The plan file is byte-identical at `aaeb62dc` (`git diff e2c75aa6 aaeb62dc -- notes/s5-plan.md` is empty).
`git diff --name-only 3dd0b2bc aaeb62dc` lists only the seven `problems/P2/` files. No commit message
carries attribution.

### A2. The certificate

- `assembly_PB.json`: `commit` = `7aa9fbd0a85f…` with no `+dirty`. `arbutil.git_commit` checks
  `git status --porcelain` in `src/p2` only, but every input is pinned by sha256. All 54 match the blobs on
  main and the files on disk.
- `python -m p2.assembly verify` (review head): `identical`. It compares `result`, `negative_controls`,
  `inputs_sha256` and `pass`, and writes nothing.
- 366 checks pass, and N1–N8b all fail as required, as recorded.

## 4. Mathematics

### B3. The cover of S_δ for every l

S_δ = U ∪ H₊ ∪ H₋, with U = {Im ω > 0} (open) and H± = {±Re ω ≥ 0, −d_δ ≤ Im ω ≤ 0} (closed). So the line
Im ω = 0 lies in H₊ ∪ H₋, not in U.
- **U, every l:** P2-0 (Im ω ≥ 0; W_l(ω) = 0 iff (l, ω) = (0, 0)).
- **H₋:** W_l(−ω̄) = conj W_l(ω) maps H₋ onto H₊ (Z4 recorded; also immediate from §1.2).
- **l ≥ 20:** P2-1 on the closed H₊.
  - h = 1/(l+½) ≤ h₀ = 2/41 iff l ≥ 20 (1/19.5 > 2/41 = 1/20.5), so l = 20 sits exactly at h₀, and P2-1's
    h-ranges are closed at h₀ (L: h² ∈ [0, h₀²]).
  - The Re-ranges [0, 8κ₋], [6κ₋, 0.9Ω_P/h], [0.8, 2]Ω_P/h and [1.5Ω_P/h, ∞) union to [0, ∞) for every
    h ≤ h₀. The only non-exact overlap is S's top against L's: 0.9Ω_P/h₀ = 90.75κ₋ > 8κ₋, and the worst
    case is h = h₀ because Ω_P/h falls with h. Independently recomputed.
  - Im ∈ [−d_δ, 0] is covered through L's design range [−11/20, 0], T's and A's ν_low balls and the S tilings
    in `cover_PB.json`.
- **l ≤ 19:** R′_l = [−a_l, a_l] × [−b, e] (P2-3B) and {abs(Re ω) ≥ a_l} × [−d_δ, 0] (P2-2).
  - The lower strip is covered because [−d_δ, 0] ⊂ [−b, e] and the Re-ranges meet at ±a_l.
  - a_l is the same rational in both certificates, for all 20 l.
  - b = 6280879/2²⁸ is the least 2⁻²⁸-dyadic ≥ d_δ, with b − d_δ = 4.54603·10⁻¹⁰. The reviewer's interval
    for d_δ has width 2.7·10⁻⁵⁸ and `params.json`'s κ₋ ball has radius 5·10⁻⁴²: rigorous by about 30
    orders of magnitude. A small margin is no fragility here; the design only needs b ≥ d_δ.
  - e = 4567911/2³⁰ = 0.09999998κ₋ is > 0 and ≤ η, so the band 0 < Im ω ≤ e lies in both R′_l and U.
  - The corners (±a_l, −d_δ), (±a_l, 0) lie in both pieces.
- **The zero at 0.** 0 ∈ int R′₀ (a₀, b, e > 0); W₀(0) = 0 (R ≡ 1); the count in R′₀ is 1, so ω = 0 is the
  only zero in R′₀ and is simple. D4 confirms it independently: W₀′(0) = 5.54913720005i.
- **Poles.** E₊ ∪ E_c lies below S_δ: d_δ < κ₊ (margin 0.01161) and d_δ < κ_c (margin 0.06998). Also
  b < κ₊, so every R′_l lies above −iκ₊. Every disc box clears E₊ ∪ E_c, with min dist²/(√2ρ)² = 2.26–2.51
  over l.

**Is the assembly check a faithful machine version?** Yes, for the bookkeeping: identical a_l, b ≥ d_δ,
e > 0, the exact traced right half and its discs, pole clearance, the re-summed windings, the h₀ ↔ l₀ seam,
the regime union, the Im ranges, the ingredient flags and the β/HV17/CF17 numerics. By design it does not
re-verify an ingredient's analytic content (the lg criterion, the increment balls, the regime inequalities);
the ingredients' own replays do that. The review re-ran the cheap ones and a 196-segment sample of S3's (§5).
- The controls: N1–N8 each fail on their named item.
- The review's 22 further mutations (`mutations.py`) are each caught. They are: a gap from the P2-2 side; a
  missing lg certificate; a removed piece; a 2⁻⁴⁰ path gap; a moved b; an increment off by π; the l = 6
  certificate in l = 7's slot; an lg criterion ≥ 2; a narrower lg strip; an uncertified winding; an S tiling
  flag; a missing regime file; L's Im range; a T box short of −d_δ; an oversized disc; a moved Krawczyk
  enclosure; inconsistent `params`; the l = 0 count; ζ_SU; h₀ = 2/39; a P2-0 step; Z4.
- So a wrong a_l, a missing l and a winding of the wrong value are all caught.
- What it cannot see (B2, B3): an increment ball changed within the ¼ tolerance, or an lg matching point
  changed. Both pass, as they should; only `p2.part3 verify` and `p2.lg_bounds verify` can see them.

### B4. From zeros to resonances; β

- **P2-D1′'s scope.** S_δ misses E₊ ∪ E_c and lies in {Im ω > −min(κ₊, κ_c)}, so both the set equivalence
  and the multiplicity clause apply on S_δ.
- **ω = 0.**
  - W_l(0) ≠ 0 for l ≥ 1 (P2-0), and m₀(0) = 1 (D4 and the R′₀ count).
  - P2-D1′(2) at σ = 0 gives pole order max_l m_l(0) = 1, one Jordan chain of length 1 in sector (0, 0), and
    algebraic multiplicity Σ(2l+1)m_l(0) = 1. So there is no Jordan block and no generalized state.
  - The state is e^{−i·0·T}R₊ = 1, a constant.
  - This is exactly HV17's "simple resonance at σ = 0, with resonant states equal to constants" (order 1,
    rank 1; Lemma 2.14(2)).
- **β.** D3's α is an infimum over all l of −Im ω over nonzero QNMs. Every nonzero QNM has Im ω < −d_δ
  (P2-4(a)), so α ≥ d_δ and β ≥ 11/20, without strictness. That is what HV17's closed half-plane hypothesis
  needs.
- **The two-sided bound.** P2-K6's certified ω_* = −0.90938097899480545552i κ₋ is a zero of W₀ off E₊ ∪ E_c
  (κ₊/κ₋ = 0.8229, 2κ₊/κ₋ = 1.6458, κ_c/κ₋ = 2.1951), so β ≤ 0.90938098 (margin 1.0·10⁻⁹). The interval
  0.55 ≤ β ≤ 0.90938098 is correct and, as the author says, not sharp.

### B5. The level of P2-4

WRITTEN_PROOF is the minimum of the ingredients: P2-0, P2-D1′ and P2-1 are WRITTEN_PROOF, and each is
needed for infinitely many (l, ω). The CERTIFIED rows (P2-par, P2-2, P2-3B) cannot replace them. The
assembly's own steps are routine: the set decomposition, the argument principle given certified counts, and
the D3 infimum. The ingredient table in `assembly.py` matches `CLAIMS.md` on main, and the test enforces
this. No level is misreported.

### C6. The resonance notion (the decisive check)

**The hypothesis.** HV17 §2.5 opens (ℓ. 1313, verbatim in `quotes_check.py`):

> The only resonance of the forward problem in Ω₂₃ in Im σ ≥ 0 is a simple resonance at σ = 0, with resonant
> states equal to constants, see the references given in the proof of Lemma 2.14, and there exists α > 0 such
> that 0 is the only resonance in Im σ ≥ −α. (This does *not* mean that the global problem for 𝒫 does not have
> other resonances in this half space!) In the notation of Proposition 2.12, we may assume −α > −γ₀ …

Prop 2.17, Cor 2.19 and Thm 2.24 take "α > 0 as above".
- The second clause continues the sentence about the Ω₂₃ forward problem.
- The parenthesis would contradict the hypothesis if the hypothesis were about the global P̂(σ)⁻¹.
- The notion is the localized one.
- P2-D1′ is stated for "the forward problem on Ω₂₃ (HV17 §2.5)". Its proof (`s1-p20-d1prime.md` §4)
  covers HV17's region [r₂ − 2δ, r₃ + 2δ]: no condition at the final boundaries, hyperbolic in r beyond the
  horizons, and smoothness above threshold selecting the analytic branch.
- The author's reading is right.

**How HV17's proof uses the global problem.**
1. *Lemma 2.15 (solvability)* uses only Lemma 2.14(1)'s Im σ > 0 part: mode stability on Ω₂₃, which is
   P2-0.
2. *Prop 2.17* expands over global resonances in −α < Im σ < −α̃ (ℓ. 1338). It removes, in r > r₁, those
   that are not Ω₂₃ resonances by Lemma 2.14(1)'s general clause (ℓ. 1342), and the one at 0 by
   Lemma 2.14(2). By hypothesis, nothing else is in the strip.
3. *The contour line* Im σ = −α: see S1.

No other use of the global notion appears in §§2.4–2.6.

**Lemma 2.14(1)'s general clause, written out** (HV17's outline, ℓ. 1266, made explicit). Let σ₀ not be
an Ω₂₃ resonance. Write Y₂₃ = [r₂ − 2δ, r₃ + 2δ] × S², and let v(σ) = P̂(σ)⁻¹f.
- 𝒬 vanishes on Y₂₃: it is supported near L₀ ∪ Γ_*, in r < r₁ (W5). So P̂(σ) = P̂₂₃(σ) there.
- The restriction v(σ)|_{Y₂₃} lies in the Ω₂₃ space: above threshold at r₂ and r₃, extendible at both ends.
  By uniqueness for the Ω₂₃ problem, v(σ)|_{Y₂₃} = P̂₂₃(σ)⁻¹(f|_{Y₂₃}) near σ₀, and this is holomorphic at
  σ₀. So every Laurent coefficient A_k f, k ≥ 1, of P̂(σ)⁻¹f at σ₀ vanishes on Y₂₃.
- The coefficients satisfy Σ_j P̂_j A_{k+j} = 0 for k ≥ 1, where P̂(σ) = Σ_j P̂_j (σ−σ₀)^j is quadratic in σ.
  Work down from the top coefficient. Each A_k f then solves, in r₁ < r < r₂ − 2δ, a wave equation for which
  r is a time function (the symbol is |μ|ξ² − r⁻²|η|² there). Its source vanishes, and so do its Cauchy data
  at r = r₂ − 2δ. Hence A_k f = 0 in r > r₁.
- This is the "more precise result": the range of the singular part is supported in r ≤ r₁.
- The same restriction identity at σ = 0 shows that the restriction to r > r₁ of the global singular part
  extends the Ω₂₃ singular part. At P_B that has order 1, rank 1 and constant range (D4, P2-0, P2-D1′),
  which is Lemma 2.14(2) at P_B.
- The identity also shows that every Ω₂₃ resonance is a global one, which S1's fix uses.

**Verdict.** The citation is sound. The general clause is outlined in HV17 but routine, so citing it under
"conversion cited" is acceptable. The row should say that it is an outline (S3).

### C7. HV17's other hypotheses

- **Setting.** Def 2.1 asks for Λ, M, Q > 0 and three simple positive roots of μ. At P_B, Q = 1.00984,
  and Δ has four simple real roots r_n < 0 < r₋ < r₊ < r_c (exact bisection: four disjoint certified sign
  changes of a quartic). The dimension is 4.
- **Trapping.** §2.2 uses that (r − r_P)(r⁻²μ)′ < 0 on (r₂, r₃) (ℓ. 927). At P_B the other critical point,
  (3 − √(9 − 8Q²))/2 = 1.04127, lies below r₊ = 1.09435, so r_P = 1.95873 is the only critical point in
  (r₊, r_c).
- **No hidden smallness.** HV17 §2 makes no smallness assumption for RNdS; "small Q" appears only in Lemma
  2.2 and Remark 2.7, both illustrations.
- **Field and data.**
  - The equation is □_g u = 0, the frozen massless scalar.
  - H_I = {t_* = 0, r₂ − 2δ ≤ r ≤ r₃ + 2δ} (ℓ. 754) is one-ended. D(H_I) = Ω ∩ {r > r₁} (ℓ. 992) reaches the
    Cauchy horizon only towards i⁺ (t₀ → ∞).
  - "A minor extension" (ℓ. 388) for a two-ended H_I up to B is not used, and nothing at B is claimed (S2).
- **α < γ₀.**
  - γ₀ = (2r_P)⁻¹(μ(r_P)(2 − 3M/r_P))^{1/2} (ℓ. 1508) is D6's λ_P/2. Note (r⁻²F)″(r_P) =
    −2(2 − 3M/r_P)/r_P⁴, so the Lyapunov exponent is √(F(r_P)(2 − 3M/r_P))/r_P.
  - d_δ/γ₀ = 0.326756 and γ₀ − d_δ = 0.0482091, recomputed. No shrinking of α is needed.
- **Exponent.** β₁ = 2/|μ′(r₁)| (ℓ. 850), so 1/2 + αβ₁ = 1/2 + 11/20 = 21/20 exactly, and
  H^{21/20−0} ⊂ H¹ locally. Near r₁ at finite t₀, b-regularity is ordinary regularity.

### C8. CF17

- **Setting.**
  - CF17's RNdS is D(Σ) for a compact Σ = S¹ × S² (ℓ. 410).
  - Thm 3.1 needs φ ∈ C^∞(M_int ∪ H⁺) and the Price law on "H⁺ = {r = r₊}, covered by" the v chart (ℓ. 771).
    That chart covers M_int ∪ H⁺_A ∪ CH⁺_B (ℓ. 496), so CF17's H⁺ is H⁺_A, the component on which v → ∞ at
    i⁺.
  - Smooth data on Σ give φ smooth there.
  - The surface gravities are ½|D′| (ℓ. 455), the same as ours.
- **Price law with p = α.**
  - HV17's H_I crosses H⁺_A. Near r₂, t₂₃ = t − F₂₃ with F₂₃′ = −(μ⁻¹ + c₂), so t₂₃ = v + (smooth in r), and
    on H⁺ t₀ − v is constant. HV17 itself says t₀ is "equivalent … to v near the event horizon" (ℓ. 355).
  - T = ∂_v on H⁺_A (ℓ. 498), and T and the Ω_i are stationary and tangent to {r = r₁}. So Thm 1.1 bounds
    T^{l+1}Ω^I u′ and Ω_iT^lΩ^I u′ by C e^{−αv} on H⁺_A ∩ D(H_I).
  - u₀ is constant, so only u′ enters. Moreover, |∇̸ψ|² = r⁻²Σ(Ω_iψ)² (ℓ. 740).
  - The integral over [v₁, v₂] × S² is ≤ (2πC/α)e^{−2αv₁}, and the compact initial piece is bounded by
    smoothness. The derivation in `s5-assembly.md` F2 is correct.
- **The condition.** 2 min{p, κ₊} = 2p = 1.1κ₋ > κ₋, because p = d_δ < κ₊. The κ₋ < 2κ₊ part is implied
  (W1).
- **What CF17 then gives.** Bounded energy for u ≥ 0 (ℓ. 810), and C⁰ ∩ H¹_loc on M_int ∪ H⁺ ∪ CH⁺ there,
  that is, towards i⁺_A. By translation, it holds near every point of CH⁺_A other than B (S2).

### C9. Overclaiming

- **FORMULATION v1.0.** Accurate apart from the "Reading of P2-5" (S2) and "complete" (W6). The scope fences
  of §7 are unchanged (the diff touches only the header and the revision log).
- **s5-report §6 rows.** Accurate apart from S1–S3; the corrected rows are in §7.
- **.tex v1.3 entry.**
  - (a) and (b) are accurate. (b) says honestly that the result is "a written proof with machine-checked
    certificates … not a fully machine-checked certificate".
  - (c) needs S2's localization. "The linear analogue of the failure of the Christodoulou formulation for
    smooth data" is right for the stated data class, away from B.
  - (d) needs S2's wording.
  - The entry claims nothing at P_A, and does not claim to be "first".
- **The CF17 bibitem.** Correct: Costa and Franzen, Ann. Henri Poincaré 18, 3371 (2017), arXiv:1607.01018;
  pages 3371–3398 and DOI 10.1007/s00023-017-0592-z per Crossref.

## 5. Independent checks, replays and tests (scope D)

All under `problems/P2/certificates/review_s5/` (README there), at the review head, `nice -n 10`.

| check | what | result |
|---|---|---|
| `indep_assembly.py` → `indep_assembly_out.json` | Shares no code with `p2`; exact rationals and `mpmath.iv` (80 digits). It recomputes P_B from FORMULATION §1.3: r_ext, Q_max², the four roots by exact bisection, κ's, r_P, Ω_P, γ₀, all overlapping `params.json`. It then checks the seams (B), all 20 lg/winding pairs (C), l₀ and the regime union (D), a 20 295-point sampled cover of S_δ × {l} with seam points and l up to 10⁶ (E), and β and D4 (F) | 280/280 pass |
| `mutations.py` → `mutations_out.json` | M1–M22 and B1–B3 against the author's `check` | M1–M22 caught on the named items; B2, B3 pass (by design); B1 caught incidentally |
| `quotes_check.py` → `quotes_check_out.json` | 38 HV17/CF17 passages the review relies on, found verbatim in the TeX, with line numbers | 38/38; source sha256 as recorded |
| `p2.part3 verify --l 0-19 --sample 0.05` → `replay_part3_PB_sample.*` | the author's replay: exact re-sum and coverage of all 20 certificates and the standard 132-segment sample | OK (780 s CPU) |
| `part3_fresh_sample.py` → `part3_fresh_sample_out.*` | 64 different segments (4 per l, l = 3–18; the S3 review replayed l = 0, 1, 2, 19 in full) recomputed with the author's code | 64/64 overlap (350 s CPU) |
| `p2.assembly verify`; `p2.params verify`; `p2.lg_bounds verify PB 19` | the author's replays | identical; OK; OK |
| `pytest problems/P2/tests` | fast | 90 passed, 17 skipped (37.5 s) |
| `P2_SLOW=1 pytest problems/P2/tests` | slow | 107 passed, 1 pre-existing warning (S2's m3 overflow), 685 s |

No run changed a tracked file: `git status --short` showed only the two new paths (this document and `review_s5/`) before the commit.

## 6. Compute and tool notes

| item | CPU |
|---|---|
| pytest, fast and `P2_SLOW=1` | 0.20 core-h (714 s) |
| `p2.part3 verify` (standard 5% sample, 2 workers) | 0.22 core-h (780 s) |
| fresh 64-segment sample (2 workers) | 0.10 core-h (350 s) |
| `assembly verify`, params, lg, `indep_assembly`, `mutations` (2 runs), `quotes_check` | < 0.01 core-h (≈ 15 s) |
| **total** | **≈ 0.52 core-h** (cap 2) |

All runs were at nice 10, with at most 3 processes of the reviewer's computing at once: the slow suite
alongside one 2-worker replay. The two replays did not overlap. No other process was touched.
Denied tool calls: none.

## 7. The ledger

**Level of P2-4: WRITTEN_PROOF.**
- Every seam is machine-checked, and independently re-checked.
- The ingredients that exclude infinitely many (l, ω) are written proofs: P2-0 (U), P2-1 (l ≥ 20) and
  P2-D1′ (the resonance form).
- The lower strip for l ≤ 19 is CERTIFIED.

**Level of P2-5: WRITTEN_PROOF, conversion cited.**
- HV17's hypotheses are exactly P2-4(b), P2-0/D4/P2-D1′ at σ = 0, and P2-par's d_δ < γ₀.
- CF17's are the Price law derived from HV17 and d_δ < κ₊.
- Neither theorem is re-proven. One cited step is an outline in HV17, checked here (S3).

**P2-4 (new row; the author's text with the review filled in):**

| P2-4 | P2 | At P_B (ΛM² = 3/50, Q = (999/1000)Q_max, δ = 1/20), for every l ≥ 0: W_l has no zero in S_δ = {Im ω ≥ −(½+δ)κ₋} other than the zero of W₀ at ω = 0, which is simple (m₀(0) = 1). Equivalently (P2-D1′), the forward problem on Ω₂₃ (HV17 §2.5) has no resonance in S_δ∖{0}, and ω = 0 is a resonance of multiplicity 1 whose resonant states are the constants. Hence D3's gap satisfies β ≥ ½ + δ = 0.55 (strict inequality not claimed); with P2-K6, 0.55 ≤ β ≤ 0.90938098 (two-sided, not sharp; the two-sided enclosure is the stretch goal P2-6). Assembled from P2-0 (Im ω > 0, every l), the conjugation symmetry (Re ω < 0), P2-1 (l ≥ 20), P2-2 with P2-3B (l ≤ 19: abs(Re ω) ≥ a_l, and R′_l) and P2-par (E₊ ∪ E_c below S_δ; S_δ inside P2-D1′'s multiplicity window). Not P_A | WRITTEN_PROOF, the minimum of the ingredients: P2-0, P2-D1′ and P2-1 are WRITTEN_PROOF; P2-par, P2-2 and P2-3B are CERTIFIED (for l ≤ 19 the lower strip −(½+δ)κ₋ ≤ Im ω ≤ 0 is CERTIFIED). Machine-checked bookkeeping (exact rationals and Arb, 256-bit; 366 checks; the 54 input certificates pinned by sha256): for each l ≤ 19, the a_l of R′_l is identical to P2-2's stored threshold, b − (½+δ)κ₋ = 4.546e−10 > 0, e > 0 and b < κ₊; the 2351 stored discs trace the right half of ∂R′_l exactly and clear E₊ ∪ E_c; the re-summed windings are 1 (l = 0) and 0. P2-1's cover replays identically; h₀ = 2/41 gives l₀ = 20, so no l is missed; the four regimes cover Re ω ≥ 0 for every h ≤ h₀. (½+δ)κ₋ < κ₊, κ_c. Negative controls N1–N8 fail as required: b just below (½+δ)κ₋; one l's certificate removed; a₅ lowered by 2⁻²⁴; h₀ = 2/43; e = 0; a wrong winding; a failed regime; γ₀ and κ₊ lowered. Written (notes/s5-assembly.md §1, routine): the decomposition of S_δ, the argument principle given the certified counts, the D3 infimum. Independent review 2026-10-26: proven as stated. Reviewer code sharing nothing with p2 (exact rationals and mpmath intervals, P_B recomputed from its definitions) re-checks every seam: the a_l identities, b − (½+δ)κ₋ (its enclosure width 3e−58), the discs and their pole clearance, the exact re-sums, l₀, the regime union and a 20 295-point sampled cover of S_δ for l up to 10⁶; 22 further mutations of the inputs each fail the check; 196 S3 segments were recomputed (the standard 5% sample and a fresh one) | problems/P2/src/p2/assembly.py; problems/P2/certificates/assembly/assembly_PB.json; problems/P2/notes/s5-plan.md, s5-assembly.md, s5-report.md; problems/P2/tests/test_assembly.py; docs/review-2026-10-26-p2-s5-adversarial.md; problems/P2/certificates/review_s5/ | 2026-10-26 |

**P2-5 (new row; S1–S3 applied):**

| P2-5 | P2 | At P_B, with α = (½+δ)κ₋ = 0.55κ₋ (< γ₀ = 1.683κ₋), by citation. (i) [HV17] (arXiv:1512.08004) Thm 1.1 in its quantitative form (Prop 2.17, Cor 2.19, Thm 2.24): every solution of □_g u = 0 with smooth initial data on HV17's one-ended initial surface H_I (it crosses the event and cosmological horizons) equals u₀ + u′ in the domain of dependence of H_I, with u₀ ∈ ℂ and abs(u′) ≤ Ce^{−αt₀} uniformly in r > r₋ (also after any number of stationary derivatives tangent to the Cauchy horizon). u is continuous up to the Cauchy horizon, and u′ ∈ e^{−αt₀}H^{1/2+α/κ₋−0} = e^{−αt₀}H^{1.05−0} in t₀ > 0; in particular u ∈ H¹_loc up to the part of the Cauchy horizon in that domain of dependence (towards i⁺). For smooth data on a complete Cauchy surface, (i) applied to Killing translates of H_I and their mirror images in the other exterior gives the same near every point of the Cauchy horizon other than its bifurcation sphere B, with constants depending on the point; nothing is claimed at B. (ii) [CF17] (arXiv:1607.01018) Thm 3.1 with p = α, for smooth data on a Cauchy surface: by (i) the Price law (Assumption 2.2) holds with p = α on CF17's event-horizon component H⁺ (the one adjacent to i⁺), and 2 min{p, κ₊} = 1.1κ₋ > κ₋; so the energy is bounded up to the adjacent Cauchy-horizon component, and the solution extends to it in C⁰ ∩ H¹_loc, near every point of it other than B. This is the linear analogue at P_B of the failure of Christodoulou's formulation of SCC for smooth data. Massless scalar on fixed RNdS only; not P_A | WRITTEN_PROOF (P2-4's level), conversion cited: neither theorem is re-proven. HV17's α hypothesis (§2.5), quoted from the TeX: "The only resonance of the forward problem in Ω₂₃ in Im σ ≥ 0 is a simple resonance at σ = 0 … and there exists α > 0 such that 0 is the only resonance in Im σ ≥ −α. (This does not mean that the global problem for P does not have other resonances in this half space!)" So P2-4 with P2-D1′ supplies the hypothesis exactly. In the proof of HV17 Prop 2.17 the global resonances that are not Ω₂₃ resonances are removed in r > r₁ by HV17 Lemma 2.14(1), whose general clause HV17 proves in outline (checked by the S5 review); global poles on the contour line are avoided by applying HV17 at some α′ ∈ (α, γ₀) that leaves 0 the only Ω₂₃ resonance in Im σ ≥ −α′ and puts no global pole on Im σ = −α′ (it exists because both kinds of resonance are finitely many in Im σ ≥ −γ for γ < γ₀); this gives the stated decay and regularity at α. 0 simple with constant states: D4 with P2-0 and P2-D1′ (HV17 Lemma 2.14(2) takes the order 1 at 0 from Schwarzschild–de Sitter references). From P2-par: α < γ₀ (margin 0.0482) and α < κ₊ (margin 0.0116), which also gives CF17's κ₋ < 2κ₊. The exponent ½ + α/κ₋ = 21/20 is exact. In CF17's condition, 2p > κ₋ holds for every δ > 0, so the CF17 route's content is the Price law with p > κ₋/2, i.e. P2-4. The numerical hypotheses are checked in certificates/assembly/assembly_PB.json (C7, C8, with negative controls). Independent review 2026-10-26: proven as stated, by citation; every HV17 and CF17 passage used was checked against the arXiv TeX (HV17 v2, CF17 v2) | problems/P2/notes/s5-assembly.md §2; problems/P2/certificates/assembly/assembly_PB.json; docs/review-2026-10-26-p2-s5-adversarial.md; problems/P2/certificates/review_s5/quotes_check.py | 2026-10-26 |

**Optional, P2-1.** The parenthesis "(the assembly P2-4 is S5)" could read "(assembled in P2-4)", as the
report suggests. Nothing else in the existing rows needs to change.

**FORMULATION v1.0: corrections.**
- "The ladder at P_B is complete up to P2-5": add "(P2-4 and P2-5 pending review and approval)" until the
  rows are entered (W6).
- Replace the bullet "Reading of P2-5" by a narrowing, recorded as such (S2):
  > **Scope of P2-5 (a narrowing of the §2 template).** (a) Data on HV17's one-ended initial surface H_I: the
  > conclusions hold in its domain of dependence, near the Cauchy horizon towards i⁺. (b) Data on a complete
  > Cauchy surface: near every point of the Cauchy horizon other than its bifurcation sphere B, with constants
  > depending on the point (Killing translates of H_I). (c) Nothing at B. For [CF17], the C⁰ ∩ H¹_loc
  > extension is to the Cauchy-horizon component adjacent to the event-horizon component carrying the Price
  > law, with the same exception.
- In the P2-5 bullet, after "with α = 0.55κ₋", add "(decay e^{−αt₀} at α exactly via some α′ ∈ (α, γ₀);
  s5-assembly H7)" (S1).

**The proposed .tex v1.3 entry: corrections** (the owner approves .tex changes).
- (c): replace "lies in $H^{1.05-0}$, in particular in $H^1_{\mathrm{loc}}$, at the Cauchy horizon." with "lies
  in $H^{1.05-0}$, in particular in $H^1_{\mathrm{loc}}$, near the part of the Cauchy horizon in its domain of
  dependence; for smooth data on a complete Cauchy surface the same holds near every point of the Cauchy
  horizon other than its bifurcation sphere."
- (d): replace "or uniform up to the bifurcation sphere of the Cauchy horizon" with "or at the bifurcation
  sphere of the Cauchy horizon".
- Optional Status sentence: "has since been established with $\delta=1/20$" → "has since been established
  with $\delta=1/20$, by a written proof with machine-checked certificates".
- The bibitem `costafranzen2017` is correct as proposed.

## 8. Before the S5 branch is merged into main

1. **S1.** Finish H7 in `s5-assembly.md` §2.2: give the α′ argument, or spell out route 1, and either drop
   route 2 or state its weaker decay.
2. **S2.** Fix the scope wording in FORMULATION v1.0 ("Reading of P2-5"), `s5-assembly.md` §2.1 (ii) and
   §2.5, and `s5-report.md` §6 (P2-5 row) and §7 (.tex (c), (d)), as in §7.
3. **S3.** Use the P2-5 row of §7, or add its two clauses on HV17's outlined Lemma 2.14(1) and on the source
   of the order-1 claim.
4. **W1–W9** as convenient. W1, W2 and W4 are one-line text fixes in the notes.
5. The owner enters P2-4 and P2-5 (§7), and approves the .tex v1.3 entry with the §7 corrections.

Under the repository's rules the merge into main must be run by the user.
[Correction 2026-10-26, §9: the user approves the merge into main, and the coordinator runs it.]

## 9. Follow-up check of the author's fixes (`3f2a17aa`), 2026-10-26

Scope.
- The author's commit `3f2a17aa` (21:05:40) sits on top of the review (`65af9e72`). The coordinator
  fast-forwarded it into `p2/s5-review`.
- It changes `FORMULATION.md` and the notes `s5-assembly.md`, `s5-plan.md` and `s5-report.md`, and
  nothing else: no code, test, certificate, `docs/` file or `CLAIMS.md` (`followup_check.py`, item 5).
- On 2026-10-26 the user approved S1–S3, W1–W9, the §7 rows, the P2-1 parenthesis "(assembled in P2-4)",
  the .tex v1.3 entry with the §7 corrections, the `costafranzen2017` bibitem, and the merge. The
  coordinator reported these approvals.

**Verdict:** every fix holds, and nothing is blocking.
- The §7 rows stand unchanged.
- One change is recommended to the .tex entry: one word, "Cauchy surface" → "initial surface" (5(i)).
  It applies the approved W3; the final text is below.
- One residual nit in FORMULATION v1.0 (item 2).

| item | verdict |
|---|---|
| 1. S1 (H7, the decay at α exactly) | DONE; complete and correct. Finiteness of both kinds of resonance, the choice of α′ then α″, the decay e^{−d_δt₀} and the regularity all check. HV17 imposes no cap beyond the radial-point threshold. The α̂ route is stated with its weaker decay and marked unused. The FORMULATION P2-5 bullet is right |
| 2. S2 (scope at B) | DONE. "Scope of P2-5" in FORMULATION v1.0 is the §7 text, recorded as a narrowing of the frozen §2 template, and §2 carries only a pointer. `s5-assembly.md` §2.1 (i)–(ii), §2.2 H2 and §2.5 match, and so do `s5-report.md` §§4, 6 and 8. No "not uniform up to B" remains. Residual nit: the stale parenthesis below |
| 3. S3 (the row names the outlines) | DONE. The §6 rows are byte-identical to the §7 rows (sha256 `624ee58e…` for P2-4 and `c3fd53dc…` for P2-5). `s5-assembly.md` §2.3 names both outlined HV17 steps correctly |
| 4. W1–W9 | DONE. The plan's pre-registered text (9976 bytes at `e2c75aa6`) is a byte prefix of the current plan; only the dated notes for W4 and W9 are appended |
| 5. The .tex entry | (c) and (d) are in verbatim, the dates are filled and the bibitem is unchanged. 5(i): change one word (below). 5(ii): accurate as it stands, with no phrase needed |
| 6. The §7 rows | They stand, with no wording change |
| 7. Tests | Fast: 90 passed, 17 skipped (36.1 s). `python -m p2.assembly verify`: identical |

**1. S1 (`s5-assembly.md` §2.2 H7).** The argument is complete and correct.
- **Finiteness.** HV17 ℓ. 1530 states that there are finitely many resonances in Im σ > −α for each
  α < γ₀. The closed half-plane Im σ ≥ −γ, γ < γ₀, lies inside such an open one.
  - These are the global resonances. The Ω₂₃ resonances are among them. If σ₀ were an Ω₂₃ pole but not
    a global one, take g ∈ C_c^∞ in the interior of Y₂₃ with P̂₂₃(σ)⁻¹g singular at σ₀; the singular
    part is finite-rank and continuous, so such g exist.
  - The restriction identity of §4.C6 gives (P̂(σ)⁻¹g)|_{Y₂₃} = P̂₂₃(σ)⁻¹g on a punctured neighbourhood
    of σ₀. The left side is holomorphic at σ₀, which is a contradiction.
  - H7 cites §4.C6 for this. The two lines above complete that reference.
- **α′ and α″.** Fix γ ∈ (d_δ, γ₀).
  - The finitely many nonzero Ω₂₃ resonances in Im σ ≥ −γ all have Im σ < −d_δ (P2-4(b)). So some
    α′ ∈ (d_δ, γ) leaves 0 as the only Ω₂₃ resonance in Im σ ≥ −α′.
  - The global poles in Im σ ≥ −α′ are finitely many, so some α″ ∈ (d_δ, α′] puts none on Im σ = −α″.
  - Then α″ satisfies §2.5 (and α″ < γ₀), and HV17's contour is pole-free.
- **The decay.** Cor 2.19 at α″ gives |u′| ≤ Ce^{−α″t₀} ≤ Ce^{−d_δt₀} in t₀ ≥ 0. The derivative bounds
  and the weighted space e^{−α″t₀}H^s ⊂ e^{−d_δt₀}H^s on t₀ > 0 follow likewise. So the decay holds at
  d_δ exactly.
- **The regularity, and caps.** HV17 at α″ gives H^{1/2+α″/κ₋−0} ⊂ H^{21/20−0}, since α″ > d_δ.
  - The only cap at r₁ is the radial-point threshold: a forward order function has s < 1/2 + β₁α″ near
    r₁ (Def 2.6, ℓ. 1029; Prop 2.17; Thm 2.24).
  - The stricter s < −1/2 + β₁α near r₁ ((EqRNdSFredholmFwOrder2), ℓ. 1232) only sets up the
    Fredholm theory. Thm 2.13 drops it ("We can strengthen the regularity at the Cauchy horizon by
    dropping …", ℓ. 1245).
  - The lower bounds s > 3/2 + max(β₂, β₃)α″ ((EqRNdSFredholmFwOrder), ℓ. 1207) sit in r > r₁ + 2δ′.
    They are compatible with s′ ≥ 0.
  - Nothing else limits the exponent. In any case the stated H^{21/20−0} is weaker than what HV17 gives
    at α″.
- **The α̂ route** is now stated with its weaker decay (e^{−α̂t₀} for each α̂ < α only) and marked
  unused. So is the through-the-pole route.
- **Strictness** (β > 0.55 at "conversion cited" level) is noted and not claimed: §1.3 Remarks and
  s5-report §8.
- **The FORMULATION P2-5 bullet** ("decay e^{−αt₀} at α exactly via some α′ ∈ (α, γ₀); s5-assembly
  H7") is right. Its α′ is H7's final exponent α″, a harmless difference of names; the P2-5 row's single
  α′ is described by both of its properties.

**2. S2.**
- FORMULATION v1.0's "Scope of P2-5" is the §7 text set as bullets (a)–(c) plus the CF17 sentence. It
  adds that it narrows the frozen §2 template, that the user approved recording it this way, and that
  §2 points to it.
- §2's P2-5 text is unchanged apart from the appended pointer "[v1.0: the scope of P2-5 at P_B is
  narrowed; see the Revision log, v1.0, "Scope of P2-5".]" (FORMULATION ℓ. 237–238), in the style of
  §2's earlier bracketed notes.
- The notes match:
  - `s5-assembly.md` §2.1 (i) has "u ∈ H¹_loc up to the part of the Cauchy horizon in that domain of
    dependence, towards i⁺", then the complete-Cauchy-surface sentence and "Nothing is claimed at B";
  - §2.1 (ii) localizes CF17 to the component adjacent to i⁺, away from B;
  - §2.2 H2 says "Nothing is claimed at B", and the "minor extension" is not used;
  - §2.5 restates (a)–(c);
  - `s5-report.md` §4 (H2), §6 and §8 match.
- `followup_check.py` (item 4) finds no "not uniform", "up to the bifurcation sphere" or "Reading of
  P2-5" in the four files. The only exception is HV17's own sentence about a "two-ended" H_I, which is
  quoted verbatim.
- **Residual nit (not blocking).** The v1.0 entry's first line still reads "assembled up to P2-5 (P2-4
  and P2-5 pending review and approval)". The same entry goes on to record the review and the user's
  approval, and the header says so too.
  - Replace the parenthesis by "(reviewed and approved 2026-10-26; not yet entered in the ledger)", or
    drop it when the rows are entered.
  - The two "Proposed, not yet in the ledger" bullets are accurate until then.

**3. S3.**
- The P2-4 and P2-5 lines in the code block of `s5-report.md` §6 are byte-identical to this review's §7
  rows (`followup_check.py`, item 1). The two other lines starting "| P2-4 |" and "| P2-5 |" in the
  report belong to its §5 levels table.
- `s5-assembly.md` §2.3 ("What the citation leans on") names both outlined steps correctly:
  - Lemma 2.14(1)'s general clause (ℓ. 1266), with the Laurent refinement covered by its "likewise for
    the more precise result" (ℓ. 1266), which refers back to the "iterative argument" of the Im σ > 0 case
    (ℓ. 1264), pointing to §4.C6;
  - the order 1 at 0, from SdS/KdS references (ℓ. 1268), with only the Im σ > 0 argument recalled
    (ℓ. 1264), and supplied at P_B by D4, P2-0 and P2-D1′.

**4. W1–W9.**
- W1: `s5-assembly.md` §2.4 F3, and `s5-report.md` §3 (N8b) and §4 (F3).
- W2: the table header now reads "among others; the JSON lists every failure", plus the N8b note.
- W3: "initial surface" in the notes and the rows; for the .tex see 5(i).
- W4: Sources, and plan note 1.
- W5: §2.3 item 3.
- W6: applied, with the residual nit above.
- W7: (c) and (d), and the level wording in the optional sentence. That sentence is kept as a separate,
  unapproved proposal.
- W8: §1.4, "Two weak checks".
- W9: plan note 2 and the report's status lines.
- The plan's pre-registered text is unchanged: the file at `e2c75aa6` is a byte prefix of the current
  one (9976 of 10723 bytes). The appended section says it "change[s] nothing in the plan", and the report
  still lists "no amendments". That is consistent: these are dated notes, not amendments.

**5. The .tex entry (`s5-report.md` §7).** The `\item[v1.3]` equals the pre-review text with exactly
the review's (c) and (d) corrections and the two dates ("approved October 26, 2026", "independently
reviewed October 26, 2026") filled in. The leading ⟨Month day⟩ is kept, and the bibitem is unchanged
(`followup_check.py`, item 2).
- **5(i): yes, change one word**: "on the one-ended Cauchy surface of \cite{hintzvasy2017}" → "on the
  one-ended initial surface of \cite{hintzvasy2017}".
  - This is W3, which the user approved, applied to the .tex.
  - After the (c) correction the same sentence sets "the one-ended Cauchy surface" against "a complete
    Cauchy surface". But H_I is not a Cauchy surface of anything; it is a partial one covering one
    exterior of two-ended RNdS.
  - In the strong-cosmic-censorship literature "one-ended" names a different class of spacetimes (one
    asymptotic end, e.g. collapse with a regular centre). "One-ended Cauchy surface" therefore invites
    the wrong reading.
  - "Smooth Cauchy data" stays: it means data for the Cauchy problem.
  - The present wording is HV17's own loose usage ("the Cauchy surface H_I", ℓ. 360; "'Cauchy' surface",
    ℓ. 473), so it is imprecise rather than wrong. The change is precision of the same kind as the §7
    corrections, not a rewrite.
- **5(ii): accurate as it stands; no phrase is needed.**
  - The CF17 sentence inherits the data classes and localization of the sentence before it, and (d)
    excludes B.
  - *Data on a complete Cauchy surface.* HV17 gives the Price law on every event-horizon component (one
    application per exterior), so Thm 3.1 applies to each black-hole interior from both sides. The
    adjacent Cauchy-horizon components together make up the whole Cauchy horizon minus B, the same set
    as in the HV17 sentence.
  - *Data on H_I.* D(H_I) ∩ M_int is a strip {u > u₁} of exactly the kind CF17's estimate covers
    (ℓ. 810; CF17's u is +∞ on H⁺_A, ℓ. 526). After any smooth extension of the data, Thm 3.1 gives the
    extension near the part of the Cauchy horizon in D(H_I), the HV17 sentence's set again.
  - The component adjacency describes a single application of Thm 3.1. It belongs in the ledger row,
    where it is.
- **Final paste-ready text.** This is `s5-report.md` §7 with the one word changed; it is also in
  `review_s5/followup_tex_item_v1.3.txt`. The bibitem stands as written in §7, and the leading ⟨Month day⟩
  is the coordinator's:

```latex
\item[v1.3] ⟨Month day⟩, 2026 (proposed October 26, 2026 by the P2 stage-S5 assembly, \texttt{problems/\allowbreak P2/\allowbreak notes/\allowbreak s5-report.md}; approved October 26, 2026). P2 only; the other problems are unchanged. \emph{P2 Status.} (a) At the primary pre-registered point $\Lambda M^2=3/50$, $Q/Q_{\max}=999/1000$, the resonance-free region of the entry point holds with $\delta=1/20$ for the linear massless scalar field: for every $\ell$, the strip $\{\operatorname{Im}\omega\ge-0.55\,\kappa_-\}\setminus\{0\}$ contains no quasinormal mode, and $\omega=0$ is a simple resonance with constant resonant states ($0.55\,\kappa_-\approx0.33\,\gamma_0$). Hence $0.55\le\beta\le0.90938098$ there; the upper bound comes from an interval-certified $\ell=0$ near-extremal mode. (b) The proof combines a written mode-stability argument on the closed upper half-plane; a large-$\ell$ exclusion, uniform in $\ell\ge20$, whose inequalities are certified in interval arithmetic and whose analytic steps are written; certified large-$|\operatorname{Re}\omega|$ exclusions and validated argument-principle counts for $\ell\le19$; and a written identification of the zeros of the radial Wronskian with the resonances used in \cite{hintzvasy2017}. The assembly is machine-checked. The result is a written proof with machine-checked certificates, independently reviewed October 26, 2026; it is not a fully machine-checked certificate. (c) By \cite[Thm.~1.1 and \S\S2.5--2.6]{hintzvasy2017}, whose spectral hypothesis concerns the forward problem near the black-hole exterior (the notion established in (a)), every solution of the linear wave equation with smooth Cauchy data on the one-ended initial surface of \cite{hintzvasy2017} lies in $H^{1.05-0}$, in particular in $H^1_{\mathrm{loc}}$, near the part of the Cauchy horizon in its domain of dependence; for smooth data on a complete Cauchy surface the same holds near every point of the Cauchy horizon other than its bifurcation sphere. \cite[Thm.~3.1]{costafranzen2017} gives bounded energy and a $C^0\cap H^1_{\mathrm{loc}}$ extension. This is the linear analogue of the failure of the Christodoulou formulation for smooth data at this point. (d) Not established: the Einstein--Maxwell problem; the second point $Q/Q_{\max}=249/250$; a two-sided enclosure of $\beta$ (the stretch goal); any statement at other parameters, or at the bifurcation sphere of the Cauchy horizon. \emph{Bibliography.} Added \cite{costafranzen2017}.
```

**6. The §7 rows** stand as written; no wording must change. They match the fixed notes:
- the P2-5 row's single α′ ("that leaves 0 the only Ω₂₃ resonance in Im σ ≥ −α′ and puts no global pole
  on Im σ = −α′") is H7's α″;
- its scope sentences are FORMULATION's (a)–(c), and its CF17 clause is the localized one;
- "one-ended initial surface" is W3;
- the evidence columns point to files that exist, including this document, which now carries the
  follow-up.

The P2-1 parenthesis "(assembled in P2-4)" is recorded in `s5-report.md` §6.

**7. Tests** at `3f2a17aa`, nice 10:
- fast P2 suite: 90 passed, 17 skipped (36.1 s; `review_s5/followup_pytest_fast.txt`);
- `python -m p2.assembly verify`: identical (2.8 s; `review_s5/followup_assembly_verify.txt`).

No slow suite was run, since no code changed.

**Correction to §8.** Its last sentence should read: the user approves the merge into main, and the
coordinator runs it. The user approved this merge on 2026-10-26.

**Compute and tools.** About 45 s of CPU (≈ 0.01 core-h of the 0.25 allowed), one process at a time,
nice 10. No tool call was denied.
