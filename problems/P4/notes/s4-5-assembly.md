# P4 S4-5 — Theorem C: the assembly (proof skeleton for the final review)

Status 2026-10-17, branch `p4/theorem-c-s4-3`, written at HEAD `c0f0effd`.  This is
a draft for the final adversarial review of the chain.  It adds no computation:
the only operations run for it were sha256 hashes, JSON reads and exact
rational arithmetic on stored endpoints (§8).  **Theorem C is not proved yet.**
The strip exclusion lacks the 30 centre pieces (x ≤ −4.5), so
`strip_exclusion_certified = False` in the current merge.  CLAIMS.md,
FORMULATION.md and every other note are unchanged; §6 proposes ledger text.
`P4/` means `problems/P4/`; other paths are relative to the repository root.

Tags used below: **[M]** a machine-checked ball/interval certificate;
**[X]** an exact symbolic identity checked by machine (fmpq_mpoly, ℚ(V₀,√3));
**[W]** written prose only; **[F]** float evidence only.

## 0. Read first: what this assembly found

**G1 (real gap; it affects only the E-form of the statement).** FORMULATION §2
(`P4/FORMULATION.md:60-61`), theorem-c-plan §1 (`P4/notes/theorem-c-plan.md:11-28`)
and the P4-8 ledger row ("hence E(κ) ≠ 0 on {Re κ ≥ 6}") state Theorem C
through the matching function E.  E is defined and certified only on R
(Theorem B′, `P4/notes/s2-theorem-b.md:525-540`: "Let E: R → ℂ").  The
equivalence E = 0 ⟺ eigenvalue is concluded "for every κ ∈ R" (`:581`).  It
rests on two box certificates of §3.5, (S) and (C) (`:455-481`), that cover
R only.  In particular the order-0 sonic solve (the 2×2 system for
(W_p, V_p)(0) with A_p(0) = 1, `P4/src/p4/validated/linsonic4.py:143-150`) must
be nonsingular for p̃, and hence E, to exist at κ.  The note asserts that it
is "singular iff σ(κ) = 0" (`s2-theorem-b.md:340-342`), but that assertion is
certified only by ball solves on the boxes of R.  The plan's "extends to Ω
verbatim" (`theorem-c-plan.md:24-28`, repeated in `P4/notes/s4-energy.md:58-60`)
has no written or machine proof outside R.  The determinant of that solve,
det M₀(κ), is affine in κ with real coefficients (only the ℓ·G row depends on
κ, `linsonic4.py:143`), so it vanishes at most once, at a real κ.  The box
solves show det M₀ ≠ 0 on R.  So the possible defect is confined to at most
one real point κ* ∈ (15, ∞) ⊂ Ω, where E (normalised by A_p(0) = 1) would not
be defined.  The note's "singular iff σ(κ) = 0", i.e. det M₀ ∝ σ(κ), would put
κ* at −0.0990 instead, but that is not proved.
*Scope.* The energy arguments exclude **eigenvalues**, not zeros of E, and the
mode count needs only that (§4.3).  So the eigenvalue form of Theorem C and
the unique-unstable-mode corollary do not depend on G1.  To close G1, either
certify det M₀(κ) = const·σ(κ) (one exact identity at the closed-form sonic
point, or a ball enclosure of the root) and write the classical Frobenius
existence of p̃, r₁, r₂ for every κ, or state Theorem C and P4-8 in
eigenvalue form only.

**B1 (bookkeeping; checked and harmless here, but not machine-checked by the
merge).** The merge's cover check (`P4/scripts/s43_parts.py:1427-1452`) is a
float check with tolerance 1e-12 on the plan's float labels (x, h).  In exact
rational arithmetic those labels do **not** tile: 371 joins between
consecutive steps of the same tube leave gaps of 8e-20 to 1.8e-16
(label rounding).  The certified steps are the arb steps of the tube caches,
and those tile exactly.  For this note I checked all four caches: every step
position has radius 0 and x_{k+1} = x_k − h_k holds exactly (construction:
`lintube.py:94,104`, `tmint.py:325`).  The layer engine uses that arb x_k
(`tc_layer.py:758`), and the driver matches the float label to the cache step
(`s44_layer.py:231`).  So the real cover has no hole.  The merge should check
this itself (§5 item 9).

**Known open items** (not new): the centre (30 pieces, no route decided); the
re-check of the review's fixes (running separately; the substitution is
pending it); the prose lemmas of the strip argument (identity, boundary terms,
endgame, conjugation), which the earlier review left for this stage
(`docs/review-2026-10-17-p4-layer-engine-adversarial.md:255-264`); and the
Lemma D condition on P4-8.  The full list, with the inconsistencies found
between notes, code and artifacts, is §5.

Apart from G1 and the centre, I found no hole in the chain.  §4 re-derives the
strip identity and the boundary terms to support that judgement, but this is
one reader's check, not a review.

## 1. The statement

### 1.1 Setting

* **Background.** The certified EC solution of Theorem A (ledger P4-5): the
  sonic velocity V₀* in the box c* ± 1e−16, c* = 0.1124394013880983 (A3 gives
  V₀* − c ∈ 6.3e−15 ± 7.5e−17, `s2-theorem-b.md:162-163`); the centre balls
  (a*, μ*) = a3 midpoint + Krawczyk offset (`P4/src/p4/validated/ec_centre.py`,
  `modecount.certified_centre`); the A1 sonic series, A2 centre series and A3
  tube.  The background is real.
* **System.** The 4D linearised KHA system, sonic-point gauge, time dependence
  e^{κs}: P(u)p′ = [DQ(u) − Ψ(u, u′) − κP_s(u)]p, p = (A_p, N_p, W_p, V_p),
  x ∈ (−∞, 0), with the sonic point at x = 0 and the centre at x → −∞, t = e^x
  (`s2-theorem-b.md:541-549`; `theorem-c-plan.md:34-42`).
* **Class (as in P4-6 and P4-8).** κ is an *eigenvalue* if a nonzero solution p
  satisfies all four conditions: (i) p is analytic in x at x = 0; (ii) p lies
  on the linearised constraint surface Σ, c(x, p) = 0 with c = (κ − A)S A_p −
  A(C̃_N N_p + C̃_W W_p + C̃_V V_p); (iii) p is regular at the centre: analytic
  in t with the scaled variables (Â_p, ñ_p, w̃_p, ṽ_p) bounded, so the
  point-mass direction ∝ t^{−3} is excluded; (iv) N_p(0) = 0 (sonic-point
  gauge).  `s2-theorem-b.md:541-549`, `s4-energy.md:30-33`.
* **E.** E(κ) = e^{3x_d} det[r₁(κ), r₂(κ), p̃(x_d; κ)]_{(Â,ñ,ṽ)}, x_d = −3, with
  p̃ the analytic sonic solution normalised by A_p(0) = 1, N_p(0) = 0 on Σ, and
  r₁, r₂ the regular centre family (`s2-theorem-b.md:525-530`).
* R = [0, 15] × [−14, 14] (closed), and Ω := {Re κ ≥ 0} \ R.  Write κ = ε + iτ.

### 1.2 Theorem C (proposed wording)

**Theorem C (eigenvalue form).** In the class of §1.1 there is no eigenvalue
in Ω.

**Theorem C (E-form), conditional on G1.** If E is defined on Ω (that is, the
order-0 sonic solve is nonsingular there), then E(κ) ≠ 0 for every κ ∈ Ω.

**Corollary (the mode count; Theorem B′ + Theorem C, eigenvalue form).** The
eigenvalues with Re κ ≥ 0 are exactly two simple zeros of E:
* κ̄ ∈ [0.35569919907, 0.35569920835], the pure-gauge mode (κ̄ = 2 − A₀ +
  2ω₀/3);
* κ₁ ∈ [2.81055254439, 2.81055255326].

So κ₁ is the unique non-gauge mode with Re κ ≥ 0, and
γ = 1/κ₁ ∈ [0.3558019218, 0.3558019231].  This is FORMULATION's target
(`FORMULATION.md:56-58`: "conditional on Theorem C") and the plan's
(`theorem-c-plan.md:19-28`).

### 1.3 Not claimed

Kink modes (non-analytic at the sonic point), non-spherical modes, nonlinear
codimension-one stability, Re κ < 0, and the gauge dependence of κ̄
(`FORMULATION.md:73-81`, `s2-theorem-b.md:546-549`).

## 2. Region decomposition

| part of Ω | covered by | constants | boundary |
|---|---|---|---|
| {Re κ > 15} | S4-2 (P4-8): no eigenvalue on the **closed** half-plane {Re κ ≥ c₀} | c₀ = 6 (`tc_energy` `c0`, `tc_commutator.C0_LEMMA_D = 6`), pole bound RBAR = 2 < c₀ (`tc_energy.py:53`) | all τ, both signs: the S4-2 argument needs no conjugation (`s4-energy.md:51-52`) |
| {6 ≤ Re κ ≤ 15, \|Im κ\| > 14} | S4-2 | as above | Re κ = 6 included |
| {0 ≤ Re κ ≤ 6, Im κ > 14} | the strip of S4-3/S4-4, certified on ε ∈ [0, 6], τ ∈ [14, ∞) | b₁ = 14, c = −1/32, m = 1, EPS_MAX = 6 (`tc_commutator.py:61-64`); 12 ε-cells with edges 0, 0.03, 0.07, 0.12, 0.2, 0.35, 0.6, 1.0, 1.6, 2.4, 3.4, 4.6, 6.0 (`tc_commutator.py:1103-1105`; every layer_t14 part records them) | ε = 0, ε = 6 and τ = 14 all included (closed cells; τ_min = 14.0 in every part) |
| {0 ≤ Re κ ≤ 6, Im κ < −14} | conjugation of the row above (§4.1 step S13) | — | — |

The union covers Ω.  A point of Ω has either Re κ > 15, or 0 ≤ Re κ ≤ 15 with
|Im κ| > 14.  The first kind and the second kind with Re κ ≥ 6 are in S4-2's
half-plane; the second kind with Re κ ≤ 6 is in the strip or its conjugate.

Boundaries, one by one:
* **Re κ = 6** is covered twice (S4-2's closed half-plane; the strip's last
  cell [4.6, 6.0]).  Nothing needs an open boundary here.
* **|Im κ| = 14** belongs to R (closed), so Ω needs only |Im κ| > 14.  The
  strip includes τ = 14, which overlaps R's top edge.  That overlap is
  consistent with Theorem B′: E has no zero on ∂R (min |E| on ∂R = 4.81e−3,
  `s2-theorem-b.md` §3.4).
* **The imaginary axis Re κ = 0 with |Im κ| > 14** is in Ω and covered by the
  strip's first cell [0, 0.03], which starts at 0.0 exactly.  The function-space
  facts used there need only Re σ(κ) < 0, and σ(κ) = −σ₀ − (1 − σ₀)κ with
  σ₀ = 0.090105970507920 ± 5e−16 gives Re σ(iτ) = −σ₀ < 0
  (`theorem-c-plan.md:51-60`; `s4-commutator.md:124-126`).
* **Re κ = 15** is R's right edge; Ω contains only Re κ > 15 there, plus the
  points with Re κ = 15 and |Im κ| > 14.  Both are inside S4-2's Re κ ≥ 6.
  S4-2 also covers [6, 15] × [−14, 14] ⊂ R.  That is consistent with B′,
  whose two zeros have Re κ < 3.
* **Near κ = 0.** Ω is bounded away from 0: every point of Ω is at distance
  more than 14.  E's pole sits at σ = 0, i.e. κ = −σ₀/(1−σ₀) = −0.0990, and
  the resonances at Re κ ≤ −1.198; all lie in Re κ < 0.
* **Apparent singularities of the reduction** (`s4-reduction.md:80-92`).
  * κ = A(x) and κ = −F_N(x) lie on real κ ∈ [κ̄, 1.89] ⊂ R.
  * The roots of the Bc-cubic (Bcr) are movable and complex.  They are
    excluded pointwise: in S4-2 by condition (d) (Re ≤ 2 < 6), and in the
    strip by Bcr ≠ 0 (§4.1 S3).
  * Ω ∩ ℝ = (15, ∞) lies in S4-2's region.
* **Im κ < 0 (conjugation).** The strip certificates cover τ ≥ 14 only, and
  the conjugate half rests on step S13 of §4.1: the background is real, so
  p ↦ p̄ maps eigenfunctions at κ to eigenfunctions at κ̄.  This is used
  but only asserted: `s4-commutator.md:101-102`, `tc_commutator.py:36-37`
  ("real coefficients pair kappa with conj kappa"), and for E,
  `s2-theorem-b.md:535`, `theorem-c-plan.md:17`.  §4.1 S13 gives a
  three-line written proof; it is routine but must be reviewed.  S4-2 does not
  use conjugation.
* **No route-(c) bridge is needed.** theorem-c-plan §6 planned a contour bridge
  in case c₀ > 15 or b₁ > 14 (`theorem-c-plan.md:360-364`).  The achieved
  constants are c₀ = 6 ≤ 15 and b₁ = 14 ≤ 14, so nothing between the stages and
  R is left over.

## 3. Inputs

Hashes are sha256 of the files at HEAD `c0f0effd` (`shasum -a 256`).  An
"aggregate" hash covers a directory of parts; its command is in §8.

**I1. Background (Theorem A; ledger P4-4, P4-5, both CERTIFIED; Theorem A
reviewed 2026-08-26).**
* *Used:* the V₀ box, the certified centre balls and the A1/A2/A3 enclosures
  (§1.1).
* *Artifacts:* `P4/results/a3_midpoint.json` `7bcb804f4956d440626505b25960a9ef4f6a5fdd34cbdb293982cff80fd43adb`
  (midpoint and `K_minus_m`; the certified balls are built from both,
  `ec_centre.py`).
* *Tube caches,* pinned by content hash in their records, with the record
  file's sha256 alongside:

  | record | record sha256 | content hash | covers |
  |---|---|---|---|
  | `P4/results/theorem_c/parts/caches/full.json` | `9a2bcc535cb2642d2d815467f8e6e1c9bab86c593fd8230a297b97d5e523059b` | `cb2531921be5ff5c387f2f1b424f351c7fa378f80f71c97dc7edd464de0b1b6e` | 344 steps, [−8, −0.05]; A3 tube to −3, then the certified centre series |
  | `.../caches/gap.json` | `d16bb0b245efa7265b3a6c1754d9637907d1d94305f01487af347b95998f0e58` | `a5428c96c3513e8f95cb0f523350b0cf19925de1be4fc4c6fc83991fe1b9e182` | 107 steps, [−0.055, −0.0085] |
  | `.../caches/deep.json` | `28c42cd72a44b29ca81fef35712d316be949f0b79c2eb411819fbed075b13d2a` | `ef544d927ef12534ed38a4f04044a4bbb6d61f5a45f0f7c5e20ec00775cfa493` | 129 steps, [−4.56, −2.95]; built from the certified centre series |
  | `.../caches/band.json` | `c4a55cd3d34a0e52bc2a90540c3336c9868cf8d1f33ceb3722ee2fc122436d1e` | `9f1e73819af48bcc63bd40e5100c807479ffc600df0522285d3fd3d1a0a4ebcb` | 60 steps, [−0.17, −0.05] |

  The cache files themselves (`P4/data/cache/tube_*.json`) are gitignored.  I
  recomputed all four content hashes from the files on disk for this note, and
  they equal the records.
* *Code:* `lintube.py` `8c034ee8…`, `tmint.py` `8dc2d5e8…`, `centre.py`
  `9f2655d5…`, `tailbound.py` `8c4b5de6…`, `ec_centre.py` `d6183108…`.  The
  centre-series evaluations raise if e^x exceeds the certified radius
  (`tailbound.py:69-71`); this guards the deep cache, which starts at
  x = −2.95.

**I2. Theorem B′ (P4-6, CERTIFIED; adversarial review 2026-08-30).**
* *Used:* E is analytic on a neighbourhood of R; it has exactly two zeros in
  R, both simple (κ̄ and κ₁); E(κ) = 0 ⟺ κ is an eigenvalue, **for κ ∈ R**
  (`s2-theorem-b.md:525-540, 560-581`).
* *Artifacts:*
  * `P4/results/theorem_b/winding_R.json` `d8977be1482ce448166afd9f27068b814bf2ddc87a48102960589fffe6cfe688`
  * `P4/results/theorem_b/analyticity_R.json` `d435c46cd1a182c425ade858e2b86d8c6c0e3274d3061daac72ee143cea478a5`
  * `P4/results/theorem_b/krawczyk.json` `7e4b02aeb3b08e8276c14b2ee42f6b434b0032230ec328b08493eaf05cf7a42a`
  * fixed-code re-runs (after the 2026-10-17 ball fix):
    `P4/results/recert_2026-10-17/{winding,analyticity,krawczyk}.json`
    (`bbd30770…`, `ee9d4237…`, `18f983e6…`), compared by certified fields as
    identical in `P4/results/recert_2026-10-17/comparison.json`
    `608a86a84b425b42bc34171fd687fc21b4f35f595a45934bec559ad21c49f7d4`.
* *Code:* `modecount.py` `356435ea…`, `analyticity.py` `34ef72d9…`,
  `krawczyk_kappa.py` `7a445915…`, `linsonic4.py` `32bfb104…`, `linsys.py`
  `fc6c8ff9…`.

**I3. The exact reduction (S4-1; supports the P4-6 caveat closure and P4-8).**
* *Used:* [X] the gauge vector g = (A′, N′ + κN, W′, V′) solves the 4D system
  for every κ (`tc_reduce.gauge_residual`).
* [X] On Σ the pair χ = ℓ(p), η = m(p) closes, with D2 = S²Δ̃²·dn·e
  (`tc_reduce.reduce_pair`); ℓ(g) = m(g) = c(g) = 0.
* [X] The atom factorizations, checked at every import of `tc_energy.Atoms`
  as exact Dp-power divisions (`tc_energy.py:79-121`): SΔ̃ = Dp·SDr,
  D2 = Dp³·D2r, Bc = Dp²·Bcr, TN1 = Dp⁵·TN1r, U1 = Dp⁶·U6, QNm = Dp¹⁰·QNr.
* Hence A0 = −TN1r/dT and B0 = −U6/dB, with dT = Dp·SDr·D2r·Bcr and
  dB = Dp³·SDr·D2r²·Bcr (`tc_commutator.build_pairs`, `tc_commutator.py:231-262`).
* [W] ker(ℓ, m)|_Σ = span g wherever dn·e ≠ 0 (`s4-reduction.md:44-48`).
* *Code:* `tc_reduce.py` `fb07a983…`, `tc_energy.py` `b5c08649…`; tests in
  `P4/tests/test_tc_reduction.py`.

**I4. P4-3 monotonicity (CERTIFIED).**
* *Used:* 𝔴 = −D > 0 on the open cone, so SΔ̃ ≠ 0 there and P is invertible
  (`s4-energy.md:55-56`).  In the strip, dT ≠ 0 gives the same per piece.
* *Artifacts:* `P4/results/theorem_c/tc_monotone.json` `2b08581f26ea609ec76f6656b2255732a99335eacc6d6d871237d9f63696f222`,
  `P4/results/theorem_c/parts/p43_centre.json` `2d128e2ae7e131d3b3d484bc05f37259d94dbf90da620937f18a2d95eaf77f22`;
  fixed-code re-runs `recert_2026-10-17/tc_monotone.json`, `p43_centre.json`,
  identical.
* *Code:* `tc_monotone.py` `def35430…`.

**I5. S4-2 direct energy, Re κ ≥ 6 (P4-8: CERTIFIED conditional on Lemma D).**
* *Used:* no eigenvalue with Re κ ≥ 6.  Conditions (a)–(d) of
  `s4-energy.md:48-49` hold pointwise on x ∈ (−∞, 0); they are certified on
  [−4.5, 0] and, as Lemma D, on x ≤ −4.5.
* *Artifacts:*
  * `P4/results/theorem_c/tc_energy.json` `7f48a49d2c89d792a97e1f2e2e17e9c072884b05fe8b19b83d1928683df01245`.
    Cite **only** its tube, sonic and closing blocks (14520 + 16666 + 1 boxes;
    min (b) 2.91e−4, min (c) 2.249 on the tube; `s4-energy.md:145-149`).
    **Its `lemma_d` block is the 2026-08-31 midpoint-ball record and must not
    be cited** (`s4-3-recovery.md:107-117, 437-439`).
  * Lemma D: `P4/results/theorem_c/parts/lemma_d_certified.json` `41fda13b765a0db01a9e3599501332c9006b29b911b1019f5d2f3d73f47de9ad`
    (ok, 1100/1100 subranges, the first [0, ·] so y → 0 is included;
    q2 ≥ 1.961836442292799, mB ≥ 3.4752428204250765, mC ≥ 15.348292102988353,
    pole gaps ≥ 0.99952, 0.99999).  Replays:
    `parts/replays/lemma_d_certified.json` `e56317674d128d3980e03c2c1db5f5863d2b7c35f34f1a31b3071f1b7dea8992`
    and `recert_2026-10-17/lemma_d.json` `d75c9dd0aa22a0cec8b0cdd2b9383a5ea60ca241ec28ac620e792e4434400762`.
  * Fixed-code re-run of the whole certificate:
    `recert_2026-10-17/tc_energy.json` `04fde3d6795c1c2c60f89c85f98588851a45e1d19bec85be9da6e1c578983727`.
    It is identical except `lemma_d.mins`, which there equal the certified
    Lemma D.
* *Code:* `tc_energy.py`; `tc_commutator.run_lemma_d` (`tc_commutator.py:1511-1624`;
  it refuses unless E2's pmin + 2 ≥ 16).
* *Status:* the prose of `s4-energy.md` §1–2 is not yet reviewed as routine,
  and the Lemma D condition is to be dropped by this review (P4-8 row on main).

**I6. Exact checks E1/E2/E3/E5 (S4-3).**
* *Used:* E1/E3/E5 justify the sonic Laurent raises; E3 also justifies the
  (tt·dT_k)′ raise; E2 (pmin = 14, so val ≥ 16) justifies Lemma D.
* *Artifact:* `P4/results/theorem_c/parts/checks.json` `1034c0ab0bb78f6b204b8ee9e5c9a95f86f763a3cbc4a754f4eb641ade14c92b`,
  valid for the current code (I checked every recorded source hash).  Results:
  E1 true, E2 [14, true], E3 and E5 all true.

**I7. Strip, pointwise certificate on 546 of the 576 pieces (S4-4 layer engine
under the substitution rule).**
* *Statement of each part* (`scope.statement`): for all ξ in the piece,
  ε ∈ [0, 6] and τ ≥ 14, dT, dB, nB ≠ 0 and Σ = S ξ²/τ² ≤ b_P < 0.
* *Coverage:*
  * sonic: 51 pieces, ξ ∈ [0, 0.0091004], i.e. x ∈ [−0.0091, 0];
  * gap: 107 pieces, [−0.055, −0.0085];
  * tube: 241 pieces, [−3.05, −0.05] without the four h = 0.0075 steps in
    [−0.1690625, −0.1390625];
  * band: 18 pieces, [−0.17, −0.134];
  * deep: 129 pieces, [−4.56, −2.95].
* *Worst bound per region* (merge detail): sonic −0.4517 (s29), gap −0.6031
  (gap-000), band −3.951 (band-000), tube −0.2094 (full-244), deep −0.02165
  (deep-128).  Every piece-norm part has relative margin ≥ 0.5000.
* *Artifacts:*
  * `P4/results/theorem_c/parts/layer_t14/` (546 files) aggregate
    `e9835a41fb2dbcc718e6beff9d6aa89b3b50b7874c996baee10fa173d0f5ff6c`.  I
    checked all 546 for this note: every recorded source hash is current,
    `git_dirty_src` is false, the driver hash is `3a86e50f…` (= current
    `s44_layer.py`), the τ-scope is [14, inf] and the ε-scope [0, 6].
  * The merge,
    `P4/results/theorem_c/layer_recompute/merge_substitute_layer_t14.json` `89c529a45a1600db25db3b63c381a74fef92e358e994db4c0f0480236b80f854`
    (head `fda58897`): substituted 546, missing the 30 centre pieces, failed 0,
    manual review 0, rejected 0, blocked 0, `checks_valid_for_code` true;
    uncovered x for τ ≤ 224 and τ > 224: (−∞, −4.56).  `git log fda58897..HEAD`
    shows no change to `P4/src`, `P4/scripts` or `P4/tests`.  The merge's
    recorded code hashes equal the current files.
  * The plan, `P4/results/theorem_c/parts/plan.json` `b3e658249452ee976b244d35c3944573980eecd4aa3aa53b7c70e28639fb10e2`.
  * Supporting files (not cited for the theorem):
    * `parts/layer/` (τ ≥ 224, 546 files) aggregate
      `974f44e1091c9504fb8a7a3917df8ea9d8d2cd3a770d45f7add09165a23cf3ce`, a
      second certificate of τ ≥ 224;
    * `layer_recompute/comparison.json` `fc4dc70c…`,
      `isolate_ballfix.json` `1f5c1c31…`, `snapshot.json` `2d8b2b94…`.
* *Code:*
  * `tc_layer.py` `29a864dd171e02af9719a762ebe531065b525192c619367b8bb7a19c6536a405`
  * `tc_layer_mv.py` `b1986b574c67c949924e9453bfaf56fc05864cad930a497beac90aa95763a5d1`
  * `tc_tm.py` `2db115ae9d729a97e0efe33b5958be2c5252c42d1fcb2c9abf6e631ebf72a422`
  * `tc_commutator.py` `11cfdc4969577218472dcc4c79673629c78b8686e5cfb758e8f9e199cf006da7`
    (`build_pairs`, `commute`, `sonic_raise_inputs`)
  * `ballutil.py` `39528e526fe96a5e5f5d9b44cb6579b4928176207b8a8a5f25270180b51072f1`
  * drivers: `P4/scripts/s44_layer.py` `3a86e50f1723e44873d230136a08f9a9ac4f6077ec71c220ac40849bcd823c1e`
    and `P4/scripts/s43_parts.py` `c323dad1d131a8e06531c5b8d9ed96c0fae8e4e5e5fa71205655045e127ed2e4`.
* *Review:*
  `docs/review-2026-10-17-p4-layer-engine-adversarial.md` `dd452d10b66748cb623e90ac1149a979f7a7b47e927ee4eebdb258b6574ac045`
  found the substitution sound and the Σ enclosure correct.  Its must-fix
  items are reported done in `s4-4-layer-engine.md` §14; **the re-check of
  those fixes is running separately, so treat the substitution as pending.**
* *Proposed level:* part of P4-9 (§6).

**I8. Strip, the 30 centre pieces: missing.**
* *Domain:* y = e^x ∈ (0, e^{−4.5}] (x ≤ −4.5).  The pieces run from c00 =
  [0, 1.337e−8] (x ≤ −18.13) to c29 = [0.006943, 0.011109]; the overlap with
  deep is [−4.56, −4.5].
* *Plan call:* `centre_piece(at, ce, lo, hi, 14, eps_cells_default(), 320)`,
  the S4-3 path (direct phase on τ ∈ [14, 224], τ-quadratic phase on
  τ ≥ 224).  No committed-code run has finished: one ran ≥ 25.3 CPU-h without
  finishing (`s4-3-recovery.md:133`).
* *Requirements for a centre certificate,* from §4.1:
  * for all y ∈ (0, e^{−4.5}], ε ∈ [0, 6] and τ ≥ 14: S < 0;
  * dT, dB and B0 (or nB) nonzero, or the S4-3 per-piece bounds that imply it.
    In the S4-3 path, the direct phase rejects boxes where |dT|, |dB| or |B0|
    is not bounded away from 0 (`tc_commutator.py:1397, 1427`), and the
    τ-quadratic phase uses thU < 1, Rbar < 14 and the real-pole bound
    (`tc_commutator.py:909-945`);
  * the closing piece must contain y = 0;
  * all of it on the current code with the certified centre balls.

**I9. The written lemmas of the strip argument**
(`s4-commutator.md:99-137`; `tc_commutator.py:20-38`; endgame
`s4-energy.md:51-60`).  Not yet reviewed as routine; §4.1 lists them as S5–S9,
S12 and S13.

**I10. The E-bridge outside R.** See G1.

## 4. The logical chain

### 4.1 The strip: fix κ = ε + iτ with ε ∈ [0, 6] and τ ≥ 14

Suppose p ≠ 0 is an eigenfunction at κ.

* **S1 [X]** g solves the system at κ (I3).  Its N-component at x is
  N(F_N + κ), which is nonzero for nonreal κ (N > 0).  So g(x) ≠ 0
  everywhere, and g_N(0) = N₀(κ − κ̄), using κ̄ = −F_N(0).
* **S2 [X+W] Reduction.** χ = ℓ(p) = Δ̃[(N′+κN)V_p − V′N_p] and η = m(p).  On Σ,
  D2·χ′ = Ac χ + Bc η and D2·η′ = Ae χ + Be η (I3).  χ is analytic at x = 0
  with χ(0) = 0 (Δ̃ ∝ D, and D(u₀) = 0 exactly).  At the centre χ, η are
  analytic in t and bounded (`s4-energy.md:64-75`).
* **S3 [M per piece] Nonvanishing.** Each of the 546 parts certifies dT ≠ 0,
  dB ≠ 0 and nB ≠ 0 on its piece × strip.  On sonic pieces this is certified
  for the Laurent-raised pivot forms, which is equivalent for ξ > 0
  (`s4-4-layer-engine.md:156-201`).  By the exact
  factorizations of I3 and the review's §1 (`review…adversarial.md:78-92`),
  dT = Dp·SDr·D2r·Bcr ≠ 0 implies all of the following:
  * D2 ≠ 0, so the pair closes;
  * Bc ≠ 0, so b = Bc/D2 ≠ 0;
  * D2r = c_κ-free·dn·(e/Dp) ≠ 0, so dn ≠ 0 and e ≠ 0;
  * SΔ̃ ≠ 0, so P is invertible.

  With nB, dB ≠ 0, B0 = nB/dB is finite and nonzero.  The factorization
  D2r = c·dn·(e/Dp) was checked in the review's scratch copy and is not a
  stored artifact; the Dp-power factorizations are re-checked at every
  import.  The κ-free tops dT_top and dB_top must also be nonzero, since
  tt = TN_top/dT_top and u2 = nB_top/dB_top.  That is enforced per piece: a
  zero-containing hull gives a non-finite Σ and the box is rejected
  (`s4-4-layer-engine.md:237-238`; review 2e).
* **S4 [X] The scalar equation.** Where b ≠ 0, eliminating η gives
  χ″ + A0χ′ + B0χ = 0 with A0 = nT/dT = −T and B0 = nB/dB = −U
  (`s4-reduction.md:52-58`; `s4-commutator.md:88-90`).  By S3 this holds on
  every certified piece.
* **S5 [W] Commutation.** Where B0 ≠ 0, χ₁ := χ′ satisfies χ₁″ + A1χ₁′ + B1χ₁ = 0
  with A1 = A0 − (log B0)′ and B1 = B0 + A0′ − A0(log B0)′.  (Differentiate
  the χ equation and eliminate χ = −(χ″ + A0χ′)/B0.)
* **S6 [W] The identity** (`s4-commutator.md:104-110`; "toy-verified numerically
  to 5e-10", `:116-117`).  Hand check for this note: put u = χ₁ and
  g = exp∫ρ with ρ = Re A1 − c·tt, where tt is the real κ-lead of T.
  Multiplying by gū and integrating over [x₁, x₂] gives
  [gūu′] − ∫g|u′|² + ∫g(c·tt + i Im A1)ūu′ + ∫gB1|u|² = 0.
  Take the real part minus (τ/c) times the imaginary part, with a = Re ūu′ =
  (|u|²)′/2, b = Im ūu′, ν = Im A1 + τ·tt and f = (c + τ²/c)tt − (τ/c)ν.
  The coefficient of a is f and that of b is −ν, so

  [BT] − ∫g|u′|² − ∫gνb + ∫gK|u|² = 0,
  K = Re B1 − (τ/c) Im B1 − (f′ + ρf)/2,
  BT = Re[(1 + iτ/c) g ū u′] + g f |u|²/2.

  This agrees with the note and with `tc_commutator.py:26-31`.  The step
  ∫gfa = [gf|u|²/2] − ∫g(f′ + ρf)|u|²/2 uses g′ = ρg.
* **S7 [W] Young, in a form that needs only the boundary limits.**
  |νb| ≤ |u′|² + ν²|u|²/4 gives ∫_{x₁}^{x₂} g S |u|² ≥ BT(x₁) − BT(x₂), with
  S = K + ν²/4.  If S < 0 on (−∞, 0), then
  ∫_{x₁}^{x₂} g|S||u|² ≤ BT(x₂) − BT(x₁).  The left side does not decrease as
  the interval grows, so if BT → 0 at both ends then u ≡ 0.  No absolute
  integrability is needed beyond that.
* **S8 [W] Boundary terms** (`s4-commutator.md:119-124` gives two sentences;
  draft lemma here):
  * *Sonic end, x → 0⁻.* u = χ₁ is analytic at 0 (class (i)), so u, u′ are
    bounded.
    * Res₀ A0 = −(1 + σ(κ)): xT → (1 − σ₀)(1 − κ), certified as
      τ₀ + τ₁κ ± 2.3e−5 in the S4-2 closing box (`s4-energy.md:66-70`).
    * Res₀(log B0)′ = −2: B0 = (1 + σ(κ))/x² + O(κ²/|x|).  The certified
      sonic u-hats are u1 = −0.909894/ξ² and u0 = +0.909894/ξ²
      (`s4-commutator.md:152-155`, resting on E1/E5), and 1 + σ(κ) ≠ 0 for
      κ ≠ 1.
    * So Res₀ A1 = 1 − σ(κ).  Also Res₀ tt = −(1 − σ₀) (sonic hat
      tt ≈ 0.91/ξ = (1 − σ₀)/ξ, `s4-commutator.md:208-211`).
    * Hence ρ = γ/x + O(1) and g ≍ |x|^γ with
      γ = 1 + σ₀ + (1 − σ₀)(ε + c) ≥ 1.0901 − 0.9099/32 = 1.0617 for ε ≥ 0.
      This is the note's exponent 2m − (1−σ₀)(1−ε+|c|).
    * f = O(1/|x|) because tt and Im A1 have simple poles, so
      BT = O(|x|^{γ−1}) → 0.  The condition γ > 1 is the note's
      (1−σ₀)(1+|c|) < 1, i.e. |c| < σ₀/(1−σ₀) = 0.0990, and it holds with
      margin 0.06 against ball widths of 2.3e−5.
  * *Centre end, t = e^x → 0.*
    * T → −3: the χ-exponents are {0, −3} (`s4-energy.md:75-77`), and
      𝒬 → 9/4 is backed exactly by E2.
    * B0 = −U → 0 and B0 is meromorphic in t, so B0 = βt^k(1 + O(t)) with
      β ≠ 0 and some integer k ≥ 1.  (The note's g ~ e^{x} is the case k = 2;
      September floats saw |B0| ~ e^{2x}, `s4-commutator.md:385-388`.)
    * Hence A1 → 3 − k and B1 → −3k.  Also tt = −(μ₊ + μ₋) → 0, because the
      κ-terms carry positive t-weights at the centre (‖C̃‖ ~ 1.40 e^x; float
      profile μ± → ±1.4011 t; `theorem-c-plan.md:62-65, 72-74`).  So ρ → 3 − k
      with O(t) corrections, and g ≍ t^{3−k}.
    * The regular χ = c₀ + O(t) gives χ₁ = θχ = O(t^k), where θ = t d/dt and
      the coefficients of t¹ … t^{k−1} vanish by the recursion j(j+3)c_j = 0.
      Likewise χ₁′ = O(t^k), and f → 0 since ν → Im(3 − k) + τ·0 = 0.
    * So BT = O(t^{3+k}) → 0.  The point-mass branch χ ~ t^{−3} is excluded
      by class (iii).
* **S9 [M per piece] The pointwise condition.** Every part certifies
  Σ = Sξ²/τ² ≤ b_P < 0, hence S < 0 at every ξ > 0 of its piece and every
  strip κ.
  * The layer form of Σ is an exact rearrangement of S (review §3 item 3:
    sympy residual 0; scratch, not stored).
  * The pivot forms, box tiling, both pseudo-boxes, zero-denominator rejection
    and the mv3 mean-value enclosure are reviewed in the review's §1 rows
    2a–3c.
  * The evaluation-ball flaw (2f) and the ι-centre (3d) are fixed: all 546
    parts were recomputed on the fixed code (`s4-4-layer-engine.md:812-863`).
* **S10 [M + exact cache check, see B1] The x-cover.**
  * The union of the piece intervals is [−4.56, 0], up to the closed
    endpoint 0, which is the singular point itself.
  * Cross-cache contacts overlap by at least 3.1e−5, far above label
    rounding.  The smallest are tube/full-002 with gap/gap-105 (3.1e−5) and
    gap/gap-006 with sonic/s50 (3.8e−5).
  * One exception: tube/full-244 meets deep/deep-003 at x = −3.0 with a
    2.7e−16 overlap, below label rounding.  It is redundant, because
    full-243 and full-244 are exactly contiguous steps of one cache.
  * Consecutive steps of one cache join exactly (B1).
  * Sonic pieces share their float endpoints exactly (50 touching joins), and
    the Laurent domain of a sonic piece is [p.lo, p.hi] itself.
  * **The rest of the line, x < −4.56, is covered by nothing yet.** It is the
    30 centre pieces' job.  **Lemma D does not cover it:** Lemma D certifies
    S4-2's conditions (a)–(d) for Re κ ≥ 6 (I5), a different inequality on a
    different κ-region.  The strip's deep tail has no substitute for the
    centre pieces.
* **S11 [M] The τ and ε ranges.**
  * For the 546 substituted pieces, one cover per part spans
    τ ∈ [14, ∞) × ε ∈ [0, 6].  τ = ∞ enters through the ι = 1/τ and w = 1/κ
    pseudo-boxes; the TSTAR coupling plays no role.
  * The merge's strict rule checks that τ_min = 14 equals the scope, that the
    ε-cells tile [0, 6], the ℓ-range [0, h] or [lo, hi], the cache hash, the
    code hashes and the driver hash (`s43_parts.py:1485-1605, 1634-1681`).
  * For the centre (S4-3 path), the direct phase covers τ ∈ [14, 224]
    (TSTAR_FAC·b₁ = 16·14, `tc_commutator.py:1308`) and the τ-quadratic phase
    covers τ ≥ 224.  Both include τ = 224, and the same 12 ε-cells are used.
* **S12 [W] The endgame** (`s4-energy.md:51-60`; S4-3 reuses it,
  `s4-commutator.md:342-352`).
  * S7–S9 give χ₁ ≡ 0, so χ is constant; B0χ = 0 then gives χ ≡ 0.
  * b ≠ 0 gives η ≡ 0, so ℓ(p) = m(p) = c(p) = 0.
  * ker(ℓ, m)|_Σ = span g (dn, e ≠ 0), so p = λ(x)g(x) with λ = p_N/g_N
    analytic.
  * P(λ′g) = 0 and P is invertible, so λ′ = 0.
  * N_p(0) = 0 gives λN₀(κ − κ̄) = 0; since κ ≠ κ̄, p = 0, a contradiction.
* **S13 [W] Conjugation, for τ ≤ −14.** P, DQ − Ψ and P_s are real on the real
  background, so p̄ solves the system at κ̄.  Each class condition is
  preserved:
  * (i) and (iii): conjugating a convergent power series in x or in t gives
    another;
  * (ii): c_κ̄(x, p̄) = conj c_κ(x, p), since A, S and C̃ are real;
  * (iv): N_p̄(0) = conj N_p(0) = 0.

  So κ̄ is an eigenvalue whenever κ is, and S1–S12 at κ̄ (τ ≥ 14) exclude
  τ ≤ −14.

### 4.2 Re κ ≥ 6 (S4-2, ledger P4-8)

* **E1 [X]** The normal form v″ = 𝒬v, v = χ/f with 2f′/f = T.  𝒬 is exactly
  rational in κ (QNm = Dp¹⁰·QNr).  The field identity QN_red(u₀) ≡ 0 holds
  in ℚ(V₀, √3).  The κ²-density is q₂ = (μ₊ − μ₋)²/4 (`s4-energy.md:11-28, 86-94`).
* **E2 [M]** Conditions (a) q₂ > 0, (b) 2c₀q₂ + q₁ − B₁ > 0,
  (c) c₀²q₂ + c₀q₁ + q₀ − B₀ > 0 and (d) all κ-poles have Re ≤ 2 hold on
  [−4.5, 0] (tube, sonic ladder and closing box; I5) and on x ≤ −4.5
  (Lemma D, I5).
* **E3 [W]** The identity ∫|v′|² + ∫𝒬|v|² = 0 with vanishing boundary terms
  (`s4-energy.md` §2): at the sonic end |v̄v′| = O(|x|^{−Re σ}); at the centre
  T∞ = −3 and v, v′ decay exponentially.
* **E4 [W]** (b) gives Im 𝒬/τ > 0, hence τ = 0; then (b) and (c) give
  𝒬(x; ε) > 0, hence v ≡ 0 and χ ≡ 0.  Condition (d) makes b ≠ 0, and the
  endgame is as in S12, with P invertible by I4 (`s4-energy.md:51-58`).

### 4.3 Assembly

* **A1 [W]** By §2, Ω ⊂ {Re κ ≥ 6} ∪ strip ∪ conj(strip).
* **A2** Theorem C (eigenvalue form) follows from §4.2 on {Re κ ≥ 6} and §4.1
  on the strip and its conjugate.  **Today it is conditional on I8 (the
  centre).**
* **A3 [M+X+W]** Theorem B′: the eigenvalues in R are exactly κ̄ and κ₁
  (E = 0 ⟺ eigenvalue on R).
* **A4** The mode-count corollary of §1.2 follows from A2 and A3.
* **A5 (E-form)** On Ω: no eigenvalue, plus E = 0 ⟹ eigenvalue, gives E ≠ 0.
  The implication needs only §3.8's κ-independent ingredients
  (`s2-theorem-b.md:560-579`: Σ invariant; p̃, r₁, r₂ ∈ Σ; e_W ∉ Σ(x_d)), but
  it needs E to be **defined** on Ω.  That is G1.

### 4.4 What is machine-checked and what is prose

| step | kind | where |
|---|---|---|
| background, Theorem B′ zeros and equivalence on R | M + X | I1, I2 |
| gauge identity, pair closure, factorizations | X | I3 |
| S4-2 pointwise conditions on (−∞, 0] | M | I5 |
| strip: S < 0 and dT, dB, nB ≠ 0 on x ∈ [−4.56, 0) | M | I7 |
| strip: the same on x ≤ −4.5 | **missing** | I8 |
| cover of the x-line by pieces | M (float, merge) + exact contiguity (this note, B1) | S10 |
| kernel lemma, elimination of η, commutation | W (routine linear algebra and ODE) | S2, S4, S5 |
| χ₁ energy identity | W (re-derived here; toy-checked to 5e-10) | S6 |
| Young / limit argument | W | S7 |
| boundary terms, both ends, both stages | W (uses certified residues and exact centre exponents) | S8, E3 |
| endgame | W | S12, E4 |
| conjugation | W | S13 |
| E-form outside R | **gap G1** | A5 |

## 5. Open items

1. **Centre, 30 pieces, missing.** See `merge_substitute_layer_t14.json`
   (`missing`, and `tau_le_tstar_uncovered_x = [(−∞, −4.56)]`) and
   `s4-4-layer-engine.md:733-741, 865-876`.  No route is decided.
   * *Requirements:* I8.
   * *Risk for the committed default (direct [14, 224] + τ-quadratic ≥ 224),*
     a prediction, not run.  The c29 diagnostic at tstar = 14
     (`P4/results/theorem_c/parts/experiments/centre-tstar_fac-1/c29_tau_quad_diag.json`;
     `s4-3-recovery.md:506-528`) shows K2 = O(y²) against a y-independent
     τ-linear slack R1 ≈ 0.441, with K0 + R0 ≈ −5.36.  On the first five
     ε-cells the vertex of K2τ² + R1τ lies at τ ≈ 374, 321, 265, 239 and 232,
     all above 224, and the certified bounds are +77.0, +65.5, +53.1, +47.4
     and +45.8.
   * If R1 does not shrink when the κ-boxes start at 224, those five cells
     fail the τ ≥ 224 phase on c29.  It gets worse toward y → 0, the same
     layer-type obstruction as the sonic band (`s4-commutator.md:277-299`).
   * `analyze_piece` runs the direct phase of each cell even after its
     τ-quadratic check fails (`tc_commutator.py:1022-1036`).  A failing centre
     launch would therefore spend the ≥ 25 CPU-h before reporting.
   * Recommended before any launch: a τ-quadratic-only diagnostic at
     tstar = 224 on c29 and c00.  The alternative is a centre layer engine
     (pivot forms in y, y² normalization; `s4-4-layer-engine.md:733-741`).
2. **The re-check of the review's fixes** is running separately.  The
   substitution, and so all of I7, is pending it.
3. **Prose: the χ₁ identity** (S6).  `s4-commutator.md:99-117` and
   `tc_commutator.py:20-38` give no derivation beyond "toy-verified".  §4.1 S6
   re-derives it by hand.  Recommend an exact symbolic check (e.g. sympy
   residual 0) stored as an artifact, like the review's Σ check.
4. **Prose: the boundary terms** (S8).  `s4-commutator.md:119-124` is two
   sentences, and its "against the L²-normalized decay" (`:120`) is not an
   argument.  §4.1 S8 is a draft lemma for review.  The S4-2 boundary terms
   (`s4-energy.md:62-79`) are of the same kind.
5. **Prose: the endgame** (S12, E4), written once (`s4-energy.md:51-60`) and
   reused for the strip (`s4-commutator.md:342-352`).  The review found its
   nonvanishing inputs sound (`review…adversarial.md:78-96`) but did not judge
   the prose.
6. **Prose: conjugation** (S13).  It is only asserted
   (`s4-commutator.md:101-102`, `tc_commutator.py:36-37`,
   `s2-theorem-b.md:535`).
7. **G1: the E-form outside R** (§0).  It also affects the E-clause of P4-8 for
   κ ∉ R ("hence E(κ) ≠ 0 on {Re κ ≥ 6}"; inside R ∩ {Re κ ≥ 6} the clause
   follows from B′).
   * `s2-theorem-b.md:560` is headed "(exact, all κ)", but its conclusion
     (`:581`) is for κ ∈ R only.
   * `theorem-c-plan.md:24-28` claims the extension "verbatim".
   * `theorem-c-plan.md:24` calls E "entire-in-κ … away from the pole", i.e.
     meromorphic, on Re κ > −1.19.  No certificate supports that outside R.
8. **P4-8's Lemma D condition.**
   * The certified artifact exists (I5), but the ledger keeps the condition
     until this review.
   * `P4/notes/s4-energy.md` §4–§5 still read "NOT yet certified" and
     "conditional" (`:122, 151-156`); the correction lives only in
     `s4-commutator.md:67-84` and `s4-3-recovery.md` §10–§11.
   * `tc_energy.json` keeps the midpoint `lemma_d` block (not replaced;
     `s4-4-layer-engine.md:920-929`).  Any citation of that file must exclude
     the block.
9. **The merge's cover check is float-level (B1).**
   * `uncovered(..., tol=1e-12)` runs on plan labels (`s43_parts.py:1435-1452`),
     and `x_interval` uses float x − h (`:1427-1432`).
   * Recommend an exact check against the caches: positions with radius 0,
     x_{k+1} = x_k − h_k, and cross-cache overlaps computed in arb.
10. **Merge output wording.**
    * `scope.tau_above_tstar` reads "OPEN on x in [(−∞, −4.56)] (direct phase
      only there)" (`s43_parts.py:1845-1846`).  Neither phase covers x < −4.56.
    * `scope.xi_open_max = Infinity` (`:1769`) is meaningless when the gap is
      the centre.
    * `scope.full_strip_tau = [14, 224]` (`:1841`) labels the direct range,
      but the substituted pieces cover [14, ∞).
    * `scope.tau_gt_tstar_by` counts only S4-3 parts and reads 0/0/0 here.
    * All of these are cosmetic; the counts and gaps are right.
11. **ε-cell count.** `s4-commutator.md:166` says "13 graded cells".  The code
    and every part use 12 cells with 13 edges (`tc_commutator.py:1103-1105`).
    `run_centre(..., ncells=6)` passes `ncells` to `eps_cells_default`, which
    ignores it (`tc_commutator.py:1103, 1130, 1141`).
12. **Superseded drafts.** The P4-9 draft at `s4-commutator.md:324-331` (piece
    counts and the "~0.048" open band) is superseded.  Do not reuse it; §6
    replaces it.
13. **Tests on the S4-3 path** (relevant if the centre runs on it).
    * `test_commutation_value_consistency` and `test_kdiv_roundtrip` fail as
      tightness assertions; their balls overlap (`s4-3-recovery.md:98-105,
      496-498`).
    * `test_centre_first_piece` has never completed (`:86-90`).
14. **Unchecked prose ingredient of S2.** "χ, η = O(1), analytic in t at the
    centre" (`s4-energy.md:72-75`) is asserted from the covector weights.  It
    is routine but should be in the written proof.
15. **Policy.** README defines CERTIFIED as a certificate "committed under
    `certificates/`" (main `README.md:1178`; layout at `:1163`).  P4 keeps all
    its artifacts under `P4/results/`, and P4 has no `certificates/`
    directory.  This is a pre-existing convention of the P4 rows, for the
    owner to note or change.
16. **P4-5's V₀ interval.** The ledger prints "± 3×10⁻¹⁵", while every tube
    uses c* ± 1e−16.  This is consistent, because A3 gives
    V₀* ∈ c* ± 7.5e−17 (`s2-theorem-b.md:162-163`), but a reader can trip on
    it.

## 6. Proposed ledger rows (text only; CLAIMS.md is not edited)

Levels below follow main `README.md:1172-1182`.  "Routine" means the final
adversarial review accepts the written lemmas S5–S8, S12 and S13 (and, for
P4-8, `s4-energy.md` §1–2) as routine.  Under that policy the row is then
CERTIFIED; otherwise it is WRITTEN_PROOF, provided the reviewer checks a
complete written proof.

**P4-9 (new): strip exclusion (Theorem C, part 2).**
* *Claim:* In the class of P4-6 there is no eigenvalue with 0 ≤ Re κ ≤ 6
  and |Im κ| ≥ 14.
  * Method: a once-commuted (m = 1) GHJS-type energy identity for χ₁ = χ′,
    weight g = exp∫(Re A1 − c·tt), c = −1/32.  The pointwise condition
    S = K + ν²/4 < 0, together with dT, dB, nB ≠ 0, is certified on
    piece × ε ∈ [0, 6] × τ ∈ [14, ∞) for every piece of an exact cover of
    x ∈ (−∞, 0).
  * Done today: the 546 pieces on [−4.56, 0) (sonic 51, gap 107, tube 241,
    band 18, deep 129), carried by S4-4 `layer_t14` parts under
    `merge --substitute-layer-t14`.  Open: the 30 centre pieces on x ≤ −4.5.
  * Im κ ≤ −14 follows by conjugation.
* *Level now:* CONJECTURED.
  * 546 of 576 pieces are certified.
  * The centre is open.  Its only evidence is a September τ-quadratic run,
    30/30 on the midpoint background with superseded code, whose artifacts
    were lost.
* *Target:* CERTIFIED when all three hold:
  * (i) the centre pieces pass on the current code with the certified centre
    balls and the merge reports `strip_exclusion_certified = true`;
  * (ii) the re-check of the 2026-10-17 review's fixes passes;
  * (iii) the S4-5 review accepts S5–S8, S12 and S13 as routine.

  WRITTEN_PROOF if (i) and (ii) hold but the reviewer finds the written
  lemmas non-routine and checks them complete.
* *Evidence:* `problems/P4/notes/s4-commutator.md` §2–§3;
  `problems/P4/notes/s4-4-layer-engine.md` §1–§5, §12–§14;
  `problems/P4/notes/s4-5-assembly.md` §4.1;
  `problems/P4/results/theorem_c/parts/layer_t14/` (aggregate `e9835a41…`);
  `problems/P4/results/theorem_c/layer_recompute/merge_substitute_layer_t14.json`;
  `problems/P4/results/theorem_c/parts/checks.json`;
  `problems/P4/src/p4/validated/{tc_layer,tc_layer_mv,tc_tm,tc_commutator}.py`;
  `docs/review-2026-10-17-p4-layer-engine-adversarial.md`.  When the centre
  lands, add its parts and the complete `tc_commutator.json`.

**P4-10 (new): Theorem C and the mode count.**
* *Claim:* In the class of P4-6 there is no eigenvalue in
  Ω = {Re κ ≥ 0} \ R, R = [0, 15] × [−14, 14]: P4-8 covers Re κ ≥ 6 and P4-9
  covers 0 ≤ Re κ ≤ 6, |Im κ| ≥ 14.
  * With Theorem B′ (P4-6), the eigenvalues with Re κ ≥ 0 are exactly the
    pure-gauge κ̄ ∈ [0.35569919907, 0.35569920835] and
    κ₁ ∈ [2.81055254439, 2.81055255326], both simple zeros of E.
  * So κ₁ is the unique non-gauge mode with Re κ ≥ 0, and
    γ = 1/κ₁ ∈ [0.3558019218, 0.3558019231].
  * Not claimed: kink modes, non-spherical modes, nonlinear stability,
    Re κ < 0.
  * The E-form, E(κ) ≠ 0 on Ω, is to be added only once G1 is closed.
* *Level now:* CONJECTURED (it inherits P4-9).
* *Target:* CERTIFIED when P4-9 is CERTIFIED and this review drops P4-8's
  Lemma D condition, accepting `s4-energy.md` §1–2 as routine (artifact
  `parts/lemma_d_certified.json`).  If P4-8 or P4-9 ends at WRITTEN_PROOF, so
  does P4-10.
* *Evidence:* `problems/P4/notes/s4-5-assembly.md`; the P4-6, P4-8 and P4-9
  evidence; `problems/P4/results/recert_2026-10-17/comparison.json`.

**Related edits for the owner, once the above holds:**
* P4-8: drop "conditional on Lemma D"; cite `parts/lemma_d_certified.json`,
  not the `tc_energy.json` `lemma_d` block; restrict "hence E(κ) ≠ 0" to R
  until G1 is closed.
* P4-2: update its scope note from "in [0,15]×[−14,14]" to Re κ ≥ 0.
* FORMULATION v1.5 (`theorem-c-plan.md:365-370`).

## 7. Checklist for the final adversarial review

1. **Statement.** Is the class of §1.1 identical in P4-6, P4-8 and both
   S4-3/S4-4 arguments?  Is the eigenvalue form what the corollary needs?
   Decide G1: prove the extension or drop the E-form.
2. **Region algebra (§2).** Check every boundary listed, and that
   Ω ⊂ {Re κ ≥ 6} ∪ strip ∪ conj(strip) with both closed sets as certified.
3. **S6.** Re-derive the χ₁ identity; ideally store a symbolic check.  Confirm
   that `direct_S_boxes` (`tc_commutator.py:1314-1325, 1427-1447`) and the
   layer Σ (`s4-4-layer-engine.md:143-147`) evaluate the same S.
4. **S7.** Check the one-sided limit argument; it needs only BT → 0.
5. **S8.** Check the sonic exponent γ > 1 for all ε ∈ [0, 6] at c = −1/32,
   using the certified residues (xT closing-box ball, tt hat, E1/E5), and the
   centre exponents for general k ≥ 1.  Check the S4-2 analogue
   (`s4-energy.md:62-79`).
6. **S3/S12.** Re-run or inspect the exact factorizations:
   * the import-time `_dpow` checks, `tc_energy.py:79-121`;
   * D2r = c·dn·(e/Dp), which is not stored;
   * that the kernel and P-invertibility steps use only polynomials dividing
     dT up to powers of Dp;
   * that dT_top and dB_top ≠ 0 is enforced per box.
7. **Layer parts.**
   * Recompute the aggregate hash (§8).
   * Spot-recompute two or three pieces bit for bit: at least deep-128, the
     thinnest at −0.02165, and sonic/s29, the thinnest sonic.
   * Confirm every part passes `layer_checks(strict=True)` on the reviewer's
     copy.
   * Confirm the re-check of the fixes.
8. **Cover (S10, B1).** Verify exact contiguity in the four caches, the
   cross-cache overlaps, the sonic touching joins, and the overlap
   [−4.56, −4.5] with the centre once it lands.
9. **Centre.** Check:
   * its route and code version, the certified centre balls, and that the
     closing piece contains y = 0;
   * that both τ phases (or a reviewed centre engine) cover τ ≥ 14 for all
     12 cells;
   * the nonvanishing requirements of I8;
   * that the merge then reports `complete` and `strip_exclusion_certified`.
10. **Conjugation (S13).** Accept the written lemma or require more.
11. **S4-2.**
    * Lemma D from `parts/lemma_d_certified.json`: E2 pmin + 2 ≥ 16, the
      certified centre balls, the y → 0 subrange.
    * `tc_energy.json`'s tube, sonic and closing blocks, and their fixed-code
      re-run.
    * The boundary terms and the logic of (a)–(d).
12. **Provenance.** No `P4/src` or `P4/scripts` change after the parts and the
    merge (`git log fda58897..HEAD` is empty for both); `checks.json` valid
    for current code; no accepted code variants in the merge.
13. **Ledger.** Apply the decisions of §6 to P4-8, P4-9, P4-10 and P4-2, and
    the README `certificates/` policy point (§5 item 15).

## 8. Reproducing the hashes and checks in this note

From the repository (worktree) root:

```
shasum -a 256 problems/P4/results/theorem_b/*.json problems/P4/results/theorem_c/tc_*.json \
  problems/P4/results/theorem_c/parts/{checks,plan,lemma_d_certified,p43_centre}.json \
  problems/P4/results/theorem_c/parts/caches/*.json problems/P4/results/theorem_c/layer_recompute/*.json
# aggregate of a part directory (paths relative to the root, C-locale sort):
find problems/P4/results/theorem_c/parts/layer_t14 -name '*.json' | LC_ALL=C sort \
  | xargs shasum -a 256 | shasum -a 256
git log --oneline fda58897..HEAD -- problems/P4/src problems/P4/scripts problems/P4/tests   # empty
```

The exact cover analysis was done with Python `fractions` on two sources:
* `parts/plan.json` endpoints: 371 same-tube label gaps, at most 1.8e−16,
  and no cross-region gap;
* the dyadic step positions of `P4/data/cache/tube_{full,gap,deep,band}.json`
  (mantissa and exponent as stored by `Tube.save`): radius 0 throughout and
  x_{k+1} = x_k − h_k exactly in all four caches.

The cache content hashes were recomputed with `s43_parts.tube_content_hash`'s
recipe.
