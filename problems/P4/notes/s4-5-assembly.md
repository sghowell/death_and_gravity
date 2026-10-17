# P4 S4-5 — Theorem C: the assembly (proof skeleton)

Status 2026-10-17, branch `p4/theorem-c-s4-3`.
* First written at HEAD `c0f0effd`, with the re-check `82728526` folded in.
* Revised at HEAD `d68b64ff` after the independent S4-5 review
  `docs/review-2026-10-17-p4-s4-5-adversarial.md` (`7d0ec1cf`, `d68b64ff`),
  which found **no flaw** in the chain.  This revision:
  * applies the review's wording fixes (its §10, items 1–8);
  * adds Lemma G1 (§4.5), which closes G1;
  * follows the review's §11 for the class, the levels and the ledger text.

Neither version ran a computation.  The only operations were sha256 hashes,
JSON reads and exact rational arithmetic on stored endpoints (§8).  The
review's own checks are committed in `P4/certificates/theorem_c/review_s45/`, and
its README says how each one was run.

**Theorem C is not proved yet.**  The strip exclusion lacks the 30 centre
pieces (x ≤ −4.5), so `strip_exclusion_certified = False`.  CLAIMS.md and
FORMULATION.md are unchanged; §6 proposes ledger text.  `P4/` means
`problems/P4/`; other paths are relative to the repository root.

*Tags:*
* **[M]** a machine-checked ball/interval certificate;
* **[X]** an exact symbolic identity checked by machine (fmpq_mpoly, sympy,
  ℚ(V₀,√3));
* **[W]** written prose;
* **[F]** float evidence only.

*Notation:* κ̄ is the gauge zero, κ₀ the real root of det M₀(κ) (Lemma G1;
σ's zero lies in the same ball), and z^c the complex conjugate of z.

## 0. Read first

* **Review verdict** (`docs/review-2026-10-17-p4-s4-5-adversarial.md` §1).
  * No flaw.  S4-2 and S4-3/S4-4 are sound for every κ they claim.
  * The review accepts the written lemmas S5–S8, S12 and S13 (§4.1) as
    routine, given its committed checks.
  * Its five gaps are all closed by those checks: the S8 sonic residues, the
    S8 and S4-2 centre limits, G1, the S4-2 join at −4.5, and the class
    wording.
* **G1 is closed by Lemma G1 (§4.5).**
  * *The gap.* The E-form of Theorem C needs E to be defined outside R.  It
    appears in FORMULATION.md:60-61, theorem-c-plan.md:11-28, and P4-8's
    "hence E(κ) ≠ 0 on {Re κ ≥ 6}".  But Theorem B′ certifies E's ingredients
    on R only (`s2-theorem-b.md:455-481, 525-540, 581`).
  * The review showed that the gap has two halves: the sonic order-0 solve,
    and the centre regular family (certificate (C)).  It settled both.
    * det M₀(κ) = d₁(κ − κ₀), with κ₀ = −0.0990290820550 ± 2.9e−14 real;
      σ's zero lies in the same ball.
    * det(nP₀ − G₀(κ)) is κ-free.
  * Lemma G1 combines them: on Re κ > −0.0990, E is defined and E = 0 ⟺
    eigenvalue.  So the E-form follows from the eigenvalue form on all of Ω.
* **B1 is resolved in code.**
  * The merge's cover check used to be a float check, with tolerance 1e-12,
    on rounded labels.
  * Since `ccc74f3a` the merge checks the cover exactly on the tube caches:
    content hash, step identity, radius 0, exact joins, Fractions.
  * The review reproduced the result with independent code
    (`P4/certificates/theorem_c/review_s45/cover_independent.txt`; review §8).
* **Open** (full list in §5):
  * the centre: 30 pieces, x ≤ −4.5, no route decided;
  * the merge-output wording, `scope.xi_open_max` and `scope.full_strip_tau`
    (review §8, "cosmetic leftovers").  This is code, to be done with the
    centre work;
  * the README policy point, under which CERTIFIED artifacts belong under
    `certificates/`.  This is an owner decision.

## 1. The statement

### 1.1 Setting

* **Background.** The certified EC solution of Theorem A (ledger P4-5):
  * the sonic velocity V₀* in the box c* ± 1e−16, c* = 0.1124394013880983
    (A3 gives V₀* − c ∈ 6.3e−15 ± 7.5e−17, `s2-theorem-b.md:162-163`);
  * the centre balls (a*, μ*) = a3 midpoint + Krawczyk offset
    (`P4/src/p4/validated/ec_centre.py`, `modecount.certified_centre`);
  * the A1 sonic series, the A2 centre series and the A3 tube.

  The background is real.
* **System.** The 4D linearised KHA system in the sonic-point gauge, with time
  dependence e^{κs}:
  P(u)p′ = [DQ(u) − Ψ(u, u′) − κP_s(u)]p, p = (A_p, N_p, W_p, V_p),
  x ∈ (−∞, 0).  The sonic point is x = 0 and the centre is x → −∞, t = e^x
  (`s2-theorem-b.md:541-549`; `theorem-c-plan.md:34-42`).
* **Class 𝒞** (explicit; review §11).  𝒞 consists of the nonzero solutions p
  of the system above, in the gauge N_p(0) = 0, that are:
  * (i) analytic at the sonic point x = 0;
  * (ii) on the linearised momentum-constraint surface Σ: c(x, p) = 0, with
    c = (κ − A)S A_p − A(C̃_N N_p + C̃_W W_p + C̃_V V_p);
  * (iii) regular at the centre: analytic in t = e^x with the scaled variables
    (Â_p, ñ_p, w̃_p, ṽ_p) bounded, so the point-mass direction ∝ t^{−3} is
    excluded.

  κ is an *eigenvalue* if 𝒞 is nonempty at κ.  No normalisation and no radius
  enter; this is the class S4-2 uses (`s4-energy.md:30-33`).
* **Relation to P4-6** (`s2-theorem-b.md:541-549`).
  * By Lemma G1 (i), every element of 𝒞 is a multiple of p̃.
  * P4-6's class (A_p(0) = 1, radius ≥ 0.065) is contained in 𝒞.  So if 𝒞 has
    no eigenvalue, neither does P4-6's class.
  * On R the two classes coincide (Theorem B′).
* **E.** E(κ) = e^{3x_d} det[r₁(κ), r₂(κ), p̃(x_d; κ)]_{(Â,ñ,ṽ)}, x_d = −3.
  Here p̃ is the analytic sonic solution normalised by A_p(0) = 1, N_p(0) = 0
  on Σ, and r₁, r₂ are the regular centre family (`s2-theorem-b.md:525-530`).
  E is defined on Re κ > −0.0990 (Lemma G1).
* R = [0, 15] × [−14, 14] (closed), and Ω := {Re κ ≥ 0} \ R.  Write κ = ε + iτ.

### 1.2 Theorem C (proposed wording)

**Theorem C (eigenvalue form).** There is no eigenvalue of 𝒞 in Ω.

**Theorem C (E-form), by Lemma G1.** E(κ) ≠ 0 for every κ ∈ Ω.

**Corollary (the mode count; Theorem B′ + Theorem C).** The eigenvalues with
Re κ ≥ 0 are exactly the following, both simple zeros of E:
* κ̄ ∈ [0.35569919907, 0.35569920835], the pure-gauge mode
  (κ̄ = 2 − A₀ + 2ω₀/3);
* κ₁ ∈ [2.81055254439, 2.81055255326].

So κ₁ is the unique non-gauge mode with Re κ ≥ 0, and
γ = 1/κ₁ ∈ [0.3558019218, 0.3558019231].  This is FORMULATION's target
(`FORMULATION.md:56-58`: "conditional on Theorem C") and the plan's
(`theorem-c-plan.md:19-28`).

### 1.3 Not claimed

* kink modes (non-analytic at the sonic point);
* non-spherical modes;
* nonlinear codimension-one stability;
* Re κ < 0;
* the gauge dependence of κ̄ (`FORMULATION.md:73-81`, `s2-theorem-b.md:546-549`).

## 2. Region decomposition

| part of Ω | covered by | constants | boundary |
|---|---|---|---|
| {Re κ > 15} | S4-2 (P4-8): no eigenvalue on the **closed** half-plane {Re κ ≥ c₀} | c₀ = 6 (`tc_energy` `c0`, `tc_commutator.C0_LEMMA_D = 6`), pole bound RBAR = 2 < c₀ (`tc_energy.py:53`) | all τ, both signs: the S4-2 argument needs no conjugation (`s4-energy.md:51-52`) |
| {6 ≤ Re κ ≤ 15, \|Im κ\| > 14} | S4-2 | as above | Re κ = 6 included |
| {0 ≤ Re κ ≤ 6, Im κ > 14} | the strip of S4-3/S4-4, certified on ε ∈ [0, 6], τ ∈ [14, ∞) | b₁ = 14, c = −1/32, m = 1, EPS_MAX = 6 (`tc_commutator.py:61-64`); 12 ε-cells with edges 0, 0.03, 0.07, 0.12, 0.2, 0.35, 0.6, 1.0, 1.6, 2.4, 3.4, 4.6, 6.0 (`tc_commutator.py:1103-1105`; every layer_t14 part records them) | ε = 0, ε = 6 and τ = 14 all included (closed cells; τ_min = 14.0 in every part) |
| {0 ≤ Re κ ≤ 6, Im κ < −14} | conjugation of the row above (§4.1 S13) | — | — |

The union covers Ω.  A point of Ω has either Re κ > 15, or 0 ≤ Re κ ≤ 15 with
|Im κ| > 14.  The first kind, and the second kind with Re κ ≥ 6, lie in S4-2's
half-plane.  The second kind with Re κ ≤ 6 lies in the strip or its conjugate.
The review confirms this region algebra (§3).

Boundaries, one by one:
* **Re κ = 6** is covered twice: by S4-2's closed half-plane and by the
  strip's last cell [4.6, 6.0].  Nothing needs an open boundary here.
* **|Im κ| = 14** belongs to R (closed), so Ω needs only |Im κ| > 14.
  * The strip includes τ = 14 and so overlaps R's top edge.
  * That is consistent with Theorem B′: E has no zero on ∂R (min |E| on
    ∂R = 4.81e−3, `s2-theorem-b.md` §3.4).
* **The imaginary axis Re κ = 0 with |Im κ| > 14** is in Ω.
  * It is covered by the strip's first cell [0, 0.03], which starts at 0.0
    exactly.
  * The function-space facts used there need only Re σ(κ) < 0.  With
    σ(κ) = −σ₀ − (1 − σ₀)κ and σ₀ = 0.090105970507920 ± 5e−16,
    Re σ(iτ) = −σ₀ < 0 (`theorem-c-plan.md:51-60`; `s4-commutator.md:124-126`).
* **Re κ = 15** is R's right edge.
  * There Ω contains only Re κ > 15, plus the points with Re κ = 15 and
    |Im κ| > 14.  Both lie inside S4-2's Re κ ≥ 6.
  * S4-2 also covers [6, 15] × [−14, 14] ⊂ R.  That is consistent with B′,
    whose two zeros have Re κ < 3.
* **Near κ = 0.** Ω is bounded away from 0: every point of Ω is at distance
  more than 14.
  * The order-0 sonic solve is singular only at κ₀ = −0.0990290820550 ±
    2.9e−14, which is real, and σ's zero lies in the same ball
    (`review_s45/g1_order0.json`).  The plan places E's pole there
    (`theorem-c-plan.md:55`).
  * The resonances σ ∈ ℕ lie at Re κ ≤ −1.198.
  * All of these are in Re κ < 0.
* **Apparent singularities of the reduction** (`s4-reduction.md:80-92`;
  review §3).
  * κ = A(x) ∈ (1, 1.89] and κ = −F_N(x) ∈ [0.3549, 1) are real and inside R.
    These ranges are float; −F_N dips below κ̄ = 0.3557 near x ≈ −0.03.
  * The roots of the Bc-cubic (Bcr) are movable and complex.  They are
    excluded pointwise: in S4-2 by condition (d) (Re ≤ 2 < 6), and in the
    strip by Bcr ≠ 0 (§4.1 S3).
  * Ω ∩ ℝ = (15, ∞) lies in S4-2's region.
* **Im κ < 0 (conjugation).** The strip certificates cover τ ≥ 14 only.
  * The conjugate half is §4.1 S13: the background is real, so p ↦ p^c maps
    𝒞 at κ onto 𝒞 at κ^c.  The review accepts this as routine (§6).
  * The other notes only assert it: `s4-commutator.md:101-102` and
    `tc_commutator.py:36-37`; for E, `s2-theorem-b.md:535`, now written
    E(κ^c) = E(κ)^c.
  * S4-2 needs no conjugation.
* **No route-(c) bridge is needed.** theorem-c-plan §6 planned a contour
  bridge in case c₀ > 15 or b₁ > 14 (`theorem-c-plan.md:360-364`).  The
  achieved constants are c₀ = 6 ≤ 15 and b₁ = 14 ≤ 14.

## 3. Inputs

Hashes are sha256 of the files (`shasum -a 256`).  An "aggregate" hash covers
a directory of parts; its command is in §8.  Except where marked, the hashes
were taken at `c0f0effd`, and no cited file changed afterwards (review §9).

**I1. Background** (Theorem A; ledger P4-4 and P4-5, both CERTIFIED; Theorem A
reviewed 2026-08-26).
* *Used:* the V₀ box, the certified centre balls and the A1/A2/A3 enclosures
  (§1.1).
* *Artifacts:* `P4/results/a3_midpoint.json`
  `7bcb804f4956d440626505b25960a9ef4f6a5fdd34cbdb293982cff80fd43adb`
  (midpoint and `K_minus_m`; `ec_centre.py` builds the certified balls from
  both).
* *Tube caches,* pinned by content hash in their records:

  | record | record sha256 | content hash | covers |
  |---|---|---|---|
  | `P4/results/theorem_c/parts/caches/full.json` | `9a2bcc535cb2642d2d815467f8e6e1c9bab86c593fd8230a297b97d5e523059b` | `cb2531921be5ff5c387f2f1b424f351c7fa378f80f71c97dc7edd464de0b1b6e` | 344 steps, [−8, −0.05]; A3 tube to −3, then the certified centre series |
  | `.../caches/gap.json` | `d16bb0b245efa7265b3a6c1754d9637907d1d94305f01487af347b95998f0e58` | `a5428c96c3513e8f95cb0f523350b0cf19925de1be4fc4c6fc83991fe1b9e182` | 107 steps, [−0.055, −0.0085] |
  | `.../caches/deep.json` | `28c42cd72a44b29ca81fef35712d316be949f0b79c2eb411819fbed075b13d2a` | `ef544d927ef12534ed38a4f04044a4bbb6d61f5a45f0f7c5e20ec00775cfa493` | 129 steps, [−4.56, −2.95]; built from the certified centre series |
  | `.../caches/band.json` | `c4a55cd3d34a0e52bc2a90540c3336c9868cf8d1f33ceb3722ee2fc122436d1e` | `9f1e73819af48bcc63bd40e5100c807479ffc600df0522285d3fd3d1a0a4ebcb` | 60 steps, [−0.17, −0.05] |

  The cache files themselves (`P4/data/cache/tube_*.json`) are gitignored.
  Their content hashes are re-derived by the merge's exact cover check.
* *Code:* `lintube.py` `8c034ee8…`, `tmint.py` `8dc2d5e8…`, `centre.py`
  `9f2655d5…`, `tailbound.py` `8c4b5de6…`, `ec_centre.py` `d6183108…`.
  Centre-series evaluations raise if e^x exceeds the certified radius
  (`tailbound.py:69-71`); this guards the deep cache, which starts at
  x = −2.95.

**I2. Theorem B′** (P4-6, CERTIFIED; adversarial review 2026-08-30).
* *Used:* E is analytic on a neighbourhood of R; it has exactly two zeros in
  R, both simple (κ̄ and κ₁); and E(κ) = 0 ⟺ κ is an eigenvalue, for κ ∈ R
  (`s2-theorem-b.md:525-540, 560-581`).  Outside R, Lemma G1 extends the
  definition and the equivalence.
* *Artifacts:*
  * `P4/certificates/theorem_b/winding_R.json` `d8977be1482ce448166afd9f27068b814bf2ddc87a48102960589fffe6cfe688`
  * `P4/certificates/theorem_b/analyticity_R.json` `d435c46cd1a182c425ade858e2b86d8c6c0e3274d3061daac72ee143cea478a5`
  * `P4/certificates/theorem_b/krawczyk.json` `7e4b02aeb3b08e8276c14b2ee42f6b434b0032230ec328b08493eaf05cf7a42a`
  * fixed-code re-runs `P4/certificates/recert_2026-10-17/{winding,analyticity,krawczyk}.json`
    (`bbd30770…`, `ee9d4237…`, `18f983e6…`).  They are compared as identical
    by certified fields in `P4/certificates/recert_2026-10-17/comparison.json`
    `608a86a84b425b42bc34171fd687fc21b4f35f595a45934bec559ad21c49f7d4`.
* *Code:* `modecount.py` `356435ea…`, `analyticity.py` `34ef72d9…`,
  `krawczyk_kappa.py` `7a445915…`, `linsonic4.py` `32bfb104…`, `linsys.py`
  `fc6c8ff9…`.

**I3. The exact reduction** (S4-1; supports the P4-6 caveat closure and P4-8).
* [X] The gauge vector g = (A′, N′ + κN, W′, V′) solves the 4D system for
  every κ (`tc_reduce.gauge_residual`).
* [X] On Σ the pair χ = ℓ(p), η = m(p) closes, with D2 = S²Δ̃²·dn·e
  (`tc_reduce.reduce_pair`); ℓ(g) = m(g) = c(g) = 0.
* [X] The atom factorizations are checked at every import of `tc_energy.Atoms`
  as exact Dp-power divisions (`tc_energy.py:79-121`):
  SΔ̃ = Dp·SDr, D2 = Dp³·D2r, Bc = Dp²·Bcr, TN1 = Dp⁵·TN1r, U1 = Dp⁶·U6,
  QNm = Dp¹⁰·QNr.
* Hence A0 = −TN1r/dT and B0 = −U6/dB, with dT = Dp·SDr·D2r·Bcr and
  dB = Dp³·SDr·D2r²·Bcr (`tc_commutator.build_pairs`, `tc_commutator.py:231-262`).
* [X] The remaining factor structure (`review_s45/factor_checks.json`):
  * D2r = 16S⁴W²·dn·(e/Dp), with e/Dp = 4SW(Q₁ + κN) and
    Q₁ + κN = N(F_N + κ);
  * SDr = 4S²W, and SΔ̃ = Dp·SDr = det P.
* [W] The kernel lemma: ker(ℓ, m)|_Σ = span g wherever dn·e ≠ 0
  (`s4-reduction.md:44-48`; re-derived in review §5).
* [X] The elimination, the normal form, the commutation and the polynomial
  forms TN1, U1, QNm are sympy-checked (`review_s45/identities_out.txt`, I1–I4).
* *Code:* `tc_reduce.py` `fb07a983…`, `tc_energy.py` `b5c08649…`; tests in
  `P4/tests/test_tc_reduction.py`.

**I4. P4-3 monotonicity** (CERTIFIED).
* *Used:* 𝔴 = −D > 0 on the open cone, so SΔ̃ ≠ 0 there and P is invertible
  (`s4-energy.md:55-56`).  In the strip, dT ≠ 0 gives the same per piece.
* *Artifacts:*
  * `P4/certificates/theorem_c/tc_monotone.json` `2b08581f26ea609ec76f6656b2255732a99335eacc6d6d871237d9f63696f222`
  * `P4/certificates/theorem_c/p43_centre.json` `2d128e2ae7e131d3b3d484bc05f37259d94dbf90da620937f18a2d95eaf77f22`
  * fixed-code re-runs `recert_2026-10-17/tc_monotone.json` and
    `p43_centre.json`, identical.
* *Code:* `tc_monotone.py` `def35430…`.

**I5. S4-2 direct energy, Re κ ≥ 6** (P4-8; proposed CERTIFIED with the
Lemma D condition dropped, §6).
* *Used:* there is no eigenvalue with Re κ ≥ 6.
  * Conditions (a)–(d) of `s4-energy.md:48-49` hold pointwise on
    x ∈ (−∞, 0).
  * They are certified on [−4.500000000000000083, 0] by tc_energy, and on
    x ≤ ln YCUT = −4.500000000000000037 by Lemma D.
  * The two regions overlap by 4.6e−17 (`review_s45/s42_join.txt`; review §7).
* *Artifacts:*
  * `P4/certificates/theorem_c/tc_energy.json`
    `7f48a49d2c89d792a97e1f2e2e17e9c072884b05fe8b19b83d1928683df01245`.
    Cite **only** its tube, sonic and closing blocks: 14520 + 16666 + 1 boxes,
    with min (b) 2.91e−4 and min (c) 2.249 on the tube (`s4-energy.md:190-194`).
    **Its `lemma_d` block is the 2026-08-31 midpoint-ball record and must not
    be cited** (`s4-3-recovery.md:107-117, 437-439`).
  * Lemma D: `P4/certificates/theorem_c/lemma_d_certified.json`
    `41fda13b765a0db01a9e3599501332c9006b29b911b1019f5d2f3d73f47de9ad`.
    * The fixed-code record; all 24 recorded source hashes are current
      (review §7).
    * ok, 1100/1100 subranges.  The first subrange is [0, ·], so y → 0 is
      included.
    * q2 ≥ 1.961836442292799, mB ≥ 3.4752428204250765,
      mC ≥ 15.348292102988353; pole gaps ≥ 0.99952 and 0.99999.
  * Fixed-code re-runs:
    * `recert_2026-10-17/lemma_d.json`
      `d75c9dd0aa22a0cec8b0cdd2b9383a5ea60ca241ec28ac620e792e4434400762`;
    * `recert_2026-10-17/tc_energy.json`
      `04fde3d6795c1c2c60f89c85f98588851a45e1d19bec85be9da6e1c578983727`,
      identical except `lemma_d.mins`, which there equal the certified Lemma D.
  * `parts/replays/lemma_d_certified.json`
    `e56317674d128d3980e03c2c1db5f5863d2b7c35f34f1a31b3071f1b7dea8992` is a
    **historical** replay.  It records pre-ball-fix source hashes for
    `tc_commutator`, `tc_energy` and `tc_tm`; cite the certified file instead.
  * The centre limit T = −3 + O(t²), exact for every κ ∉ {0, 1}
    (`review_s45/centre_layers.json`), backs the S4-2 centre boundary term.
* *Code:* `tc_energy.py`; `tc_commutator.run_lemma_d`
  (`tc_commutator.py:1511-1624`; it refuses unless E2's pmin + 2 ≥ 16).

**I6. Exact checks E1/E2/E3/E5** (S4-3).
* *Used:* E1/E3/E5 justify the sonic Laurent raises; E3 also justifies the
  (tt·dT_k)′ raise; E2 (pmin = 14, so val ≥ 16) justifies Lemma D.
* *Artifact:* `P4/certificates/theorem_c/checks.json`
  `1034c0ab0bb78f6b204b8ee9e5c9a95f86f763a3cbc4a754f4eb641ade14c92b`, valid
  for the current code.  Results: E1 true, E2 [14, true], E3 and E5 all true.

**I7. Strip, pointwise certificate on 546 of the 576 pieces** (the S4-4 layer
engine under the substitution rule).
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
  * `P4/results/theorem_c/parts/layer_t14/` (546 files), aggregate
    `e9835a41fb2dbcc718e6beff9d6aa89b3b50b7874c996baee10fa173d0f5ff6c`.
    * Every recorded source hash is current, `git_dirty_src` is false and the
      driver is `3a86e50f…` (= `s44_layer.py`).
    * τ-scope [14, inf], ε-scope [0, 6].
    * All 546 also pass the strict rule after `a0a1b5ff`: each records a git
      head and a clean tree (review §8).
  * The merge, `P4/results/theorem_c/layer_recompute/merge_substitute_layer_t14.json`
    `89c529a45a1600db25db3b63c381a74fef92e358e994db4c0f0480236b80f854`
    (head `fda58897`).
    * Substituted 546, missing the 30 centre pieces, failed 0, manual
      review 0, rejected 0, blocked 0; `checks_valid_for_code` true.
    * Uncovered x, for τ ≤ 224 and for τ > 224: (−∞, −4.56).
    * The review re-ran the merge with the exact cover code (`75869b91`) and
      got identical per-piece detail (review §8).
  * The plan, `P4/results/theorem_c/parts/plan.json`
    `b3e658249452ee976b244d35c3944573980eecd4aa3aa53b7c70e28639fb10e2`.
  * Supporting files, not cited for the theorem:
    * `parts/layer/` (τ ≥ 224, 546 files), aggregate
      `974f44e1091c9504fb8a7a3917df8ea9d8d2cd3a770d45f7add09165a23cf3ce`;
    * `layer_recompute/comparison.json` `fc4dc70c…`,
      `isolate_ballfix.json` `1f5c1c31…`, `snapshot.json` `2d8b2b94…`.
* *Code:*
  * `tc_layer.py` `29a864dd171e02af9719a762ebe531065b525192c619367b8bb7a19c6536a405`
  * `tc_layer_mv.py` `b1986b574c67c949924e9453bfaf56fc05864cad930a497beac90aa95763a5d1`
  * `tc_tm.py` `2db115ae9d729a97e0efe33b5958be2c5252c42d1fcb2c9abf6e631ebf72a422`
  * `tc_commutator.py` `11cfdc4969577218472dcc4c79673629c78b8686e5cfb758e8f9e199cf006da7`
    (`build_pairs`, `commute`, `sonic_raise_inputs`)
  * `ballutil.py` `39528e526fe96a5e5f5d9b44cb6579b4928176207b8a8a5f25270180b51072f1`
  * drivers: `P4/scripts/s44_layer.py`
    `3a86e50f1723e44873d230136a08f9a9ac4f6077ec71c220ac40849bcd823c1e`, and
    `P4/scripts/s43_parts.py`
    `6e95af680036f8b1d38de807b89a97f3709304aa9fab9248f505ad4632637fa0`
    (at `d68b64ff`, after `a0a1b5ff` and `ccc74f3a`).
* *Reviews:* `docs/review-2026-10-17-p4-layer-engine-adversarial.md`
  `9155516367185c17a35e2f82fb3138e1527767271dd362daf3230b9b61e81bd8`.
  * It found the substitution sound and the Σ enclosure correct.
  * Its fixes are recorded in `s4-4-layer-engine.md` §14.
  * The re-check (§6 of that document, `82728526`) confirmed the fixes.  It
    recomputed parts bit for bit, sampled 11,445 independent points with none
    above a bound, and judged that the 546-piece statement "merits CERTIFIED".
* *Level:* the pointwise statement on the 546 pieces is CERTIFIED by the
  reviews.  As an input it is part of P4-9 (§6).

**I8. Strip, the 30 centre pieces: missing.**
* *Domain:* y = e^x ∈ (0, YCUT], i.e. x ≤ −4.5.
  * The pieces run from c00 = [0, 1.337e−8] (x ≤ −18.13) to
    c29 = [0.006943, 0.011109].
  * c29 reaches x = −4.5 − 3.7e−17, so the overlap with the deep end at
    −4.56 is 0.06 (`review_s45/cover_independent.txt`).
* *Plan call:* `centre_piece(at, ce, lo, hi, 14, eps_cells_default(), 320)`,
  the S4-3 path: the direct phase on τ ∈ [14, 224] and the τ-quadratic phase
  on τ ≥ 224.  No committed-code run has finished; one ran ≥ 25.3 CPU-h
  without finishing (`s4-3-recovery.md:133`).
* *Requirements for a centre certificate,* from §4.1:
  * S < 0 for all y ∈ (0, YCUT], ε ∈ [0, 6] and τ ≥ 14;
  * dT, dB and B0 (or nB) nonzero, or the S4-3 per-piece bounds that imply it.
    In the S4-3 path, the direct phase rejects boxes where |dT|, |dB| or |B0|
    is not bounded away from 0 (`tc_commutator.py:1397, 1427`), and the
    τ-quadratic phase uses thU < 1, Rbar < 14 and the real-pole bound
    (`tc_commutator.py:909-945`);
  * the closing piece must contain y = 0;
  * everything on the current code, with the certified centre balls, written
    through the driver with the provenance of `a0a1b5ff`.
* *Pre-flight* (review §8): nothing structural stops a closing piece at y = 0
  from certifying dT, dB, nB ≠ 0.
  * The lowest centre layers of dT and dB vanish only at κ ∈ {0, 1}
    (exact, `review_s45/centre_layers.json`).
  * U6's t¹⁰ layer has no κ-root in the strip
    ([F], `review_s45/centre_U6_roots.json`).

**I9. The written lemmas of the strip argument** (`s4-commutator.md:99-137`;
`tc_commutator.py:20-38`; endgame `s4-energy.md:51-60`).  These are §4.1 S5–S8,
S12 and S13.  The S4-5 review accepts them as routine, given its committed
checks (review §1, rows 3a–6).

**I10. E outside R:** Lemma G1 (§4.5).

**I11. The S4-5 review and its checks.**
* `docs/review-2026-10-17-p4-s4-5-adversarial.md`
  `939e15912a32ae961316295aa8ff40f93c45de6ef5b7208ddd57e8b5f0145f8f`.
* Outputs in `P4/certificates/theorem_c/review_s45/` (the README there gives each
  command):

  | file | sha256 | what |
  |---|---|---|
  | `identities_out.txt` | `c14549fe98da0e7c1c7c1eb1ef76312f23b9e9a64dc7ff9dd7089456253e8b65` | sympy I1–I7: elimination, normal form, T/U/𝒬 polynomial forms, commutation, the χ₁ identity, the layer normalization; residuals 0, controls nonzero |
  | `g1_order0.json` | `904b16d5b7d478d1b458221d3af863bc71936a3708ed3fa8dc7d5e2f317354dc` | det M₀ = d₀ + d₁κ and its root κ₀, σ and its zero, E rank one, det D |
  | `centre_pattern.json` | `5ab1fdaec034a0398d82f71eba2b0da816ca30dfab775e06f323c2dce669265a` | exact zero pattern of P₀, G₀(κ) at t = 0 |
  | `g1_centre.txt` | `1085cd9a75883da985ed59b604f4883825e2cf1165ef2383bb92b197ee7f4ec4` | P₀, G₀(κ) balls; det(nP₀ − G₀(κ)), n = 1…6, at five κ |
  | `centre_layers.json` | `fe07ce151a0e0b9ab33d17a824379aa26a1dd4081fb14a11f49dcd735931c627` | exact centre layers of TN1r, denT, U6, denU, tt |
  | `centre_U6_roots.json` | `5b74e9d9e500e4ea04a87d00fffae690cc826ec086a54e04aee06c3774e556e9` | κ-roots of U6's lowest centre layer [F] |
  | `sonic_residues.json` | `b9dfc3d0318ecb0637b7e11ccd0b1925132c52f5b81d4c697e9df61e240ff875` | R_T, R_B identities, Res tt and γ at x = 0 (balls) |
  | `factor_checks.json` | `9c38452915da96d91bc7116e167a80560901eeb6f0a0306bd052b49eded1607c` | exact factor structure |
  | `s42_join.txt` | `199907362f0edb79a0ed7d5a5b8fbc378e0af7ddf8c8c1b072f68af46bf22d9e` | the S4-2 join at −4.5 |
  | `cover_independent.txt` | `06c34efe7afb9814aa6c87cf2b8f4870d608e00d26c8b3ac791437ca59242ab2` | the strip x-cover, independent code |

## 4. The logical chain

### 4.1 The strip: fix κ = ε + iτ with ε ∈ [0, 6] and τ ≥ 14

Suppose p ∈ 𝒞 at κ.

* **S1 [X]** g solves the system at κ (I3).
  * Its N-component at x is N(F_N + κ), which is nonzero for nonreal κ
    (N > 0).  So g(x) ≠ 0 everywhere.
  * g_N(0) = N₀(κ − κ̄), using κ̄ = −F_N(0).
* **S2 [X+W] Reduction.**
  * χ = ℓ(p) = Δ̃[(N′+κN)V_p − V′N_p] and η = m(p).
  * On Σ, D2·χ′ = Ac χ + Bc η and D2·η′ = Ae χ + Be η (I3).
  * χ is analytic at x = 0 with χ(0) = 0, since Δ̃ ∝ D and D(u₀) = 0
    exactly.
  * At the centre χ and η are analytic in t and bounded (`s4-energy.md:67-82`;
    review §2 (d)).
* **S3 [M per piece + X] Nonvanishing.**
  * Each of the 546 parts certifies dT ≠ 0, dB ≠ 0 and nB ≠ 0 on its
    piece × strip.  On sonic pieces this is certified for the Laurent-raised
    pivot forms, which is equivalent for ξ > 0 (`s4-4-layer-engine.md:156-201`).
  * By the exact factor structure (I3; `review_s45/factor_checks.json`),
    dT = Dp·SDr·D2r·Bcr ≠ 0 implies:
    * D2 ≠ 0, so the pair closes;
    * Bc ≠ 0, so b = Bc/D2 ≠ 0;
    * D2r = 16S⁴W²·dn·(e/Dp) ≠ 0, so dn ≠ 0 and e ≠ 0, and with them
      g_N = N(F_N + κ) ≠ 0;
    * SΔ̃ = det P ≠ 0, so P is invertible.
  * With nB, dB ≠ 0, B0 = nB/dB is finite and nonzero.
  * The κ-free tops dT_top and dB_top are inverted per piece by `TM.inv`.  It
    raises unless the top is nonzero on the whole piece (review §5;
    `tc_layer.py:355, 379`).
* **S4 [X] The scalar equation.** Where b ≠ 0, eliminating η gives
  χ″ + A0χ′ + B0χ = 0, with A0 = nT/dT = −T and B0 = nB/dB = −U
  (`s4-reduction.md:52-58`; `s4-commutator.md:88-90`; sympy I1/I3 in
  `review_s45/identities_out.txt`).  By S3 this holds on every certified piece.
* **S5 [X] Commutation.** Where B0 ≠ 0, χ₁ := χ′ satisfies
  χ₁″ + A1χ₁′ + B1χ₁ = 0, with A1 = A0 − (log B0)′ and
  B1 = B0 + A0′ − A0(log B0)′.  (Differentiate the χ equation and eliminate
  χ = −(χ″ + A0χ′)/B0; sympy I4.)
* **S6 [X] The identity** (`s4-commutator.md:104-110`; `tc_commutator.py:26-31`).
  Put u = χ₁ and g = exp∫ρ, with ρ = Re A1 − c·tt and tt the real κ-lead of T.
  Multiplying by gū, integrating over [x₁, x₂], and taking the real part minus
  (τ/c) times the imaginary part gives

  [BT] − ∫g|u′|² − ∫gν Im(ūu′) + ∫gK|u|² = 0,

  with
  * ν = Im A1 + τ·tt;
  * f = (c + τ²/c)tt − (τ/c)ν;
  * K = Re B1 − (τ/c) Im B1 − (f′ + ρf)/2;
  * BT = Re[(1 + iτ/c) g ū u′] + g f |u|²/2.

  The pointwise form BT′ = g|u′|² + gν Im(ūu′) − gK|u|² is checked in sympy
  with residual 0, and flipping the sign of the (τ/c) Im B1 term gives a
  nonzero residual (I5 in `review_s45/identities_out.txt`; review §4).
* **S7 [W] Young, in a form that needs only the boundary limits.**
  * |ν Im(ūu′)| ≤ |u′|² + ν²|u|²/4 gives
    ∫_{x₁}^{x₂} g S |u|² ≥ BT(x₁) − BT(x₂), with S = K + ν²/4.
  * If S < 0 on (−∞, 0), then ∫_{x₁}^{x₂} g|S||u|² ≤ BT(x₂) − BT(x₁).
  * The left side does not decrease as the interval grows.  So if BT → 0 at
    both ends, then u ≡ 0.
* **S8 [M + X + W] Boundary terms** (written here and reviewed, review §4;
  `s4-commutator.md:119-124` gives two sentences).
  * *Sonic end, x → 0⁻.* u = χ₁ is analytic at 0 (class (i)), so u and u′
    are bounded.  Write R_T(κ) := lim xT and R_B(κ) := lim x²B0.
    * Balls at the sonic point, on the V₀ box, with
      E_D = D′(0) = 7.56057255908 ± 4.7e−12 (`review_s45/sonic_residues.json`),
      verify two κ-polynomial identities coefficientwise (every coefficient
      ball contains 0; U6(u₀) ≡ 0 by E5):
      * TN1r(u₀; κ) − (1 + σ(κ))·E_D·SDr·D2r·Bcr(u₀; κ) ≡ 0;
      * −(dU6/dx)(0; κ) − (1 + σ(κ))·E_D³·SDr·D2r²·Bcr(u₀; κ) ≡ 0.
    * The same file gives Res₀ tt = TN1r₆/(E_D SDr D2r₂ Bcr₃) =
      −0.909894029492 ± 6e−14 = −(1 − σ₀).
    * *Integer rigidity makes the two identities exact* (review §4).  For
      σ(κ) ∉ ℤ the χ-equation is Fuchsian at 0, with one exponent m_a ≥ 1
      (the analytic branch, since Δ̃(0) = 0) and one 1 + σ + j (j ≥ 0).
      Both branches are ≢ 0: otherwise η ≡ 0, and the endgame puts the
      solution in span g.  So R_T − σ = m_a + j is an integer-valued rational
      function of κ on an open set, hence constant, and the balls pin it to 1.
    * So m_a = 1 and j = 0.  The product of the exponents then gives
      R_B = 1 + σ(κ).
    * Hence R_T = R_B = 1 + σ(κ) exactly, Res₀ A1 = 2 − R_T = 1 − σ(κ), and
      Res₀ ρ = γ with
      **γ = 1 + σ₀ + (1 − σ₀)(ε + c) ≥ 1.061671782 ± 9e−11** for ε ≥ 0 and
      c = −1/32.  This is the note's exponent 2m − (1−σ₀)(1−ε+|c|).
    * So g ≍ |x|^γ.  ν is bounded, because the residues of Im A1 and τ·tt
      cancel, and f = O(1/|x|).  Hence BT = O(|x|^{γ−1}) → 0.
    * Two older citations are *not* used here.  The S4-2 closing box, with
      remainder 1.75e−5, is valid only on Re κ ≥ 6.  The u-hats of
      `s4-commutator.md:152-155` are large-κ hats, not fixed-κ residues.
  * *Centre end, t = e^x → 0.*
    * Exact centre layers (`review_s45/centre_layers.json`; fmpq_mpoly in E2's
      grading p = −e_N + 2e_W + e_V with A → 1):
      * val(TN1r) = val(denT) = 7 and val(TN1r + 3·denT) = 9;
      * val(denU) = 8 and val(U6) = 10;
      * val(TN1r₆) = 9 > val(denT₅) = 7.
    * The lowest layers of denT and denU are −12288 n⁵w⁶ κ(κ−1)⁴ and
      −786432 n¹⁰w⁹ κ(κ−1)⁶.  Both are nonzero at the centre for κ ∉ {0, 1}.
    * Hence, for every κ ∈ Ω: T = −3 + O(t²), A0 → 3, B0 = O(t²) and
      tt = O(t²).
    * Write B0 = βt^k(1 + O(t)) with β ≠ 0 and k ≥ 2.  (k = 2 in the strip by
      the float diagnostic `review_s45/centre_U6_roots.json`; any k ≥ 1
      suffices.)  Then A1 → 3 − k and B1 → −3k, and ρ → 3 − k with O(t)
      corrections, so g ≍ t^{3−k}.
    * The regular χ = c₀ + O(t^k): the coefficients of t¹ … t^{k−1} vanish by
      the recursion j(j+3)c_j = 0.  So χ₁ and χ₁′ are O(t^k), and f → 0
      since ν → 0.
    * Hence BT = O(t^{3+k}) → 0.  The point-mass branch χ ~ t^{−3} is
      excluded by class (iii).
* **S9 [M per piece] The pointwise condition.** Every part certifies
  Σ = Sξ²/τ² ≤ b_P < 0, hence S < 0 at every ξ > 0 of its piece and every
  strip κ.
  * The layer form of Σ is an exact rearrangement of S (sympy I6, committed).
  * `direct_S_boxes` and the layer engine bound the same S (review §4, "The
    code paths").
  * The pivot forms, box tiling, both pseudo-boxes, zero-denominator rejection
    and the mv3 mean-value enclosure are covered by the layer-engine review's
    §1, rows 2a–3c.
  * The evaluation-ball flaw (2f) and the ι-centre (3d) are fixed: all 546
    parts were recomputed on the fixed code (`s4-4-layer-engine.md:812-863`).
* **S10 [X] The x-cover.**
  * Since `ccc74f3a`, the merge checks the cover exactly on the tube caches
    (`exact_x_cover`, `exact_gaps`).  It checks the content hash and the step
    identity, requires radius 0 and exact joins x_{k+1} = x_k − h_k, sets
    sonic pieces to [−hi, −lo] and centre pieces to exact logs, and sweeps
    with no tolerance.
  * The review reproduced it independently (`review_s45/cover_independent.txt`).
    The union is exactly [−5134103575202365/2⁵⁰, 0] = [−4.56, 0], with no
    gap.
  * The smallest cross-region overlaps are:
    * 1.80e−16 for full-244 / deep-003.  This one is redundant, because
      full-243 and full-244 are exactly contiguous steps of one cache;
    * 3.12e−5 for full-002 / gap-105;
    * 3.79e−5 for gap-006 / s50.
  * **The rest of the line, x < −4.56, is covered by nothing yet.** It is the
    30 centre pieces' job.
  * **Lemma D does not cover it.** Lemma D certifies S4-2's conditions
    (a)–(d) for Re κ ≥ 6 (I5), a different inequality on a different κ-region.
* **S11 [M] The τ and ε ranges.**
  * For the 546 substituted pieces, one cover per part spans
    τ ∈ [14, ∞) × ε ∈ [0, 6].  τ = ∞ enters through the ι = 1/τ and w = 1/κ
    pseudo-boxes.
  * The strict rule checks that τ_min = 14 equals the scope, that the ε-cells
    tile [0, 6], and the ℓ-range [0, h] or [lo, hi], the cache hash, the code
    hashes, the driver hash and the git provenance (review §8).
  * For the centre (S4-3 path), the direct phase covers τ ∈ [14, 224]
    (TSTAR_FAC·b₁ = 16·14, `tc_commutator.py:1308`) and the τ-quadratic phase
    covers τ ≥ 224.  Both include τ = 224 and use the same 12 ε-cells.
* **S12 [W + X] The endgame** (`s4-energy.md:51-60`; reused for the strip,
  `s4-commutator.md:345-355`; re-derived in review §5).
  * S7–S9 give χ₁ ≡ 0, so χ is constant; B0χ = 0 then gives χ ≡ 0.
  * b ≠ 0 gives η ≡ 0, so ℓ(p) = m(p) = c(p) = 0.
  * ker(ℓ, m)|_Σ = span g (dn, e ≠ 0), so p = λ(x)g(x) with λ = p_N/g_N
    analytic.
  * P(λ′g) = 0 and P is invertible, so λ′ = 0.
  * N_p(0) = 0 gives λN₀(κ − κ̄) = 0.  Since κ ≠ κ̄, p = 0: a contradiction.
* **S13 [W] Conjugation, for τ ≤ −14.** P, DQ, Ψ and P_s have rational
  coefficients and the background is real, so G(κ)^c = G(κ^c), and p^c solves
  the system at κ^c.  Each class condition is preserved:
  * (i) and (iii): conjugate power series in x, or in t with real scalings;
  * (ii): c_κ(x, p)^c = c_{κ^c}(x, p^c), since A, S and C̃ are real;
  * the gauge: N_{p^c}(0) = N_p(0)^c = 0.

  So κ^c is an eigenvalue whenever κ is, and S1–S12 at κ^c (τ ≥ 14) exclude
  τ ≤ −14 (review §6).

### 4.2 Re κ ≥ 6 (S4-2, ledger P4-8)

* **E1 [X]** The normal form v″ = 𝒬v, v = χ/f with 2f′/f = T.
  * 𝒬 is exactly rational in κ (QNm = Dp¹⁰·QNr; sympy I2/I3c).
  * The field identity QN_red(u₀) ≡ 0 holds in ℚ(V₀, √3).
  * The κ²-density is q₂ = (μ₊ − μ₋)²/4 (`s4-energy.md:11-28, 98-106`).
* **E2 [M + X]** Conditions (a) q₂ > 0, (b) 2c₀q₂ + q₁ − B₁ > 0,
  (c) c₀²q₂ + c₀q₁ + q₀ − B₀ > 0, and (d) all κ-poles have Re ≤ 2.
  * They hold on [−4.500000000000000083, 0] (tube, sonic ladder and closing
    box; I5) and on x ≤ −4.500000000000000037 (Lemma D, I5).
  * The join overlaps by 4.6e−17, checked exactly (`review_s45/s42_join.txt`).
  * The other joins are exact: the double −0.05, and the double −1e−8.
* **E3 [W + M + X]** The identity ∫|v′|² + ∫𝒬|v|² = 0 with vanishing boundary
  terms (`s4-energy.md` §2; review §4).
  * At the sonic end, xT → 1 + σ(κ): the closing-box enclosure (remainder
    1.75e−5, valid on Re κ ≥ 6), made exact by the rigidity lemma of S8.
    Then |v̄v′| = O(|x|^{−Re σ}), with −Re σ ≥ 5.55.
  * At the centre, T = −3 + O(t²) exactly (S8, `centre_layers.json`).  Then
    |f| ~ e^{−3x/2}, v and v′ are O(e^{3x/2}), and [v̄v′] → 0.
* **E4 [W]**
  * (b) gives Im 𝒬/τ > 0, hence τ = 0.  Then (b) and (c) give
    𝒬(x; ε) > 0, hence v ≡ 0 and χ ≡ 0.
  * Condition (d) makes b ≠ 0 and dn·e ≠ 0.
  * The endgame is as in S12, with P invertible by I4 (`s4-energy.md:51-58`;
    review §7).

### 4.3 Assembly

* **A1 [W]** By §2, Ω ⊂ {Re κ ≥ 6} ∪ strip ∪ conj(strip).
* **A2** Theorem C (eigenvalue form) follows from §4.2 on {Re κ ≥ 6} and from
  §4.1 on the strip and its conjugate.  **Today it is conditional on I8, the
  centre.**
* **A3 [M+X+W]** Theorem B′: the eigenvalues in R are exactly κ̄ and κ₁
  (E = 0 ⟺ eigenvalue on R).
* **A4** The mode-count corollary of §1.2 follows from A2 and A3.
* **A5 (E-form)** Ω ⊂ {Re κ > −0.0990}, so by Lemma G1 (§4.5), E is defined on Ω
  and E = 0 ⟺ eigenvalue there.  A2 then gives E ≠ 0 on Ω.  This uses both
  halves of G1: the sonic order-0 solve, and the centre family (certificate
  (C) of `s2-theorem-b.md` §3.5).

### 4.4 What is machine-checked and what is prose

| step | kind | where |
|---|---|---|
| background; Theorem B′ zeros and the equivalence on R | M + X | I1, I2 |
| gauge identity, pair closure, factor structure | X | I3, I11 |
| elimination, normal form, commutation, χ₁ identity, layer normalization | X (sympy) | S4–S6, S9, I11 |
| S4-2 conditions on (−∞, 0], with the join at −4.5 | M + X | I5, E2 |
| strip: S < 0 and dT, dB, nB ≠ 0 on x ∈ [−4.56, 0) | M | I7 |
| strip: the same on x ≤ −4.5 | **missing** | I8 |
| cover of the x-line by pieces | X (merge exact cover; independent reproduction) | S10 |
| boundary terms, both ends, both stages | M + X + W (residue identities, rigidity; centre layers) | S8, E3 |
| Young / limit argument | W (routine) | S7 |
| kernel lemma and endgame | W + X (routine) | S12, E4 |
| conjugation | W (routine) | S13 |
| E outside R | M + X + W | §4.5 |

### 4.5 Lemma G1: E outside R

**Lemma G1.** For every κ with Re κ > −0.0990:
* (i) the admissible sonic family {p analytic at x = 0 : c(0) = 0, N_p(0) = 0}
  is one-dimensional, and its nonzero elements have A_p(0) ≠ 0.  So p̃ exists
  and is unique;
* (ii) the solutions regular at the centre are exactly span(r₁, r₂), a
  two-dimensional family, normalised by (w̃, ṽ)(0) = (1, 0) and (0, 1);
* (iii) E(κ) is defined, and E(κ) = 0 ⟺ κ is an eigenvalue of 𝒞.

Ω ⊂ {Re κ ≥ 0} ⊂ {Re κ > −0.0990}, so Theorem C (eigenvalue form) gives
E ≠ 0 on Ω.

By (i), every element of 𝒞 is a multiple of p̃.  P4-6's class (A_p(0) = 1,
radius ≥ 0.065) is contained in 𝒞, so if 𝒞 has no eigenvalue, neither does
P4-6's class.  On R the two classes coincide (Theorem B′).

*Proof.*

(a) **Sonic, order 0** [M] (`P4/certificates/theorem_c/review_s45/g1_order0.json`).
* The check builds the 2×2 order-0 solve M₀ for (W_p, V_p)(0) exactly as
  `linsonic4.py:142-150` does, on the certified V₀ box.
* Only its ℓ·G row depends on κ, so
  det M₀(κ) = d₀ + d₁κ = d₁(κ − κ₀), with d₁ = 11.8434668208 ± 2.8e−11 and
  κ₀ = −0.0990290820550 ± 2.9e−14 real.
* The zero of σ(κ) = σ_a + σ_bκ lies in the same ball:
  −σ_a/σ_b ∈ −0.0990290820550 ± 3.1e−14, and d₀σ_b − d₁σ_a ∈ [± 1.9e−13].
  The proportionality det M₀ ∝ σ is established only in balls.  Its
  structural reason is that C̃ annihilates the null vector of P_fl(u₀)
  (C̃·r_null ∈ [± 1.0e−14]).
* Both zeros lie below −0.0990.  Hence det M₀ ≠ 0 on Re κ > −0.0990.

(b) **Sonic, higher orders** [X + M] (same file).
* Rows 0–2 of E₀ and E₁ are exactly zero, so E(κ) has rank ≤ 1.  With
  det D = −9.60210331438 ± 3.5e−12, this gives
  det(nD + E(κ)) = det D · n³(n − σ(κ)) identically.
* σ_b = −0.909894029492 ± 6e−14 < 0, and σ's zero lies below −0.0990.  So
  Re σ(κ) < 0 on Re κ > −0.0990, and nD + E(κ) is invertible for every n ≥ 1:
  no resonance.
* σ(κ) ≠ 0 there, so E(κ) ≠ 0 and its kernel, {p₀ : ℓ·G₀(κ)p₀ = 0}, is
  3-dimensional.
* On that kernel, the conditions N_p(0) = 0 and c(0) = 0 leave the 2×2
  system M₀(κ)·(W_p, V_p)(0) = A_p(0)·b(κ) of (a), with
  b = (−ℓG₀(κ)_A, (κ − A₀)S₀).
* Since det M₀ ≠ 0, the order-0 data form a one-dimensional space,
  parametrised by A_p(0).  A_p(0) = 0 forces p₀ = 0.
* An analytic solution whose order-0 term vanishes is ≡ 0, because every
  higher level is invertible.  This proves (i).

(c) **Centre** [X + M].
* `P4/certificates/theorem_c/review_s45/centre_pattern.json` gives the exact zero
  pattern of `linscaled.full_system` at t = 0:
  * rows 1 and 3 of P₀ are (0, 1, 0, 0) and (0, 0, n∞, 0);
  * rows 1 and 3 of G₀ vanish;
  * the κ-part of G₀ is the single entry G₀[2][2] = −3κ.
* Expanding det(nP₀ − G₀(κ)) along rows 1 and 3 removes column 2, so it is
  κ-free:

  det(nP₀ − G₀(κ)) = −n·(n∞n)·(n + 3)·(P₂₃n − G₂₃),

  with n∞ = 1.2366, P₂₃ = 28.7929 and G₂₃ = −86.3787.  This is nonzero for
  every n ≥ 1.
* `P4/certificates/theorem_c/review_s45/g1_centre.txt` gives the balls on the
  certified (a*, μ*) series: −569.685, −3560.53, … for n = 1 … 6, identical at
  κ ∈ {0, 1, 7, 3 + 50i, 1000i}.
* G₀(κ) has rank 2 for every κ, through the κ-free (Â, ñ)-minor
  diag(−3, 28.2435).
* So the exponents are {0, 0, −3, G₂₃/P₂₃} for every κ, with
  G₂₃/P₂₃ < 0.  It equals −3.0000 only within balls.
* The exponent-0 solutions form a two-dimensional analytic family,
  parametrised by (w̃, ṽ)(0), whose recursion is solvable at every order
  n ≥ 1.
* The solutions with the two negative exponents are unbounded in the scaled
  variables.  This proves (ii).

(d) **Convergence** [W, classical].
* At x = 0 the system reads p′ = x⁻¹B(x)p with B = x·adj P·G/det P analytic:
  a regular singular point of the first kind.  The same holds at t = 0 in the
  scaled variables.
* Formal series with non-resonant recursions converge on the disc where B is
  analytic.  The solutions continue along the real line, where P is invertible
  (P4-3).

(e) **κ-independent ingredients** [X + M], verbatim from `s2-theorem-b.md`
§3.8:
* Σ is invariant (16SD·c′ = Λ_lin c, with Λ_lin κ-free).
* The constraint exponents at both ends, γ_lin ∉ ℕ⁺ and ρ = −3 ∉ ℕ₀, put p̃,
  r₁ and r₂ in Σ.
* e_W ∉ Σ(x_d).
* r₁(x_d) and r₂(x_d) are independent.

With (i) and (ii), the argument of `s2-theorem-b.md` §3.5 (ii) applies at every
κ with Re κ > −0.0990:

E(κ) = 0 ⟺ p̃(x_d) ∈ span(r₁(x_d), r₂(x_d)) ⟺ p̃ is regular at the centre
⟺ κ is an eigenvalue.

This proves (iii). ∎

## 5. Open items

1. **Centre, 30 pieces, missing.** See `merge_substitute_layer_t14.json`
   (`missing`; `tau_le_tstar_uncovered_x = [(−∞, −4.56)]`) and
   `s4-4-layer-engine.md:733-741, 865-876`.  No route is decided.
   * *Requirements and pre-flight:* I8.
   * *Risk for the committed default* (direct phase on [14, 224] plus the
     τ-quadratic phase on τ ≥ 224).  This is a prediction, not a run.
     * The c29 diagnostic at tstar = 14
       (`P4/results/theorem_c/parts/experiments/centre-tstar_fac-1/c29_tau_quad_diag.json`;
       `s4-3-recovery.md:506-528`) shows K2 = O(y²) against a y-independent
       τ-linear slack R1 ≈ 0.441, with K0 + R0 ≈ −5.36.
     * On the first five ε-cells, the vertex of K2τ² + R1τ lies at
       τ ≈ 374, 321, 265, 239 and 232, all above 224.  The certified bounds
       there are +77.0, +65.5, +53.1, +47.4 and +45.8.
   * If R1 does not shrink when the κ-boxes start at 224, those cells fail
     the τ ≥ 224 phase on c29.  It gets worse toward y → 0: the same
     layer-type obstruction as the sonic band (`s4-commutator.md:280-302`).
   * `analyze_piece` runs the direct phase of each cell even after its
     τ-quadratic check fails (`tc_commutator.py:1022-1036`).  So a failing
     centre launch would spend ≥ 25 CPU-h before reporting.
   * Recommended before any launch: a τ-quadratic-only diagnostic at
     tstar = 224 on c29 and c00.  The alternative is a centre layer engine,
     with pivot forms in y and a y² normalization
     (`s4-4-layer-engine.md:733-741`).
2. **Merge output wording** (code; review §8 "cosmetic leftovers", §10 item 9).
   To be done with the centre work.
   * `scope.xi_open_max = Infinity` and `scope.full_strip_tau = [14, 224]` still
     mislabel the open range and the τ range of the substituted pieces, which
     cover [14, ∞).
   * `run_centre(..., ncells=6)` still passes an ignored argument
     (`tc_commutator.py:1103, 1130, 1141`).
3. **Owner decision: the README policy.  Resolved 2026-10-17: moved.**
   * README defines CERTIFIED as a certificate "committed under
     `certificates/`" (main `README.md:1178`; layout at `:1163`).
   * The CERTIFIED rows' artifacts (P4-3, P4-6, P4-8, review checks) moved
     unchanged to `P4/certificates/` (its `README.md` maps old paths to new).
   * P4-9's artifacts stay under `results/` while the row is CONJECTURED.
4. **Optional hardening of the S4-2 join** (review §7).  The join at −4.5
   overlaps by 4.6e−17 "by rounding luck".  Running Lemma D to a ycut slightly
   above e^{−4.5}, or keeping one more tube step, would give a margin.  This is
   not a gap: the exact check is committed.
5. **Superseded drafts.** The P4-9 draft at `s4-commutator.md:327-334` (piece
   counts and the "~0.048" open band) is superseded by §6; do not reuse it.
6. **Notes outside this revision's scope.**
   * `theorem-c-plan.md:17` writes E(κ̄) = E(κ)̄ for conjugation, the same
     clash as fixed in `s2-theorem-b.md`.
   * `theorem-c-plan.md:24-28` asserts the extension to Ω, now proved by
     Lemma G1.
   * Pointers can be added when the plan is next edited.
7. **P4-5's V₀ interval.** The ledger prints "± 3×10⁻¹⁵", while every tube
   uses c* ± 1e−16.  This is consistent, since A3 gives V₀* ∈ c* ± 7.5e−17
   (`s2-theorem-b.md:162-163`).  It is informational only.

Resolved since the first version of this note:
* G1, by Lemma G1 (§4.5), with the centre half added (review §2).
* B1, the float cover check.  `ccc74f3a` made the cover exact, and the review
  reproduced it independently (§8).
* The two provenance gaps (git head outside a checkout; late module loads),
  by `a0a1b5ff`, verified end to end in review §8.
* The two S4-3 tightness tests, by `967ce2b0`: they now evaluate at one point
  ξ, and they pass (review §9).
* The ε-cell count in `s4-commutator.md`, by `75869b91`.
* The written lemmas S5–S8, S12 and S13, accepted as routine given the
  committed checks (review §1).
* The S8 citations, now replaced by the residue identities, the rigidity
  lemma and the centre layers (review §10 item 1).
* P4-8's Lemma D condition.  `s4-energy.md` §4–§5 now carry a status banner,
  and the S4-2 join is checked (`s42_join.txt`).
* The class wording, now explicit (§1.1, §6).
* The complex-conjugate notation, now z^c (S13, `s2-theorem-b.md:446, 535`).

## 6. Proposed ledger rows (text only; CLAIMS.md is not edited)

Levels follow main `README.md:1172-1182`.  The S4-5 review accepts the written
lemmas as routine, given its committed checks, so each row's target is
CERTIFIED.

**Class in all three rows** (review §11).  The nonzero solutions p of the 4D
linearised KHA system about the certified EC background, in the sonic-point
gauge N_p(0) = 0 with time dependence e^{κs}, that are:
* analytic at the sonic point;
* on the linearised momentum-constraint surface;
* regular at the centre (analytic in t = eˣ with bounded scaled variables).

By Lemma G1 (i), on Re κ > −0.0990 every element of this class is a multiple
of p̃.  P4-6's class (A_p(0) = 1, radius ≥ 0.065) is contained in this class,
so if this class has no eigenvalue, neither does P4-6's class.  On R the two
coincide (P4-6, Theorem B′).

**P4-8 (update): direct-energy exclusion (Theorem C, part 1).**
* *Claim:* in the class above there is no eigenvalue with Re κ ≥ 6.  Hence,
  by Lemma G1 (`s4-5-assembly.md` §4.5), E(κ) ≠ 0 on {Re κ ≥ 6}.
* *Method:*
  * the GHJS-transfer Liouville normal form v″ = 𝒬v, with 𝒬 exactly rational
    in κ and κ²-density (μ₊ − μ₋)²/4 > 0;
  * boundary terms that vanish in the class.  At the sonic point
    xT → 1 + σ(κ) with Re σ < 0.  At the centre T = −3 + O(t²) exactly.
  * conditions (a)–(d), certified pointwise on all of x ∈ (−∞, 0):
    * order-2 jet arithmetic on [−4.5 − 8.3×10⁻¹⁷, 0]: tube 14520 + sonic
      ladder 16666 + the closing box, the sonic point included.  The exact
      identities are QNm = D¹⁰·QN_red and QN_red(u₀) ≡ 0 in ℚ(V₀, √3);
    * Lemma D on y = eˣ ∈ (0, YCUT], with YCUT = 0.011108996538242306 <
      e^{−4.5}: Taylor models in y on the certified (a*, μ*) balls, 1100/1100
      subranges, and the exact layer identity E2;
    * the two regions overlap by 4.6×10⁻¹⁷.
* *Level:* **CERTIFIED**, with the "conditional on Lemma D" qualifier
  dropped.  The review's checks `centre_layers` and `s42_join` are committed
  (`7d0ec1cf`, `d68b64ff`).
* *E-clause:* it is kept only because of Lemma G1.  Without the lemma it would
  be restricted to [6, 15] × [−14, 14], where it follows from P4-6.
* *Evidence:*
  * `problems/P4/notes/s4-energy.md` §1–§3, §6;
  * `problems/P4/certificates/theorem_c/lemma_d_certified.json`;
  * `.../checks.json`;
  * `.../tc_energy.json`: the tube, sonic and closing blocks only.  Its
    `lemma_d` block is the superseded midpoint record;
  * `problems/P4/certificates/recert_2026-10-17/{tc_energy,comparison}.json`;
  * `problems/P4/certificates/theorem_c/review_s45/{centre_layers.json,s42_join.txt,sonic_residues.json,identities_out.txt}`;
  * `docs/review-2026-10-17-p4-s4-5-adversarial.md`;
  * for the E-clause: `problems/P4/notes/s4-5-assembly.md` §4.5 and
    `review_s45/{g1_order0.json,centre_pattern.json,g1_centre.txt}`.

**P4-9 (new): strip exclusion (Theorem C, part 2).**
* *Claim:* in the class of P4-8 there is no eigenvalue with 0 ≤ Re κ ≤ 6 and
  |Im κ| ≥ 14.
* *Method:*
  * a once-commuted (m = 1) energy identity for χ₁ = χ′, with weight
    g = exp∫(Re A1 − c·tt) and c = −1/32;
  * the pointwise condition S = K + ν²/4 < 0, together with dT, dB, nB ≠ 0,
    certified on piece × ε ∈ [0, 6] × τ ∈ [14, ∞) for every piece of an exact
    cover of x ∈ (−∞, 0);
  * boundary terms: γ = 1 + σ₀ + (1 − σ₀)(ε + c) ≥ 1.0617 at the sonic
    point, and BT = O(t^{3+k}) → 0 at the centre, with B0 = O(t^k), k ≥ 2;
  * Im κ ≤ −14 by conjugation.
* *Level now:* **CONJECTURED.**  546 of 576 pieces are certified on
  x ∈ [−4.56, 0] (exact cover; S4-4 `layer_t14` parts under
  `merge --substitute-layer-t14`).  The 30 centre pieces (x ≤ −4.5) are
  missing.
* *Target:* CERTIFIED when the centre pieces pass on the current code with
  the certified centre balls and the merge reports
  `strip_exclusion_certified = true`.
* *Evidence:*
  * notes: `s4-commutator.md` §2–§3, `s4-4-layer-engine.md` §1–§5 and
    §12–§14, `s4-5-assembly.md` §4.1;
  * `parts/layer_t14/` (aggregate `e9835a41…`);
  * `layer_recompute/merge_substitute_layer_t14.json`;
  * `certificates/theorem_c/checks.json`;
  * the two 2026-10-17 review documents and `review_s45/`.

  When the centre lands, add its parts and the complete `tc_commutator.json`.

**P4-10 (new): Theorem C and the mode count.**
* *Claim:* in the class of P4-8 there is no eigenvalue in
  Ω = {Re κ ≥ 0} \ [0, 15] × [−14, 14].  P4-8 covers Re κ ≥ 6; P4-9 and
  conjugation cover the rest.
  * With P4-6, the eigenvalues with Re κ ≥ 0 are exactly
    κ̄ ∈ [0.35569919907, 0.35569920835] (pure gauge) and
    κ₁ ∈ [2.81055254439, 2.81055255326].  Both are simple zeros of E.
  * So κ₁ is the unique non-gauge mode with Re κ ≥ 0, and
    γ = 1/κ₁ ∈ [0.3558019218, 0.3558019231].
  * Corollary (Lemma G1): E(κ) ≠ 0 on Ω.
  * Not claimed: kink modes, non-spherical modes, nonlinear stability,
    Re κ < 0.
* *Level now:* **CONJECTURED** (it inherits P4-9).
* *Target:* CERTIFIED with P4-9.  If P4-8 or P4-9 ends at WRITTEN_PROOF, so
  does P4-10.
* *Evidence:* the P4-6, P4-8 and P4-9 evidence; `s4-5-assembly.md`;
  `docs/review-2026-10-17-p4-s4-5-adversarial.md`.

**Related edits for the owner:**
* P4-2: widen its scope note from "in [0,15]×[−14,14]" to Re κ ≥ 0 once
  P4-10 is CERTIFIED.
* FORMULATION v1.5: record Theorem C, with Lemma G1 for the E-form
  (`theorem-c-plan.md:365-370`).

## 7. Review status and what remains to check

The S4-5 review (`docs/review-2026-10-17-p4-s4-5-adversarial.md`) covered this
note's original checklist, except the centre:

| item | status |
|---|---|
| statement, class, eigenvalue form vs E-form | done (review §2, §11); E-form by Lemma G1 |
| region algebra and boundaries | done (review §3) |
| S6 identity; `direct_S_boxes` vs the layer Σ | done (review §4; sympy I5, I6) |
| S7 limit argument | done (review §4) |
| S8 boundary terms, both ends; the S4-2 analogue | done (review §4; `sonic_residues`, `centre_layers`) |
| factorizations, nonvanishing, endgame | done (review §5; `factor_checks`) |
| layer parts, aggregate hashes, spot recompute | done (review §8, §9; the re-check of the layer-engine review) |
| exact cover, joins, the centre overlap | done (review §8; `cover_independent`) |
| conjugation | done (review §6) |
| S4-2: Lemma D, tc_energy blocks, boundary terms, the join | done (review §7; `s42_join`) |
| provenance | done (review §8, after `a0a1b5ff`) |
| **centre** | **open**: see below |

When the centre lands, check:
* its route and code version, the certified centre balls, and that the
  closing piece contains y = 0;
* that both τ phases, or a reviewed centre engine, cover τ ≥ 14 for all
  12 cells;
* the nonvanishing requirements of I8;
* that the merge reports `complete` and `strip_exclusion_certified`;
* the merge-output wording fix (§5 item 2);
* then apply §6 to P4-9 and P4-10.

## 8. Reproducing the hashes and checks in this note

From the repository (worktree) root:

```
shasum -a 256 problems/P4/certificates/theorem_b/*.json problems/P4/certificates/theorem_c/*.json \
  problems/P4/results/theorem_c/tc_leading.json problems/P4/results/theorem_c/parts/plan.json \
  problems/P4/results/theorem_c/parts/caches/*.json problems/P4/results/theorem_c/layer_recompute/*.json \
  problems/P4/certificates/theorem_c/review_s45/*.json problems/P4/certificates/theorem_c/review_s45/*.txt
# aggregate of a part directory (paths relative to the root, C-locale sort):
find problems/P4/results/theorem_c/parts/layer_t14 -name '*.json' | LC_ALL=C sort \
  | xargs shasum -a 256 | shasum -a 256
git log --oneline fda58897..HEAD -- problems/P4/src   # empty
```

The first version's own exact cover analysis used Python `fractions` on two
sources:
* the `parts/plan.json` endpoints: 371 same-tube label gaps of at most
  1.8e−16, which are label rounding, and no cross-region gap;
* the dyadic step positions of `P4/data/cache/tube_{full,gap,deep,band}.json`:
  radius 0 and exact joins.

The merge's exact cover (`ccc74f3a`) and the review's
`review_s45/cover_independent.txt` now supersede it.  The review's checks are
re-run as `P4/certificates/theorem_c/review_s45/README.md` describes.

---

*Paths (2026-10-17): the files behind CERTIFIED rows moved unchanged from `results/` to `certificates/` (`problems/P4/certificates/README.md`). The paths in this note are updated; the module mains still write their replays to `results/`. The note gained no lines above this one, so line references into it still hold.*
