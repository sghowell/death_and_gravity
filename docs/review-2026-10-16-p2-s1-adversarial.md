# Independent adversarial review: P2 stage S1

Date: 2026-10-16. Reviewer: a separate agent with no part in authoring
P2. It read the code, certificates and proofs before the S1 report. The
review was read-only: repository `git status` was unchanged, and commands
that rewrite files were run in a scratch copy.

**Overall verdict:** no false mathematical claim. S1 is a sound foundation
for S2 and S3 at both points, once the KAT-6 Krawczyk bound is fixed and
D1′ is restated for the localized resonance notion. The exception is the
stretch target P2-6, which needs new foundations (§3).

## 1. Verdicts

| Claim | Verdict | Level |
|---|---|---|
| Parameters (params.json) | PROVEN AS STATED | CERTIFIED |
| P2-0, mode stability | PROVEN AS STATED | WRITTEN_PROOF |
| D4, simple zero mode | PROVEN AS STATED | WRITTEN_PROOF |
| D1′, Wronskian zeros = resonances | GAP (presentation for P2's use; claim-level read literally) | WRITTEN_PROOF after the wording fix |
| P2-2, high-frequency thresholds | PROVEN AS STATED | CERTIFIED on −(½+δ)κ₋ ≤ Im ω ≤ 0 |
| KAT-6, mode enclosures | FLAW (implementation-level, non-fatal) | winding CERTIFIED; Krawczyk CERTIFIED after the fix |

**Parameters.**
- Recomputed independently with mpmath at 60 digits: polynomial roots
  with sign brackets, and γ₀ from the geodesic potential F/r². This
  matches HV17 §2.7, where γ₀ = μ(r_P)ν_min/4.
- Every ball and margin matches. In κ₋ units: P_B 0.27289 / 1.64507 /
  1.13321 / 0.64578; P_A 0.16684 / 0.48009 / 0.25430 / 0.35369.
- Presentation problems:
  - s1-report.md:35 says the stored radii are below 10⁻⁷⁰, but the
    params.json strings carry radii up to 3.6e-39; only the in-memory
    balls are ≤ 2e-72.
  - params.py:219 `CF17_2alpha` (κ₋ < (1+2δ)κ₋) is a tautology.
  - `vieta_sum` and `Qmax2_two_forms` are overlap checks labelled
    "certified".

**P2-0.**
- The reviewer's own sympy derivation from the (t*, r) metric (block
  determinant −1) reproduces: E_l[u] and the pointwise identity;
  b·Re − a·Im = −[b(D+LN+|ω|²𝒲) + |ω|²B]; J(r_c) − J(r₊) = iωB; and the
  weight 12r⁴/(Λ(r−r₋)(r−r_n)(r_c−r₊)²).
- The written step is correct: r²Fu′ū → 0 at both horizons because u is
  analytic there; b = 0 gives a²B = 0, which contradicts u₊(r₊) ≠ 0; and
  ω = 0 gives D + LN = 0.
- The machine-checked part is algebra only, and several "identities"
  substitute into hand-typed formulas: C3 (p20_identity.py:94–106), C7a/b
  (horizon p₁, p₂ typed in, :134–153), C8 (limit after G(r)→G(r_p),
  :174–175), C10 (:220–228). Regularity at the horizons, the integration
  step, and "positivity ⇒ u ≡ 0" are prose. That is the WRITTEN_PROOF case.

**D4.**
- Re-derived by hand: Δẏ′ ≡ −ir₊² and Δż′ ≡ ir_c², so
  W₀′(0) = i(r₊²+r_c²)/(r₊r_c).
- An independent solver (Frobenius with the unfactored indicial
  polynomial, plus mpmath odefun) with central differences gives
  5.54913720005i (P_B) and 5.32785787624i (P_A). These agree with the
  formula to 3e-14 along both step directions.
- "Limits exact" overstates: C8 is a one-line substitution, and
  analyticity in ω and constancy of Δẏ′ are prose.

**D1′.**
- The citations are accurate against the arXiv TeX of all three papers
  (Vasy Thm 2.6, Rem 2.8, §2.6, Lemma 3.1; HV17 Lemma 2.14, the §2.5
  hypothesis, t₂₃ with s₂ = −1, s₃ = +1; Hintz EqIBoxSol, the T̃ typo,
  the remark after Thm 5.5).
- Correct for HV17's Ω₂₃ forward problem and Hintz's QNM set: branch
  selection (|t|^{iσ/κ} ∈ H^s iff s < ½ − Im σ/κ), exclusion of E₊ and
  E_c, and the Abel-integral upper bound and cutoff lower bound for
  multiplicity.
- The gap is what "resonance" means. notes:87–90 and :97–99 equate
  resonances with poles of Vasy's unlocalized (P_σ−iQ_σ)⁻¹ and of HV17's
  global P̂(σ)⁻¹ (Thm 2.13). HV17 itself warns that the global problem can
  have other resonances (§2.5; Lemma 2.14(1): states supported in
  r₀ ≤ r ≤ r₁). Vasy proves independence from Q only for the localized
  ψ(P_σ−iQ_σ)⁻¹ (§3). So statement (1)'s "any of the three sources"
  (line 132), and the lower-bound step, which ignores the −iQ_σa term,
  need that localization.

**P2-2.**
- Checked, all correct: θ/p; both Volterra kernels; the Grönwall and h′/p
  bounds; the Wronskian expansion; the Jost identification (needs only
  β < κ₊, κ_c); |Im p| ≤ a|b|/√A and the G majorant; J dr = |θ/p| dx; cell
  coverage; the cumulative weights; the end-cell bound (|A_h| = 1/(2κ_h)
  exactly); outward rounding; and "error < 1 ⇒ W ≠ 0".
- lg_bounds uses no endpoint series; the end cells are analytic power-law
  bounds. The Frobenius and Taylor tail majorants in arbw.py are correct.
- Independent spot checks:
  - a float re-evaluation of the same majorant integrals for l = 0, 1, 2,
    10, 25, 50 at both points gives criterion values 1.80–1.91, against
    certified upper values 1.98–2.00;
  - the actual |W e^{iΘ}/(2iω) − 1| from an independent solver at
    Re ω = W_l (l = 0, 1, 50) is 0.004–0.015;
  - across all 102 certificates the criterion lies in [1.9479, 1.99997],
    s ≤ 0.767, and the thresholds are monotone in l.
- The Im ω > 0 part rests on P2-0, so the full row takes P2-0's level.
- Presentation: FORMULATION.md:197–199 promises Olver control "along the
  real progressive path". For horizon-referenced dominant solutions the
  real path is anti-progressive (|kernel| = e^{2∫|Im p|} ≥ 1). The
  implemented weighted bound is valid but different. The threshold symbol
  W_l also clashes with the Wronskian W_l(ω).

**KAT-6.**
- winding.py:149 bounded |W′(z) − W′(z₀)| by 2Mr/(R−r)². That holds only
  for |z − z₀| ≤ r, but the Krawczyk set Z is a square with corners at
  √2·r. The correct bound is 2√2·Mr/(R−√2r)²: 1.79e-6 against the coded
  1.27e-6 at P_A. So the committed `certified_unique_zero` rested on an
  invalid enclosure of W′(Z).
- Re-run with the corrected bound: still certifies, with contraction
  3.96e-7 (P_A) and 4.2e-10 (P_B) and the same enclosures. The winding
  part is valid, and independent 45-digit Newton zeros lie inside both
  enclosures.
- **Fixed on 2026-10-16** (winding.py now uses the corner distance √2·r
  and checks R > √2·r; certificate regenerated, both cases certify).

## 2. Replays

Read-only commands were run in the repository; commands that rewrite
files were run in a scratch copy.

- `p2.params verify`: OK (0.24 s).
- `p2.p20_identity verify`: OK (0.6 s).
- `p2.lg_bounds verify PB 50` and `PA 50`: both OK (2.3 s each).
- `p2.winding`: KAT-6 pass (11.6 s).
- `p2.kats`: KAT-1 (72 s), KAT-2 (32 s), KAT-3 (86 s), KAT-4 (54 s), KAT-5
  (11 s) and KAT-7 (41 s) all pass; every output is identical to the
  committed files apart from commit and timing fields.
- `pytest`: 40 passed, 3 skipped (9.5 s). `P2_SLOW=1 pytest`: 43 passed
  (95 s).
- Replay hygiene: `test_kat6_rerun`, `p2.winding` and `p2.kats` overwrite
  the committed certificate and result files; KAT-6 has no verify mode;
  test_kats.py:27–31 only reads the stored `pass` flags.

## 3. The P2-6 obstacle

P2-6's strip edge at −0.92κ₋ at P_B passes E₊'s first point,
−iκ₊ = −0.823iκ₋. W₀ has a simple pole there (|W|·ε ≈ 0.035 at both
ε = 1e-4 and 1e-6), so the planned "l = 0 winding 2" would count zeros
minus poles. That edge is also outside D1′'s multiplicity window and
outside the Liouville–Green condition β < κ₊. P2-6 therefore needs new
foundations. The main target needs the strip only down to −(½+δ)κ₋
(−0.55κ₋ at P_B), well above that pole, and is unaffected.

## 4. Wording fixes recommended for the ledger rows

- **P2-0:** level WRITTEN_PROOF. Evidence: "sympy algebra C1–C6 (C3, C7,
  C8 and C10 are substitutions); regularity, integration and sign argument
  written (notes §2); D4 written (§3), Arb finite-difference cross-check
  1e-30."
- **P2-2:** rename the threshold (e.g. a_l). Level: "CERTIFIED on
  −(½+δ)κ₋ ≤ Im ω ≤ 0 (Arb, weighted non-progressive Volterra bounds);
  Im ω > 0 by P2-0; Re ω < 0 by exact symmetry."
- **P2-D1′:** "Resonances of the forward problem on Ω₂₃ [HV17 §2.5] =
  poles of ψ(P_σ−iQ_σ)⁻¹ [Vas13 §3] = QNM set [Hin25 §1]; multiplicity
  for Im σ > −min(κ₊, κ_c); WRITTEN_PROOF."
- **P2-par:** drop the "radii < 1e-70" claim or store full-precision
  balls; relabel the consistency checks; drop `CF17_2alpha`.
- **P2-KAT:** KAT-6 is certified only after the Krawczyk fix.
- **s1-report §7.1:** "shallowest mode" is float evidence, not certified.
