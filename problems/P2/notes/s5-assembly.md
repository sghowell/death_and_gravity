# P2 S5: the assembly P2-4 and the conversion P2-5 at P_B

Status: written 2026-10-26 on branch `p2/s5-assembly`, from main at `3dd0b2bc`. It has not yet been
independently reviewed.
- Plan: `notes/s5-plan.md` (`e2c75aa6`).
- Machine check: `src/p2/assembly.py` (`7aa9fbd0`).
- Certificate: `certificates/assembly/assembly_PB.json` (`942f89c1`).
- Test: `tests/test_assembly.py`.

**Sources** (TeX read):
- **[HV17]** P. Hintz, A. Vasy, arXiv:1512.08004, file `cauchyhor-r.tex`, dated "December 25,
  2015. Final revision: July 12, 2017" (J. Math. Phys. 58, 081509).
  - Theorem numbers are the ones this file produces: one counter per section, shared by theorems,
    propositions, lemmas, corollaries, definitions and remarks. They were checked by counting the
    environments.
  - Labels are given in brackets, and ℓ. is a line of the TeX file.
- **[CF17]** J. L. Costa, A. T. Franzen, arXiv:1607.01018, file `RN_deSitterRevised6.tex` (Ann.
  Henri Poincaré 18, 3371).

**Dictionary.**
- [HV17]: r₁, r₂, r₃ are our r₋, r₊, r_c; κ₁ is κ₋; μ is F; λ = Λ/3. Its σ is our ω; both carry the
  time dependence e^{−iσt*}.
- [CF17]: D is F, e is Q, and κ_* is ours.

**Notation.**
- d_δ := (½+δ)κ₋, which is (11/20)κ₋ at P_B (δ = 1/20).
- S_δ = {Im ω ≥ −d_δ}.
- E₊ = {−inκ₊ : n ≥ 1} and E_c = {−inκ_c : n ≥ 1}.

## 1. P2-4 at P_B

### 1.1 Statement

**P2-4** (P_B: ΛM² = 3/50, Q = (999/1000)Q_max, δ = 1/20).

- **(a)** For every l ≥ 0, W_l has no zero in S_δ other than ω = 0 at l = 0, and that zero is
  simple: m₀(0) = 1.
- **(b)** Equivalently, by P2-D1′:
  - the forward problem for □_g on Ω₂₃ ([HV17] §2.5) has no resonance in S_δ∖{0};
  - σ = 0 is a resonance of multiplicity 1 (a simple pole of rank 1), and its resonant states are
    the constants.
- **(c)** Hence β ≥ ½ + δ = 0.55 (D3). With P2-K6, 0.55 ≤ β ≤ 0.90938098.

(a) is FORMULATION's P2-4 in D1's terms: "S_δ∖{0} contains no QNM of any degree, and 0 is a simple
QNM". (b) is the same statement in the resonance terms that [HV17] uses, and it is what P2-5
consumes.

### 1.2 Ingredients

| row | level on main `3dd0b2bc` | used for |
|---|---|---|
| P2-par | CERTIFIED | d_δ < κ₊ and d_δ < κ_c (E₊ ∪ E_c lies below S_δ; S_δ lies in D1′'s window); d_δ < γ₀ and κ₋ < 2κ₊ (P2-5) |
| P2-0 | WRITTEN_PROOF | Im ω > 0 for every l; W_l(0) ≠ 0 for l ≥ 1; D4: W₀(0) = 0 and W₀′(0) = i(r₊² + r_c²)/(r₊r_c) |
| P2-D1′ | WRITTEN_PROOF | (b): zeros of W_l are the resonances of the Ω₂₃ forward problem, with multiplicity (2l+1)·m_l on Im ω > −min(κ₊, κ_c) |
| P2-1 | WRITTEN_PROOF | l ≥ 20 on H₊ := {Re ω ≥ 0, −d_δ ≤ Im ω ≤ 0} |
| P2-2 | CERTIFIED on the lower strip | l ≤ 19 on {abs(Re ω) ≥ a_l, −d_δ ≤ Im ω ≤ 0} |
| P2-3B | CERTIFIED (the counts) | l ≤ 19 on R′_l = [−a_l, a_l] × [−b, e] |
| P2-K6 | CERTIFIED zero; the β bound at P2-D1′'s level | the upper bound in (c) |

The conjugation symmetry W_l(−ω̄) = conj W_l(ω) is also used, as P2-1 and P2-3B use it. It is
checked in `certificates/part3/zeq_identity.json` (Z4), and it is also immediate from FORMULATION
§1.2:
- (R) has real coefficients and depends on ω only through ω²;
- the exponents ∓iω/(2κ) go to their conjugates under ω ↦ −ω̄;
- c₊ and c_c are real.

Nothing else is used.

### 1.3 Proof

Write U = {Im ω > 0} and H± = {±Re ω ≥ 0, −d_δ ≤ Im ω ≤ 0}, so that S_δ = U ∪ H₊ ∪ H₋.

**Step 0 (W_l on S_δ).** By P2-par, d_δ < κ₊ and d_δ < κ_c, with certified margins 0.01161 and
0.06998. Every point of E₊ ∪ E_c has Im ω ≤ −min(κ₊, κ_c) < −d_δ. So:
- W_l is analytic on a neighbourhood of S_δ (FORMULATION §1.2), and D1 applies at every point of
  S_δ;
- S_δ ⊂ {Im ω > −min(κ₊, κ_c)}, the window of P2-D1′'s multiplicity clause.

**Step 1 (U).** By P2-0, for Im ω ≥ 0, W_l(ω) = 0 iff (l, ω) = (0, 0). So no W_l vanishes on U, and
W_l(0) ≠ 0 for l ≥ 1.

**Step 2 (H₋).** The map ω ↦ −ω̄ takes H₋ onto H₊, and W_l(−ω̄) = conj W_l(ω). So W_l has a zero
in H₋ iff it has one in H₊.

**Step 3 (l ≥ 20).** By P2-1, for every l ≥ l₀ = 20, W_l(ω) ≠ 0 on the closed half-strip H₊. With
Steps 1 and 2, W_l has no zero in S_δ for any l ≥ 20.

**Step 4 (l ≤ 19).** Fix l ≤ 19.
- P2-3B gives three facts, with b = 6280879/2²⁸ and e = 4567911/2³⁰:
  - W_l is analytic on a neighbourhood of R′_l = [−a_l, a_l] × [−b, e];
  - W_l does not vanish on ∂R′_l;
  - its winding number around ∂R′_l is 1 for l = 0 and 0 for 1 ≤ l ≤ 19.

  By the argument principle, counting multiplicity, W₀ has one zero in R′₀, and W_l has none in R′_l
  for l ≥ 1.
- P2-2: W_l(ω) ≠ 0 for abs(Re ω) ≥ a_l and −d_δ ≤ Im ω ≤ 0, with the same a_l. R′_l was built from
  P2-2's stored threshold (`W_l_exact`), and the assembly check confirms that the two numbers are
  identical.
- b ≥ d_δ (b − d_δ = 4.546·10⁻¹⁰ > 0, certified) and e > 0. Hence

      {−d_δ ≤ Im ω ≤ 0} ⊂ R′_l ∪ ({abs(Re ω) ≥ a_l} × [−d_δ, 0]).

So for 1 ≤ l ≤ 19, W_l has no zero in H₊ ∪ H₋, and with Step 1 none in S_δ.

For l = 0:
- W₀(0) = 0 by D4: R ≡ 1 solves (R) at l = 0, ω = 0 and is analytic at both horizons, so
  R₊ = R_c = 1.
- 0 is an interior point of R′₀.

The single zero counted in R′₀ is therefore ω = 0, with multiplicity m₀(0) = 1. W₀ has no other zero
in R′₀, none in the P2-2 region, and none in U.

**Step 5 (every l).** {0, …, 19} ∪ {l ≥ 20} = ℕ₀, so Steps 3 and 4 give (a). Given W₀(0) = 0, Step 4
also re-proves D4's simplicity at the CERTIFIED level.

**Step 6, (b).** By Step 0, every σ ∈ S_δ satisfies σ ∉ E₊ ∪ E_c and Im σ > −min(κ₊, κ_c).
- P2-D1′(1): σ is a resonance of the Ω₂₃ forward problem iff W_l(σ) = 0 for some l.
- By (a), the only such σ in S_δ is σ = 0, and only l = 0 contributes there (Step 1).
- P2-D1′(2):
  - the pole order at 0 is max_l m_l(0) = 1;
  - there is one Jordan chain, of length m₀(0) = 1, in the single sector (l, m) = (0, 0);
  - the multiplicity is Σ_l (2l+1)m_l(0) = 1.
- The resonant state is u₀(r)Y₀₀ with u₀ = e^{−i·0·T}R₊ = 1, a constant.

**Step 7, (c).** D3 defines α = inf{−Im ω : ω a QNM of some degree, ω ≠ 0}.
- By (a), every nonzero QNM lies outside S_δ, i.e. has −Im ω > d_δ. So α ≥ d_δ, and
  β = α/κ₋ ≥ 11/20.
- P2-K6: W₀ has a certified simple zero ω_*, off E₊ ∪ E_c, with
  ω_*/κ₋ ∈ −0.909380978994805455524933762622i ± 3.2·10⁻³¹
  (`certificates/part3/krawczyk_PB_l0.json`; KAT-6 agrees). So α ≤ −Im ω_* and
  β ≤ 0.90938098. ∎

**Remarks.**
- β ≥ 0.55 is stated for the closed strip, which is what [HV17]'s hypothesis asks for (§2.3).
  Strict inequality β > 0.55 is not claimed. It would follow from discreteness of the QNMs below the
  strip, uniformly in l, and no P2 row states that.
- The two-sided bound is not sharp. The S0 float value is β ≈ 0.909381 (the l = 0 NE mode), and
  nothing between 0.55κ₋ and 0.909κ₋ is explored. P2-6, the two-sided enclosure, stays the stretch
  goal and needs the new foundations recorded in FORMULATION v0.2.

### 1.4 The machine check

`python -m p2.assembly` reads only committed certificates: 54 files, each pinned by its sha256 in the
output. It checks with exact rationals and 256-bit Arb balls. At `7aa9fbd0`, all 366 checks pass
(`certificates/assembly/assembly_PB.json`), in 2.7 s of CPU.

What it checks (names as in the JSON):
- **C1, the point.** qfrac, δ, Λ and M. `Point("PB")` is recomputed, and it overlaps every stored
  ball.
- **C2, each l ≤ 19.**
  - The lg certificate: point, l and δ; strip depth = d_δ; criterion upper bound < 2.
  - a_l in R′_l is identical to P2-2's `W_l_exact`.
  - b ≥ d_δ (margin 4.546·10⁻¹⁰), e > 0, and b < κ₊ (analyticity on R′_l).
  - The stored discs, 2351 pieces in all:
    - they trace the right half of ∂R′_l exactly and contiguously, from (0, −b) through
      (a_l, −b) and (a_l, e) to (0, e);
    - each piece lies within ρ/2 of its centre;
    - each disc box clears E₊ ∪ E_c.
  - The increments, re-summed in Arb and divided by π, lie within 1/4 of the stored winding: 1 at
    l = 0, 0 otherwise.
  - All 20 values of l are present.
- **C3, l ≥ 20.**
  - `part1.cover()`, replayed from the regime certificates, reproduces the committed `cover_PB.json`.
  - Every regime certificate passes, is at P_B with h₀ = 2/41, and carries the cover's design
    constants. The S4 full replay (stride 1) passed.
  - l₀ = ⌈1/h₀ − ½⌉ = 20 = 19 + 1, so no l is missed.
  - For every h ∈ (0, h₀], the Re ω ranges of L, S, T and A union to [0, ∞). S's top at h₀ is
    90.75κ₋, above L's top of 8κ₋. The other overlaps are exact.
  - The Im ω range is [−d_δ, 0]:
    - L's design range is −11/20;
    - T's and A's lower ends enclose −d_δ;
    - S is covered through the cover's tiling.
  - The symmetry identity Z4 holds.
- **C4.** d_δ < κ₊ and d_δ < κ_c.
- **C5.**
  - 0 ∈ int R′₀.
  - W₀′(0) = 5.5491372 i ≠ 0.
  - P2-0's identity steps are recorded ok.
- **C6.** β ≤ 0.90938098 (margin 1.005·10⁻⁹), and ω_* lies outside S_δ.
- **C7, C8.** See §2 below.
- **C9.** The ingredient levels and their minimum. The test cross-checks the levels against
  `CLAIMS.md`.

**Negative controls.** Each mutates a copy of the inputs and must fail on its named item. All nine
do.

| id | mutation | fails on |
|---|---|---|
| N1 | b → 6280878/2²⁸ at l = 0, just below d_δ | C2.l0.b_ge_d (and the containment that uses it); nothing else |
| N2 | `winding_PB_l7.json` removed | C2.l7.present, C2.all_l_0_to_19_present |
| N3 | a₅ → a₅ − 2⁻²⁴ in R′₅ | C2.l5.a_identical |
| N4 | h₀ → 2/43 | C3.no_l_gap (l₀ = 21), and the h₀ consistency checks |
| N5 | e → 0 at l = 0 | C2.l0.e_pos, C5.zero_interior_R0 |
| N6 | stored winding 1 at l = 3 | C2.l3.winding_value, C2.l3.resum |
| N7 | T_PB fails | C3.regime_pass.T_PB, C3.cover_pass |
| N8a | γ₀ → κ₋/2 | C7.alpha_lt_gamma0 |
| N8b | κ₊ → 0.45κ₋ | C8.km_lt_2kp, C8.condition |

**What the check does not do.** It does not re-verify the analytic content of the regime, lg or
winding certificates. Their own replays do that: `p2.part1 verify`, `p2.lg_bounds verify` and
`p2.part3 verify`. The check confirms that the certified statements fit together. In S5,
`p2.params verify` and `p2.lg_bounds verify PB 19` were replayed, and both pass.

### 1.5 Level

The level of P2-4 is the minimum of its ingredients' levels: **WRITTEN_PROOF**.

Why the minimum, and why the certified pieces cannot raise it:
- Every point of S_δ∖{0} must be excluded for every l. The pieces that do this are all
  WRITTEN_PROOF:
  - P2-0 on U, for every l;
  - P2-1 for l ≥ 20, infinitely many l;
  - P2-D1′, the translation into resonances, for every σ.
- The CERTIFIED rows (P2-par, P2-2, P2-3B) cover only the lower strip for l ≤ 19 and the numerical
  inequalities. None of them replaces P2-0 on U, P2-1 for l ≥ 20, or D1′.
- The assembly adds a machine-checked bookkeeping certificate and only routine written steps: the
  set decomposition, the argument principle given the certified count, and the D3 infimum. It does
  not lower the level either.

P2-4 would be CERTIFIED only if P2-0's written steps, P2-1's W1 and W2, and D1′ were machine-checked.

Finer statement: for l ≤ 19, part (a) on the lower strip −d_δ ≤ Im ω ≤ 0 is CERTIFIED (P2-2 with
P2-3B and P2-par). Everything else is WRITTEN_PROOF. The README requires a check by an independent
reviewer for WRITTEN_PROOF, and this note has not had one yet.

### 1.6 What P2-4 does not say

- Nothing at P_A (P2-7), and nothing at other parameters.
- No sharp value of β (§1.3, Remarks).
- Nothing below the strip.

## 2. P2-5 at P_B: the conversion by citation

### 2.1 Statement

**P2-5** (P_B; conversion cited). Let α := d_δ = (11/20)κ₋.

**(i)** [HV17] Thm 1.1 applies with this α, in its quantitative form (Prop 2.17, Cor 2.19,
Thm 2.24). Let u solve □_g u = 0 with smooth initial data on [HV17]'s Cauchy surface H_I. Then, in
the domain of dependence of H_I, u = u₀ + u′ with:
- u₀ ∈ ℂ;
- |u′| ≤ Ce^{−αt₀} uniformly in r > r₋, and the same for any number of derivatives along stationary
  vector fields tangent to the Cauchy horizon;
- u continuous up to the Cauchy horizon;
- u′ ∈ e^{−αt₀}H^{1/2+α/κ₋−0} = e^{−αt₀}H^{21/20−0} in t₀ > 0.

In particular, u ∈ H¹_loc up to the Cauchy horizon.

**(ii)** [CF17] Thm 3.1 applies with p = α to solutions with smooth data on a Cauchy surface of RNdS
at P_B. They have bounded energy up to the Cauchy horizon, and they extend to it in C⁰ ∩ H¹_loc.

Neither theorem is re-proven. The hypotheses they take from P2 are exactly P2-4(b), D4 and P2-par.

This is the linear analogue, at P_B, of the failure of Christodoulou's formulation of SCC for smooth
data (FORMULATION P2-5). It concerns only the linear massless scalar on fixed RNdS (FORMULATION
§7).

### 2.2 [HV17]: hypotheses quoted and matched

**H1, the setting.**
- Thm 1.1 [ThmIntroMain] (ℓ. 368): "Let g be a non-degenerate Reissner–Nordström–de Sitter metric
  with non-zero charge Q, or a non-degenerate Kerr–de Sitter metric with small non-zero angular
  momentum a, with spacetime dimension ≥ 4."
- Def 2.1 [DefRNdSNonDegenerate] (ℓ. 512): "We say that the Reissner–Nordström–de Sitter spacetime
  with parameters Λ > 0, M > 0, Q > 0 is non-degenerate if μ has 3 simple positive roots
  0 < r₁ < r₂ < r₃."
- (EqRNdSMetric) (ℓ. 507): g = μ dt² − μ⁻¹dr² − r²dω², with μ = 1 − 2M/r + Q²/r² − λr² and
  λ = Λ/3.

Match:
- μ is FORMULATION's F (§1.1).
- At P_B, Λ = 3/50, M = 1 and Q = 1.00984… > 0.
- Δ = r²F has four simple real roots r_n < 0 < r₋ < r₊ < r_c (P2-par: four disjoint unique-root
  enclosures). So F has exactly three positive roots, all simple.
- The dimension is 4.

[C7.Lambda_M_Q_positive, C7.three_simple_positive_roots]

**H2, the field and the data.**
- Thm 1.1: "if u is the solution of the Cauchy problem □_g u = 0 with smooth initial data".
- §1 (ℓ. 355): "We consider the Cauchy problem for the linear wave equation with Cauchy data posed
  on a surface H_I as indicated in Figure 1."
- Prop 2.4 (EqRNdSMfdHyper, ℓ. 749): H_I = {t* = 0, r₂ − 2δ ≤ r ≤ r₃ + 2δ}. This is a spacelike
  hypersurface that crosses the event horizon and the cosmological horizon.
- §2.3 (ℓ. 992): "Forward solutions to the wave equation □_g u = f in the domain of dependence of
  H_I, i.e. in Ω ∩ {r > r₁}, are not affected by any modifications of the operator □_g outside, i.e.
  in r ≤ r₁."
- §1 (ℓ. 388): "a minor extension of our arguments yields analogous boundedness, decay and
  regularity results for the Cauchy problem with a 'two-ended' Cauchy surface H_I up to the
  bifurcation sphere B".

Match:
- The field is FORMULATION's: a massless, minimally coupled scalar.
- The data live on H_I, a one-ended partial Cauchy surface. The conclusions hold in its domain of
  dependence, up to the part of the Cauchy horizon to the future of H_I, towards i⁺. So P2-5(i) is
  stated for data on H_I, and FORMULATION's "smooth Cauchy data on RNdS" is read that way.
- Take smooth data on a complete Cauchy surface. Apply (i) to the restrictions of the data to the
  translates of H_I in Killing time, and to their mirror images in the other exterior. This gives
  u ∈ H^{21/20−0}_loc near every point of the Cauchy horizon except the bifurcation sphere B, with
  constants that depend on the translate.
- Up to B, [HV17] only sketches the two-ended extension. That case is not claimed.

**H3, the resonance notion.** See §2.3 below. Result: it is the localized notion.

**H4, 0 simple with constant resonant states.**
- §2.5 (ℓ. 1313): "The only resonance of the forward problem in Ω₂₃ in Im σ ≥ 0 is a simple
  resonance at σ = 0, with resonant states equal to constants, see the references given in the
  proof of Lemma 2.14".
- Lemma 2.14(2) [LemmaRNdSResonanceSupp] (ℓ. 1257): "If R denotes the restriction of distributions on
  Y to r > r₁, then the only pole of R∘P̂(σ)⁻¹ with Im σ ≥ 0 is at σ = 0, has rank 1, and the space
  of resonant states consists of constant functions."
- Its proof (ℓ. 1268) reduces the order of that pole to "the order of the 0-resonance of the forward
  problem for □_g in Ω₂₃, which is known to be equal to 1, see the references above".
  - Those references are for Schwarzschild–de Sitter and Kerr–de Sitter.
  - [HV17] says that the integration-by-parts argument for Im σ > 0 "applies without change in the
    present setting" (ℓ. 1264).
  - For the order at 0, it gives no argument specific to RNdS.

Match: at P_B, P2-4(b) (Step 6) supports [HV17]'s own claim independently:
- D4: m₀(0) = 1, with W₀′(0) = 5.549i;
- P2-0: no contribution from l ≥ 1, and no other resonance in Im σ ≥ 0;
- P2-D1′: multiplicity 1, pole order 1, constant state.

**H5, α < γ₀.**
- §2.5 (ℓ. 1313): "In the notation of Proposition 2.12, we may assume −α > −γ₀ so that we have high
  energy estimates in Im σ ≥ −α."
- Prop 2.12 [PropRNdSSemiclassical2] (ℓ. 1186): the estimates hold "provided h⁻¹Im z ∈ [−γ, c₂],
  where −γ₀ < −γ < c₂".
- §2.7 (ℓ. 1508): "α < γ₀ := μ(r_P)ν_min/4 = (2r_P)⁻¹(μ(r_P)(2 − 3r_P⁻¹M))^{1/2}", where r_P is
  the critical point of r⁻²μ.

Match:
- This γ₀ is FORMULATION D6's: λ_P/2, with r_P = (3M + √(9M² − 8Q²))/2, the critical point of
  r⁻²F.
- P2-par: d_δ < γ₀, with margin 0.04821 (d_δ/γ₀ = 0.3268).
- §2.5 allows [HV17] to shrink α to get below γ₀. Here no shrinking is needed.

[C7.alpha_lt_gamma0]

**H6, the exponent.**
- (EqRNdSRadialPointBeta, ℓ. 849): β₁ = 2/|μ*′(r₁)|. (EqRNdSSurfGrav, ℓ. 854): "κ_j := β_j⁻¹ is
  equal to the surface gravity of the horizon at r = r_j". Here μ* = μ near r₁.
- Prop 2.17 [PropRNdSPartialAsymp] (ℓ. 1315): "u′ ∈ H_b^{1/2+αβ₁−0,α} near r = r₁".
- Thm 2.24 [ThmRNdSPartialAsympConormal] (ℓ. 1452): "V₁⋯V_N u′ ∈ H_b^{1/2+αβ₁−0,α} for all
  N = 0, 1, … and all vector fields V_j ∈ V_b(Ω) which are tangent to the Cauchy horizon r = r₁".
- Thm 1.1: "u′ … lie[s] in the weighted spacetime Sobolev space e^{−αt₀}H^{1/2+α/κ₁−0} in t₀ > 0".

Match:
- κ₁ = κ₋ = ½|F′(r₋)| (FORMULATION §1.1).
- ½ + α/κ₋ = ½ + 11/20 = 21/20 > 1 exactly, so H^{21/20−0} ⊂ H¹ [C7.exponent_gt_1].
- This is FORMULATION's H^{1+δ−0}.

**H7, implicit points.**
- **The specific α.**
  - Thm 1.1 says "there exists α > 0".
  - The quantitative statements take "α > 0 as above", meaning any α that satisfies the
    hypotheses of §2.5. The introduction (ℓ. 403) says that α can be any number below the spectral
    gap.
  - P2-5 uses Prop 2.17, Cor 2.19 and Thm 2.24 with α = d_δ, which satisfies the hypotheses of
    §2.5 as stated, on the closed half-plane Im σ ≥ −α (see §2.3).
- **The contour.** This is a remark on [HV17]'s proof, not a gap in the hypotheses.
  - The proof of Prop 2.17 shifts the inverse Mellin transform to Im σ = −α. That needs the
    global family P̂(σ)⁻¹ to have no pole on the line. [HV17] does not discuss this.
  - A global pole on that line is not an Ω₂₃ resonance, by P2-4(b). So its Laurent singular part
    vanishes in r > r₁ (Lemma 2.14(1)).
  - Alternatively, the argument can be run on the line Im σ = −α̂, for α̂ < α arbitrarily close.
    The global poles are discrete: by the high-energy estimates of Prop 2.12 there are finitely
    many in the strip. This gives u′ ∈ H^{1/2+α̂/κ₋−0} with decay e^{−α̂t₀}. As α̂ ↑ α, the Sobolev
    conclusion H^{21/20−0} is unchanged.
- **Time function and signs.**
  - P̂(σ) = e^{iσt*}Pe^{−iσt*}, so resonances with Im σ < 0 decay, as for our ω.
  - D1′ shows that resonances and their multiplicities do not depend on the choice of t*
    (`s1-p20-d1prime.md` §4.1).

### 2.3 The resonance notion (the decisive check)

The question, raised by the S1 review and the S5 brief: is [HV17]'s "0 is the only resonance in
Im σ ≥ −α" about the forward problem on Ω₂₃, or about the global P̂(σ)⁻¹?
- The forward problem on Ω₂₃ is the localized notion, which P2-D1′ identifies with the zeros of
  W_l.
- The global P̂(σ)⁻¹ can have other resonances (Lemma 2.14(1)).

Here is the hypothesis verbatim. It is the paragraph that opens §2.5 [SubsecRNdSAsymp] (ℓ. 1313),
and Prop 2.17 refers to it as "α > 0 as above":

> The only resonance of the forward problem in Ω₂₃ in Im σ ≥ 0 is a simple resonance at σ = 0, with
> resonant states equal to constants, see the references given in the proof of Lemma 2.14, and there
> exists α > 0 such that 0 is the only resonance in Im σ ≥ −α. (This does *not* mean that the global
> problem for P does not have other resonances in this half space!) In the notation of
> Proposition 2.12, we may assume −α > −γ₀ so that we have high energy estimates in Im σ ≥ −α.

**Reading.**
- "0 is the only resonance in Im σ ≥ −α" is the second clause of a sentence whose subject is
  "resonance of the forward problem in Ω₂₃".
- The next sentence says explicitly that the global problem may have other resonances in that
  half-space. That parenthesis would contradict the hypothesis if the hypothesis were about the
  global P̂(σ)⁻¹.

**The notion is the localized one.**

**How [HV17] handles the global resonances.** The proof of Prop 2.17 (ℓ. 1332–1343) shifts the
contour past "the resonances of P in −α < Im σ < −α̃". These are the global ones. It then removes
them in r > r₁ (ℓ. 1342):

> Now, by Lemma 2.14 (1), all resonant states of P which are not resonant states of the forward
> problem in Ω₂₃ must in fact vanish in r > r₁, and by part (2) of Lemma 2.14, the only term in
> (EqRNdSPartialAsympProof) that survives upon restriction to r > r₁ is the constant term.

Lemma 2.14(1) (ℓ. 1257) reads:

> Every resonant state u ∈ ker_{H^s_fw} P̂(σ) corresponding to a resonance σ with Im σ > 0 is
> supported in the artificial exterior region {r₀ ≤ r ≤ r₁}; more precisely, every element in the
> range of the singular part of the Laurent series expansion of P̂(σ)⁻¹ at such a resonance σ is
> supported in {r₀ ≤ r ≤ r₁}. In fact, this holds more generally for any σ ∈ ℂ which is not a
> resonance of the forward problem for the wave equation in a neighborhood Ω₂₃ of the black hole
> exterior.

So a global resonance σ ≠ 0 in the strip contributes nothing in r > r₁ exactly when σ is not an Ω₂₃
resonance. The localized hypothesis guarantees this, and P2-4(b) supplies it at α = d_δ. At σ = 0,
Lemma 2.14(2) applies (H4).

**Match.** P2-D1′ identifies the zeros of W_l, with multiplicities, with the resonances of "the
forward problem on Ω₂₃ (HV17 §2.5)" on Im σ > −min(κ₊, κ_c) ⊃ S_δ. So P2-4(b) is literally
[HV17]'s hypothesis with α = d_δ. **No extra argument is needed, and P2-5 is not conditional.** The
plan's stop rule is not triggered.

**What the citation leans on.** These are for the reviewer.
1. **Lemma 2.14(1)'s "more general" clause.**
   - Its proof is given only in outline (ℓ. 1266): "follows along the same lines …; likewise for
     the more precise result".
   - The Laurent-coefficient refinement is likewise outlined: "An iterative argument, similar to
     [BVW, Proof of Lemma 8.3]".
   - This is [HV17]'s lemma. It is cited, not re-proven.
2. **The order-1 claim at 0 for RNdS.** [HV17] takes it from Schwarzschild–de Sitter references.
   At P_B it is also P2-0/D4 with D1′.

### 2.4 [CF17]: hypotheses quoted and matched

**F1, the setting.**
- §2.1 (ℓ. 410–438):
  - "(M, g) … the maximal domain of dependence D(Σ) = M of a compact Cauchy hypersurface
    Σ = S¹ × S²";
  - D(r) = 1 − 2M/r + e²/r² − Λr²/3, with "the mass M > 0, the charge e ≠ 0 and the cosmological
    constant Λ > 0";
  - "The subextremality assumption follows by simply requiring the existence of 3 positive and
    distinct roots 0 < r₋ < r₊ < r_c of the potential D".
- Thm 3.1 [mainThm] (ℓ. 798): "Let (M_int, g) be a black hole interior region of a subextremal
  Reissner–Nordström–de Sitter spacetime. Let φ ∈ C^∞(M_int ∪ H⁺) be a solution of the wave
  equation satisfying Price's Law (Assumption 2.2)."

Match:
- D = F and e = Q.
- P_B is subextremal (H1) [C8.subextremal].
- A solution with smooth data on Σ is smooth on M ⊃ M_int ∪ H⁺. So (ii) is stated for data on a
  Cauchy surface Σ, which is [CF17]'s setting.

**F4, the surface gravities** (eq. surfGrav, ℓ. 454): κ_* = ½|D′(r_*)|. These are FORMULATION's κ's.

**F2, the Price law.** Assumption 2.2 [pricesAss] (ℓ. 760):

> We will say that a function φ ∈ C^∞(M_int ∪ H⁺) satisfies Price's Law, on Reissner–Nordström–de
> Sitter, provided there exists p > 0 such that for any function of the form ψ = T^lΩ^Iφ … we have,
> in H⁺ = {r = r₊}, covered by the Eddington–Finkelstein coordinates (v, r, ω),
> ∫_{v₁}^{v₂}∫_{S²}[(∂_vψ)² + |∇̸ψ|²](v, r₊, ω) dω dv ≲_{l,I} e^{−2pv₁}, for all 0 ≤ v₁ ≤ v₂.

[HV17] gives this with p = α. Restrict the data to a translate of H_I that meets H⁺. Here [CF17]'s
H⁺ is H_A⁺, the component on which v runs to +∞, adjacent to i⁺; it is the component that [HV17]'s
H_I crosses. On H⁺ ∩ D(H_I) = H⁺ ∩ {v ≥ v_I}:
- **Derivative bounds.** T = ∂_t and the rotation fields Ω_i are stationary and tangent to every
  level set of r, in particular to the Cauchy horizon. Thm 1.1 says "The same bound, with a
  different constant C, holds for derivatives of u′ along any finite number of stationary vector
  fields which are tangent to the Cauchy horizon". So |T^kΩ^J u′| ≤ C_{k,J}e^{−αt₀} uniformly in
  r > r₋, which includes H⁺.
- **Only u′ enters.** u₀ is constant, so ψ = T^lΩ^Iφ and its T- and Ω-derivatives see only u′.
  - In [CF17]'s coordinates T = ∂_v ((v_metric), ℓ. 498), so ∂_vψ = T^{l+1}Ω^Iu′.
  - |∇̸ψ|² = r⁻²Σ_i(Ω_iψ)² (roundGrad, ℓ. 739).
- **t₀ against v.** Near H⁺, t₀ − v is a bounded smooth function of r. [HV17] §1 says that t₀ "is
  equivalent … to the Eddington–Finkelstein coordinate v near the event horizon". Concretely,
  t₂₃ = t − F₂₃ with F₂₃′ = −(μ⁻¹ + c₂) near r₂, so t₂₃ = v + (smooth).

Hence the integrand is at most C′e^{−2αv}, and the integral over [v₁, v₂] × S² is at most
(2πC′/α)e^{−2αv₁}. On the compact piece {0 ≤ v ≤ v_I}, the integral is bounded by smoothness, hence
by a constant times e^{−2αv₁}. So Assumption 2.2 holds with **p = α = (11/20)κ₋**.

**F3, the condition of Thm 3.1**, (mainCond): "2 min{p, κ₊} > κ₋".
- 2 min{p, κ₊} > κ₋ ⟺ (2p > κ₋ and 2κ₊ > κ₋).
- **The tautological part.** With p = (½+δ)κ₋, 2p > κ₋ ⟺ δ > 0. This is the tautology that the S1
  review found in `params.json`'s CF17_2alpha [C8.two_p_gt_km: 2p/κ₋ = 11/10]. It carries no
  information beyond the choice p = (½+δ)κ₋.
- **The substantive part.** 2κ₊ > κ₋ is P2-par's CF17_2kp, with margin 0.02747
  (2κ₊/κ₋ = 1.6458) [C8.km_lt_2kp]. At P_B also p < κ₊ (margin 0.01161), so
  min{p, κ₊} = p and 2 min{p, κ₊} = 1.1κ₋.
- **What the CF17 route rests on.** The Price law with p > κ₋/2, i.e. β > ½. That is P2-4, through
  [HV17]'s expansion. So here [CF17] is not independent of [HV17]. It gives a second conclusion by
  the energy method (bounded energy, C⁰ ∩ H¹_loc) from the same spectral input.
- [CF17]'s Assumption 4.3 [energyAss] is internal to the proof: it is derived from Price's law by
  Prop 4.1, with κ̄ = min{κ₊, 2p}. It is not a further hypothesis of Thm 3.1.

### 2.5 Verdict and level

Every hypothesis of [HV17] Thm 1.1 (§§2.5–2.6) and [CF17] Thm 3.1 is matched at P_B. The resonance
notion in [HV17]'s α hypothesis is the localized one, which P2-D1′ provides. P2-5 closes by citation,
at P2-4's level: **WRITTEN_PROOF, conversion cited**.

Two qualifications are part of its statement:
- the data live on [HV17]'s H_I, and for (ii) on a Cauchy surface Σ;
- the regularity is stated near the Cauchy horizon in the domain of dependence of H_I, not
  uniformly up to the bifurcation sphere.

## 3. β at P_B

0.55 ≤ β ≤ 0.90938098. The bound is two-sided but not sharp.
- The lower bound is P2-4 (WRITTEN_PROOF).
- The upper bound is P2-K6: the certified l = 0 zero, which is a resonance at P2-D1′'s level.

P2-6 stays the stretch goal.
