# P2 S5 plan: the assembly P2-4, the conversion P2-5 by citation, FORMULATION v1.0

Status: written 2026-10-26 on branch `p2/s5-assembly` (from main at `3dd0b2bc`), before any S5 code.
The user approved S5 on 2026-10-26. Any later change to this plan goes in a dated **Amendments**
section at the end, never by silent edit.

Notation as in S4. d_δ := (½+δ)κ₋ = (11/20)κ₋ at P_B. β is D3's gap α/κ₋. The decay rate given to
[HV17] is α := d_δ.

## 1. Deliverables

1. `src/p2/assembly.py`, run as `python -m p2.assembly` (write) or `python -m p2.assembly verify`
   (recompute and compare). It writes `certificates/assembly/assembly_PB.json`.
2. `tests/test_assembly.py`: a fast test that replays the check and every negative control.
3. `notes/s5-assembly.md`:
   - the written proof of P2-4 at P_B;
   - the P2-5 citation analysis, quoting the TeX sources with their labels: [HV17]
     arXiv:1512.08004 (final revision of 2017-07-12) and [CF17] arXiv:1607.01018.
4. `notes/s5-report.md`:
   - the verdict;
   - paste-ready CLAIMS rows P2-4 and P2-5 (not applied);
   - a proposed P2 status entry for `docs/problems/open-problems-theoretical-cosmology-2026.tex`,
     as text only;
   - what S5 does not establish.
5. FORMULATION v1.0: the header and a revision-log entry for S5.

Out of scope: P_A, P2-6, and any compute beyond replays and the assembly check.

## 2. P2-4: what is assembled

S_δ = U ∪ H₊ ∪ H₋, where U = {Im ω > 0} and H± = {±Re ω ≥ 0, −d_δ ≤ Im ω ≤ 0}.

| piece | degrees | source (level on main `3dd0b2bc`) |
|---|---|---|
| U | every l | P2-0: W_l(ω) = 0 with Im ω ≥ 0 iff (l, ω) = (0, 0) (WRITTEN_PROOF) |
| H₋ | every l | the symmetry W_l(−ω̄) = conj W_l(ω) (FORMULATION §1.2; `zeq_identity` Z4) maps it onto H₊ |
| H₊ | l ≥ 20 | P2-1: `cover_PB.json`, h = 1/(l+½) ≤ h₀ = 2/41 (WRITTEN_PROOF) |
| H₊ ∪ H₋ | l ≤ 19 | R′_l = [−a_l, a_l] × [−b, e] (P2-3B: count 1 at l = 0, else 0; CERTIFIED), plus {abs(Re ω) ≥ a_l} × [−d_δ, 0] (P2-2; CERTIFIED) |

The other ingredients:
- **Analyticity.** W_l is analytic on S_δ, and D1's notion of QNM applies there, because
  E₊ ∪ E_c lies below S_δ (P2-par: d_δ < κ₊, κ_c; CERTIFIED).
- **The zero at 0.** W₀(0) = 0 by D4 (R ≡ 1). The zero is simple by D4, and independently by the
  R′₀ count of 1. For l ≥ 1 there is no zero at 0 (P2-0).
- **Resonance form.** P2-D1′ applies, since S_δ ⊂ {Im ω > −min(κ₊, κ_c)} (WRITTEN_PROOF). So the
  Ω₂₃ forward problem has no resonance in S_δ∖{0}, and 0 is a resonance of multiplicity 1 with a
  constant resonant state.
- **β.** D3 gives β ≥ ½ + δ = 0.55. With P2-K6, 0.55 ≤ β ≤ 0.90938098.

**Level rule.** A union of exclusions is no stronger than its weakest member. The minimum over the
ingredients is WRITTEN_PROOF (P2-0, P2-D1′ and P2-1). Each of them is needed for every l (U and
the resonance form) or for infinitely many l (P2-1). The expected level is therefore WRITTEN_PROOF,
pending S5's independent review.

## 3. The assembly check

**Inputs.** Committed certificates only:
- `params.json`;
- `lg/lg_PB_l{0..19}.json` and `part3/winding_PB_l{0..19}.json`;
- `part1/cover_PB.json`, plus the regime certificates through a replay of `part1.cover()`;
- `p2_0_identity.json`, `part3/zeq_identity.json` and `part3/krawczyk_PB_l0.json`.

**Arithmetic.** Dyadics and rationals are exact Fractions. Arb balls at 256 bits are parsed from
`params.json`. `Point("PB")` is recomputed, and its balls must overlap the stored ones.

**Checks.** Each check is a named boolean, with a margin where one applies.
- **C1, point.** qfrac 999/1000, δ = 1/20 and Λ = 3/50. The recomputed balls overlap
  `params.json`.
- **C2, l ≤ 19, for each l.**
  - Inputs:
    - a certificate exists for every l in 0..19;
    - the lg claim is certified (criterion upper bound < 2), and its point, l and strip depth
      match.
  - The box R′_l:
    - its a equals lg's `W_l_exact` exactly;
    - b ≥ d_δ (Arb), e > 0 (exact), and b < κ₊ (analyticity on the box).
  - The stored path:
    - it traces the right half of ∂R′_l exactly, contiguous from (0, −b) through (a, −b) and
      (a, e) to (0, e);
    - each piece lies within ρ/2 of its centre;
    - each disc box clears E₊ ∪ E_c.
  - The count:
    - the sum of the increments, divided by π, lies within 1/4 of the stored winding;
    - the winding is 1 at l = 0 and 0 for l ≥ 1.
  - Containment: [−d_δ, 0] ⊂ [−b, e], and the box and the P2-2 region meet at abs(Re ω) = a_l with
    no gap.
- **C3, l ≥ 20.**
  - The replayed cover equals the committed `cover_PB.json`.
  - h₀ = 2/41 gives l₀ = ⌈1/h₀ − ½⌉ = 20 = 19 + 1, so no l is missed.
  - For every h ∈ (0, h₀], the regime ranges in Re ω union to [0, ∞) (exact and Arb):
    - L is [0, 8κ₋];
    - S is [6κ₋, ζ_SU·Ω_P/h];
    - T is [ζ_TL, ζ_TU]·Ω_P/h;
    - A is [ζ_A·Ω_P/h, ∞).
  - The Im ω range is [−d_δ, 0], with d_δ/κ₋ = 11/20.
- **C4, poles.** d_δ < κ₊ and d_δ < κ_c. This also puts S_δ inside D1′'s multiplicity window.
- **C5, the zero mode.**
  - 0 lies in the interior of R′₀ (exact).
  - D4's W₀′(0) = i(r₊² + r_c²)/(r₊r_c) is nonzero (Arb).
  - Every step of the P2-0 identity chain is recorded ok.
- **C6, β.**
  - The lower bound is 11/20 (exact).
  - The upper bound comes from `krawczyk_PB_l0.json`: abs(Im ω_*)/κ₋ ≤ 0.90938098, and ω_* lies
    outside S_δ (a consistency check).
- **C7, HV17 numerics.**
  - Λ, M and Q are positive.
  - F has exactly three positive roots, all simple (κ's > 0, r_n < 0).
  - 0 < α < γ₀.
  - The exponent ½ + α/κ₋ = 21/20 exceeds 1.
- **C8, CF17 numerics.** The same subextremality holds. With p := α:
  - min{p, κ₊} = p, because p < κ₊;
  - 2p > κ₋ holds exactly (11/10 > 1). This is the tautology noted in P2-par's review;
  - κ₋ < 2κ₊ (Arb) is the non-tautological part.
- **C9, levels.** The ingredient levels as on main `3dd0b2bc`, and their minimum. The test
  cross-checks this table against `CLAIMS.md`.

**Output.** `certificates/assembly/assembly_PB.json`, with the sha256 of every input, the commit,
the precision and the timing.

**Negative controls.** Each must make the check fail on the named item. They are recorded in the
JSON and re-run by the test.
- **N1, b just below d_δ.** At l = 0, b → ⌊d_δ·2²⁸⌋/2²⁸, with the path coordinates moved with it.
- **N2, one l missing.** `winding_PB_l7.json` is dropped.
- **N3, a gap in Re ω.** In R′₅, a₅ → a₅ − 2⁻²⁴, with the path moved with it.
- **N4, an l gap.** h₀ → 2/43, so l₀ = 21.
- **N5, e = 0** at l = 0.
- **N6, a wrong count.** Winding 1 is stored at l = 3.
- **N7, a failed regime.** T_PB's pass flag is false in the cover.
- **N8, the hypotheses are not vacuous.** γ₀ → 0.5κ₋ must fail HV17's α < γ₀, and κ₊ → 0.45κ₋
  must fail CF17's κ₋ < 2κ₊.

## 4. P2-5: hypotheses to match

The quotes, with their labels, go in `s5-assembly.md`.

[HV17]:
- **H1, setting.** Thm 1.1 [ThmIntroMain], Def 2.1 [DefRNdSNonDegenerate] and (EqRNdSMetric),
  against FORMULATION §1.1.
- **H2, field and data.**
  - The equation □_g u = 0, with smooth Cauchy data on H_I (EqRNdSMfdHyper).
  - The conclusions hold in r > r₁, the domain of dependence of H_I.
  - How this relates to "smooth Cauchy data on RNdS" in FORMULATION P2-5.
- **H3, the resonance notion in the α hypothesis (§2.5 [SubsecRNdSAsymp]).** Is it the forward
  problem on Ω₂₃ (localized), or the global P̂(σ)⁻¹? To be checked against Lemma 2.14
  [LemmaRNdSResonanceSupp] (1)–(2) and the proof of Prop 2.17 [PropRNdSPartialAsymp].
- **H4, 0 simple with constant resonant states.** HV17's own claim (Lemma 2.14(2) and its
  references), against P2-0/D4 together with D1′.
- **H5, α < γ₀.** Prop 2.12 [PropRNdSSemiclassical2] and the γ₀ formula of §2.7, against P2-par
  and D6.
- **H6, the exponent.**
  - β₁ = 1/κ₁ (EqRNdSRadialPointBeta, EqRNdSSurfGrav).
  - u′ ∈ e^{−αt₀}H^{1/2+α/κ₁−0} = H^{1.05−0}, by Thm 1.1, Thm 2.24
    [ThmRNdSPartialAsympConormal] and Cor 2.19 [CorRNdSBoundedness].
- **H7, implicit hypotheses**, for example the contour line Im σ = −α in the proof of Prop 2.17.

[CF17]:
- **F1, setting.** Subextremal RNdS (§2), φ ∈ C^∞(M_int ∪ H⁺) and □_gφ = 0.
- **F2, Assumption 2.2 [pricesAss].** Does HV17's expansion give it, and with which rate p?
- **F3, Thm 3.1 [mainThm].** In the condition 2 min{p, κ₊} > κ₋, which part is a tautology, and
  which part is checked?
- **F4, κ conventions** (eq. surfGrav).

## 5. Acceptance criteria

- `python -m p2.assembly` passes on the committed certificates.
- Every negative control fails on its named item.
- The fast test passes, and the P2 fast suite stays green.
- `s5-assembly.md` proves P2-4 from the listed rows only, and states its level and why.
- Every item in §4 is:
  - quoted with its label;
  - matched to a P2 row, a `params.json` inequality, or an explicit argument.

  The resonance notion is settled by quotation.
- The FORMULATION v1.0 entry changes no frozen convention.
- The ledger and .tex texts are proposed, not applied.
- Only `problems/P2/` changes. Commits carry no AI attribution.

## 6. Stop rule

- Suppose HV17's "0 is the only resonance in Im σ ≥ −α" refers to the global P̂(σ)⁻¹, and P2-4
  together with Lemma 2.14 does not close the gap. Then stop and report. P2-5 may then be stated
  only as explicitly conditional. It is never weakened silently.
- Suppose the assembly check finds a gap in the committed certificates, such as an a_l mismatch, a
  missing l, or b < d_δ. Then stop and report. S5 regenerates no certificate.

## 7. Budget

The assembly check and its tests take seconds to a minute of CPU, at nice 10 with at most 4
workers. That is far below the 2 core-h allowed.
