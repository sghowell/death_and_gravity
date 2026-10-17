# P2 S1 report: certified parameters, D1′, P2-0, KAT-1…KAT-7 and P2-2

Status: S1 complete, 2026-10-16. Branch `p2/s1-foundations`, from `p2/s0-formulation`.
- Code commit for every certificate and result: `a78a3057`.
- Conventions: `../FORMULATION.md` v0 (unchanged). Nothing in `docs/problems/*.tex` or
  `CLAIMS.md` was touched.
- S1 serves both outcomes of the pending v1.2 decision: every item was done at both points,
  P_B (Q/Q_max = 999/1000, δ = 1/20) and P_A (249/250, δ = 1/100), with ΛM² = 3/50.

## 0. Summary

| item | outcome | level |
|---|---|---|
| 1. Arb parameters | all frozen constants as 256-bit balls; every needed inequality certified, both points | CERTIFIED (`certificates/params.json`) |
| 2. Lemma D1′ | written proof checked against [Vas13], [HV17], [Hin25] definitions; machine-checked parts C0, C7, C9, C10 | proof note + exact identities |
| 3. P2-0 | exact sympy chain C1–C6 plus a written integration step. **Bonus: D4 proved exactly**, W₀′(0) = i(r₊² + r_c²)/(r₊r_c) | CERTIFIED (`certificates/p2_0_identity.json`) |
| 4. KAT-1…KAT-7 | **all seven pass** their pre-registered criteria | floats (KAT-6 is validated) |
| 5. P2-2 | W_l certified for all 0 ≤ l ≤ 50 at both points (102 certificates) | CERTIFIED (`certificates/lg/`) |
| 6. Tests | 43 tests: 40 fast (~10 s) and 3 gated by `P2_SLOW=1` (full run ~100 s) | — |

The S1 exit criterion of s0-gate §7 ("all KATs pass") is met.

## 1. Certified parameters (item 1)

`src/p2/params.py`, `certificates/params.json` (replay identical: `python -m p2.params verify`).

Method:
- Interval Newton on the integer cubic r³ − 25r + 25 (giving r_ext), then on Δ with ball
  coefficients. Four disjoint enclosures, each containing a unique root, account for all roots of
  the quartic.
- Independent cross-check: each horizon ball meets exactly one certified root of the irreducible
  degree-12 resultant Res_s(s³ − 25s + 25, Δ(r; Q² = (c²/2)s(3−s))). The reduction
  Q_max² = r_ext(3 − r_ext)/2 is exact.

Midpoints are printed below to 13 digits; the stored ball radii are below 10⁻⁷⁰.

| | P_B | P_A |
|---|---|---|
| Q | 1.009843167712 | 1.006810605647 |
| r_n, r₋, r₊, r_c | −7.960335971229, 0.9973665935631, 1.094348347310, 5.868621030356 | −7.960063333143, 0.9492759837530, 1.143117525527, 5.867669823862 |
| κ₋, κ₊, κ_c | 0.04254199481802, 0.03500736272070, 0.09338246410782 | 0.09426070199554, 0.06379974180852, 0.09332620731341 |
| r_P | 1.958730370969 | 1.971873721151 |
| γ₀/κ₋ | 1.683211531268 | 0.7643026659779 |
| Ω_P/κ₋ | 4.918835254757 | 2.209566241242 |
| strip depth (½+δ)κ₋ / γ₀ | 0.3267563165906 | 0.6672749196125 |
| ρ = κ₋/κ₊ | 1.215229926271 | 1.477446449211 |
| ϰ (Hintz), λ₁⁺, λ₂⁺ | 0.8587686195, 2.105900378, 3.190120458 | 0.8485391579, 2.114618219, 3.205730211 |

Certified inequalities, with margins in κ₋ units:

| inequality | P_B | P_A |
|---|---|---|
| (½+δ)κ₋ < κ₊ (Jost / E₊ margin) | 0.2729 | 0.1668 |
| (½+δ)κ₋ < κ_c | 1.645 | 0.4801 |
| (½+δ)κ₋ < γ₀ (HV17's α < γ₀) | 1.133 | 0.2543 |
| κ₋ < 2κ₊ (CF17) | 0.6458 | 0.3537 |

Also certified: the ordering r_n < 0 < r₋ < r₊ < r_c; r₊ < r_P < r_c; F(r_P) > 0; and ρ < 2
(Ros25). [CGNS18]'s ρ < 9/7 is certified true at P_B and certified **false** at P_A (context only).
The S0 float table of FORMULATION §1.3 rounds correctly at every printed digit
(`tests/test_params.py`). Cost: 0.2 s CPU.

## 2. Lemma D1′ (item 2)

`notes/s1-p20-d1prime.md` §4. Definitions were read from the TeX sources:
- [Vas13]: Thm 2.6, Remark 2.8, §2.6, Lemma 3.1.
- [HV17]: Thm 2.13, Lemma 2.14, Prop 2.17, Thm 2.24, and the §2.5 hypothesis sentence. "Simple"
  means pole order 1 and rank 1.
- [Hin25]: the QNM definition in §1, and the multiplicity remark after Thm 5.5.

Result: for σ ∉ E₊ ∪ E_c, σ is a resonance iff some W_l(σ) = 0. The resonant states are
u_l(r)Y_lm, and the geometric multiplicity is Σ(2l+1). For Im σ > −min(κ₊, κ_c):
- the pole order is max_l m_l(σ);
- each (l, m) sector has exactly one Jordan chain, of length m_l(σ);
- the space of generalized resonant states has dimension Σ(2l+1)m_l(σ).

The time functions of the three sources differ from ours by smooth f(r), so all three notions
coincide.

Machine-checked parts (exact sympy):
- C0: separation of □_g;
- C1: hyperboloidal transform;
- C7: horizon residues of T′ and exponents {0, iσ/κ};
- C9: Abel/Wronskian identity;
- C10: residue of p₁/p₂, i.e. Abel weight ~ |r − r_h|^{−iσ/κ_h}.

**Scope note.** The multiplicity part is proved on Im σ > −min(κ₊, κ_c), where the Abel integrals
converge absolutely. That half-plane contains S_δ at both points (certified margins above), so this
is all P2 uses. A Hadamard finite-part extension to all σ ∉ E is sketched but not claimed.
§4.4 of the note checks that P2-4 together with D4 and the certified inequalities is exactly
[HV17]'s quoted hypothesis.

## 3. P2-0 and D4 (item 3)

`src/p2/p20_identity.py` and `certificates/p2_0_identity.json`. The chain is C0–C10, 12 exact
identities. Each has a **negative control**: a deliberately perturbed version must fail, and all
12 do. Written steps: `notes/s1-p20-d1prime.md` §2–3. Cost: 0.5 s CPU.

- **P2-0.** For every l and every ω with Im ω ≥ 0: W_l(ω) = 0 iff (l, ω) = (0, 0).
  - The proof uses C2–C5 to reach the identity b(D + LN + |ω|²𝒲) + |ω|²B = 0. This is s0-gate's
    (E) multiplied by b, so there is no division by b.
  - The weight 4r⁴/((Λ/3)(r_c − r₊)²(r − r_n)(r − r₋)) > 0 comes from C6, and u₊(r₊) ≠ 0 from C7.
  - The argument uses only that Δ has four simple real roots with r₊ < r_c and Λ > 0. So it holds
    at every non-degenerate RNdS point, not only at P_A and P_B.
- **D4 (stronger than FORMULATION asked at S1).** W₀′(0) = i(r₊² + r_c²)/(r₊ r_c) ≠ 0, so the zero
  mode is simple. This is exact: the C8 limits are machine-checked.
  - Cross-check: a central difference of validated Arb values agrees to < 10⁻³⁰ at both points
    (`tests/test_arbw.py`).
  - FORMULATION planned D4 as a by-product of the S3 l = 0 winding number; that is now redundant
    for D4 itself, though still needed for P2-3.

## 4. Known-answer tests (item 4)

Criteria were written in the `src/p2/kats.py` docstring **before** any run and were not changed.
All seven pass. Results are in `results/kat_float_{12,345,7}.json` and
`certificates/kat6_winding_krawczyk.json`, all at commit `a78a3057`.

| KAT | criterion (abridged) | outcome | key numbers | CPU |
|---|---|---|---|---|
| 1 | CCDHJ18 Table I, all 9 values, each component rounds to the printed 10 decimals | **PASS 9/9** | max component deviation 4.5×10⁻¹¹ (rounding); solver self-consistency < 10⁻⁴⁶; l = 100: 222.06029002491898 − 0.76430944457316i | 81 s |
| 2 | eikonal: D_l = (l+½)²(−Im ω/κ₋ − γ₀/κ₋) and E_l = (l+½)(Re ω/κ₋ − (l+½)Ω_P/κ₋) within 5% of their l = 100 values; l = 100 digits | **PASS** | D_l = 0.068525, 0.068481, 0.068467, 0.068466 and E_l = −0.11220 … −0.11228 for l = 10, 20, 50, 100; predicted 222.0614 against 222.0603 | 33 s |
| 3 | CCDHJ18 Fig. 3 at ΛM² = 0.06, l = 10 proxy: β = ½ crossing rounds to 0.9920, NE takeover rounds to 0.9967 | **PASS** | q₁ = 0.99202246, q₂ = 0.99673082; at q₁: PS10 0.5, NE0 0.826, dS1 0.974 (PS dominant) | 86 s |
| 4 | Hin25 Thm 1.1 at q ∈ {0.999, 0.9999, 0.99999}: errors strictly decrease and e(0.99999) ≤ ½e(0.9999) | **PASS** | l=0 n=0: 0.0906, 0.0324, 0.0106; l=0 n=1: 0.145, 0.061, 0.021; l=1: 0.194, 0.065, 0.021; l=2: 0.299, 0.100, 0.032. Each ratio is ≈ 0.32–0.34, i.e. e ∝ (1−q)^{1/2} ∝ r₊ − r₋ | 56 s |
| 5 | near-Nariai at ΛM² = 0.14 (supplement formula), ε = (r_c − r₊)/r_N ∈ {0.04, …, 0.005}: monotone, last ratio ≤ 0.65 | **PASS** | l=1: 0.0148 → 0.0020; l=2: 0.0259 → 0.0035; ratios 0.506, 0.505 (first order). The printed closed form of Υ agrees with r_N/(r_N − 3/2); it passes through complex intermediates | 11 s |
| 6 | validated winding = 1 on small boxes; Krawczyk boxes contain the float values | **PASS** | P_A l=1 PS: winding [1.00000 ± 1.5×10⁻⁷]; unique zero in 3.2426164126110423036862081 ± 5×10⁻²⁶ − (0.7958326323112821365072551 ± 4×10⁻²⁶)i (κ₋ units), within 1.1×10⁻¹¹ of Table I. P_B l=0 NE: winding [1.000000000 ± 2.7×10⁻¹¹]; unique zero −0.90938097899480545552493376262i ± 6×10⁻³⁰ | 11 s |
| 7 | independent second method agrees to ≥ 10 digits on the 8 P_B modes of §1.3 (and converged in N) | **PASS** | hyperboloidal Chebyshev collocation (log-stretched variable, N = 60–120, acb Newton) vs Frobenius/Taylor: relative differences 4×10⁻⁴¹ to 2×10⁻⁴⁰ (≈ 40 digits) on all 8 | 42 s |

KAT-6 method note:
- The validated Wronskian (`arbw.py`) costs about 0.1 s per point at 256 bits, with relative
  radius about 10⁻⁶¹.
- Interval evaluation over an ω-box overestimates by a roughly constant factor of **about 4×10⁷**
  (measured at P_A, l = 1). That is the wrapping effect of ball recurrences, and it makes naive
  segment enclosures useless.
- The winding therefore uses a rigorous Cauchy-interpolation Taylor model. One box evaluation gives
  a crude sup bound M; 64 validated point values on a circle give DFT coefficients with an explicit
  aliasing bound; there is an explicit truncation bound. The 4×10⁷ is crushed by (ρ′/ρ)^64.

KAT-7 by-product (float evidence): the converged spectral eigenvalues at P_B with
Im ω > −0.6κ₋, for l = 0, 1, 2, 10, consist of ω = 0 (l = 0) only.

## 5. P2-2 thresholds (item 5)

Route N2 (s0-gate §3): `src/p2/lg_bounds.py`, proof in `notes/s1-p2-2-lg.md`. Replay:
`python -m p2.lg_bounds verify PB 50` (≈ 3 s per point).

**Statement certified for each X ∈ {P_B, P_A} and 0 ≤ l ≤ 50:** W_l(ω) ≠ 0 for all ω with
|Re ω| ≥ W_l and Im ω ≥ −(½+δ)κ₋, with W_l the exact dyadic in `lg_X_l{l}.json`.
- Im ω ≥ 0 is P2-0, and Re ω < 0 follows by the exact symmetry W_l(−ω̄) = conj W_l(ω).

Method:
- The Liouville–Green Volterra equations are referenced at the horizons. The QNM-side Jost
  solutions are dominant there for Im ω < 0. The kernel weights exp(2∫|Im p|) are integrable
  because the strip depth is below κ₊ and κ_c; the end-cell exponents s = βG/κ_h are certified < 1
  (0.67 and 0.25 at P_B; 0.75–0.77 and 0.52 at P_A).
- Majorants valid uniformly for a ≥ W_l, b ∈ [−β, 0].
- Upper Riemann sums on about 2.2k r-cells per side, using the exact tortoise antiderivative
  r* = ΣA_i ln|r − r_i| for the weights, plus analytic end cells.
- Criterion (*): e^{Ψ₊+Ψ_c}(1 + (Ψ̃₊+Ψ̃_c)/2) < 2. It holds in all 102 certificates, with values
  in [1.947, 2).

| l | 0 | 1 | 2 | 5 | 10 | 20 | 30 | 50 |
|---|---|---|---|---|---|---|---|---|
| P_B: W_l/κ₋ | 5.401 | 10.31 | 15.61 | 31.19 | 56.57 | 106.6 | 156.4 | 255.5 |
| P_B: W_l/((l+½)Ω_P) | 2.196 | 1.398 | 1.269 | 1.153 | 1.095 | 1.058 | 1.042 | 1.029 |
| P_A: W_l/κ₋ | 2.695 | 4.951 | 7.411 | 14.49 | 25.96 | 48.56 | 70.99 | 115.6 |
| P_A: W_l/((l+½)Ω_P) | 2.439 | 1.494 | 1.342 | 1.192 | 1.119 | 1.072 | 1.053 | 1.036 |

- W_l is below the planning value 1.25(l+½)Ω_P + 2κ₋ of s0-gate §3 for every l, except P_B, l = 0
  (5.40 against 5.07 κ₋). The part-(iii) boxes R_l are therefore up to 18% narrower than planned
  at large l.
- Empirically W_l/((l+½)Ω_P) − 1 ∝ (l+½)^{−0.75} at both points.

Float check of the whole chain (`tests/test_lg.py`): at Re ω = W_l and Im ω = −β, the actual
|W_l e^{iΘ}/(2iω) − 1| is about 0.01, against the certified bound of about 0.98. The bound is
conservative but correct.

**Cost:** 1.0 s CPU per (point, l), threshold search included. That is 104.6 s CPU for all 102
certificates, and about 0.05 s per l to replay. Projection: l ≤ 500 at both points would take
about 17 min CPU. The cost is flat in l because the cell count is fixed.

## 6. Cost summary

All runs used at most 2 cores, under `nice -n 10`. Production CPU (final regeneration at
`a78a3057`):

| item | CPU |
|---|---|
| params | 0.2 s |
| P2-0/D4/D1′ identities | 0.5 s |
| P2-2 (102 certificates) | 105 s |
| KAT-1 | 81 s |
| KAT-2 | 33 s |
| KAT-3 | 86 s |
| KAT-4 | 56 s |
| KAT-5 | 11 s |
| KAT-6 | 11 s |
| KAT-7 | 42 s |
| **total** | **≈ 7 min** |

The tests take about 10 s (fast) and about 100 s (`P2_SLOW=1`). Development runs added roughly 30
min of CPU. All of this is far below the 6-hour budget.

## 7. Findings bearing on the v1.2 decision and on S2/S3

1. **Nothing found in S1 argues against v1.2.**
   - At P_B, the two independent float methods agree to ~40 digits on the whole §1.3 spectrum.
   - The spectral strip scan finds nothing above −0.6κ₋ except ω = 0.
   - The dominant l = 0 NE mode now has a *certified* enclosure,
     −0.90938097899480545552493376262i·κ₋ (± 6×10⁻³⁰). That is the Krawczyk half of P2-6 at the
     shallowest mode.
   - At P_A, all nine Table I values and the eikonal trend reproduce, with O(l⁻²) in Im and
     O(l⁻¹) in Re.
2. **D4 is now exact.** The HV17 hypothesis "0 is simple with constant resonant states" reduces to
   the strip exclusion alone (P2-1, P2-2, P2-3), via P2-0 + D4 + D1′.
3. **For S3 (part iii):**
   - Naive interval evaluation of the D-finite recurrences over ω-boxes overestimates by about
     4×10⁷. Validated winding therefore needs ω-Taylor models (as FORMULATION §4 plans) or the
     Cauchy-interpolation scheme of `winding.py`. The latter needs only point evaluations
     (≈ 0.1 s each at small l) plus one crude box bound per disc. With K = 64, a disc costs about
     6.5 s at l ≤ 2.
   - A rough estimate for l = 10 at P_B is about 160 discs, or ≈ 17 min per l on one core. This is
     comparable to s0-gate §5's tile plan, not better. The disc radius could likely be increased;
     S3 should measure.
4. **For S2 (part i):** the weighted-Volterra LG machinery here is already the above-barrier regime
   (regime 1 of s0-gate §4) with explicit constants. The certified threshold approaches the barrier
   top like (l+½)^{−0.75}, close to the Airy h^{2/3} scaling. The genuinely new work is still the
   barrier-top layer and the sub-barrier/low-frequency regimes, as the gate memo says.
5. **D1′ scope.** The multiplicity statement is proved on Im σ > −min(κ₊, κ_c), which contains S_δ
   at both points. It is not proved on all of ℂ∖E (FORMULATION D1′'s wording). Nothing in P2 needs
   more; a revision-log note may record the precise scope.
6. **Hintz's printed T̃** lacks the −1 (already noted in the digest). It does not affect anything,
   because resonances do not depend on the smooth choice of t*.

## 8. Draft CLAIMS rows (for the coordinator; not applied)

| id | problem | statement | level | evidence | date |
|---|---|---|---|---|---|
| P2-0 | P2 | Massless minimally coupled scalar on any non-degenerate RNdS (Λ > 0; Δ with four simple real roots r_n < 0 < r₋ < r₊ < r_c), with the radial Wronskian W_l of FORMULATION v0 §1.2: for every l ≥ 0 and every ω with Im ω ≥ 0, W_l(ω) = 0 iff (l, ω) = (0, 0). Moreover the zero of W₀ at 0 is simple: W₀′(0) = i(r₊² + r_c²)/(r₊r_c) (D4) | CERTIFIED (exact sympy identity chain C1–C8, each with a failing negative control, plus the written integration step; D4 limits exact) | problems/P2/src/p2/p20_identity.py; certificates/p2_0_identity.json; notes/s1-p20-d1prime.md §2–3; tests/test_p20_identity.py, tests/test_arbw.py | 2026-10-16 |
| P2-2 | P2 | At P_B (ΛM² = 3/50, Q = (999/1000)Q_max, δ = 1/20) and at P_A (Q = (249/250)Q_max, δ = 1/100), for every 0 ≤ l ≤ 50, with W_l the exact dyadic stored per (point, l): W_l(ω) ≠ 0 for all ω with \|Re ω\| ≥ W_l and Im ω ≥ −(½+δ)κ₋. W_l/((l+½)Ω_P) decreases from 2.196 (l=0) to 1.029 (l=50) at P_B and from 2.439 to 1.036 at P_A; e.g. W_10 = 56.57κ₋ (P_B), 25.96κ₋ (P_A) | CERTIFIED (Arb 128-bit: LG Volterra bounds with horizon-referenced weights e^{2∫\|Im p\|}, uniform majorants in a ≥ W_l, b ∈ [−β, 0], upper Riemann sums with the exact tortoise antiderivative and analytic end cells; criterion e^{Ψ₊+Ψ_c}(1+(Ψ̃₊+Ψ̃_c)/2) < 2 in all 102 certificates; Im ω ≥ 0 by P2-0; Re ω < 0 by the exact symmetry). Parameter balls from params.json | problems/P2/src/p2/lg_bounds.py; certificates/lg/lg_{PB,PA}_l{0..50}.json; notes/s1-p2-2-lg.md; tests/test_lg.py | 2026-10-16 |

Suggested companion rows, if the coordinator wants them:
- **P2-par**: the params.json statements (certified inequalities, CERTIFIED).
- **P2-KAT**: KAT-1…KAT-7 pass. CONJECTURED/evidence level, except KAT-6, which is a certified
  enclosure of two individual QNMs.
- **P2-D1′**: the lemma, at "proof note + machine-checked parts". It needs the independent review
  before any level is assigned.

## 9. Reproduction (cwd `problems/P2/src`, python = repo `.venv`)

```
python -m p2.params verify          # params.json replay                         (~0.3 s)
python -m p2.p20_identity verify    # P2-0 / D4 / D1' identity chain replay        (~1 s)
python -m p2.lg_bounds verify PB 50 # P2-2 certificates replay (also PA)          (~3 s each)
python -m p2.winding                # KAT-6 (rewrites the certificate)            (~12 s)
python -m p2.kats 12 | 345 | 7      # float KATs (rewrite results/)               (1.9 / 2.5 / 0.7 min)
python -m pytest problems/P2/tests -q                    # from the repo root, ~10 s
P2_SLOW=1 python -m pytest problems/P2/tests -q          # ~100 s
```
