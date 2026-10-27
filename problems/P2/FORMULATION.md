# P2 — A certified resonance-free strip for subextremal RNdS (β > 1/2, linear SCC failure for smooth data)

Status: formulation **v1.1** (S6 start), 2026-10-27; frozen at v0 (S0, 2026-10-16),
with every later change recorded in the Revision log.
- The ledger has rows P2-par, P2-0, P2-D1′, P2-1w, P2-1, P2-2, P2-3B, P2-4, P2-5 and
  P2-K6 (`CLAIMS.md`).
- S5 assembled P2-4 and P2-5 at P_B (`notes/s5-assembly.md`).
- The independent review (`docs/review-2026-10-26-p2-s5-adversarial.md`, 2026-10-26)
  found both proven as stated, with nothing blocking. Its fixes S1–S3 and W1–W9 are
  applied, and its follow-up check (§9) passed (2026-10-26).
- P2-4 and P2-5 at P_B (`notes/s5-report.md` §6, the review's §7 rows) were entered
  on 2026-10-26, with the user's approval of that day.
- P2-6 (the stretch goal) and P2-7 (P_A) are still *targets*. S6 (P2-6 at P_B) is in
  progress (`notes/s6-preregistration.md`; Revision log v1.1).

Changes after S1 begins go in the Revision log, never by silent edit (house
discipline of `problems/P9/b/FORMULATION.md`).

Pre-registered statement: `docs/problems/open-problems-theoretical-cosmology-2026.tex`
§P2, **v1.2** (approved 2026-10-17). The S0 gate (`notes/s0-gate.md`)
recommended **GO with a modified entry point** and proposed the v1.2 entry
(`s0-gate.md` §7). It was approved with one change: the two-sided enclosure of
β (P2-6) is a stretch goal, not part of the entry point. So P_B is the primary
target and P_A the second (§1.3).

Literature basis: `notes/literature-digest.md`, whose bracketed keys are used
below. Float numbers quoted here come from `src/s0_sanity.py` (mpmath, 40
digits). They are not certified.

## 1. Frozen conventions

### 1.1 Units, metric, field model

Units are G = c = 1 and M = 1. The RNdS metric is
g = −F dt² + F⁻¹dr² + r²dΩ², with
F(r) = 1 − 2M/r + Q²/r² − Λr²/3 = Δ(r)/r² and
Δ(r) = −(Λ/3)r⁴ + r² − 2Mr + Q².

Δ has four real roots r_n < 0 < r₋ < r₊ < r_c (Cauchy, event and cosmological
horizons), with r_n = −(r₋+r₊+r_c). Surface gravities are
κ_* = ½|F′(r_*)| = |Δ′(r_*)|/(2r_*²) for * ∈ {−, +, c}, as in [CCDHJ18] eq. (4).

**Field model (frozen).** A massless, minimally coupled scalar field, □_gφ = 0,
on the fixed subextremal RNdS exterior r₊ < r < r_c. The exterior is extended
across the future event and cosmological horizons by a hyperboloidal time
t* = t − T(r), with T′ = T̃/F, T̃(r₊) = −1, T̃(r_c) = +1 ([Hin25] §1.1, with the
−1 normalization restored).

Justification:
1. The entry point is the *linear*, smooth-data statement ("the first
   theorem-grade statement of linear strong-cosmic-censorship violation for
   smooth data").
2. The two linear conversion theorems, [HV17] and [CF17], are proven for exactly
   this field on fixed RNdS (§3).
3. The only high-precision published known answers are for this field:
   [CCDHJ18] Table I and Fig. 3.
4. After separation, the radial equation is a linear ODE with polynomial
   coefficients. It is Fuchsian, with regular singular points r_n, r₋, r₊, r_c
   and ∞. Every certificate is therefore a D-finite computation on exactly
   defined algebraic coefficients.

What other work adds:
- **[DRS18]** (arXiv:1808.02895) is the paper that treats coupled
  gravito-electromagnetic perturbations of RNdS. It uses the Kodama–Ishibashi
  master variables Φ_±, which are isospectral across vector and scalar type,
  and finds a *stronger* violation: arbitrarily smooth extensions for large
  near-extremal black holes. That model is the natural successor of P2 and is
  out of scope here (§7).
- **[DERS18]** (arXiv:1801.09694) treats Kerr–de Sitter, not RNdS. It is cited
  by the .tex for "Kerr–de Sitter appears safe".

### 1.2 Separation, exponents, the Wronskian

Write φ = e^{−iωt}R(r)Y_lm and L = l(l+1). Then

    (Δ R′)′ + (ω² r⁴/Δ − L) R = 0,          (R)

which is equivalent to Φ″ + (ω² − V_l)Φ = 0 in r* (dr* = dr/F), with Φ = rR and
V_l = F(L/r² + F′/r) ([CCDHJ18] eqs. 6–7).

The indicial exponents of (R) are ±iω/(2κ₊) at r₊ and ±iω/(2κ_c) at r_c. The
QNM boundary conditions of [CCDHJ18] eq. (8) select:
- **R₊(r; ω, l) = (r − r₊)^{−iω/(2κ₊)} g₊(r)**, with g₊ analytic at r₊ and
  g₊(r₊) = 1;
- **R_c(r; ω, l) = (r_c − r)^{−iω/(2κ_c)} g_c(r)**, with g_c analytic at r_c and
  g_c(r_c) = 1.

R₊ is defined and analytic in ω off E₊ = {−inκ₊ : n ≥ 1}, and R_c off
E_c = {−inκ_c : n ≥ 1}.

**Wronskian (frozen):**

    W_l(ω) := Δ(r)·(R₊ R_c′ − R₊′ R_c)·exp(−iω c₊ + iω c_c)/(r₊ r_c).

- W_l is independent of r and analytic on ℂ∖(E₊ ∪ E_c).
- c₊ = Σ_{i≠+} A_i ln|r₊ − r_i| and c_c = Σ_{i≠c} A_i ln|r_c − r_i|, where
  A_i = r_i²/Δ′(r_i) are the coefficients in r* = Σ A_i ln|r − r_i|.
- The exponential factor is entire and never zero. It normalizes R₊ and R_c to
  Jost solutions, Φ ~ e^{∓iωr*}, without moving any zero. A certificate may
  drop it, provided it says so.

### 1.3 Frozen parameter points

Set M = 1 and Λ = 3/50 (ΛM² = 0.06). Let r_ext be the smallest positive root of
r³ − 25r + 25 = 0, which is the double-root condition 2Λr³/3 − r + 1 = 0 at
Λ = 3/50. Set Q_max² = r_ext − r_ext⁴/50. Both points are exact algebraic
numbers, and certificates carry them as Arb balls produced from these
definitions.

- **P_B (primary target; v1.2 approved 2026-10-17): Q = (999/1000)·Q_max.**
  - Published support, [CCDHJ18] Fig. 3, ΛM² = 0.06 panel (tick labels read
    from the figure): β > ½ for Q/Q_max > 0.9920, and the l = 0 near-extremal
    mode is dominant for Q/Q_max > 0.9967. Its curve lies between −0.8κ₋ and
    −1.0κ₋ up to extremality.
  - S0 float spectrum, in units of κ₋:
    - l=0: NE −0.909381i and −1.855200i.
    - l=1: NE −1.911925i; dS −3.324477i; PS 7.214864 − 1.753306i.
    - l=2: NE −2.891324i; PS 12.196453 − 1.707868i.
    - l=10: PS 51.62353 − 1.68459i.
  - Eikonal limit: γ₀/κ₋ = 1.683212.
  - So β ≈ 0.90938, attained by an isolated l = 0 mode.
- **P_A (second target): Q = (249/250)·Q_max.** This is the representative point of [CCDHJ18] Table I.
  - All nine tabulated values are known answers (KAT-1). S0 reproduced six of
    them to |Δ| ≤ 4.6e−11.
  - β = γ₀/κ₋ = 0.764303. It is attained **only as l → ∞**, by the photon-sphere
    family: the Table I values decrease from 0.7958 (l=1) to 0.76431 (l=100).

Float constants (S0, from `src/s0_sanity.py geom`):

| quantity | P_B (Q/Q_max = 0.999) | P_A (Q/Q_max = 0.996) |
|---|---|---|
| Q/M | 1.00984316771 | 1.00681060565 |
| r_n, r₋, r₊, r_c | −7.96033597, 0.99736659, 1.09434835, 5.86862103 | −7.96006333, 0.94927598, 1.14311753, 5.86766982 |
| κ₋, κ₊, κ_c | 0.04254199, 0.03500736, 0.09338246 | 0.09426070, 0.06379974, 0.09332621 |
| κ₊/κ₋, κ_c/κ₋ | 0.8229, 2.1951 | 0.6768, 0.9901 |
| r_P, Ω_P/κ₋ | 1.95873037, 4.91883525 | 1.97187372, 2.20956624 |
| γ₀ = λ_P/2, γ₀/κ₋ | 0.07160718, 1.68321153 | 0.07204371, 0.76430267 |
| δ (frozen, D2) | 1/20 | 1/100 |
| strip edge (½+δ)κ₋ as a fraction of γ₀ | **0.327** | **0.667** |
| ρ = κ₋/κ₊ | 1.215 | 1.478 |

The last-but-one row is the decisive one for part (i) (`s0-gate.md` §1, §4).
Both points satisfy (½+δ)κ₋ < min(κ₊, κ_c), so
**E₊ ∪ E_c does not meet the strips S_δ** (D2), and (½+δ)κ₋ < γ₀.

### 1.4 Definitions (frozen)

- **D1 (QNM).** Let l ∈ ℕ₀ and ω ∈ ℂ∖(E₊ ∪ E_c). Then ω is a QNM of degree l if
  W_l(ω) = 0. Its radial multiplicity m_l(ω) is the order of the zero, and its
  multiplicity as a resonance of □_g is (2l+1)·m_l(ω).
- **D1′ (equivalence; an S1 lemma, written proof with citations).** For
  ω ∉ E₊ ∪ E_c, D1 is equivalent to σ = ω being a resonance (a QNM) in the
  sense of [Vas13], [HV17] and [Hin25] [v0.2: read as the resonances of the
  forward problem on Ω₂₃ (HV17 §2.5) = the poles of the localized
  ψ(P_σ−iQ_σ)⁻¹ (Vas13 §3) = the QNM set of Hin25 §1]: a smooth resonant state in (t*, r, θ, φ)
  coordinates, with the same multiplicity [v0.1: the multiplicity clause is
  proved for Im ω > −min(κ₊, κ_c) only; see the Revision log]. The ingredients are separation of
  variables, and the fact that in the hyperboloidal variable u = e^{−iωT}R
  (so that φ = e^{−iωt*}u Y_lm) the exponents at each horizon are {0, iω/κ}.
  Smoothness selects the analytic branch whenever iω/κ ∉ ℕ, which holds on S_δ
  at both points.
- **D2 (strip).** S_δ := {ω : Im ω ≥ −(½+δ)κ₋}, with δ_B = 1/20 and
  δ_A = 1/100.
- **D3 (gap, β).** α := inf{−Im ω : ω a QNM of some degree l, ω ≠ 0}, and
  β := α/κ₋. This is [CCDHJ18] eq. (9) and the α₀ of [HV17].
- **D4 (zero mode).** For l = 0, R ≡ 1 solves (R) at ω = 0 and is analytic at
  both horizons, so W₀(0) = 0 exactly. *Requirement:* the zero is simple,
  m₀(0) = 1. For l ≥ 1, ω = 0 is not a QNM (P2-0).
- **D5 (families, descriptive only).**
  - PS: Re ω ≈ (l+½)Ω_P, and −Im ω → γ₀ as l → ∞ ([SBZ97], [Dya16]).
  - dS: ≈ −ilκ_c^{dS}, with κ_c^{dS} = √(Λ/3) ([CCDHJ18] eqs. 11–12).
  - NE: ≈ −iκ₋(λ_l⁺ + n) ([Hin25] Thm 1.1).

  Every theorem below counts all zeros of W_l. A family label never enters a
  proof.
- **D6 (trapping constants).** Photon sphere r_P = (3M + √(9M² − 8Q²))/2, which
  does not depend on Λ. Then Ω_P = √F(r_P)/r_P,
  λ_P = r_P⁻¹√(F(r_P)(2 − 3M/r_P)) and γ₀ = λ_P/2 ([HV17] §2.7).
- **D7 (boxes).** R_l := [−W_l, W_l] × [−(½+δ)κ₋, η], with η = κ₋/10 and W_l
  the part-(ii) threshold.

## 2. Theorem ladder (claim templates for CLAIMS.md)

Every statement is at a frozen point X ∈ {P_B, P_A} with δ = δ_X, unless it
says otherwise.

### P2-0 Mode stability on the closed upper half-plane [target: CERTIFIED]
For every l and every ω with Im ω ≥ 0, W_l(ω) = 0 if and only if
(l, ω) = (0, 0).

Proof route: the exact hyperboloidal energy identity
D + L·N + |ω|²·𝒲 = −(|ω|²/Im ω)·B for Im ω ≠ 0 (`s0-gate.md` §2), plus the
boundary-flux argument on the real axis. Evidence: a sympy identity chain plus a
three-line written argument. Part of P2-0 is that the identity is checked as an
exact identity, not quoted.

### P2-1 Part (i): large-l exclusion [target: CERTIFIED]
There is an explicit l₀ such that for every l ≥ l₀, W_l has no zero in S_δ.
"Every" means uniformly in Re ω ∈ ℝ: above-barrier, barrier-top, sub-barrier
and low-frequency regimes.

The required width as a fraction of the trapping gap is 0.327·γ₀ at P_B and
0.667·γ₀ at P_A. The method is chosen at the S2 gate (`s0-gate.md` §4). Its
certificate format is fixed in §4 below.

### P2-2 Part (ii): high-frequency exclusion per l [target: CERTIFIED]
For each l < l₀ there is an explicit W_l such that W_l(ω) ≠ 0 on
S_δ ∩ {|Re ω| ≥ W_l}. The route is Liouville–Green above the barrier, with
Olver-type Volterra error control ([Olv74] Ch. 6) [v0.2: implemented as
weighted bounds along the real line, which is anti-progressive for the
horizon-referenced solutions; see the Revision log]. The error-control
integrals are enclosed in Arb.

### P2-3 Part (iii): validated count on the compact remainder [target: CERTIFIED]
For each l < l₀:
- W_l is analytic on a neighbourhood of R_l;
- W_l does not vanish on ∂R_l;
- the winding number of W_l around ∂R_l is 1 for l = 0 (the simple zero at 0)
  and 0 for 1 ≤ l < l₀.

The l = 0 case also certifies D4. The implementation is ω-Taylor-model tiles
over a D-finite connection, following the P4 Theorem B design
(`problems/P4/notes/s2-theorem-b.md` §3).

### P2-4 Assembly [target: CERTIFIED, level = min of P2-0…P2-3]
At X, S_δ∖{0} contains no QNM of any degree, and 0 is a simple QNM with a
constant resonant state. Hence **β ≥ ½ + δ** (β ≥ 0.55 at P_B, β ≥ 0.51 at
P_A).

### P2-5 Conversion [theorem by citation; level of P2-4, marked "conversion cited"]
Take α = (½+δ)κ₋ < γ₀, which holds at both points (§1.3). Then [HV17] Thm 1.1
(quantitative form §§2.5–2.6) applies. Every solution of □_gφ = 0 with smooth
Cauchy data on RNdS at X equals u₀ + u′, with u′ ∈ e^{−αt₀}H^{1+δ−0} near the
Cauchy horizon. In particular φ ∈ H¹_loc up to the Cauchy horizon.

This is the linear analogue of the failure of Christodoulou's SCC for smooth
data, at X. [CF17] Thm 3.1 gives the same H¹ conclusion with a C⁰ extension,
since 2 min{α, κ₊} > κ₋ at both points. Neither theorem is re-proven. [v1.0: the
scope of P2-5 at P_B is narrowed; see the Revision log, v1.0, "Scope of P2-5".]

### P2-6 Two-sided β at P_B [target: CERTIFIED; stretch goal, not part of the entry point (v1.2)]
A Krawczyk enclosure of the simple l = 0 zero ω_* ≈ −0.909381iκ₋, together
with P2-1…P2-3 run with the strip edge at −0.92κ₋ (where the l = 0 winding is
2), gives **β ∈ [β_lo, β_hi]** with β_hi − β_lo ≤ 10⁻⁸. Part (i) is then
needed up to 0.547·γ₀. [v0.2: this edge passes the first point of E₊,
−iκ₊ ≈ −0.823iκ₋, where W₀ has a simple pole; P2-6 needs new foundations,
see the Revision log.] [v1.1: the winding parenthesis is corrected, and the route changes; see the
Revision log, v1.1.]

### P2-7 The P_A statement [target: CERTIFIED; second target]
P2-0…P2-5 at P_A with δ = 1/100. Part (i) is needed up to 0.667·γ₀.

**Not targeted:** a two-sided β at P_A. It would need part (i) up to γ₀ − ε,
the sharp essential gap with explicit constants.

## 3. Published conversion theorems: hypotheses at the frozen points

| result | conclusion | hypotheses | P_B | P_A |
|---|---|---|---|---|
| [HV17] Thm 1.1, §§2.5–2.6 (linear, fixed RNdS) | u′ ∈ H^{1/2+α/κ₋−0} at the Cauchy horizon | 0 is the only resonance in Im σ ≥ −α, 0 simple with constant states, α < γ₀ | α = 0.55κ₋ < 1.683κ₋ ✓ (gives H^{1.05−0}) | α = 0.51κ₋ < 0.764κ₋ ✓ (gives H^{1.01−0}) |
| [CF17] Thm 3.1 (linear, fixed RNdS) | bounded energy; C⁰ ∩ H¹_loc extension | Price law with rate p on H⁺ for all T^lΩ^I derivatives (Assumption 2.2); 2 min{p, κ₊} > κ₋ | 2κ₊ = 1.646κ₋ ✓ | 2κ₊ = 1.354κ₋ ✓ |
| [CGNS18] Thm 1.1 (nonlinear, spherical, real scalar) | no mass inflation; C⁰ metric, Γ ∈ L²_loc, φ ∈ H¹_loc | two-sided exponential Price law on the dynamical H⁺ (assumed); ρ < 9/7 and β > 7/9 ([CCDHJ18] fn. 1) | ρ = 1.215 ✓, β ≈ 0.909 ✓ (float) | ρ = 1.478 ✗ |
| [Ros25] (nonlinear, spherical, charged/massive scalar) | H¹ extension with Γ ∈ L²_loc | Price-law upper bound with s > K₋ (β > ½); ρ < 2 | ✓ | ✓ |

The P2 ledger targets are the **linear** rows only. In P2-5 the rate
p = α < α₀ needed by [CF17] comes from the resonance expansion of [HV17]. The
nonlinear rows are context: their Price-law hypotheses concern the dynamical
spacetime, and P2 does not establish them.

## 4. Certificate formats and deliverables

- `src/p2/` (from S1): `params.py`, which turns the exact definitions of §1.3
  into Arb balls, with the defining polynomials recorded; `frobenius.py`, the
  series at r₊ and r_c as ω-Taylor models with tail bounds in the style of
  [Mez19] majorants, following P4's sonic-series pattern; `connect.py`, the
  Taylor-series continuation with rigorous truncation; `wronskian.py`;
  `modecount.py` (winding), `analyticity.py` and `krawczyk.py`, all following
  P4's Theorem B; `lg_bounds.py` (part ii); `large_l/` (part i).
- `certificates/`, as JSON:
  - `params.json`;
  - `winding_{X}_l{l}.json`, with tile endpoints as exact dyadics, argument
    increments as exact dyadic balls, model remainders and the precision used;
  - `analyticity_{X}_l{l}.json`;
  - `lg_{X}_l{l}.json`, the error-control enclosures that produce W_l;
  - `large_l_{X}.json`, the pointwise-inequality boxes with their margins,
    including closing boxes at h = 1/(l+½) → 0, or the method-specific
    equivalent fixed at S2;
  - `krawczyk_PB.json`.

  Every file records the point's exact definition, δ, the l range, the code
  commit hash, precision and timing. Re-check tests live in `tests/`.
- Ledger levels follow `problems/P9/b/FORMULATION.md` §3.3:
  - CONJECTURED: float evidence.
  - VERIFIED_N: an independent second implementation agrees (KAT-7).
  - CERTIFIED: Arb certificates committed and re-checked by tests.

## 5. Known-answer tests (gate: nothing above HEURISTIC until these pass)

1. **KAT-1.** [CCDHJ18] Table I at P_A: all nine tabulated values (l = 0, 1,
   2, 3, 4, 5, 10 and 100, plus dS l = 1) to at least 10 decimals. S0
   reproduced 6 of 9 in floats with |Δ| ≤ 4.6e−11.
2. **KAT-2.** Eikonal consistency at P_A. γ₀/κ₋ = 0.7643026660 lies below the
   Table I l = 100 value 0.7643094446, and the difference shrinks like O(l⁻²).
   Also Re ω/(l+½) → Ω_P (at l = 100: 222.0614 predicted against 222.0603).
3. **KAT-3.** [CCDHJ18] Fig. 3 at ΛM² = 0.06: the crossings β = ½ at
   Q/Q_max ≈ 0.9920 and the NE takeover at ≈ 0.9967 must be reproduced to the
   printed digits, using an l = 10 PS proxy as the figure does.
4. **KAT-4.** The [Hin25] Thm 1.1 limit. At Q/Q_max ∈ {0.9999, 0.99999}, the
   NE modes for l = 0, 1, 2 divided by κ₋ approach −i(λ_l⁺ + n), with ϰ taken
   at the current r₊ and r_c.
5. **KAT-5.** The [CCDHJ18] supplementary Nariai formula near minimal charge
   for one ΛM² > 1/9. This exercises the code at r₊ ≈ r_c.
6. **KAT-6.** A validated winding on small boxes around known float modes
   returns 1: Table I l = 1 PS, and the P_B l = 0 NE mode. The Krawczyk boxes
   contain the float values.
7. **KAT-7.** Two independent float methods agree to at least 10 digits on the
   P_B spectrum of §1.3 before any P_B number is promoted. The methods are the
   Frobenius/Taylor connection of `src/s0_sanity.py` and a hyperboloidal
   Chebyshev collocation in the manner of [Jan17].

## 6. Stages

Sizing and kill criteria are in `s0-gate.md` §7.
- **S0** (this document, the digest, the gate memo).
- **S1**: exact parameters in Arb; lemma D1′; the P2-0 identity in sympy;
  KAT-1…KAT-7; the P2-2 thresholds W_l. About 1–2 sessions.
- **S2**: a feasibility probe for part (i) on the barrier-top model, then the
  method choice. This is the kill gate.
- **S3**: the part-(iii) pipeline and runs for l < l₀.
- **S4**: certification of part (i).
- **S5**: assembly, adversarial review, ledger rows, and a v1 of this
  formulation.

## 7. Scope fences (what is NOT claimed)

- No nonlinear Einstein–Maxwell statement. [CGNS18] and [Ros25] are context,
  and their Price-law hypotheses are not established here.
- None of the following: Kerr–dS or Kerr–Newman–dS; charged, massive or
  conformally coupled scalars (the [AD26] potential carries an extra
  +(2Λ/3)F); gravito-electromagnetic Kodama–Ishibashi perturbations ([DRS18]).
- No sharpness claim. P2-4 is a lower bound on β. The only upper bound
  targeted is P2-6 at P_B. Nothing is claimed about generic solutions being no
  smoother than H^{1/2+β}.
- No statement at parameters other than P_B and P_A, and no open-neighbourhood
  claim. A certified neighbourhood may be proposed at S5 by revision-log entry;
  it is not pre-registered.
- Nothing about the near-extremal limit, and in particular nothing on [Hin25]
  Conj. 1.5. Rough data ([DSR18]) and quantum effects ([HWZ20]) are out of
  scope.
- Family labels (D5) are descriptive. No claim is made that the PS, dS and NE
  families exhaust the spectrum, because the certificates are family-blind.

## Revision log

- **v0 (2026-10-16, S0).** Initial formulation. Frozen:
  - the field model of §1.1;
  - the radial normalization and Wronskian of §1.2;
  - the points P_B (Q/Q_max = 999/1000) and P_A (249/250) at ΛM² = 3/50, with
    exact algebraic definitions;
  - δ_B = 1/20 and δ_A = 1/100;
  - definitions D1–D7, the ladder P2-0…P2-7, and KAT-1…KAT-7.

  The roles of P_B (primary) and P_A (second) depend on approval of the
  proposed .tex v1.2 entry (`s0-gate.md` §7). If it is not approved, P_A is the
  only target and P_B, P2-6 are dropped. The literature was checked on
  2026-10-16 (`notes/literature-digest.md`).

- **v0.1 (2026-10-16, S1 close-out).** D1′ narrowed to what is proved. The
  equivalence "zero of W_l ⟺ resonance" holds as stated on ω ∉ E₊ ∪ E_c,
  but the clause "with the same multiplicity" is proved only for
  Im ω > −min(κ₊, κ_c) (`notes/s1-p20-d1prime.md`, part 2). That region
  contains S_δ at both frozen points (certified in `certificates/params.json`;
  margins in `notes/s1-report.md`), so no target changes. No other frozen
  item changes.

- **v0.2 (2026-10-16, after the independent review of S1;
  `docs/review-2026-10-16-p2-s1-adversarial.md`).**
  - D1′ is read in the localized sense: the resonances of the forward
    problem on Ω₂₃ (HV17 §2.5), which are the poles of the localized
    ψ(P_σ−iQ_σ)⁻¹ (Vas13 §3) and the QNM set of Hin25 §1. HV17's global
    P̂(σ)⁻¹ can have other resonances (HV17 Lemma 2.14(1)). Its level is
    WRITTEN_PROOF.
  - P2-2's error control is a weighted Volterra bound along the real line,
    which is anti-progressive for the horizon-referenced solutions, not
    Olver's progressive-path control. The bound is valid. The threshold is
    written a_l in new text, to avoid a clash with the Wronskian W_l(ω).
  - P2-6 needs new foundations. At P_B its strip edge (−0.92κ₋) passes the
    first point of E₊ (−iκ₊ ≈ −0.823iκ₋), where W₀ has a simple pole, so a
    winding count there counts zeros minus poles. That edge is also outside
    D1′'s multiplicity window and outside the Liouville–Green condition
    β < κ₊. The main target needs the strip only down to −(½+δ)κ₋
    (−0.55κ₋ at P_B) and is unaffected.
  - KAT-6's Krawczyk step used the inscribed-disc radius r instead of the
    square's corner distance √2·r. It is fixed, and the certificate is
    regenerated; both modes still certify.

- **v0.3 (2026-10-17, S2 start).** The .tex v1.2 entry proposed in
  `s0-gate.md` §7 is approved (2026-10-17) and recorded in the .tex revision
  log, with one change: the two-sided certified enclosure of β is a stretch goal
  and not part of the entry point.
  - P_B is the primary target and P_A the second (§1.3 role labels updated;
    the conditional "if v1.2 is not approved" fallback of v0 is void).
  - P2-6 is a stretch goal. It also needs the new foundations noted in v0.2
    (the −0.92κ₋ edge passes E₊'s first point, lies outside D1′'s multiplicity
    window, and violates the Liouville–Green condition β < κ₊).
  - The header status is updated to match. No frozen convention, definition,
    point, δ or KAT changes.

- **v0.4 (2026-10-17, S2 gate; `notes/s2-report.md`).** The P2-1 route is
  chosen at the S2 gate, as §2 provides. It is **M3**: rotate the tortoise
  coordinate by π/4 through the photon sphere, then prove a variable-phase
  coercivity bound certified by a ground-state (Allegretto–Piepenbrink)
  pointwise inequality (`notes/s2-m3.md`). M2 was not implemented
  (`notes/s2-preregistration.md` §2).
  - The part-(i) certificate format is `certificates/m3/m3_{X}_geometry.json`
    (contour, sector analyticity, tails), `m3_{X}_l{l}.json` (window, ground
    state, margins, negative control, P2-2 overlap) and `m3_{X}_summary.json`.
    It replaces the placeholder `large_l_{X}.json` of §4. The h → 0 closing
    boxes of S4 will use the same inequality in the scaled variable.
  - S2 result at P_B: the barrier-top window
    [ω_l − 8γ₀, ω_l + 6γ₀] × [−(½+δ)κ₋, 0] is QNM-free for l = 20, 50, 100
    (machine-checked inequalities plus written lemmas, so WRITTEN_PROOF; the
    independent review, `docs/review-2026-10-17-p2-s2-adversarial.md`, found it
    proven as stated). Margins μ_cert/μ_HO are 0.528, 0.568, 0.560. Verdict: GO
    for S3 under the pre-registered rule; the margin at l = 100 is 0.56 μ_HO.
  - The projected l₀ = 2 (fit C = 0.94) is only the formal output of the
    pre-registered fit. The fit is mis-specified (s2-report §3). m(l) is a sum of
    design losses with opposite l-dependence:
    - the 5% design slack;
    - one ground state for the whole ω² hull, a loss of about 1/l;
    - the spline sawtooth at the fixed node spacing 0.01, a loss growing about
      like √l;
    - the mean-value cells, which cost < 1e-4 μ_HO.
  - The three points certify l = 20, 50 and 100 only. They say nothing about
    other l, and nothing below l = 20. With the design held fixed, the √l term
    would eventually exhaust the margin, so a fixed design cannot be uniform in l.
    S4 must scale the design: node spacing ∝ (l+½)^{−1/2} (fixed in
    σ = s(l+½)^{1/2}) and window width ∝ the oscillator quantum.
  - Still open for P2-1: uniformity in h, and the sub-barrier and low-frequency
    regimes (S4). P_A has not been probed.

- **v0.5 (2026-10-17, S3 exit; `notes/s3-preregistration.md`, `notes/s3-report.md`).** Part (iii)
  is implemented as follows; the deviations from §4 were declared in the pre-registration.
  - **W_l.** It is computed from the phase-factored Jost solutions
    Z₊ = e^{iω(r*−c₊)}R₊ and Z_c = e^{−iω(r*−c_c)}R_c (`src/p2/zsol.py`), via the exact identity
    W_l = [Δ(Z₊Z_c′ − Z₊′Z_c) + 2iωr²Z₊Z_c]/(r₊r_c), which is checked in sympy. The S1 code,
    which evaluates the free phase, gives useless box bounds at high |ω|.
  - **Counting.** Cauchy-interpolation discs (S1 `winding.py`) replace §4's ω-Taylor-model tiles.
    Only the right half of ∂R′_l is computed, by the conjugation symmetry. The box is
    R′_l = [−a_l, a_l] × [−b, e], with dyadic b ≥ (½+δ)κ₋ and 0 < e ≤ η.
  - **Certificates.** `certificates/part3/winding_{X}_l{l}.json`.
    - Each certificate stores every disc (c, ρ). The pole clearance √2ρ < dist(c, E₊ ∪ E_c) is
      enforced at compute time, and can be re-checked from (c, ρ); the S3 review's
      `check_certs.py` does so. There is no separate `analyticity_*.json`.
    - The l = 0 Krawczyk enclosure is `certificates/part3/krawczyk_PB_l0.json`.
  - **S3 result at P_B, l = 0–19.** The counts on R′_l are CERTIFIED: 1 at l = 0 (the simple zero
    at 0) and 0 for 1 ≤ l ≤ 19. They use the analyticity of W_l off E₊ ∪ E_c (§1.2). The
    independent review (`docs/review-2026-10-17-p2-s3-adversarial.md`) found them proven as stated.
  - **The strip corollary is at P2-0's level.** With P2-0 (Im ω > 0) and P2-2 (|Re ω| ≥ a_l),
    S_δ contains no QNM of degree l ≤ 19 other than ω = 0. Its Im ω > 0 part, including the band
    0 < Im ω ≤ e with |Re ω| ≥ a_l, rests on P2-0's written argument.
  - **Open.** P2-3 for 20 ≤ l < l₀, which waits for S4; P_A; and the two-sided β (P2-6).
  - The wording of this entry was corrected after the S3 review (2026-10-17): the corollary's
    level, the reference for analyticity, ∂R′_l, the strip depth written as (½+δ)κ₋ (β is D3's
    gap), and the description of the pole checks.

- **v0.6 (2026-10-18, S4 exit; `notes/s4-preregistration.md`, `notes/s4-report.md`, `notes/s4-proof.md`).**
  Part (i) is implemented at P_B, uniformly in l ≥ l₀ = 20 (h = 1/(l+½) ≤ h₀ = 2/41). The
  deviations from §4 were declared in the pre-registration and its dated Amendments 1–6.
  - **Four regimes** cover {l ≥ 20} × {Re ω ≥ 0, −(½+δ)κ₋ ≤ Im ω ≤ 0}. ζ := h Re ω.
    - **L** (Re ω ≤ 8κ₋): a τ-scaled Frobenius series at each horizon (r = r_h ± h²τ), uniform in
      h² ∈ [0, h₀²], in affine Arb arithmetic. Then Riccati invariance through the forbidden bulk.
    - **S** (6κ₋ ≤ Re ω ≤ 0.9Ω_P/h): half-line coercivity on C³ contours passing below each
      turning point, in the normalised form where h enters only through ζ and 1/Re ω. Below
      ζ = 0.02Ω_P, a near-horizon Taylor model closes ζ → 0.
    - **T** (0.8–2.0 Ω_P/h): S2's rotated contour with the exact splitting
      h²E = E₀ + hE₁ + h²E₂ and a linear ground state. This replaces S2's per-l spline and window.
    - **A** (≥ 1.5Ω_P/h): monotone in Re ω.
  - **The part-(i) certificate format** is `certificates/part1/`:
    - `{T,A}_PB.json`, `S_PB_{plus,c}.json` and `L_PB_{plus,c}.json`;
    - `bulk_PB.json`, `cover_PB.json`, `controls_PB.json` and `replay_PB.json`.
    It replaces the `m3_{X}_*` format of v0.4 for part (i). The S2 geometry certificate is
    replayed as K1.
  - **S4 result at P_B.** W_l(ω) ≠ 0 for every l ≥ 20 on that half-strip. With the conjugation
    symmetry and P2-0, this gives **P2-1 at P_B with l₀ = 20**.
    - Level: WRITTEN_PROOF. The inequalities are certified in Arb (128-bit). W1 (vertical-strip
      analyticity of r(x)) and W2 (Jost continuation and decay on the curved contours) are
      non-routine written arguments. The Im ω > 0 part is at P2-0's level.
    - The target level of P2-1 (CERTIFIED) would need W1 and W2 machine-checked.
  - **With P2-3B, every l < l₀ = 20 is counted at P_B, on P2-3B's boxes R′_l.** The S3 fallback
    for 20 ≤ l < l₀ was not needed. P2-1 contains P2-1w's statement.
  - **The independent review** (`docs/review-2026-10-18-p2-s4-adversarial.md`, 2026-10-18/19) found
    the theorem proven as stated, with nothing blocking.
    - Every certificate was replayed in full, and reviewer code sharing nothing with `p2` confirms
      every certified bound.
    - S3's counting pipeline independently finds no zero of W₂₀ in R′₂₀.
  - **The review's fixes (2026-10-19).**
    - S1: the S closing box's truncation term bounded P_N on the real segment, while the contour
      visits complex ρ. It now uses rigorous complex bounds, and the regenerated margins agree with
      the first ones to 8 digits.
    - The rest: a dated correction to Amendment 6 (the launch-time K4 rule, and the secant
      mechanism), the simply connected union of strips and windows in W1/W2, and wording.
  - **Open.** The assembly P2-4/P2-5 (S5); P_A; P2-6.

- **v1.0 (2026-10-26, S5 exit; `notes/s5-plan.md`, `notes/s5-assembly.md`, `notes/s5-report.md`).**
  The ladder at P_B is assembled up to P2-5 (P2-4 and P2-5 entered in the ledger on 2026-10-26). The
  status at P_B, with levels:
  - **Already in the ledger.**
    - P2-0: WRITTEN_PROOF.
    - P2-1: WRITTEN_PROOF, with l₀ = 20.
    - P2-2: CERTIFIED on the lower strip, for l ≤ 50 and l = 100.
    - P2-3: P2-3B counts every l < l₀ = 20 on R′_l. The counts are CERTIFIED; the strip corollary
      is at P2-0's level.
  - **P2-4: WRITTEN_PROOF**, the minimum of its ingredients (P2-0, P2-D1′ and P2-1). Entered in
    the ledger on 2026-10-26.
    - For every l, S_δ∖{0} contains no QNM, and 0 is a simple QNM with a constant resonant state.
    - So β ≥ 0.55. With P2-K6, 0.55 ≤ β ≤ 0.90938098 (two-sided, not sharp).
    - The bookkeeping is machine-checked, with negative controls:
      `certificates/assembly/assembly_PB.json` (`python -m p2.assembly`).
  - **P2-5: WRITTEN_PROOF, conversion cited** (P2-4's level). Entered in the ledger on 2026-10-26.
    - [HV17]'s α hypothesis concerns the forward problem on Ω₂₃ (quoted from its §2.5). P2-4 with
      P2-D1′ supplies it exactly, so the conversion is not conditional.
    - For smooth data on [HV17]'s H_I, u′ ∈ e^{−αt₀}H^{1.05−0} with α = 0.55κ₋ (decay e^{−αt₀} at α
      exactly via some α′ ∈ (α, γ₀); s5-assembly H7). [CF17] applies with p = α.
  - **Achieved against target levels.** §2 sets CERTIFIED as the target for P2-0, P2-1 and P2-4;
    the achieved level is WRITTEN_PROOF. CERTIFIED would need P2-0's written steps, P2-1's W1–W2
    and D1′ to be machine-checked.
  - **Scope of P2-5 (a narrowing of the §2 template).**
    - (a) Data on HV17's one-ended initial surface H_I: the conclusions hold in its domain of
      dependence, near the Cauchy horizon towards i⁺.
    - (b) Data on a complete Cauchy surface: near every point of the Cauchy horizon other than its
      bifurcation sphere B, with constants depending on the point (Killing translates of H_I).
    - (c) Nothing at B.

    For [CF17], the C⁰ ∩ H¹_loc extension is to the Cauchy-horizon component adjacent to the
    event-horizon component carrying the Price law, with the same exception.

    This narrows the frozen §2 P2-5 template ("Every solution of □_gφ = 0 with smooth Cauchy data
    on RNdS at X … φ ∈ H¹_loc up to the Cauchy horizon"). The user approved recording it this way
    on 2026-10-26. The §2 text is not rewritten; it points here.
  - **Open.** P2-6 (the stretch goal, which needs the new foundations of v0.2) and P2-7 (P_A, where
    no part-(i) run exists).
  - The scope fences of §7 are unchanged. No frozen convention, definition, point, δ or KAT
    changes.
  - **The independent review** (`docs/review-2026-10-26-p2-s5-adversarial.md`, 2026-10-26) found
    P2-4 and P2-5 at P_B proven as stated, with nothing blocking.
    - It confirms from the arXiv TeX that [HV17]'s α hypothesis concerns the forward problem on
      Ω₂₃.
    - It writes out the general clause of [HV17] Lemma 2.14(1), which [HV17] proves only in
      outline.
    - Reviewer code sharing nothing with `p2` re-checks every seam of the cover of S_δ.
    - 22 further mutations of the inputs each make the assembly check fail.
  - **The review's fixes (2026-10-26), approved by the user.**
    - S1: the decay e^{−αt₀} at α = d_δ exactly, via some α′ ∈ (α, γ₀) (`notes/s5-assembly.md` §2.2
      H7).
    - S2: the scope at the bifurcation sphere, recorded above as a narrowing of the §2 template.
    - S3: the P2-5 row names the step of [HV17] Lemma 2.14(1) that is proved only in outline, and
      the source of the order-1 claim at 0.
    - W1–W9: wording, the source versions ([HV17] arXiv v2, [CF17] arXiv v2), a note on an [HV17]
      typo, and dated notes appended to the plan.
  - **The ledger (2026-10-26).** The rows P2-4 and P2-5 are entered verbatim from the review's
    §7, which its §9 confirms, and P2-1's parenthesis now reads "(assembled in P2-4)".
  - **The problems document (2026-10-26).** Revision 1.3 adds the v1.3 entry (`notes/s5-report.md`
    §7, with the review's §9 one-word change to "initial surface") and the bibitem
    `costafranzen2017`, as the user approved on 2026-10-26. The PDF is rebuilt with tectonic, and
    the v1.2 PDF is archived in `docs/problems/revisions/`. The P2 statement and entry point are
    unchanged from v1.2. The optional Status sentence (`s5-report.md` §7) is not applied.

- **v1.1 (2026-10-27, S6 start; `notes/s6-preregistration.md` at `c98ca16d`, with its Amendments
  A1–A6).** The user approved S6 on 2026-10-27, and the coordinator approved Phase 2 the same day.
  The coordinator approved recording this entry as coordinator: it changes no frozen definition,
  no target and no edge, and he informs the user. Five items:
  - **1. A correction to P2-6's parenthesis "(where the l = 0 winding is 2)".**
    - At depth 0.92κ₋ the full l = 0 box holds two zeros (0 and ω_*) and one simple pole, at
      ω₀ = −iκ₊ (the first point of E₊; v0.2).
    - The winding of W₀ is therefore 1 on the full box. On S6's band box
      [−a′₀, a′₀] × [−0.92κ₋, −κ₋/2] it is 0.
  - **2. P2-6's route changes method.** The frozen text reads "P2-1…P2-3 run with the strip edge at
    −0.92κ₋". S6 instead:
    - covers S_δ by P2-4;
    - covers the band −0.92κ₋ ≤ Im ω ≤ −κ₋/2 by part (i) at depth 0.92 for l ≥ 20, with T
      redesigned, S's design grid extended, and L using the pole-free series R̃ = (1 + 2ρ)R;
    - replaces part (ii) for l ≤ 19 by regime A at fixed h = 1/(l+½). Liouville–Green on the real
      line cannot reach Im ω < −κ₊;
    - uses pole-aware band counts for part (iii), l ≤ 19;
    - adds Lemma E: ω₀ is not a resonance for any l, and W_l has a simple pole there.
  - **3. Additive.** The band boxes B_l (D7 is unchanged), Lemma E, and the S6 KATs B1–B7.
  - **4. Unchanged.**
    - The edge, 0.92κ₋.
    - The target, β ∈ [β_lo, β_hi] with width ≤ 10⁻⁸. It is to be realised by P2-K6's Krawczyk
      ball, of width 6.4·10⁻³¹, reused unchanged.
    - D1, D1′ and D3. Lemma E shows that D1's QNMs and HV17's resonances agree in the strip.
    - The scope fences of §7.
  - **5. The level.** §2 states CERTIFIED as P2-6's target. The level expected from S6 is
    WRITTEN_PROOF, the minimum of its ingredients (`notes/s6-preregistration.md` §6).
  - The §2 text of P2-6 is not rewritten. It carries a bracketed pointer here.
