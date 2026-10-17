# P2 — A certified resonance-free strip for subextremal RNdS (β > 1/2, linear SCC failure for smooth data)

Status: formulation **v0** (stage S0), 2026-10-16. Nothing here is certified.
Every claim id is a ledger *target*: S0 adds no `CLAIMS.md` rows.

Changes after S1 begins go in the Revision log, never by silent edit (house
discipline of `problems/P9/b/FORMULATION.md`).

Pre-registered statement: `docs/problems/open-problems-theoretical-cosmology-2026.tex`
§P2, v1.1. That document is frozen. The S0 gate (`notes/s0-gate.md`) recommends
**GO with a modified entry point** and proposes a v1.2 revision-log entry
(`s0-gate.md` §7). The proposal has **not been approved**. This formulation is
written so that both outcomes are well-defined: see §1.3, and the roles of
P_B and P_A there.

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

- **P_B (primary target if v1.2 is approved): Q = (999/1000)·Q_max.**
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
- **P_A (second target; the sole target if v1.2 is not approved):
  Q = (249/250)·Q_max.** This is the representative point of [CCDHJ18] Table I.
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
  sense of [Vas13], [HV17] and [Hin25]: a smooth resonant state in (t*, r, θ, φ)
  coordinates, with the same multiplicity. The ingredients are separation of
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
Olver complex-variable error control along the real progressive path
([Olv74] Ch. 6). The error-control integrals are enclosed in Arb.

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
since 2 min{α, κ₊} > κ₋ at both points. Neither theorem is re-proven.

### P2-6 Two-sided β at P_B [target: CERTIFIED; stretch of the primary target]
A Krawczyk enclosure of the simple l = 0 zero ω_* ≈ −0.909381iκ₋, together
with P2-1…P2-3 run with the strip edge at −0.92κ₋ (where the l = 0 winding is
2), gives **β ∈ [β_lo, β_hi]** with β_hi − β_lo ≤ 10⁻⁸. Part (i) is then
needed up to 0.547·γ₀.

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
