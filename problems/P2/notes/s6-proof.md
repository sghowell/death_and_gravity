# P2 S6: the two-sided β at P_B (P2-6), proof note

Status: written 2026-10-27 on branch `p2/s6-beta` (from main at `300c22a7`), at the S6 exit. Not
yet reviewed.
- Plan: `notes/s6-preregistration.md` (`c98ca16d`), with its Amendments A1–A6 (`f82114bf`).
- Code: `src/p2/beta_*.py` (`7168b2e4`); the assembly check `src/p2/beta_assembly.py` (`5ee4796f`).
- Certificates: `certificates/beta/` (each records the commit it was computed at).
- Tests: `tests/test_beta.py`.
- Report: `notes/s6-report.md`.

**Notation.**
- e := 23/25, the strip edge of FORMULATION P2-6; S_e := {Im ω ≥ −eκ₋}.
- The band B := {−eκ₋ ≤ Im ω ≤ −κ₋/2}. It overlaps S_δ = {Im ω ≥ −(11/20)κ₋}.
- E₊ = {−inκ₊ : n ≥ 1}, E_c = {−inκ_c : n ≥ 1}; ω₀ := −iκ₊ = −0.8228895441i κ₋.
- ω_* is P2-K6's zero of W₀: ω_*/κ₋ ∈ −0.909380978994805455524933762622i ± 3.2·10⁻³¹ (real part
  ± 6.6·10⁻³⁵), `certificates/part3/krawczyk_PB_l0.json`.
- h = 1/(l+½), h₀ = 2/41, l₀ = 20. x := r* = Σ_i A_i ln|r − r_i|, A_i = r_i²/Δ′(r_i), so A₊ = 1/(2κ₊)
  and A_c = −1/(2κ_c); c₊, c_c are FORMULATION 1.2's constants.
- Float values at P_B (FORMULATION 1.3): r₊ = 1.0943483, r_c = 5.8686210, κ₋ = 0.0425420,
  κ₊/κ₋ = 0.8228895, κ_c/κ₋ = 2.1950905, Ω_P/κ₋ = 4.9188353.

## 0. The result

**P2-6 at P_B.** The forward problem for □_g on Ω₂₃ ([HV17] §2.5) has exactly two resonances σ with
Im σ ≥ −eκ₋, namely σ = 0 and σ = ω_*. Both are of degree l = 0, and each is a simple zero of W₀.
Hence D3's α is attained at ω_*, and

  **β = −Im ω_*/κ₋ ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹]**,

two-sided, of width 6.4·10⁻³¹ ≤ 10⁻⁸. In D1's terms: the QNMs of any degree in S_e are 0 and ω_*.
D1's QNMs and HV17's resonances agree in S_e.

The enclosure is P2-K6's Krawczyk ball, reused unchanged (A4). What S6 adds is the other side.
Apart from 0 and ω_*, there is no QNM and no resonance anywhere in S_e, whose edge −0.92κ₋ lies
0.0106κ₋ below ω_* (§6).

## 1. Reduction

1. **S_δ, every l: P2-4** (WRITTEN_PROOF). W_l has no zero in S_δ∖{0}; 0 is the simple l = 0 zero.
2. **One E-point.** S_e ∩ (E₊ ∪ E_c) = {ω₀}, since eκ₋ < 2κ₊ (margin 0.7258κ₋) and eκ₋ < κ_c
   (margin 1.2751κ₋); and ω₀ ∈ S_e (κ₊ < eκ₋, margin 0.0971κ₋).
3. **P2-D1′'s set clause** (WRITTEN_PROOF, any depth): for σ ∉ E₊ ∪ E_c, σ is a resonance iff
   W_l(σ) = 0 for some l.
4. So P2-6 follows from:
   - **(a) the band:** for every l, W_l has no zero in B, except ω_* when l = 0 (§§3–5);
   - **(b) Lemma E:** ω₀ is not a resonance, for any l (§2);
   - **(c)** ω_* is a resonance: it is a zero of W₀ (P2-K6) and ω_* ∉ E₊ ∪ E_c, so (3) applies.
5. **Multiplicity.** The radial order m₀(ω_*) = 1 is certified (Krawczyk; the band count §5). The
   resonance multiplicity at ω_* lies outside D1′'s window Im σ > −min(κ₊, κ_c) and is **not
   claimed**. β is an infimum over a set, so it does not need it.
6. **What each part of Lemma E is needed for.**
   - D3 defines β through D1's QNMs, which exclude E₊ ∪ E_c by definition. For β itself, (a) and
     (c) therefore suffice, together with the bound "W_l has at most a simple pole at ω₀" (§2.3).
     That bound is what turns the band windings into zero counts (§5).
   - The non-resonance half of Lemma E (W̃_l ≠ 0) is needed for the resonance statement: the
     forward problem's resonances in S_e are exactly 0 and ω_*, and D1's QNMs and HV17's
     resonances agree in S_e.
   - For l ≤ 19 the band counts force the pole, and with it W̃_l ≠ 0, independently of §§2.4–2.5.

## 2. Lemma E

**Lemma E.** For every l ≥ 0, ω₀ = −iκ₊ is not a resonance of the forward problem on Ω₂₃, and W_l
has a simple pole at ω₀ with residue

  Res_{ω₀} W_l = [i(2κ₊r₊ + L)/(2r₊²)]·W̃_l ≠ 0,

where W̃_l is defined in §2.1. Moreover, for l ≥ 1,

  **W̃_l = −(1/r₊)∫_{r₊}^{r_c} (L + rF′(r))·Z2(r) dr/r < 0**,   (E-int)

and W̃₀ = +0.038285581814409 ± 4.4·10⁻¹⁶.

### 2.1 Normalisation (A1)

These are zsol's conventions (`src/p2/zsol.py`), and they are the ones B5 tests.
- The Z-equations: with s = +1 at r₊ and s = −1 at r_c,
  ΔZ″ + (Δ′ − 2isωr²)Z′ − (2isωr + L)Z = 0 (′ = d/dr).
  They come from R = e^{−isωx}Z in (R).
- **Z₊ := e^{iω(x−c₊)}R₊ = Σ_{m≥0} a_m t^m, with a₀ = 1 and t = r − r₊** (`zsol.frob_z`). Z_c :=
  e^{−iω(x−c_c)}R_c is analytic at r_c with Z_c(r_c) = 1. R₊, R_c are FORMULATION 1.2's.
- W_l(ω) = [Δ(Z₊Z_c′ − Z₊′Z_c) + 2iωr²Z₊Z_c]/(r₊r_c), at any r ∈ (r₊, r_c). This is exactly
  FORMULATION 1.2's W_l: the factors e^{∓iω(x−c)} cancel against exp(−iωc₊ + iωc_c).
- At ω₀: **Z2 := t + Σ_{m≥2} b_m t^m** is the exponent-1 solution of the s = +1 equation (§2.2), and
  **W̃_l := [Δ(Z2·Z_c′ − Z2′·Z_c) + 2iω₀r²·Z2·Z_c]/(r₊r_c)**. Being Δ times the Wronskian of two
  solutions of (R), up to a constant factor, it does not depend on r.

### 2.2 Step 1: exponents, the log term, and the criterion

Write Δ = Σ_{k≥1} Δ_k t^k, with Δ₁ = Δ′(r₊) = 2κ₊r₊² > 0. Write also B(t) = Δ′ − 2iωr² = Σ B_k t^k
and C(t) = −(2iωr + L) = C₀ + C₁t.
- The coefficient of t^{m−1} gives D_m a_m + (terms in a_{m−1}, a_{m−2}, a_{m−3}) = 0, with
  D_m = m(Δ₁(m − 1) + B₀) = Δ₁m(m − iω/κ₊), since B₀ = Δ₁ − 2iωr₊² = Δ₁(1 − iω/κ₊).
- The m = 1 step is **D₁a₁ + C₀a₀ = 0**, with C₀ = −(2iωr₊ + L). sympy confirms the t⁰ coefficient
  (`epole_PB.json`, `sympy_first_step`).
- At ω₀ (iω₀ = κ₊): B₀ = 0 exactly, so **D₁ = 0**, while **C₀ = −(2κ₊r₊ + L) < 0**. Also
  D_m = Δ₁m(m − 1) ≠ 0 for m ≥ 2. The indicial exponents are {0, iω₀/κ₊} = {0, 1}.
  - The exponent-1 solution Z2 exists, and its series converges for |t| < r₊ − r₋ (the nearest
    other singular point).
  - The exponent-0 solution has a logarithm. It is 1 + c·Z2·ln t + Σ_{m≥1} d_m t^m, and the t⁰
    coefficient of the equation gives Δ₁c + C₀ = 0. So c = (2κ₊r₊ + L)/Δ₁ ≠ 0, and the solution
    contains t ln t.
- **Smoothness and the threshold.** By D1′'s set-up (`notes/s1-p20-d1prime.md` §4.3, with
  [Vas13] §2.6 and the proof of [HV17] Lemma 2.14), a resonant state at σ has smooth
  spherical-harmonic components. Each is u(r)Y_lm, where u = e^{−iσT}R solves the hyperboloidal
  radial equation on [r₊ − ε, r_c + ε] and lies in H^s above the threshold s > ½ − Im σ/κ at each
  horizon.
  - Near r₊, T = −x + (smooth), because T′ = T̃/F with T̃(r₊) = −1. So u = (smooth, nowhere
    zero)·Z, where Z = e^{iσ(x−c₊)}R solves the s = +1 equation. The exponents of u at r₊ are those
    of Z.
  - At σ = ω₀ the threshold at r₊ is ½ + 1 = 3/2. t ln t lies in H^{3/2−ε} but not in H^{3/2}. So
    the only admissible solutions at r₊ are the multiples of Z2, which is analytic through r₊.
  - At r_c the exponents of the s = −1 equation are {0, iω₀/κ_c} = {0, κ₊/κ_c}, with
    κ₊/κ_c = 0.3749. That is not an integer, so there are two Frobenius solutions and no logarithm:
    - Z_c (analytic, value 1);
    - Z̃_c = (r_c − r)^{κ₊/κ_c}·g̃(r), with g̃ analytic and g̃(r_c) = 1.

    Z̃_c lies in H^{0.875−ε} only, below the threshold ½ + κ₊/κ_c = 0.875. So the admissible
    solutions at r_c are the multiples of Z_c.
- **Criterion.** As in D1′'s proof (1), ω₀ is a resonance of degree l iff Z2 and Z_c, continued along
  (r₊, r_c) where the equation is regular, are linearly dependent. That holds **iff W̃_l = 0**. If
  W̃_l = 0, the common solution, continued analytically past both horizons, is a smooth resonant
  state. The continuation is unique, as in D1′.

### 2.3 Step 2: the residue

Take |ω − ω₀| < κ₊. Then D_m(ω) ≠ 0 for every m ≥ 2, since D_m vanishes only at −imκ₊.
- The recursion is linear in (a₀, a₁). So Z₊ = U(t; ω) + a₁(ω)·V(t; ω), where:
  - U and V are generated from (a₀, a₁) = (1, 0) and (0, 1);
  - both are analytic in ω near ω₀;
  - V(·; ω₀) = Z2.
- a₁ = −C₀(ω)/D₁(ω), and D₁(ω) = Δ₁(κ₊ − iω)/κ₊ = −2ir₊²(ω − ω₀). So
  (ω − ω₀)a₁(ω) → C₀(ω₀)/(2ir₊²) = i(2κ₊r₊ + L)/(2r₊²).
- Hence **Z₊ = [i(2κ₊r₊ + L)/(2r₊²)]·Z2/(ω − ω₀) + O(1)**. The bound holds near r₊, and along
  (r₊, r_c) by analytic dependence of the continuation on ω and on the initial data.
- Z_c is analytic in ω near ω₀, since ω₀ ∉ E_c.
- Therefore W_l has at most a simple pole at ω₀, and **Res_{ω₀}W_l = i(2κ₊r₊ + L)/(2r₊²)·W̃_l**. The
  prefactor is nonzero. So W_l has a simple pole at ω₀ iff W̃_l ≠ 0 iff ω₀ is not a resonance.
  Otherwise ω₀ is a removable point of W_l.

**B5 tests exactly this formula.** It compares ε·W_l(ω₀ + ε), where W_l is computed both by zsol and
by the independent S1 code `arbw` (both in FORMULATION 1.2's normalisation), with the residue
above.
- Settings: ε = 10⁻⁵κ₋ and 10⁻⁷κ₋, along (1+i)/√2, for l = 0, 1, 5, 19.
- Relative differences at ε = 10⁻⁵κ₋: 9.16·10⁻⁵, 1.42·10⁻⁵, 4.49·10⁻⁵, 6.76·10⁻⁵. They are 100
  times smaller at 10⁻⁷, so they are O(ε), and the two codes agree to every printed digit.
- Criterion: below 100ε. Pass.

### 2.4 Step 3: l ≥ 1 (written)

At ω₀, ω₀² = −κ₊². Everything below is real: the s = +1 equation has real coefficients, so Z2 is
real for r ∈ (r₊, r_c), and so are Z_c and Z̃_c.

**(a) V_l > 0 for l ≥ 1.** Write V_l = F(L + rF′)/r², with F > 0 on (r₊, r_c), and
rF′ = (2/r)(M − Q²/r) − (2Λ/3)r².
- The first term is positive for r > r₊ > Q²/M. Arb: r₊ − Q²/M = 0.0745651239349 ± 4.6·10⁻¹⁴.
- The second term is ≥ −(2Λ/3)r_c² > −2. Arb: 2 − (2Λ/3)r_c² = 0.622371488082 ± 4.9·10⁻¹³.
- So L + rF′ > L − 2 ≥ 0 for L ≥ 2.

These are in `epole_PB.json`, `V_l_positive_for_l_ge_1`.

**(b) The solutions** (Φ = rR solves Φ_xx = (V_l + κ₊²)Φ):
- Ψ := e^{κ₊x}, which solves Ψ_xx = κ₊²Ψ.
- **Φ₂ := r·R^{(2)} = r·e^{−κ₊(x−c₊)}Z2** with R^{(2)} := e^{−iω₀(x−c₊)}Z2, since −iω₀ = −κ₊.
  - Its hyperboloidal form is u^{(2)} := e^{−iω₀T}R^{(2)} = e^{−κ₊(T+x−c₊)}·Z2.
  - T + x is smooth near r₊, because (T + x)′ = (T̃ + 1)/F and T̃ + 1 vanishes at r₊. So the
    hyperboloidal factor e^{−κ₊(T+x−c₊)} is smooth and positive there.
  - So Φ₂ = rR^{(2)}, R^{(2)} = e^{iω₀T}u^{(2)}, and u^{(2)} is Z2 times a smooth positive factor:
    Φ₂ is the solution that is smooth at r₊ in the sense of §2.2.
- **Φ_c := r·R_c = r·e^{κ₊(x−c_c)}Z_c** (R_c = e^{iω₀(x−c_c)}Z_c). Its hyperboloidal form is
  u_c = e^{−iω₀T}R_c = e^{−κ₊(T−x+c_c)}·Z_c, and T − x is smooth near r_c. So Φ_c is the solution
  smooth at r_c.
- **Φ̃_c := r·e^{κ₊(x−c_c)}Z̃_c**, the other branch at r_c. Its hyperboloidal form is
  e^{−κ₊(T−x+c_c)}·Z̃_c, which is not smooth (§2.2).

**(c) Asymptotic expansions, with error terms.**
- **At r₊ (x → −∞).** x = (1/(2κ₊))ln t + c₊ + O(t), so τ := e^{2κ₊(x−c₊)} = t·h(t), with h analytic
  and h(0) = 1. Inverting, t = τ(1 + O(τ)), analytic in τ. Hence
  - Φ₂ = e^{κ₊(x−c₊)}G(τ), where G is analytic at 0 and G(0) = r₊;
  - Φ₂,x = κ₊e^{κ₊(x−c₊)}[G(τ) + 2τG′(τ)] (dτ/dx = 2κ₊τ).

  That is, Φ₂ = r₊e^{κ₊(x−c₊)}(1 + O(e^{2κ₊x})), and likewise Φ₂,x with the factor κ₊.
- **At r_c (x → +∞).** s := r_c − r. Then x = −(1/(2κ_c))ln s + c_c + O(s), so σ_c := e^{−2κ_c(x−c_c)}
  = s(1 + O(s)) and s = σ_c(1 + O(σ_c)). Hence:
  - Φ_c = e^{κ₊(x−c_c)}H(σ_c), with H analytic, H(0) = r_c, and
    Φ_c,x = e^{κ₊(x−c_c)}[κ₊H − 2κ_cσ_cH′];
  - Φ̃_c = e^{κ₊(x−c_c)}σ_c^{κ₊/κ_c}H̃(σ_c) = e^{−κ₊(x−c_c)}H̃(σ_c), with H̃ analytic and
    H̃(0) = r_c·(lim s/σ_c)^{κ₊/κ_c} = r_c. Also Φ̃_c,x = e^{−κ₊(x−c_c)}[−κ₊H̃ − 2κ_cσ_cH̃′].

  The relative error terms are O(σ_c) = O(e^{−2κ_cx}). Every solution is O(e^{κ₊x}) at +∞.

**(d) Φ₂ > 0 on ℝ, i.e. Z2 > 0 on (r₊, r_c).** By (c), Φ₂ > 0 and Φ₂,x > 0 for x below some x₀. Let
x₁ be the supremum of the x with Φ₂ > 0 and Φ₂,x > 0 on (−∞, x). Suppose x₁ < ∞.
- On (−∞, x₁), Φ₂,xx = (V_l + κ₊²)Φ₂ > 0. So Φ₂,x is increasing there, and Φ₂,x(x₁) > 0; also
  Φ₂(x₁) > 0.
- By continuity both stay positive past x₁, which contradicts the definition of x₁.

**(e) The Wronskian identity.** Write W[f, g] := f g_x − f_x g.
- W[Φ₂, Ψ]_x = Φ₂Ψ_xx − Φ₂,xxΨ = −V_lΦ₂Ψ.
- At −∞, by (c), W[Φ₂, Ψ] = −2κ₊τG′(τ)e^{κ₊(2x−c₊)} = O(e^{4κ₊x}) → 0. Hence
  W[Φ₂, Ψ](x) = −∫_{−∞}^{x} V_lΦ₂Ψ dy.
- **Convergence at both ends.**
  - At −∞: V_l = O(F) = O(τ) = O(e^{2κ₊y}), and Φ₂Ψ = O(e^{2κ₊y}). The integrand is O(e^{4κ₊y}).
  - At +∞: V_l = O(σ_c) = O(e^{−2κ_cy}), and Φ₂Ψ = O(e^{2κ₊y}). The integrand is
    O(e^{−2(κ_c−κ₊)y}). It is integrable because **κ₊ < κ_c**. Arb: κ_c − κ₊ = 0.0583751013871 ±
    2.4·10⁻¹⁴.
- Write Φ₂ = AΦ_c + BΦ̃_c (A, B real). By (c):
  - W[Φ_c, Ψ] = 2κ_cσ_cH′(σ_c)e^{κ₊(2x−c_c)} = O(e^{−2(κ_c−κ₊)x}) → 0;
  - W[Φ̃_c, Ψ] = e^{κ₊c_c}[2κ₊H̃ + 2κ_cσ_cH̃′] → 2κ₊r_c e^{κ₊c_c}.
- Letting x → +∞: **B·2κ₊r_c e^{κ₊c_c} = −∫_ℝ V_lΦ₂Ψ dx.**
- For l ≥ 1 the integrand is positive (V_l > 0, Φ₂ > 0, Ψ > 0). So B < 0. Φ₂ is not a multiple of
  Φ_c, Z2 and Z_c are independent, and **W̃_l ≠ 0**.

**(f) The closed form (E-int).**
- By (b), Δ(R^{(2)}R_c′ − R^{(2)′}R_c) = e^{κ₊(c₊−c_c)}r₊r_cW̃_l. For Φ = rR,
  W[rR₁, rR₂] = r²F(R₁R₂′ − R₁′R₂) = Δ(R₁R₂′ − R₁′R₂). Hence
  W[Φ₂, Φ_c] = e^{κ₊(c₊−c_c)}r₊r_cW̃_l.
- On the other hand, W[Φ₂, Φ_c] = B·W[Φ̃_c, Φ_c], a constant equal to its limit, 2κ₊r_c²·B, by (c).
- With (e) and V_lΦ₂Ψ dx = e^{κ₊c₊}(L + rF′)Z2 dr/r (since V_l/F = (L + rF′)/r²), this gives (E-int).
- The integrand is O(t) at r₊. At r_c it is O(s^{−κ₊/κ_c}), integrable since κ₊/κ_c < 1, because
  Z2 = e^{κ₊(x−c₊)}Φ₂/r = O(e^{2κ₊x}) = O(s^{−κ₊/κ_c}).
- For l ≥ 1, (a) and (d) give W̃_l < 0, as certified for every l ≤ 19 and at the samples
  (§2.6).
- **Float check, not used by the proof** (`tests/test_beta.py::test_lemma_E_closed_form`; scipy
  DOP853 with rtol 10⁻¹³ and the end-point singularity integrated analytically): (E-int) reproduces
  the stored Arb values of W̃_l for l = 0, 1, 2 and 5 within 10⁻⁸ relative, with their signs.
  At l = 0, for instance, it gives 0.03828558181442 against 0.038285581814409 ± 4.4·10⁻¹⁶. This
  confirms the normalisations and signs of (b)–(f) independently of zsol's Wronskian evaluation.

**Why l = 0 is excluded.** At L = 0 the integrand of (E-int) is F′(r)Z2(r), and V₀ = F·F′/r changes
sign. F′ > 0 on (r₊, r_F) and F′ < 0 on (r_F, r_c), where r_F is the maximum of F (F(r₊) = F(r_c) = 0).
So neither (e) nor (E-int) fixes the sign of B or of W̃₀. In fact W̃₀ > 0: the negative part near r_c,
where Z2 is large, dominates. l = 0 is therefore done in Arb (§2.5).

### 2.5 Step 4: l = 0, and l ≤ 19, in Arb (`certificates/beta/epole_PB.json`, `src/p2/beta_epole.py`)

The computation, at 128 bits:
- Z2's Frobenius series at r₊ + t₀, where t₀ is half of zsol's majorant radius (exact rational), with
  a geometric tail bound;
- zsol's Taylor continuation to r_m ≈ r_P;
- zsol's Z_c (`zsol.propagate("c")`), propagated to r_m;
- the formula of §2.1.

| l | W̃_l(ω₀) (real part; the imaginary parts are balls around 0) |
|---|---|
| 0 | +0.038285581814409 ± 4.4·10⁻¹⁶ |
| 1 | −22.800094349715 ± 4.3·10⁻¹³ |
| 2 | −20790.693939592 ± 2.5·10⁻¹⁰ |
| 3 | −1.6146662780386·10⁷ |
| 5 | −9.3233687215158·10¹² |
| 10 | −2.5336088082825·10²⁷ |
| 19 | −2.9600527340171·10⁵³ ± 3.3·10³⁹ |
| 20, 25, 30, 40, 50, 60, 100 (samples) | all negative, certified nonzero |

All 20 values for l ≤ 19 are certified nonzero.

### 2.6 Status of Lemma E, per l range

| l | ω₀ not a resonance; simple pole | level |
|---|---|---|
| 0 | Arb: W̃₀ ≠ 0 (§2.5). Also the band count (§5): winding 0 with ω_* inside forces the pole | **CERTIFIED** (written steps §§2.2–2.3 routine; first step sympy) |
| 1–19 | written: §2.4 (positivity), with (E-int). Arb: W̃_l ≠ 0 (§2.5). The band counts (§5): winding −1 forces the pole | **WRITTEN_PROOF and CERTIFIED** |
| ≥ 20 | written: §2.4. Also part (i)'s L certificate at depth e (§3): in the grid box containing ω₀, the solution regular at r₊ (R̃, which at ω₀ is the exponent-½ = Z2-type solution) and R_c have log-derivatives with Re y₊ ≥ 0.1617 > 0 > −0.1201 ≥ Re y_c, for every h ≤ h₀. So they are not proportional (W5, W6′). Arb samples to l = 100 | **WRITTEN_PROOF**; the L certificate is a second, independent route at the same level |

**Consequence.** D1's QNMs and HV17's resonances agree in S_e:
- off E₊ ∪ E_c by D1′'s set clause;
- at ω₀ neither notion applies (W_l has a pole there, and ω₀ is not a resonance).

So D3's α is HV17's α₀ for strips down to −eκ₋. No definitional change is needed.

**Literature.** No source is cited for Lemma E; the argument above is self-contained.

## 3. Part (i) at depth e for l ≥ 20 (`certificates/beta/part1/`)

S4's four regimes are run with the depth d = eκ₋ (`beta_common.deep_point`) and h₀ = 2/41, so
l ≥ 20. Everything is uniform in h ≤ h₀, as in P2-1.

| regime | region (Re ω ≥ 0) | design | result |
|---|---|---|---|
| L | [0, 8κ₋] × [−19/20, −11/20]κ₋ (the band rows; S_δ above is P2-4) | S4's L with R̃ = (1+2ρ)R in every box (W6′); 1280 grid boxes of side κ₋/20, none refined | all pass; min Re y₊ ≥ 0.1284, max Re y_c ≤ −0.0770 |
| S (both sides) | [6κ₋, (9/10)Ω_P/h] × [−eκ₋, 0] | S4's procedure, design grid extended (pre-registration O4) | 45/45 jobs per side; the S+ box [0.06, 0.08]Ω_P needs the extension |
| T | [(4/5), (13/10)]Ω_P/h × [−eκ₋, 0] | redesigned T6: s1 = 12, φ_∞ = −1, taper (3, 6) | 10/10 ζ-boxes, μ/a* = 0.334 … 0.239; no fallback design used |
| A | Re ω ≥ (6/5)Ω_P/h, [−eκ₋, 0] | ζ_A = 6/5, W7 | μ/a* = 1.972; W7: φ₀ = −1/5 < −0.0076 |

**W6′ (E₊ inside L).** S4's W6 used n + 2ρ ≠ 0, where ρ = −iω/(2κ₊) is the exponent of R₊.
- In S_e, Re(n + 2ρ) = n + Im ω/κ₊ ≥ n − e/(κ₊/κ₋) = n − 1.118 > 0 for n ≥ 2.
- n + 2ρ = 0 only for n = 1 and ω = ω₀.
- R̃ := (1 + 2ρ)R has ã₀ = 1 + 2ρ and ã₁ = −(A₃ρ(ρ−1) + B₂ρ + C₁)/A₂ (`beta_common.frob_tau_pf`). The
  factor 1 + 2ρ is cancelled exactly. R̃ is the solution regular at r₊ for every ω in the box: off
  ω₀ it is a nonzero multiple of R₊, and at ω₀ it is the exponent-(ρ+1) solution.
- y₊ = R̃′/R̃ does not depend on the normalisation, so W5 (Riccati invariance) is unchanged.
- y₊ has a finite certified enclosure in every box (Re y₊ > 0), which also shows that R̃ ≢ 0
  there.

**The cover** (`cover_PB.json`, replayed identically by the assembly). 11 inequalities and
6 tilings:
- the overlaps 6 < 8, 0.8 < 0.9 and 1.2 < 1.3;
- (9/10)Ω_P/h₀ > 8κ₋ (by 82.75κ₋);
- Re ω + Im ω > 0 in S, T and A;
- eκ₋ < 2κ₊, eκ₋ < κ_c and κ₊ < eκ₋;
- L's rows reach −eκ₋ and overlap S_δ;
- h₀ = 2/41 gives l₀ = 20;
- the tilings of T, A, S (both sides) and L.

The negative Jost margin κ₊ − eκ₋ < 0 is used nowhere. S, T and A use rotated contours (tails need
only Re ω + Im ω > 0), and L uses Frobenius series.

**Controls** (`controls_PB.json`):
- NB3: L's ω₀ box with S4's series fails ("affine form may vanish");
- NB4: S4's T design at depth e on [0.8, 0.85] gives μ < 0 (−0.0475a*);
- NB5: A at h₀ = 2 with ζ_A = 5 gives μ < 0 (−0.392a*).

**Level.** WRITTEN_PROOF, as P2-1: the lemmas W1, W2, W5 and W6′ are written, and the numerics are
certified.

## 4. Part (ii′): regime A at fixed h_l, l ≤ 19 (`certificates/beta/part2/A_fixed_PB.json`)

For each l ≤ 19, run S4's regime A (design A: φ₀ = −1/5) with h₀ := h_l = 1/(l+½).
- Its lower bound F = E₀/h₀ + E₁ + h₀min(0, E₂) is valid for every h ≤ h_l, in particular at h = h_l.
- At ζ = ζ_A(l) and ν ∈ [−eκ₋, 0] it gives E ≥ μ/h > 0.
- Above ζ_A(l) it is monotone in Re ω, by S4's Lemma W7 (φ ≤ φ₀ < −arctan(eκ₋h_l/(ζ_AΩ_P))).
- So **W_l(ω) ≠ 0 for Re ω ≥ a′_l := ζ_A(l)Ω_P/h_l (rounded up to 2⁻³⁰), −eκ₋ ≤ Im ω ≤ 0**. On the
  rotated contour the Jost solutions are subdominant whenever Re ω + Im ω > 0 (S4 W2). Re ω < 0
  follows by the symmetry.

The ζ_A rule (pre-registration §3): take the smallest multiple of 1/20 in [1.1, 12]Ω_P with float
margin ≥ 0.02a* and W7, then certify it in Arb. Every l passed at its float value, with no Arb retry.

| l | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7–9 | 10–19 |
|---|---|---|---|---|---|---|---|---|---|
| ζ_A/Ω_P | 21/4 | 37/20 | 3/2 | 27/20 | 5/4 | 6/5 | 6/5 | 23/20 | 11/10 |
| a′_l/κ₋ | 12.912 | 13.650 | 18.446 | 23.241 | 27.668 | 32.464 | 38.367 | 42.425, 48.082, 53.738 | 56.813 … 105.509 |
| μ/a* | 0.0348 | 0.0386 | 0.0892 | 0.1014 | 0.0573 | 0.0533 | 0.1812 | 0.081, 0.175, 0.269 | 0.047 … 0.596 |

**Level.** WRITTEN_PROOF, through S4's written W1, W2 and W4, with certified numerics.

## 5. Part (iii): pole-aware band counts, l ≤ 19 (`certificates/beta/part3/band_PB_l{l}.json`)

**The boxes.** B_l = [−a′_l, a′_l] × [−b′, −b₀], where:
- a′_l is part (ii′)'s threshold (identical by the assembly check);
- b′ = eκ₋ rounded up to 2⁻³⁰, which is 21012395/2²⁹ (b′ − eκ₋ = 3.9·10⁻¹⁰, i.e. 9.1·10⁻⁹κ₋);
- b₀ = κ₋/2 rounded down to 2⁻³⁰, which is 22839559/2³⁰.

**Meromorphy.** E₊ ∪ E_c meets a neighbourhood of B_l only in ω₀, which is interior (b₀ < κ₊ < b′).
W_l is analytic elsewhere there, and has at most a simple pole at ω₀ (§2.3).

**The argument principle.** With no zero or pole on ∂B_l, winding(W_l, ∂B_l) = Z − P, where:
- Z is the number of zeros in B_l with multiplicity;
- P ∈ {0, 1} is the order of the pole at ω₀.

**The computation.**
- S3's disc machinery (`part3.plan_edge`, `part3.segment_with_split`) gives certified increments of
  arg W_l along the right half H of ∂B_l, (0, −b′) → (a′, −b′) → (a′, −b₀) → (0, −b₀).
- Each piece lies in a disc whose box clears E₊ ∪ E_c. The planner is pole-aware, and the assembly
  re-checks the clearance exactly.
- W_l(−ω̄) = conj W_l(ω) (zeq_identity Z4) gives the left half's increment equal to H's. So the
  winding is (sum over H)/π, and the re-summed balls must lie within 1/4 of an integer.

**Results.**

All 20 certificates certify their winding. The pieces are exactly the 2562 planned segments, each
certified at its planned disc: there were no splits and no retries, and K ≤ 104 against K_MAX = 160.
Computed at `7168b2e4` (launch record `results/beta/part3/launch_20261027T163622Z.json`; 3756 s wall
with 4 workers at nice 10, 14 963 s of segment CPU) and assembled at `5ee4796f`.

| l | winding | (sum over H)/π | Z, P | pieces | a′_l/κ₋ | min abs(W_l) on H | K_max | CPU (s) |
|---|---|---|---|---|---|---|---|---|
| 0 | 0 | [+/- 8.69e-5] | 1, 1 | 41 | 12.912 | 0.0321 | 64 | 69 |
| 1 | −1 | [-1.000 +/- 1.09e-4] | 0, 1 | 43 | 13.650 | 0.351 | 80 | 84 |
| 2 | −1 | [-1.00 +/- 2.28e-4] | 0, 1 | 53 | 18.446 | 0.477 | 88 | 118 |
| 3 | −1 | [-1.000 +/- 1.81e-4] | 0, 1 | 63 | 23.241 | 0.6 | 88 | 160 |
| 4 | −1 | [-1.00 +/- 2.68e-4] | 0, 1 | 71 | 27.668 | 0.71 | 88 | 201 |
| 5 | −1 | [-1.00 +/- 2.05e-4] | 0, 1 | 81 | 32.464 | 0.817 | 88 | 254 |
| 6 | −1 | [-1.000 +/- 2.62e-4] | 0, 1 | 92 | 38.367 | 0.913 | 88 | 324 |
| 7 | −1 | [-1.00 +/- 2.90e-4] | 0, 1 | 101 | 42.425 | 0.986 | 88 | 381 |
| 8 | −1 | [-1.000 +/- 3.32e-4] | 0, 1 | 112 | 48.082 | 1.1 | 88 | 466 |
| 9 | −1 | [-1.000 +/- 4.59e-4] | 0, 1 | 123 | 53.738 | 1.13 | 96 | 554 |
| 10 | −1 | [-1.00 +/- 5.02e-4] | 0, 1 | 128 | 56.813 | 1.27 | 96 | 616 |
| 11 | −1 | [-1.000 +/- 5.37e-4] | 0, 1 | 141 | 62.223 | 1.31 | 96 | 729 |
| 12 | −1 | [-1.00 +/- 4.83e-4] | 0, 1 | 151 | 67.634 | 1.43 | 96 | 842 |
| 13 | −1 | [-1.000 +/- 5.16e-4] | 0, 1 | 163 | 73.045 | 1.5 | 96 | 979 |
| 14 | −1 | [-1.00 +/- 5.34e-4] | 0, 1 | 173 | 78.455 | 1.61 | 96 | 1110 |
| 15 | −1 | [-1.000 +/- 6.59e-4] | 0, 1 | 182 | 83.866 | 1.68 | 96 | 1253 |
| 16 | −1 | [-1.000 +/- 6.03e-4] | 0, 1 | 195 | 89.277 | 1.78 | 96 | 1429 |
| 17 | −1 | [-1.000 +/- 5.97e-4] | 0, 1 | 205 | 94.688 | 1.82 | 104 | 1597 |
| 18 | −1 | [-1.00 +/- 6.84e-4] | 0, 1 | 216 | 100.098 | 1.93 | 104 | 1792 |
| 19 | −1 | [-1.00 +/- 7.01e-4] | 0, 1 | 228 | 105.509 | 1.96 | 104 | 2004 |

**Accounting.**
- **l = 0:** winding 0 = Z − P. ω_* ∈ B₀ (§6), so Z ≥ 1, hence P ≥ 1, so P = 1 and Z = 1. Thus ω_* is
  the only zero of W₀ in B₀, it is simple, and W₀ has a simple pole at ω₀ (W̃₀ ≠ 0, again).
- **1 ≤ l ≤ 19:** winding −1 = Z − P, so P = Z + 1 ≥ 1. Hence P = 1 and Z = 0: no zero in B_l, and a
  simple pole at ω₀ (W̃_l ≠ 0, again).

**KATs and controls** (`kat_controls_PB.json`). Full-boundary windings at the KAT boxes (κ₋ units),
all as expected:

| KAT | l | box | contains | winding |
|---|---|---|---|---|
| B2 | 0 | [−0.03, 0.03] × [−0.93, −0.89] | ω_* | +1 |
| B3 | 0 | [−0.03, 0.03] × [−0.86, −0.79] | ω₀ (the pole) | −1 |
| B3 | 1 | the same | ω₀ (the pole) | −1 |
| B4 | 0 | [−0.03, 0.03] × [−0.93, −0.79] | ω₀ and ω_* | 0 |
| B6 | 0 | [−0.05, 0.05] × [−1.90, −1.80] | the overtone ≈ −1.855iκ₋ (below the edge) | +1 |
| B6 | 1 | [6.9, 7.5] × [−1.80, −1.70] | the PS mode ≈ 7.2149 − 1.7533iκ₋ (below the edge) | +1 |

The controls also behave as required:
- **NB1.** The l = 0 band box with its bottom edge at −0.905κ₋ (above ω_*) gives −1.
- **NB2.** An edge at Im ω = −0.9093809790κ₋ (Im ω_* to 10 digits) does not certify: "piece enclosure
  contains 0 at the maximal depth".
- **The symmetry.** The left-half sums at l = 0 and l = 2 overlap the stored right-half sums.
- **B1.** P2-K6's Krawczyk ball replays to the identical enclosure (contraction 5.3·10⁻¹⁵), and
  `krawczyk_PB_l0.json` is not regenerated (A4).
- **The 5% replay** (`replay_PB.json`). 141 stored pieces recomputed with their stored K all
  overlap.

**Level.** CERTIFIED counts. The written steps are routine: the argument principle for a meromorphic
function, the symmetry, and the at-most-simple pole of §2.3.

## 6. Assembly: the cover of S_e for every l

| region | l | covered by |
|---|---|---|
| S_δ = {Im ω ≥ −(11/20)κ₋} | all | P2-4 (`certificates/assembly/assembly_PB.json`, pinned) |
| Re ω < 0 | all | the symmetry W_l(−ω̄) = conj W_l(ω) |
| B ∩ {Re ω ≥ 0} | ≥ 20 | part (i) at depth e (§3) |
| B ∩ {abs(Re ω) ≤ a′_l} | ≤ 19 | the band counts (§5) |
| B ∩ {Re ω ≥ a′_l} | ≤ 19 | part (ii′) (§4) |
| ω₀ | all | Lemma E (§2) |

h₀ = 2/41 gives l₀ = ⌈1/h₀ − ½⌉ = 20, so l ≤ 19 and l ≥ 20 cover every degree.

**Conclusion.**
- The zeros of the W_l in S_e are exactly 0 and ω_*, both at l = 0 and both simple. ω_* lies in B₀
  (0.0106κ₋ above its bottom edge, abs(Re ω_*) ≤ 6.6·10⁻³⁵κ₋), off E₊ ∪ E_c, and below S_δ.
- By D1′'s set clause and Lemma E, the resonances in S_e are exactly 0 and ω_*.
- Hence α = −Im ω_* (attained), and β = −Im ω_*/κ₋ ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹].

**The machine check** (`src/p2/beta_assembly.py`, `certificates/beta/assembly_beta_PB.json`):
- 332 checks pass, at 256 bits with exact rationals, at `5a78bbcb` (clean).
- The six S2–S5 inputs are pinned by sha256, and each equals its blob on main `300c22a7` (A2).
  P2-4's 54 pinned inputs are unchanged. 38 inputs are pinned in all.
- What is checked:
  - C1: only ω₀ of E₊ ∪ E_c in S_e, and the band bounds;
  - C2: P2-4;
  - C3: Lemma E for l ≤ 19 and B5;
  - C4: part (ii′) for every l ≤ 19, at depth e;
  - C5: for every l ≤ 19, a′_l is identical to C4's; b′ ≥ eκ₋ (by 3.9·10⁻¹⁰); the 2562 stored
    pieces trace the right half of ∂B_l exactly, lie in their discs, and clear E₊ ∪ E_c; the
    re-summed windings are 0 and −1, consistent with C3;
  - C6: part (i) at depth e, and the cover replays identically;
  - C7: l₀ = 20;
  - C8: P2-K6's ball, ω_* ∈ B₀ (0.0106κ₋ above its bottom edge), off E₊ ∪ E_c and below S_δ,
    width 6.3·10⁻³¹;
  - C9: the KATs, the controls and the replay;
  - C10: the levels, and no `+dirty` commit tag in any S6 certificate.
- **Negative controls NB6a–NB6i** each fail on their named items:
  - a band certificate removed;
  - a′₇ lowered by 2⁻²⁴;
  - a winding set to 0 at l = 3;
  - b′ the dyadic just below eκ₋ (A5);
  - W̃₅ replaced by a ball containing 0;
  - part (ii′) at l = 9 only to depth 0.705κ₋;
  - h₀ = 2/43;
  - an S5 input off its main blob;
  - a β ball of width 2·10⁻⁷.
- `python -m p2.beta_assembly verify` recomputes the check identically.

## 7. Levels

| ingredient | level |
|---|---|
| P2-4 (S_δ) | WRITTEN_PROOF |
| P2-D1′, set clause | WRITTEN_PROOF |
| part (i) at depth e, l ≥ 20 (W1, W2, W5, W6′) | WRITTEN_PROOF |
| part (ii′), A at fixed h, l ≤ 19 (W1, W2, W4) | WRITTEN_PROOF |
| Lemma E, l ≥ 1 (§2.4) | WRITTEN_PROOF |
| Lemma E, l = 0 and l ≤ 19 (Arb) | CERTIFIED |
| part (iii) band counts, l ≤ 19 | CERTIFIED |
| ω_*'s Krawczyk ball (P2-K6) | CERTIFIED |

**P2-6 at P_B: WRITTEN_PROOF**, the minimum. FORMULATION §2 targets CERTIFIED. That would need P2-0's
written steps, W1/W2/W4/W5/W6′, D1′ and Lemma E's argument for l ≥ 1 to be machine-checked, as for
P2-4.

## 8. What S6 does not establish

- **No resonance multiplicity at ω_*.** Only the radial order m₀(ω_*) = 1 is certified. ω_* lies
  below −κ₊, outside D1′'s multiplicity window.
- **Nothing below −eκ₋.** The B6 KAT boxes at −1.855iκ₋ (l = 0) and 7.2149 − 1.7533iκ₋ (l = 1) are
  method tests, not claims.
- **Not P_A.** FORMULATION does not target a two-sided β at P_A.
- **No new conversion.** With β known, [HV17] could be applied with any α < −Im ω_* (and α < γ₀),
  which would sharpen P2-5's H^{1.05−0}. S6 does not claim this; P2-5 is unchanged.
- **Float evidence is not used.** The float check of (E-int) and KAT B7 (|W_l| at 100 random band
  points for l ∈ [20, 1000]) support the proof but are not part of it.
- **The level is below target:** WRITTEN_PROOF, not CERTIFIED.
- **Massless scalar on fixed RNdS only**, as in FORMULATION §7.
