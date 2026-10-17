# P2 S1: P2-0, D4 and Lemma D1′ (written proofs; machine-checked parts)

Status: S1, 2026-10-16. Conventions are those of `../FORMULATION.md` v0 §1. The exact
identities cited as C0–C10 are checked in `src/p2/p20_identity.py` and stored in
`certificates/p2_0_identity.json`; each comes with a negative control (the identity with a
deliberately wrong term must fail, and does). Replay: `python -m p2.p20_identity verify`.
Source checks (TeX sources fetched 2026-10-16): [Vas13] arXiv:1012.4391, [HV17]
arXiv:1512.08004, [Hin25] arXiv:2504.01734. Theorem numbers below are those of the arXiv
sources (counted from the TeX; labels given in brackets).

Notation. Δ = r²F, κ₊ = Δ′(r₊)/(2r₊²), κ_c = −Δ′(r_c)/(2r_c²). The hyperboloidal time is
t* = t − T(r), T′ = η/F with η(r) = (2r − r₊ − r_c)/(r_c − r₊) (so η(r₊) = −1, η(r_c) = 1,
|η| < 1 between). For φ = e^{−iωt}R(r)Y_lm = e^{−iωt*}u(r)Y_lm, u = e^{−iωT}R. L = l(l+1).

## 1. The radial operator and its hyperboloidal form

- **C0.** For g = −F dt² + F⁻¹dr² + r²dΩ²,
  □_g(e^{−iωt}R Y) = e^{−iωt}r⁻²[(ΔR′)′Y + ω²r⁴Δ⁻¹RY + R Δ_{S²}Y]. With Δ_{S²}Y_lm = −L Y_lm
  this is the radial equation (R) of FORMULATION §1.2, (ΔR′)′ + (ω²r⁴/Δ − L)R = 0.
- **C1.** With R = e^{iωT}u,
  e^{−iωT}[(ΔR′)′ + (ω²r⁴/Δ − L)R] = E_l[u] :=
  (Δu′)′ + iω(2r²ηu′ + (r²η)′u) + ω²r²(1−η²)F⁻¹u − Lu.
  So e^{iσt*}□_g e^{−iσt*}(uY_lm) = r⁻²E_l[u]Y_lm, for any σ = ω.
- **C6.** r²(1−η²)/F = 4r⁴/((Λ/3)(r_c−r₊)²(r−r_n)(r−r₋)), which is positive and bounded on
  [r₊, r_c]; E_l therefore has polynomial-over-(r−r_n)(r−r₋) coefficients, regular except at
  the zeros of Δ.
- **C7.** Res_{r₊} T′ = −1/(2κ₊) and Res_{r_c} T′ = −1/(2κ_c), so
  e^{−iωT} = (r−r₊)^{iω/(2κ₊)}e^{−iωT_reg} near r₊ (T_reg analytic) and likewise at r_c. The
  indicial exponents of E_l are {0, iω/κ₊} at r₊ and {0, iω/κ_c} at r_c.

Consequence (Lemma A). Let ω ∉ E₊ = {−inκ₊ : n ≥ 1}. Near r₊ (on both sides) the solutions of
E_l[u] = 0 are spanned by an analytic one, u₊ = e^{−iωT}R₊ = e^{−iωT_reg}g₊ with u₊(r₊) ≠ 0, and
a second one which is |r−r₊|^{iω/κ₊}·(analytic) on each side if iω/κ₊ ∉ ℤ, or contains
log|r−r₊| (ω = 0). The second solution is not C^∞ at r₊, and it lies in H^s near r₊ only for
s < ½ − Im ω/κ₊. (For iω/κ₊ ∈ ℕ, i.e. ω ∈ E₊, both solutions can be smooth: this is why D1 and
D1′ exclude E₊ ∪ E_c.) The same holds at r_c with E_c and κ_c.

## 2. P2-0 (closed upper half-plane)

**Claim.** For every l and every ω with Im ω ≥ 0: W_l(ω) = 0 ⟺ (l, ω) = (0, 0).

Proof. Write ω = a + ib. If W_l(ω) = 0, then R₊ = cR_c, so u := e^{−iωT}R₊ is analytic on
[r₊, r_c] (Lemma A at both ends; E₊ ∪ E_c lies in Im ω < 0). Then D = ∫r²F|u′|², N = ∫|u|²,
X = ∫r²η Im(u′ū), 𝒲 = ∫r²(1−η²)F⁻¹|u|² (finite by C6) and B = r_c²|u(r_c)|² + r₊²|u(r₊)|²
are finite. Integrate the pointwise identity **C2**,

  ū E_l[u] = J′ − r²F|u′|² − 2ωr²η Im(u′ū) + ω²r²(1−η²)F⁻¹|u|² − L|u|²,
  J = r²F u′ū + iωr²η|u|²,

over [r₊, r_c] (written step: fundamental theorem of calculus; u ∈ C¹[r₊, r_c]). By **C3**,
J(r_c) − J(r₊) = iωB. So 0 = iωB − D − 2ωX + ω²𝒲 − LN. By **C4** this splits into
−D − bB − 2aX + (a²−b²)𝒲 − LN = 0 and aB − 2bX + 2ab𝒲 = 0, and by **C5**

  b(D + LN + |ω|²𝒲) + |ω|²B = 0.        (E′)

- b > 0: every term of (E′) is ≥ 0, so 𝒲 = 0 (ω ≠ 0). The weight of 𝒲 is positive on
  (r₊, r_c) (C6), so u ≡ 0, a contradiction.
- b = 0, a ≠ 0: (E′) gives a²B = 0, so u(r₊) = 0, contradicting u₊(r₊) ≠ 0 (C7).
- ω = 0: the identity reads D + LN = 0, so u′ ≡ 0. For l ≥ 1 this forces u ≡ 0; for l = 0,
  u is constant. Conversely R ≡ 1 solves the l = 0, ω = 0 equation and is analytic at both
  horizons with value 1, so R₊ = R_c = 1 and W₀(0) = 0. ∎

This is the s0-gate §2 argument with every algebraic step machine-checked. (E′) is (E) of
s0-gate §2 multiplied by b, so no division by b is needed.

## 3. D4 (the zero mode is simple)

**Claim.** W₀′(0) = i(r₊² + r_c²)/(r₊r_c) ≠ 0, so m₀(0) = 1.

Proof. Let y = R₊(·; ω) and z = R_c(·; ω) (l = 0). They are analytic in ω near 0, since
0 ∉ E₊ ∪ E_c and the Frobenius coefficients are rational in ω with poles only on E₊ (resp. E_c).
At ω = 0, y = z = 1. Differentiate (ΔR′)′ + (ω²r⁴/Δ)R = 0 in ω at ω = 0 (**C8**): (Δẏ′)′ = 0, so
Δẏ′ is a constant, and it equals its limit at r₊. Write ẏ = −(i/(2κ₊))ln(r−r₊) + (analytic at
r₊), using R₊ = (r−r₊)^{−iω/(2κ₊)}g₊. Then Δẏ′ → −(i/(2κ₊))Δ′(r₊) = −ir₊² (C8, exact limit).
Likewise Δż′ = ir_c². Hence d/dω[Δ(yz′ − y′z)] at ω = 0 is ir_c² − (−ir₊²). The Jost factor
exp(−iωc₊ + iωc_c)/(r₊r_c) is entire and nonzero, and W₀(0) = 0, so
W₀′(0) = i(r₊² + r_c²)/(r₊r_c). ∎

Numerical cross-check: `tests/test_arbw.py::test_D4_derivative_matches_exact_formula` (a
central difference of validated Arb values agrees to < 10⁻³⁰ at both points).

## 4. Lemma D1′ (radial Frobenius QNMs = resonances of [Vas13], [HV17], [Hin25])

### 4.1 What the sources define (checked against the TeX sources)

- **[Vas13]** §2 (abstract framework), applied in §6 to Kerr–de Sitter spaces.
  - P_σ − iQ_σ : 𝒳^s → 𝒴^s = H^{s−1} is an analytic family of Fredholm operators on
    ℂ_s = {Im σ > β⁻¹(1 − 2s)} (Thm 2.6 [thm:classical-absorb]; for the wave operator the order
    is 2). Its inverse is meromorphic, and the resonances are its poles. The poles do not depend
    on s above threshold (Remark 2.8 [rem:inv-indep-of-s]).
  - §2.6: the nullspace on H^s above threshold is finite-dimensional and its elements are in
    C^∞(X).
  - In the expansion of Lemma 3.1 [lemma:Mellin-expand], u = Σ_j Σ_{κ≤m_j} τ^{iσ_j}(log τ)^κ a_{jκ}
    + u′ with a_{jκ} ∈ C^∞, "simple poles corresponding to m_j = 0". So m_j + 1 is the order of
    the pole.
- **[HV17]** §2.4–2.5.
  - The resonances are the poles of P̂(σ)⁻¹, where P̂(σ) = e^{iσt*}𝒫e^{−iσt*} is the
    Mellin-transformed normal operator family, on extendible/supported variable-order spaces
    (§2.4, Thm 2.13 [ThmRNdSFredholm]).
  - The relevant set is "the forward problem for the wave equation in a neighborhood Ω₂₃ of the
    black hole exterior", Ω₂₃ = {t₂₃ ≥ 0, r₂ − 2δ ≤ r ≤ r₃ + 2δ}, whose resonant states are
    smooth at r₂, r₃ by radial-point propagation above threshold (proof of Lemma 2.14
    [LemmaRNdSResonanceSupp]).
  - Their time function is t₂₃ = t − F₂₃(r), with F₂₃′ = s_j(μ⁻¹ + c_j) near r_j, s₂ = −1,
    s₃ = +1, and c_j smooth.
  - The hypothesis used in Prop 2.17 [PropRNdSPartialAsymp] and Thm 2.24
    [ThmRNdSPartialAsympConormal]: "The only resonance of the forward problem in Ω₂₃ in
    Im σ ≥ 0 is a simple resonance at σ = 0, with resonant states equal to constants, and there
    exists α > 0 such that 0 is the only resonance in Im σ ≥ −α … we may assume −α > −γ₀"
    (§2.5).
  - "Simple" is made precise in Lemma 2.14: order of the pole 1, rank 1.
- **[Hin25]** §1 (eq. [EqIBoxSol]).
  - σ ∈ QNM iff there is 0 ≠ u ∈ C^∞(X), X = [r₋', r₊'] × S² with r₋' ∈ (r_C, r_e) and
    r₊' > r_c, such that □_g(e^{−iσt*}u) = 0.
  - His t* = t − T, T′ = T̃/F. T̃ is printed as 2(r−r_e)/(r_c−r_e), which lacks the −1 needed
    for the normalization T̃(r_e) = −1, T̃(r_c) = +1 that he states. Our η is the normalized
    version.
  - Multiplicity: with order-1 poles, "the multiplicity of a QNM ς is equal to the dimension
    of the nullspace" (§5, after Thm 5.5 [ThmPfPN]).

Choice of time function. If t*₁ and t*₂ are two such time functions, then t*₁ − t*₂ = f(r) with f
smooth on a neighborhood of [r₊, r_c]. Near r₊ both have derivative −1/F + (smooth), and similarly
at r_c. So u₂ = e^{iσf}u₁ preserves smoothness, H^s spaces and polynomial-in-t* structure. The
resonances and their multiplicities therefore do not depend on the choice (ours, HV17's t₂₃, or
Hintz's). In particular the three sources define the same set and the same pole orders for this
operator.

### 4.2 Statement

Let σ ∉ E₊ ∪ E_c. Then:

1. **(set)** σ is a resonance (in the sense of any of the three sources) iff W_l(σ) = 0 for some l.
   The resonant states are spanned by u_l(r)Y_lm, m = −l…l, over the l with W_l(σ) = 0, where
   u_l = e^{−iσT}R₊ extended analytically across both horizons. The geometric multiplicity is
   Σ_{l: W_l(σ)=0}(2l+1).
2. **(order and multiplicity)** If moreover Im σ > −min(κ₊, κ_c), then:
   - the order of the pole of P̂(σ)⁻¹ at σ equals max_l m_l(σ);
   - in each (l, m) sector there is exactly one Jordan chain, of length m_l(σ);
   - the space of generalized resonant states (solutions Σ_{k≤K} t*^k e^{−iσt*}a_k, a_k smooth)
     has dimension Σ_l (2l+1)m_l(σ).
   This is the multiplicity (2l+1)m_l(σ) of D1.

The strip S_δ lies in {Im σ > −min(κ₊, κ_c)} at both frozen points. This is certified:
`params.json`, `strip_vs_kp` and `strip_vs_kc`.

Scope note (relative to FORMULATION D1′). The set equality (1) holds for all σ ∉ E₊ ∪ E_c. The
multiplicity statement (2) is proved here on the half-plane Im σ > −min(κ₊, κ_c), where the
horizon Abel integrals converge absolutely (C10). This is the only region P2 uses. (Hadamard
finite parts would extend (2) to all σ ∉ E₊ ∪ E_c; that extension is not needed and is not claimed.)

### 4.3 Proof

**Separation.** P̂(σ) commutes with the rotations, and the function spaces of all three sources
are rotation invariant. Projection onto Y_lm commutes with P̂(σ), and by C0/C1
P̂(σ)(uY_lm) = r⁻²E_l[u]Y_lm. A resonant state, which is smooth by [Vas13] §2.6 and the proof of
[HV17] Lemma 2.14, has smooth spherical-harmonic components. Each component solves E_l[u] = 0 on
[r₊ − ε, r_c + ε] and lies in H^s above threshold.

**(1)** By Lemma A at r₊, a solution of E_l[u] = 0 on a two-sided neighborhood of r₊ that is in
H^s with s > ½ − Im σ/κ₊ has no component on the singular branch, on either side. So it is a
multiple of u₊ (extended analytically through r₊; the Frobenius series converges for
|r − r₊| < r₊ − r₋). The same holds at r_c with u_c := e^{−iσT}R_c. A nonzero component therefore
exists iff u₊ ∝ u_c, i.e. iff W_l(σ) = 0. The FORMULATION Wronskian is a nonzero multiple of
μ(u₊u_c′ − u₊′u_c) by C1/C9. Conversely, if W_l(σ) = 0, then u_l(r)Y_lm, extended analytically
past both horizons, is a smooth resonant state.

Beyond the horizons the extension is unique for HV17's forward problem. There the region is
hyperbolic in r and the artificial final boundaries carry no condition. Hintz's X needs only
r₋' > r_C, and the series at r₊ converges up to r₋ = r_C. Only finitely many l can contribute,
by finite dimensionality ([Vas13] §2.6). This proves the set statement and the geometric
multiplicity, since each sector kernel is one-dimensional.

**(2), Jordan chains ↔ generalized resonant states.** Let e^{−iσ₀t*}Σ_{k≤K} t*^k a_k ∈ ker □ with
a_K ≠ 0. This holds iff there is a polynomial family a(σ) with a(σ₀) ≠ 0 and
P̂(σ)a(σ) = O((σ−σ₀)^{K+1}). Write U = ∮ e^{−iσt*}a(σ)(σ−σ₀)^{−K−1}dσ/2πi; then □U is the
corresponding contour integral of e^{−iσt*}P̂(σ)a(σ). Both directions are this computation.

**Lower bound: there is a chain of length m = m_l(σ₀).** Take a second solution v(σ) with
D(σ) := μW(u_c, v) ≠ 0 near σ₀. Then u₊ = αu_c + β̃v on the interior, with β̃ = −W̃/D,
W̃ := μW(u₊, u_c), and α holomorphic. Since W̃ ∝ W_l, β̃ vanishes to order exactly m. With a
cutoff χ (= 1 near r ≤ r₁, = 0 near r ≥ r₂, r₊ < r₁ < r₂ < r_c), set a(σ) = χu₊ + (1−χ)αu_c.
This is smooth on the extended domain, a(σ₀) = u₊(σ₀) ≠ 0, and
P̂(σ)a(σ) = [P̂(σ), χ](β̃v) = O((σ−σ₀)^m). Since this right-hand side is smooth and compactly
supported in the interior, P̂(σ)⁻¹ has a pole of order ≥ m at σ₀.

**Upper bound: pole order ≤ m.** Let f be smooth. For σ ≠ σ₀ near σ₀, a(σ) := P̂(σ)⁻¹f is smooth
([Vas13] §2.6 regularity). By C9 and C10 the Abel weight satisfies μ/p₂ ~ |r − r_h|^{−iσ/κ_h}, which
is integrable at r_h when Im σ > −κ_h. The boundary terms μW(a, u_h) vanish at the horizon r_h, so

  μW(a, u_c)(r) = −∫_r^{r_c} μu_c f/p₂,   μW(a, u₊)(r) = ∫_{r₊}^{r} μu₊ f/p₂.

Both are holomorphic in σ. Solving for the coefficients of a in the basis {u₊, u_c} on the
interior gives a(r; σ) = [holomorphic]/W̃(σ). So the Laurent coefficients of P̂(σ)⁻¹f of order
below −m vanish on (r₊, r_c).

These coefficients are generalized resonant states: the leading one lies in ker P̂(σ₀), and the
next ones satisfy the chain equations. Each one that vanishes on the exterior vanishes identically,
by uniqueness across the horizons (Lemma A and H^s above threshold). By induction they all vanish.
Smooth f are dense in 𝒴^s and the Laurent coefficients are bounded operators, so the pole order
is ≤ m.

Together: in each sector the geometric multiplicity is 1 and the pole order is m_l(σ₀). So there is
exactly one chain, of length m_l. The generalized resonant states of the sector are spanned by the
m_l states built from it (t*-degrees 0, …, m_l − 1). Summing over (l, m) gives the dimension
Σ(2l+1)m_l(σ₀), and the pole order of the full family is max_l m_l. ∎

### 4.4 Consequence for the conversion P2-5 (HV17 hypotheses)

At a frozen point X, suppose P2-4 holds: no zero of any W_l in S_δ∖{0}, together with P2-0, D4
and the certified inequalities (½+δ)κ₋ < min(κ₊, κ_c, γ₀). Then by D1′:

- 0 is the only resonance with Im σ ≥ −α, α = (½+δ)κ₋;
- the pole at 0 has order m₀(0) = 1 (D4), the l ≥ 1 sectors do not contribute (P2-0), and the
  resonant states are the constants (l = 0 sector, u = 1);
- α < γ₀.

These are exactly the hypotheses quoted from [HV17] §2.5 in §4.1. Nothing here re-proves [HV17];
P2-5 remains a theorem by citation.
