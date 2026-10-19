# P2 S4: part (i) at P_B, uniformly in l ≥ 20 — written proof of the regime certificates

Status: S4, 2026-10-18. This is the proof behind `src/p2/part1.py`, `part1_T.py`, `throat.py`,
`lowfreq.py` and `affine.py`, and behind the certificates `certificates/part1/*.json`.
Pre-registration: `notes/s4-preregistration.md` (with Amendments 1–6).
Replay: `python -m p2.part1 verify`. Conventions: FORMULATION v0.6. Point: P_B (ΛM² = 3/50,
Q = 0.999 Q_max, δ = 1/20).

Notation.
- h = 1/(l+½), so l ≥ 20 ⟺ h ∈ (0, h₀] at the discrete values, with h₀ = 2/41.
- The strip depth is d_δ = (½+δ)κ₋ = 0.55κ₋.
- ω = Re ω + iν, with ν ∈ [−d_δ, 0].
- ζ = h·Re ω, and ε = 1/Re ω.
- x = r* is the tortoise coordinate; Φ = rR solves −Φ_xx + (V_l − ω²)Φ = 0, with y := Φ_x/Φ.
- U = F/r² and W̃ = FF′/r − U/4, where F = Δ/r².

## Theorem (S4)

  For every l ≥ 20 and every ω with Re ω ≥ 0 and −d_δ ≤ Im ω ≤ 0, W_l(ω) ≠ 0.

**Corollary (P2-1 at P_B, l₀ = 20).** For l ≥ 20, W_l has no zero in S_δ. The parts are:
- Im ω ≥ 0: P2-0, at P2-0's level.
- Re ω < 0: the conjugation symmetry W_l(−ω̄) = conj W_l(ω) (FORMULATION §1.2; `zeq_identity`
  Z4).

**Level: WRITTEN_PROOF.**
- The inequalities are machine-checked in Arb.
- W1 and W2 below are non-routine written arguments, as in P2-1w.
- The Im ω > 0 part is at P2-0's level.

## 0. Exact facts used throughout (`bulk_PB.json`, W8)

- **The splitting.** V_l = h⁻²U + W̃ exactly, since L = l(l+1) = h⁻² − ¼; checked in sympy.
- **U is monotone.** Δ′r − 4Δ = −2(r² − 3r + 2Q²) (sympy), so U_r = −2(r − r_P)(r − r_P′)/r⁵ with
  r_P′ = 3 − r_P.
  - r_P′ = 1.04127 < r₊ = 1.09435 (Arb). Hence U is strictly increasing on (r₊, r_P) and strictly
    decreasing on (r_P, r_c).
- **The W̃ bound.** W̃ ≥ −W_bound on [r₊, r_c], with W_bound = 99/12500 = 0.00792. The certified
  lower bound over 6000 cells is −0.0079101. The float minimum is −0.0078910, at r = 3.9461 (review
  W2).
  - C_W := 4·max(κ₊², W_bound) = 0.03168.
- **Bulk positivity.** On a real interval where h⁻²U ≥ M and M ≥ max(4|ω|², C_W):
  Re(V − ω²) = h⁻²U + W̃ − (Re ω)² + ν² ≥ M − M/4 − M/4 > 0.

## 1. The cover (W9; `cover_PB.json`)

| regime | domain (every h ∈ (0, h₀]) | certificate |
|---|---|---|
| L | 0 ≤ Re ω ≤ 8κ₋ | `L_PB_plus.json`, `L_PB_c.json` |
| S | 6κ₋ ≤ Re ω ≤ 0.9Ω_P/h | `S_PB_plus.json`, `S_PB_c.json` |
| T | 0.8Ω_P/h ≤ Re ω ≤ 2.0Ω_P/h | `T_PB.json` |
| A | Re ω ≥ 1.5Ω_P/h | `A_PB.json` |

The certified inequalities are:
- 6 < 8;
- 0.8 < 0.9;
- 1.5 < 2.0;
- 0.9Ω_P/h₀ > 8κ₋;
- d_δ < κ₊ < κ_c (the poles E₊ ∪ E_c lie below the strip);
- 6κ₋ > d_δ, 0.8Ω_P/h₀ > d_δ and 1.5Ω_P/h₀ > d_δ (Re ω + Im ω > 0 in S, T and A).

So every (h, ω) with h ≤ h₀ and ω in the closed half-strip lies in at least one regime. Each
certificate's parameter boxes tile its domain; `cover_PB.json` checks the tiling.
- T's 24 ζ-boxes tile [0.8, 2.0]Ω_P.
- S's ζ-leaves tile [0, 0.9]Ω_P. In each leaf the ε-boxes tile [0, ε_max(ζ_lo)], where
  ε_max(ζ_lo) ≥ min(1/(6κ₋), h₀/ζ) ≥ ε for every (h, Re ω) in S.
- L's leaves tile [0, 8] × [−0.55, 0] (κ₋ units). Since `279878d4` the cover checks this box by
  box: in every grid box of side 1/20, the leaves are disjoint and fill it (review W3).

## 2. W1: vertical-strip analyticity of r(x) (S2 Lemma 1, extended by Amendment 1)

**Lemma W1.** The inverse r(x) of r* continues analytically, by vertical continuation from the
real axis, to each open vertical strip, at every height:
- Π₋ = {Re x < X₀},
- Π₀ = {X₀ < Re x < 0},
- Π₊ = {Re x > 0},

where X₀ = Re r*(0) = −2.29664. The values avoid {0, ∞}. Continuation across the line Re x = X₀
(resp. 0) is possible through the segments |Im x| ≤ 0.1413 (resp. 2.438), and r is analytic
there.

*Proof.* This is S2 Lemma 1's argument with the sector replaced by the strip.
- Along any continuation Re r*(r(x)) = Re x, and Σ A_i = 0. So the continuation can fail only
  where r → 0 (Re x = X₀) or r → ∞ (Re x = 0).
- Off those lines the vertical continuation exists at every height. Its x-derivative solves
  ∂_X R = F(R) (S2 Lemma 1, last bullet), so it is analytic in x.
- On the two lines, S2 integrated r along the segments up to heights 0.1413 and 2.438
  (`certificates/m3/m3_PB_geometry.json`, replayed as K1), with |r| bounded away from 0 and ∞.
  Conjugation (r(x̄) = conj r(x)) gives the segments below the axis.
- Each strip is simply connected, so the continuation is single-valued there.
- **The union is simply connected** (added after the review, finding S3).
  - Let N₀ and N₁ be open neighbourhoods of the checked windows
    {X₀ + iY : |Y| ≤ 0.1413} and {iY : |Y| ≤ 2.438} on which r continues analytically. They exist
    because the vertical continuation depends continuously on X near X₀ and 0, the windows stay in a
    compact subset of ℂ∖{0}, and Cauchy–Riemann holds as above.
  - Then D := Π₋ ∪ Π₀ ∪ Π₊ ∪ N₀ ∪ N₁ is open and connected. Its complement in the Riemann sphere is
    {∞} together with the four rays {X₀ + iY : |Y| > 0.1413} and {iY : |Y| > 2.438}, minus the parts
    inside N₀, N₁. These pieces are joined at ∞, so the complement is connected and D is simply
    connected.
  - Hence r, and every solution of the radial equation in x, is single-valued on D. Equivalently,
    the region between the real axis and any S contour is a curvilinear trapezoid inside D.
  - So Φ_c continued along a c-side contour, which may cross Re x = X₀ and Re x = 0 inside the
    windows, is the same function as Φ_c continued vertically. On the tail half-plane of W2 that
    function is the Frobenius representation (∗). ∎

**Where each contour lies** (`w1_check` in each S leaf; S2's geometry for T and A):
- T and A use S2's contour, which lies in S2's sectors.
- **S, side +.** Re x(t) ≤ X_a < x_P < X₀, so the contour lies in Π₋. X_a < X₀ is certified per
  leaf.
- **S, side c.** Re x(t) = X_a + (t_a − t)/k grows as t decreases.
  - If X_a < X₀, the contour crosses Re x = X₀ at t* = t_a − k(X₀ − X_a). The crossing height
    |Im x(t*)| = |η(t*)|/k is enclosed and certified < 0.1413.
  - It crosses Re x = 0 at a height certified < 2.438.
  - Both crossings occur within k·(X − X_a) ≤ 0.22 of the exit in t, where the contour is real or
    barely bent. The certified heights are at most 0 at Re x = X₀ and 0.054 at Re x = 0; they are
    recorded per leaf.
- **Closing box.** The contour stays in the half-plane {|e^{±2κ_h(x − c_h)}| < R} (§8 below),
  which lies in Π₋ (side +) resp. Π₊ (side c). This is certified through r*(r_h ± d_r).

## 3. W2: the Jost solutions on the contours, and their decay

**Lemma W2.** Let Φ₊ and Φ_c be the FORMULATION §1.2 Jost solutions (Φ = rR). Continued along
the contours of §§4–6, they are analytic. Along the tails they decay with their t- (or s-)
derivatives, provided Re ω + Im ω > 0:
- T and A: S2's contour x_P + e^{iπ/4}s, as s → ±∞;
- S, side +: the ray x(t) = X_a + (t(1+i) − t_a)/k, as t → −∞;
- S, side c: x(t) = X_a − (t(1+i) − t_a)/k, as t → −∞;
- the closing boxes: the same rays from x_T(ζ), with k = κ_h.

*Proof.* This is S2 Lemma 3's argument, which uses only Re x ≤ X₋ (resp. ≥ X₊) at any height.
- **The disc around r₊.** For Re x ≤ X_L, the identity Re r* = Re x and a first-exit argument give
  |r − r₊| ≤ d. The disc bound is computed per leaf (`tail_disc`): d ≤ 10⁻⁴ in all S leaves.
  - Inside the disc, Φ₊ = c·r·g₊(r)·exp(−iω(x − c₊ − R̃₊(r))), with g₊ and R̃₊ analytic.
  - The right side is analytic in x on the strip region {Re x ≤ X_L} ⊂ Π₋, which is connected
    and contains a real half-line. So the identity theorem gives (∗) on the whole region.
- **Decay.** On the + ray, Im(ωx(t)) = const + (Re ω + Im ω)t/k. So
  |Φ₊(x(t))| ≤ C e^{(Re ω + Im ω)t/k} → 0. The same holds for u′ = x′Φ_x.
  - In S, Re ω ≥ 6κ₋ > d_δ ≥ −Im ω.
  - In T and A, Re ω ≥ 0.8Ω_P/h₀ > d_δ.
- **The c side** is the same at r_c. The tail region {Re x ≥ X_L} lies in Π₊, since
  X_L = X_a + (t_a − t_L)/k > 0 in every leaf; this is checked at assembly.
  - The contour reaches the tail from a real exit in Π₋ or Π₀ through the windows. By the simply
    connected union D of W1, the function continued along it is the one represented by (∗) on the
    tail.
- **The closing box.** r(x) = r_h + s·ψ(e^{±2κ_h(x − c_h)}) on that half-plane, by the identity
  theorem (both sides are analytic and agree on the real axis). (∗) holds there directly. ∎

## 4. W3: half-line coercivity on a curved contour (regime S)

Let x(t), t ≤ t_a, be C², with x(t_a) real and x′(t_a) real, and let φ be C¹ with φ(t_a) = 0. Let
u(t) := Φ(x(t)), with u and u′ decaying as t → −∞ (W2). Then u solves

  −u″ + (x″/x′)u′ + Q u = 0,   Q = x′²(V − ω²).

Multiply by e^{−iφ}ū, integrate over (−∞, t_a], integrate the first term by parts and take real
parts:

  Re(e^{−iφ(t_a)} u′ū)(t_a) = ∫ cos φ |u′|² + Re ∫ e^{−iφ}(x″/x′ − iφ′)u′ū + ∫ Re(e^{−iφ}Q)|u|².

**Lemma W3.** Let κ_tot := |x″/x′ − iφ′| and cos φ > 0, and use Young's inequality with the
pointwise weight cos φ/κ_tot. Then

  Re(e^{−iφ(t_a)} u′ū)(t_a) ≥ ∫ q |u|²,   q = Re(e^{−iφ}Q) − κ_tot²/(4 cos φ).

Where κ_tot = 0, the cross term vanishes. Suppose q ≥ μ > 0 on (−∞, t_a] and u ≢ 0. Then the
left side is positive. With φ(t_a) = 0 and u′ = x′Φ_x, this gives
  s·Re y(x_a)·|Φ(x_a)|²/k > 0,
where x′(t_a) = s/k with s = +1 on the + side and s = −1 on the c side. Hence Φ(x_a) ≠ 0, and
  Re y₊(x_a) > 0 (side +),   Re y_c(x_a) < 0 (side c). ∎

**The normalised form.** x = X_a + s(Y − t_a)/k, so x′² = Y′²/k² and x″/x′ = Y″/Y′. Also
V − ω² = (Re ω)²[U/ζ² − 1 + 2iε|ν| + ε²(ν² + W̃)] exactly. Hence q = (Re ω/k)² q̂, with

  q̂ = Re(e^{−iφ}Y′²[U/ζ² − 1 + 2iε|ν| + ε²(ν² + W̃)]) − (kε)²|Y″/Y′ − iφ′|²/(4 cos φ).

h enters only through ζ and ε (Amendment 3).

## 5. W4: global coercivity with a C¹ ground state (regimes T and A)

S2 §§3–4 gives, for S2's contour and a C¹ phase φ with c := cos φ − η|φ′| ≥ 0, η = ½,

  J(u) := Re ∫ e^{−iφ}(−u″ + Qu)ū ≥ ∫ E|u|²,   E = q + c′G′ − cG′² + cG″,

for every g = e^{−G}. Here G′ must be C¹ (S2 had piecewise linear G′; the integration by parts
needs only absolute continuity of G′), q = Re(e^{−iφ}Q) − |φ′|/(4η) and Q = i(V − ω²).
- A QNM gives a nonzero decaying u with J(u) = 0 (W2). So E ≥ μ > 0 on ℝ excludes it.

**The semiclassical splitting (W8).** Put G′ = Γ/h and ω = ζ/h + iν. Then exactly

  h²E = E₀ + hE₁ + h²E₂,
  E₀ = Re(ie^{−iφ}(U − ζ²)) − cΓ²,
  E₁ = 2ζν cos φ + (cΓ)′,
  E₂ = Re(ie^{−iφ}(W̃ + ν²)) − |φ′|/(4η).

This is checked against S2's E at sample points to 10⁻²⁰ (`test_part1.py`).
- If E₀ ≥ 0 and F := E₀/h₀ + E₁ + h₀·min(0, E₂) ≥ μ, then E ≥ μ/h for every h ∈ (0, h₀].
- **Monotonicity in the parameters.** sin φ ≤ 0 and cos φ > 0, so:
  - E₀ is nondecreasing in ζ²;
  - E₁ is nonincreasing in ζ and nondecreasing in ν (ν ≤ 0);
  - E₂ is nonincreasing in ν².
  - Hence on a box [ζ_lo, ζ_hi] × [ν_lo, 0], F ≥ E₀(ζ_lo)/h₀ + E₁(ζ_hi, ν_lo) + h₀·min(0, E₂(ν_lo)).
    No hull is lost.
- **T's design (Amendment 4).**
  - φ = −1.4375 s⁴/(7⁴ + s⁴).
  - Γ = (1 − 0.05)a_d·s·θ(s), with a_d = 32178552/2³⁰. θ is a C¹ cosine taper: 1 on |s| ≤ 2 and
    0 on |s| ≥ 4.
  - Near s = 0, E₀ = s²Ẽ₀, where U − U_P ∈ (s²/2)·conv U_ss([0, s]) because U_s(0) = 0. So
    Ẽ₀ ≥ 0 on a cell with an endpoint at 0 gives E₀ ≥ 0 there.
- **A (W7).** φ = −0.2 + (−1.4375 + 0.2)s⁴/(7⁴ + s⁴) and Γ = 0. Then ∂q/∂Re ω = 2|ω| sin(arg ω − φ).
  - This is ≥ 0, because φ ≤ φ₀ = −0.2 < −arctan(d_δh₀/ζ_A) ≤ arg ω ≤ 0 < φ + π (all certified).
  - So E ≥ μ_A/h at Re ω = ζ_A/h implies E ≥ μ_A/h for every larger Re ω.
- **The tails** |s| ≥ 64, 32 use S2's horizon discs (Γ = 0 there). E₀ ≥ ζ_lo²|sin φ(S)| − sup|U|
  there.

## 6. W5: Riccati invariance in the forbidden bulk (`s0-gate.md` §4)

**Lemma W5.** Suppose Re(V − ω²) > 0 on a real interval [x₁, x₂].
1. If Re y₊(x₁) > 0, then Φ₊ has no zero on [x₁, x₂], and Re y₊ > 0 there.
2. If Re y_c(x₂) < 0, then Φ_c has no zero on [x₁, x₂], and Re y_c < 0 there.

*Proof.*
- y′ = V − ω² − y², so (Re y)′ = Re(V − ω²) − (Re y)² + (Im y)².
- **Part 1.** At a first point where Re y₊ = 0, (Re y₊)′ > 0; this contradicts reaching 0 from
  above.
  - Suppose Φ₊ has a first zero x* in (x₁, x₂], with Re y₊ > 0 on [x₁, x*).
  - The zero is simple (else Φ₊ ≡ 0), so y₊ ~ 1/(x − x*) → −∞ as x → x*⁻. Contradiction.
- **Part 2** is the mirror argument (leftward). ∎

**Use.** If both hold on [x_a⁺, x_a^c] ∋ x_P, then
  Re y₊(x_P) > 0 > Re y_c(x_P).
So y₊ ≠ y_c and Φ₊, Φ_c ≠ 0 at x_P. Hence W_l(ω) ∝ Φ₊Φ_c(y_c − y₊)(x_P) ≠ 0.

## 7. W6: the scaled Frobenius series is the Jost solution, uniformly in h (regime L)

At r_h, with s = ±1, put r = r_h + s h²τ. Then Δ = h²D(τ), with
D = Σ_k d_k s^k h^{2(k−1)} τ^k.

The radial equation Δ²R″ + ΔΔ′R′ + (ω²r⁴ − LΔ)R = 0 becomes exactly

  D²R_ττ + sDE R_τ + (ω²r⁴ − (1 − h²/4)D) R = 0,   E = Δ′(r).

This is polynomial in τ with coefficients polynomial in h².
- **The solution.** The indicial roots are ±iω/(2κ_h). The FORMULATION Jost solution is
  R = (s(r − r_h))^ρ S(r), with ρ = −iω/(2κ_h), S analytic and S(r_h) = 1. So
  R = (h²)^ρ τ^ρ Σ a_n τ^n, with a₀ = 1 and a_n from the recurrence with denominator
  d₁²n(n + 2ρ).
  - n + 2ρ ≠ 0, since Im ω > −κ_h.
  - The series converges for τ < (distance from r_h to the nearest other root)/h². It is
    certified at τ_a by arbw's majorant bound, whose derivation uses only the degrees of
    D², DE and the C-polynomial.
- **The certified quantity.**
  y = Φ_x/Φ = h²D/r³ + s(D/r²)·R_τ/R. The factor (h²)^ρτ^ρ cancels in R_τ/R.
- **Uniformity.** The enclosure holds for every ω in the box and every h² in each of the 16
  sub-intervals of [0, h₀²]. These contain every h = 1/(l+½), l ≥ 20, and the Bessel limit h = 0.
- **Affine arithmetic (`affine.py`).** Each quantity is c₀ + c₁δω + c₂δh + E, with |δω| ≤ R_ω,
  |δh| ≤ R_h and |E| ≤ e.
  - Sums are exact.
  - A product keeps the linear terms exactly. All second-order terms are bounded by
    dev(x)·dev(y) into e.
  - The inverse is 1/c₀ − L/c₀² + E′, with |E′| ≤ e/|c₀|² + s²/(|c₀|²(|c₀| − s)) and s ≥ |L + E|.
  - All coefficients are Arb balls, so this is a rigorous enclosure (tested against sampled
    points). ∎

## 8. The S closing box ζ ∈ (0, 0.02Ω_P] (the exponential throat)

- **The contour.** Take k = κ_h exactly and x(t) = x_T(ζ) + sY(t)/κ_h. Then
  ρ := e^{±2κ_h(x − c_h)} = ρ_T e^{2Y}, where ρ_T = ρ(x_T(ζ)) ∈ (0, ρ₁].
- **The near-horizon map.**
  - ρ(r) = s(r − r_h)f(r), with f analytic and f(r_h) = 1.
  - On |r − r_h| ≤ d_r, |f| ≥ m_f (term-wise distance bounds).
  - By Rouché, for |ρ| < R := m_f·d_r the inverse ψ exists, with ψ(ρ) in that disc.
- **G.** G(ρ) := U(r_h + sψ(ρ)) is analytic there, with |G| ≤ M_G (the sup of |U| over 24×24
  sub-boxes).
  - G = Σ_{n≤14} g_nρⁿ + E₁₄, with g_n from series reversion and composition (Arb),
    and |E₁₄| ≤ M_G(|ρ|/R)¹⁵/(1 − |ρ|/R).
- **The turning point.** G is increasing on [0, ρ₁] and G(ρ₁) ≥ ζ₁², so ρ_T(ζ) ≤ ρ₁ for every
  ζ ≤ ζ₁.
- **Ũ.** Ũ = G(ρ_T e^{2Y})/G(ρ_T) = e^{2Y}P(ρ_T e^{2Y})/P(ρ_T), with P = G/ρ.
  - On the radial segment, P(λz) ∈ g₁ + Σ g_n zⁿ⁻¹·[0, 1] for λ ∈ [0, 1], for complex z.
  - **The error terms** need |P_N(ρ)| bounded on the complex set the contour visits,
    ρ = λρ₁e^{2Y(t)}. arg ρ = 2η(t) covers [−π, 0] on the bend, and the g_n have mixed signs.
    - Lower bound: the larger of |g₁| − Σ_{n≥2}|g_n||ρ|^{n−1} and the complex-segment enclosure's lower
      end. Upper bound: the smaller of |g₁| + Σ_{n≥2}|g_n||ρ|^{n−1} and its upper end.
    - The tail uses the disc bound.
    - The first version bounded |P_N| on the real segment only, which is not a bound at complex ρ
      (review S1). It was corrected at `68de8184`; the regenerated margins agree with the first ones to
      8 digits.
  - W̃(r_h) = 0 gives |W̃| ≤ M_W(|ρ|/R)/(1 − |ρ|/R).
  - ρ₁e^{2t_a} < R is certified.
- **The check.** q̂ is then an explicit function of t, with ρ_T as an interval parameter. The
  mean-value form is applied to the polynomial model, and the remainder bounds are taken as sups
  over each cell. No contour integration is needed.

## 9. Proof of the theorem

Fix l ≥ 20, so h ≤ h₀, and ω with Re ω ≥ 0 and ν ∈ [−d_δ, 0]. By §1, (h, ω) lies in one of the
regimes.

- **L** (Re ω ≤ 8κ₋). Take an L-leaf box containing ω.
  - The certificates give Re y₊ > 0 at r₊ + h²τ_a and Re y_c < 0 at r_c − h²τ_a^c. Both hold for
    all h² ∈ [0, h₀²] ∋ h² (W6, §7).
  - The exit rule gives h⁻²U ≥ M ≥ max(4|ω|², C_W) at both exits.
  - By U's monotonicity this holds on the whole interval between them, so §0 gives
    Re(V − ω²) > 0 there.
  - W5 then gives W_l(ω) ≠ 0.
- **S** (6κ₋ ≤ Re ω ≤ 0.9Ω_P/h). Put ζ = h Re ω ∈ (0, 0.9Ω_P] and ε = 1/Re ω.
  - Then ε ≤ min(1/(6κ₋), h₀/ζ), so (ζ, ε) lies in an S-leaf (ζ-box × ε-box) on each side.
  - On each side, W2 and W3 with q̂ ≥ μ > 0 (certified on every t-cell and on the tail) give
    Re y₊(x_a⁺) > 0 and Re y_c(x_a^c) < 0.
  - The certified exit conditions give x_a⁺ < x_P < x_a^c and Ũ(x_a) − 1 > ε_max²W_bound.
  - By U's monotonicity, Ũ ≥ min(Ũ(x_a⁺), Ũ(x_a^c)) on [x_a⁺, x_a^c]. So there
    Re(V − ω²) = (Re ω)²[Ũ − 1 + ε²(ν² + W̃)] > 0.
  - W5 gives W_l(ω) ≠ 0.
- **T** (0.8 ≤ ζ/Ω_P ≤ 2.0). The ζ-box containing ζ gives E ≥ μ_T/h > 0 on S2's contour (W4, §5),
  so J(u) > 0 for every nonzero decaying u. A zero of W_l would give one with J(u) = 0 (W2).
- **A** (ζ ≥ 1.5Ω_P). Likewise, through W7.

Hence W_l(ω) ≠ 0. ∎

## 10. What the computer checks

*(Final numbers in `notes/s4-report.md` §2; certificate files in `certificates/part1/`.)*

- **T.** 24 ζ-boxes in [0.8, 2.0]Ω_P, with ν ∈ [−d_δ, 0] and h ∈ (0, h₀].
  - Cells are mean-value forms of E₀, E₁ and E₂ (E₀ also in second order), with the s²-form at 0.
  - Cells are bisected to a*/10, depth ≤ 12.
  - Tails as in S2. c ≥ 0 on every cell.
- **A.** The box ζ = ζ_A, the same checks, and W7.
- **S.** For each side: the closing box plus 44 ζ-boxes of width 0.02Ω_P, halved where needed.
  - Each leaf has its design, (X_a, k), the validated contour, mean-value cells, the tail disc,
    the exit and order checks, and the W1 check.
- **L.** 1760 boxes of side κ₋/20 on [0, 8] × [−0.55, 0] (κ₋ units), quadrisected on failure,
  with 16 h²-sub-intervals each.
- **Bulk, cover, controls (N1–N5, K1–K4) and replay**: as listed in the pre-registration.
