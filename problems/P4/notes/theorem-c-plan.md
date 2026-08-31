# P4 S4 — Theorem C (exclusion outside the rectangle): scoping plan

Status: scoping, 2026-08-31. No new certified claims; every number here is float
unless it carries a ball radius. Sizing code: `src/p4/validated/tc_sizing.py`
(float only, ~10 min; `PYTHONPATH=problems/P4/src uv run python -m p4.validated.tc_sizing`).
Inputs: FORMULATION §2 (Theorem C template), `s2-theorem-b.md` (E, Theorem B′),
the S1 monotonicity diagnostics (ledger P4-3), and the Stage-1–3 validated code.

## 1. Statement and what it adds

**Theorem C (target).** Let E(κ) be the matching function of Theorem B′
(`s2-theorem-b.md` §3.8: E(κ) = e^{3x_d} det[r₁, r₂, p̃(x_d; κ)], x_d = −3, on the
certified EC background, KHA sonic-point gauge, class = analytic at the sonic
point, on the linearised constraint surface, regular at the centre). Then E has
no zeros in
    Ω := {Re κ ≥ 0} \ R = {Re κ > 15} ∪ {0 ≤ Re κ ≤ 15, |Im κ| > 14},
R = [0, 15] × [−14, 14]. By E(κ̄) = E(κ)̄ it suffices to treat Im κ ≥ 0.

With Theorem B′ (exactly two zeros in R: the gauge zero κ̄ and κ₁) this closes
the linear part of P4(b): **κ₁ ∈ [2.81055254439, 2.81055255326] is the unique
non-gauge mode of the EC solution with Re κ ≥ 0 in the stated class**, and
γ = 1/κ₁ ∈ [0.3558019218, 0.3558019231] unconditionally at the linear level
(still excluded from the claim: kink modes, non-spherical modes, nonlinear
codimension-one stability — FORMULATION §3). E is entire-in-κ on Re κ > −1.19
away from the pole κ = −0.0990 (§2.2 below), so "no zeros in Ω" is well-posed;
the eigenvalue ⟺ E = 0 equivalence of §3.8 extends to Ω verbatim (its three
ingredients are κ-independent or hold for all κ with Re σ(κ) < 0, which is
automatic on Re κ ≥ 0).

## 2. Large-|κ| structure of our system

### 2.1 κ-dependence of the coefficient matrices

The 4D linearised system (`linsys.py`, `linsonic4.py`; p = (A_p, N_p, W_p, V_p)):
    P(u) p′ = [DQ(u) − Ψ(u, u′) − κ P_s(u)] p,
P κ-free and block-diagonal — row 1: (S, 0, 0, 0); row 2: (0, 1, 0, 0); fluid
block P_fl = [[3(1+NV)S, 4W(N+V)], [(4V+N+3NV²)S, 4W(1+V²+2NV)]] on (W_p, V_p).
κ enters ONLY through −κP_s with P_s supported on the fluid rows:
    P_s,fl = [[3S, 4VW], [4VS, 4(1+V²)W]]      (S = 1 − V², from `systems.sonic_system`).
So the coefficient matrix is **affine in κ**:
    𝒜(x; κ) = P⁻¹(DQ − Ψ) − κ P⁻¹P_s =: B(x) − κ C(x),
    C(x) = 0₂ ⊕ P_fl⁻¹ P_s,fl   (rank 2: rows/cols A_p, N_p are κ-free).
The centre-scaled form (`linscaled.full_system`) is the Λ₄-conjugate; the
eigenvalues of C are unchanged. (The reduced 3D system has the spurious factor
(κ − A(x)) and κ² terms — as in Stages 1–3, all asymptotics are done in 4D.)
The constraint slaves A_p algebraically: A_p = A(C̃_N N_p + C̃_W W_p + C̃_V V_p)/((κ−A)S)
= O(1/κ) · (fluid data) as |κ| → ∞ — the κ-affine 4D form needs no such division.

### 2.2 Singular exponents at the two ends

**Sonic point x = 0** (regular singular; `linsonic4.py`): M_n(κ) = nD₄ + E₄(κ),
E₄ affine in κ of rank one, exponents {0, 0, 0, σ(κ)} with the single κ-affine
exponent (balls from the A1 background, `linsonic4.sigma`):
    σ(κ) = −σ₀ − (1 − σ₀) κ,   σ₀ = 0.090105970507920 ± 5e−16,
so σ(1) = −1 exactly, the pole of E sits at σ = 0 ⟺ κ = −σ₀/(1−σ₀) = −0.09903,
and the resonances σ(κ) = n ∈ ℕ are at κ = −0.09903 − 1.09903 n (S1's list).
On Re κ ≥ 0: Re σ(κ) ≤ −σ₀ < 0 uniformly — the non-analytic Frobenius solution
x^{σ(κ)} blows up (algebraically at least like |x|^{−σ₀} even on the imaginary
axis), and the admissible (analytic) family stays 3-dimensional for ALL κ ∈ Ω;
the normalisation A_p(0) = 1 loses nothing (σ(κ) ≠ 0 on Ω).

**Centre t = e^x → 0** (Fuchsian; `lincentre.py`): exponents {−3, −3, 0, 0} in t,
κ-INDEPENDENT; the κ-terms of the scaled system carry positive t-weights
(‖C̃(x)‖ ~ 1.40 e^x, §2.3), so they are a regular perturbation of size O(|κ| t).
The regular (exponent-0) family r₁, r₂ is analytic in κ for every κ (§2.4 of
`s2-theorem-b.md`, rank G₀ = 2 with a κ-free minor).

### 2.3 WKB/Levinson form as |κ| → ∞ in the closed right half-plane

Eigenvalues of C(x): {0, 0, μ₋(x), μ₊(x)}, the μ's = eigenvalues of
P_fl⁻¹P_s,fl = the two inverse characteristic slopes ds/dx of sound rays in
(x, s) — real on the real background. Float profile (`tc_sizing.py mu`):
    x·μ₋(x) → 1 − σ₀ = 0.909894 as x → 0⁻  (μ₋ → −∞; the sonic characteristic),
    μ₊(0) = +0.434254 (finite);  μ±(x) → ±1.4011 e^x as x → −∞.
Both ∫μ converge at the centre end;  **Φ₊ := ∫_{−∞}^0 μ₊ dx = 0.8617**
(pieces 0.8177 [−8, −0.1] + 0.0435 [−0.1, 0] + 0.00047 tail).

Diagonalise C: T(x)⁻¹CT(x) = M(x) = diag(0, 0, μ₋, μ₊) (T real-analytic,
κ-free, from the background tube). y = T⁻¹q̃ satisfies y′ = (−κM + B̃)y,
B̃ = T⁻¹BT − T⁻¹T′. On any [x₁, −δ] with δ > 0 fixed, Levinson's theorem
(Eastham Ch. 1; Coddington–Levinson Ch. 3 Thm 8.1) gives four fundamental
solutions
    f_i(x; κ) = (e_i + r_i(x; κ)) exp(∫^x (−κ μ_i + B̃_ii)),   ‖r_i‖ = O(1/|κ|)
(the O(1/κ) rate needs one integration by parts of the oscillatory couplings;
B̃ ∈ C^ω suffices). **Dichotomy**: the phase differences are κ∫(μ_i − μ_j) with
∫(μ_i − μ_j) REAL and of fixed sign, so Re[κ∫Δμ] = Re κ · ∫Δμ:
- Re κ > 0: strict exponential dichotomy, f₊ (the μ₊ mode) is uniformly
  DOMINANT toward the centre, f₋ recessive, the two κ-free modes in between.
- Re κ = 0: all phases unimodular — the **Stokes set of E in the κ-plane is
  exactly the imaginary axis** (and the real axis for the conjugate pairing);
  no other Stokes lines enter the closed right half-plane. Levinson still
  applies there (bounded, not decaying, dichotomy), so the asymptotic form
  holds UNIFORMLY on Re κ ≥ 0 — consistent with the observed single-term
  behaviour of E on the axis (§5: no dips, steady phase rotation).
Two κ-dependent layers bound the Levinson zone:
- **sonic layer** |x| ≲ 1/|κ|: ξ = κx is O(1); μ₋ ~ (1−σ₀)/x merges with the
  Frobenius exponent σ(κ) (x^{σ(κ)} = e^{−κ(1−σ₀)ln|x|}·|x|^{−σ₀}); the model
  problem is a 2×2 Fuchsian-with-large-exponent layer (Bessel/Whittaker of
  large order ~ (1−σ₀)κ and argument ~ κx: Olver's uniform large-order theory);
- **centre layer** |κ|e^x ≲ 1, i.e. x ≲ x_t(κ) = −ln|κ|: the κ-terms fall
  below O(1) and the system relaxes to the κ-free Fuchsian one; matching to
  r₁, r₂ is a regular perturbation of size O(|κ| e^{x_d}) at x_d < x_t.

Predicted form (to be pinned in S4-1, all factors from the three regions):
    E(κ) = c(κ) e^{Φ₊ κ} (1 + O(1/κ)),   c(κ) ~ C κ^{−p} slowly varying,
with the measured p ≈ 1.00 on the imaginary axis (τ|E(iτ)| = 0.63–0.64 over
τ ∈ [14, 60]) and effective p ≈ 0.7 along the real ray with Φ₊ = 0.8617 fixed
(a slowly-varying factor, e.g. κ^{−1+O(σ₀)} or a Γ-ratio from the sonic layer,
is not yet resolved by the fits); d(arg E)/dτ on the axis → 0.859 ≈ Φ₊. The
1/κ law is exactly what the constraint slaving A_p = O(1/κ)·(fluid) suggests
for E ∝ A_p-amplitude. Nonvanishing of the limit c(κ)κ^p is the crux of
route (a).

## 3. Literature digest (exact citations; web-checked 2026-08-31)

**Koike–Hara–Adachi.** PRL 74 (1995) 5170 (gr-qc/9503007); PRD 59 (1999)
104008 (preprint gr-qc/9607010). The 1999 paper, §V.G.3: "a thorough search
… in the region 0 ≤ Re κ ≤ 15, |Im κ| ≤ 14" found only the relevant mode and
the gauge mode; they state explicitly that "it is theoretically impossible to
cover the whole values of κ ∈ ℂ" by search and that they "have not found a
beautiful argument which can restrict possible eigenvalues". The gap is
covered only heuristically by their Lyapunov/PDE-evolution analysis (§VI,
App. D: Gram–Schmidt extraction of the top exponents, valid "even if the
eigenmodes do not form a complete set"). → Our box R is literally theirs;
KHA provide corroboration, no exclusion technique.

**Gundlach.** PRD 55 (1997) 695 (gr-qc/9604019): matching determinant
det A(λ), holomorphy checked numerically, zero counting by ∮ A′/A dλ
(App. E); no large-|λ| exclusion, and §V.C wishes for one. Living Rev. Rel.
10 (2007) 5 (arXiv:0711.4620): mode searches reviewed as purely numerical;
completeness flagged open (§3.3). → Argument-principle-in-a-box transfers
(= our Theorem B); nothing for Ω.

**Guo–Hadžić–Jang(–Schrecker).** LP existence + monotonicity ρ̂′ < 0, ω̂′ > 0:
CMP 386 (2021) 1551 (arXiv:2011.01013, Remark 4.16). Relativistic LP: Ann.
PDE 9 (2023) art. 4 (arXiv:2112.10826). **Mode stability of LP: GHJS,
arXiv:2509.12435, Thm 3.4** — exactly one eigenvalue with Re λ ≥ 0 (the
trivial λ = 1). Anatomy: Lemmas 3.8–3.9 (energy identity ⇒ no eigenvalues
Re λ > 1); Prop. 3.10 + 3.14 (high frequency: |Im λ| ≥ 8, Re λ ∈ [0, 1]
excluded by commuting the 2nd-order eigenvalue ODE m times with
∂_y(∂_y + 2/y), shifting the first-order coefficient by 2m w′/w < 0, w ≃
distance to the sonic point; real+imag energy identities combine to a
coercive inequality, constants b₀ = 1/5, b₁ = 8 certified by interval
arithmetic); Prop. 3.15 (|Im λ| ≤ 1/5 via the quotient by the trivial
eigenfunction); Prop. 3.16 (intermediate box by VNODE-LP validated Frobenius
matching). They note the commutator damping "appears to have a universal
character" for fluid-implosion problems. → Closest structural match; our
P4-3 diagnostics (v_rel ↑, ρ̂ ↓, w = 1/3 − v_rel² ↓ on the cone) were run
precisely to enable this transfer. Their profile is Newtonian-isothermal 2D;
ours is relativistic 4D with a constraint — the reduction to a single
2nd-order equation with the right sonic weight is the transfer work.

**Computer-assisted exclusion outside a compact set.**
- Barker–Zumbrun, Math. Models Methods Appl. Sci. 26 (2016) 2451
  (arXiv:1601.00837): interval-arithmetic Evans function winding number on a
  compact contour + analytic high-frequency energy estimates outside — the
  exact procedural precedent for Theorem B + C. Also Barker, JDE 257 (2014)
  2950. Humpherys–Lyng–Zumbrun, ARMA 194 (2009) 1029: dedicated
  "high-frequency bounds" sections confining Evans zeros to an explicit
  compact set; Humpherys–Sandstede–Zumbrun, Numer. Math. 103 (2006) 631
  (analytic bases / well-defined E at large |λ|).
- Chen–Hou, arXiv:2210.07191 + 2305.05660: no spectral localisation at all —
  frequency-uniform weighted energy estimates with computer-verified
  constants. Alternative philosophy; does not yield a box statement.
- Buckmaster–Cao-Labora–Gómez-Serrano, Forum Math. Pi 13 (2025)
  (arXiv:2208.09445) §7: L = A₀ − δ_g + K (maximally dissipative + compact),
  so σ(L) ∩ {Re λ > −δ_g/2} is finite — compact confinement "for free", no
  explicit box; non-radial: Cao-Labora–Gómez-Serrano–Shi–Staffilani,
  arXiv:2310.05325. → the soft half only.
- Costin–Donninger–Glogić(–Huang/Xia), CMP 343 (2016) 299, CMP 351 (2017)
  959, Nonlinearity 29 (2016): global right-half-plane mode stability via
  quasi-solutions / hypergeometric connection coefficients at a singular
  point with λ-affine exponent — needs an explicit/rational profile (GHJS
  say so explicitly), which we lack.
- Warnick, CMP 333 (2015) 959: QNMs at a horizon with exponent affine in the
  spectral parameter as honest H^k eigenvalues — the function-space framework
  for κ-affine-exponent singular points; no large-frequency exclusion there.
- Rigorous asymptotics toolbox: Olver, *Asymptotics and Special Functions*
  (1974), Ch. 6 (Liouville–Green with computable error bounds) and Chs.
  10–12 (uniform large-order Bessel); Eastham, *Asymptotic Solution of
  Linear Differential Systems* (1989), Ch. 1 (Levinson); Coddington–Levinson
  (1955) Ch. 3 Thm 8.1.

What does NOT transfer: CDG-type explicit-profile hypergeometrics (no closed
form for EC); BCG sectoriality alone (gives finiteness, not a box); KHA/
Gundlach numerics (not rigorous). What transfers: GHJS's three-regime plan
and the Barker–Zumbrun pipeline shape (energy at infinity + validated
winding on the compact leftover), with our sonic point playing the role of
their sonic point/shock end, and Re σ(κ) ≤ −σ₀ < 0 supplying the boundary
terms' sign at x = 0.

## 4. Candidate routes

### 4.1 Route (a): Levinson/integral-equation asymptotics with explicit K₀

**Claim shape.** For |κ| ≥ K₀, Re κ ≥ 0: |E(κ)| ≥ (c₀/|κ|) e^{Φ₊ Re κ} > 0,
with explicit c₀, K₀. Then Ω is covered by {|κ| ≥ K₀} plus the compact
leftover R′ \ R, R′ = [0, K₀] × [−K₀, K₀] (route (c) machinery, §4.3).

**Formulas for our system.** Split (−∞, −δ] at x_t(κ) = −ln|κ| + ln ε₀:
1. Levinson zone [x_t, −δ]: y = T⁻¹q̃ as in §2.3; the four solutions f_i with
   remainders bounded by the affine-contraction estimate
   ‖r_i‖ ≤ ρ(κ) = (C₁/|κ|) e^{C₂},  C₁, C₂ = weighted integrals of |B̃_ij|,
   |B̃′_ij|, |μ_i − μ_j|⁻¹ over the tube (all data already enclosed per step in
   `Tube.system_data`; the 1/|κ| comes from one integration by parts, kernel
   K(x,y) = exp(−κ∫_y^x Δμ) with ∂_y-antiderivative 1/(κΔμ)).
2. Sonic connection (Lemma L1, the research item): the analytic solution p̃
   normalised by A_p(0) = 1 satisfies
   p̃(−δ; κ) = ĉ_s(κ) f₊(−δ) + d₁(κ)f₁ + d₂(κ)f₂ + O(f₋-recessive),
   ĉ_s(κ) = κ^{−q} c_s (1 + O(1/κ)), c_s ≠ 0 explicit — from the ξ = κx layer
   (2×2 fluid block ⇒ scalar 2nd-order model with a regular singular point of
   exponents {0, σ(κ)}; uniform large-order Bessel/Whittaker asymptotics,
   Olver Chs. 10–12, give the layer's connection coefficients as Γ-ratios).
3. Centre connection (Lemma L3): at x_d(κ) = x_t(κ): r₁, r₂ from the §2.4
   Fuchsian family, certified for a κ-box covering |κ| ≤ K₀ and evaluated at
   t_d = ε₀/|κ| (the series improves as t_d shrinks); E's determinant is then
   E(κ) = ĉ_s(κ) e^{−κ∫_{−δ}^{x_d}μ₊ + ∫B̃₄₄} D(κ) (1 + O(ρ + |κ|t_d)),
   D(κ) = e^{3x_d} det[r₁, r₂, (T e₊)](x_d) → D_∞ ≠ 0 (κ-free Fuchsian data).
**Finite computation.** Interval enclosures of μ±, T, T⁻¹, B̃ and their
integrals over the 244-step tube + centre extension (reuse `lintube`); the
Γ-ratio lower bound (Arb's rigorous acb Γ); D_∞ by the existing centre
machinery; assemble c₀ and the threshold K₀ = min{|κ|: total relative error
< 1/2}. **Analytic lemmas.** L2 (Levinson with explicit constants, marginal
dichotomy on the axis included) — standard but laborious; L3 (κt-perturbation
of a Fuchsian system) — fits the repo's Banach-contraction style; L1 (sonic
layer, uniform in arg κ ∈ [−π/2, π/2]) — genuine research, though 2×2 and
with a century of special-function technology behind it. **Difficulty:**
high (L1), medium (L2), low (L3). **Certification:** every constant is a
ball; no new integrator needed. Estimated K₀ ~ 30–100 (from §5 the
asymptotic regime is visibly established by |κ| ≈ 20).

### 4.2 Route (b): GHJS-style energy/commutator estimates (physical variables)

**Claim shape.** (b1) No eigenvalues with Re κ > c₀ (direct energy identity);
(b2) no eigenvalues with Re κ ∈ [0, c₀], |Im κ| ≥ b₁ (commutator damping).
If c₀ ≤ 15 and b₁ ≤ 14, Theorem B already covers the rest and NO new contour
run is needed; else top up by route (c) on [0, max(c₀,15)] × [14, b₁].

**Formulas for our system.** Reduce the 4D system to one 2nd-order scalar
equation: rows A, N slave A_p (constraint, §2.1) and N_p (quadrature of
N_p′ = N_p F_N + (κ-free source)); the fluid pair gives, for a suitable
combination χ (S1's diagnostics suggest the v_rel-perturbation),
    (𝔴(x) χ′)′ + a₁(x; κ) χ′ + a₀(x; κ) χ = 0,
with 𝔴 vanishing simply at x = 0 (𝔴 ∝ Δ̃/(product of positive factors), the
sonic weight; 𝔴 ∝ det P_fl, D₁ = 7.5606 ≠ 0 from §3.8 of `s2-theorem-b.md`)
and a₀ quadratic in κ. GHJS Prop. 3.10's mechanism: commuting m times with
the analogue of ∂_y(∂_y + 2/y) shifts the first-order coefficient by
2m 𝔴′/𝔴 < 0 near the sonic point; the P4-3 monotonicity (v_rel ↑ to 1/√3,
ρ̂ ↓, w = 1/3 − v_rel² ↓) supplies the interior sign conditions. Real and
imaginary parts of ∫ χ̄·(eq) with weights give
    |Im κ|² ∫ χH|χ|² + ∫ χ|D χ|² ≤ ∫ χ H̃ |χ|²,  H ≤ −c < 0 for |Im κ| ≥ b₁.
**Finite computation.** Pointwise interval sign-checks of the coefficient
combinations (H, H̃, 𝔴′/𝔴, weights) along the certified tube + endpoint
asymptotic sign lemmas — cheap (one pass over 344 steps); certify P4-3
(currently CONJECTURED) as a tube sign certificate on the way. **Analytic
lemmas.** The 2nd-order reduction (exact, `fmpq_mpoly`, like
`linear_constraint_propagation`); the function-space setting: exclusion must
hold in a class CONTAINING ours — energy methods exclude H^m-eigenfunctions,
and our class (analytic at sonic, regular centre) embeds once Re σ(κ) < 0
(§2.2) makes the boundary terms at x = 0 vanish; centre boundary terms
vanish by the {−3,−3,0,0} exponents. **Difficulty:** medium — GHJS closed
the exact analogue (their Re λ ∈ [0,1], |Im λ| ≥ 8, b's certified by IA);
the risks are (i) the relativistic 4D→2nd-order reduction may not produce a
single sign-definite weight, (ii) the gauge freedom: our κ̄ sits at
Re κ < 1, harmless if c₀, b₁ leave it inside R. **Certification:** trivially
Arb-friendly (interval inequalities); no oscillatory integration at all —
the κ-independence of the a-priori bound is the whole point.

### 4.3 Route (c): growing contours with uniform κ-Taylor-model control

**What it is.** The Theorem-B machinery scales to bigger rectangles R′:
tile half-width can stay 0.25 (E varies on the κ-scale ≈ 2.5, oscillation
scale 2π/Φ₊ ≈ 7.3, both κ-independent); the only growing cost is the sonic
series order K ≥ max Re |σ(κ)| ≈ 0.91|κ| (+ margin) and the per-tile Taylor
degree m. Estimate for R′ = [0, 60] × [−60, 60]: contour length 360 (vs 86),
~1450 tiles at ~40–80 s ⇒ ~4–8 h on 8 workers; min |E| on the new contour
≈ 0.64/60 ≈ 1.1e−2 (§5) — larger than the 4.8e−3 already handled on ∂R.
**What it cannot do.** Close Ω alone: a finite computation certifies finitely
many tiles; with no decay-at-∞ input there is always leftover contour. Its
role is exactly to bridge [14, K₀] (route a) or [14, b₁] (route b), and it is
the fallback if the analytic thresholds come out larger than expected.
**Certification.** Existing `modecount`/`analyticity` code, parameters only
(K, box widths near the axis where |E| ~ 0.6/|κ|); winding 2 on ∂R′ +
Theorem B′'s two zeros ⇒ no zeros in R′ \ R.

## 5. Numerical sizing

Float, `tc_sizing.py` (DOP853 rtol 1e−12, x_end = −8, E = E_fin =
A_p(x_end)e^{x_end} — same zeros as E₃ by §3.3 of `s2-theorem-b.md`).

Rays (46–61 samples each; "dips" = interior local minima of |E|):

| path | |E| range | min |E| (where) | dips |
|---|---|---|---|
| A: κ = 15 + t, t ∈ [0, 45] | 9.9e3 → 2.6e20 | 9.9e3 (t = 0) | none |
| B: κ = s + 14i, s ∈ [0, 15] | 4.5e−2 → 9.4e3 | 4.5e−2 (s = 0) | none |
| C: κ = iτ, τ ∈ [14, 60] | 4.5e−2 → 1.1e−2 | 1.1e−2 (τ = 60) | none |
| strip 0.5 + iτ, τ ∈ [14, 60] | 6.7e−2 → 1.6e−2 | at τ = 60 | none |
| strip 1.3 + iτ (the ∂R minimum's abscissa) | 1.3e−1 → 3.3e−2 | at τ = 60 | none |
| strip 2.0 + iτ | 2.3e−1 → 5.9e−2 | at τ = 60 | none |

Arcs |κ| = r, first quadrant (θ ∈ [0°, 90°], 61 samples; conjugation covers
the fourth quadrant):

| r | |E| at θ=0 | |E| at θ=90° | r·|E|(ir) | dips |
|---|---|---|---|---|
| 20 | 6.4e5 | 3.2e−2 | 0.633 | none |
| 30 | 2.7e9 | 2.1e−2 | 0.636 | none |
| 50 | 5.5e16 | 1.3e−2 | 0.639 | none |

Readings. (i) |E| is monotone along every ray/arc probed: **no near-zeros
anywhere in the first quadrant out to |κ| = 60**, and none appear on the way.
(ii) On the imaginary axis τ|E(iτ)| = 0.627 → 0.639 (τ = 14 → 60): a clean
1/|κ| law with unit-modulus constant ≈ 0.64; the phase rotates at
d(arg E)/dτ = 0.822 → 0.859, converging to Φ₊ = 0.8617 like O(1/τ) — E
behaves as a SINGLE WKB term, no two-term interference. (iii) Off the axis
|E| grows like e^{Φ₊ Re κ} (ray A: ln-slope 0.84 over t ∈ [0, 45]).
(iv) The asymptotic regime is established by |κ| ≈ 20 already (the 1/τ law
holds at τ = 14 to 2 %): K₀ in the 30–60 range looks realistic, and the
compact leftover for route (a) is modest. (v) The weakest |E| on any ∂R′ is
the axis value ≈ 0.64/|κ|, which stays above Theorem B's certified contour
minimum 4.8e−3 for |κ| ≲ 130 — route-(c) tiles remain tractable there.

## 6. Recommendation and staged plan

**Recommendation.** Primary route (b) — the GHJS commutator/energy template —
with route (c) as the bridging computation, and route (a) kept warm as the
fallback (its cheap S4-1 float validation doubles as a diagnostic for both).
Reasons: (b) is the only route whose analytic core has already been carried
out end-to-end on a directly analogous sonic-point problem (GHJS 2509.12435,
Thm 3.4), our P4-3 diagnostics verify its structural hypotheses on the EC
profile, its certification layer is pointwise interval inequalities (no
oscillatory quadrature, κ-independent), and it degrades gracefully: whatever
(c₀, b₁) come out, route (c) bridges the gap. Route (a)'s L1 (sonic layer,
uniform in arg κ) is genuinely novel analysis; keep it as the fallback since
§5 shows the asymptotics are real and one-termed. Decide (b) vs (a) at the
S4-2 gate below.

**Stages** (sized like S1/S2; each ends with a note + ledger update):
- **S4-1 (reduction + asymptotics, ~1 session).** Derive the exact 2nd-order
  reduction of the 4D linearised system (fmpq_mpoly identity, the sonic
  weight 𝔴, the slaved rows) and float-check that its regular eigenvalues
  reproduce κ₁, κ̄ to 1e−9. Derive the route-(a) leading term (c_s, q, D_∞)
  and float-validate E(κ) = c(κ)e^{Φ₊κ}(1+O(1/κ)) on |κ| ∈ [30, 300]
  (explains the p ≈ 1 vs 0.7 puzzle of §2.3). Output:
  `notes/s4-reduction.md`, `tc_` scripts. Ledger: new row P4-7 (structure of
  E at large |κ|: 1/|κ| law, Φ₊, single-term behaviour) at CONJECTURED;
  P4-3 unchanged.
- **S4-2 (energy identities, ~1–2 sessions).** (b1): Re κ > c₀ exclusion
  with c₀ certified (tube interval pass + endpoint lemmas). GATE: if the
  reduction of S4-1 yields no sign-definite weight (risk (i) of §4.2),
  switch primary to route (a) and re-scope. Ledger: P4-8 "no eigenvalues
  with Re κ > c₀" at CERTIFIED; P4-3 promoted to CERTIFIED (tube sign
  certificate).
- **S4-3 (high frequency, ~2–3 sessions, research core).** (b2): commutator
  damping for Re κ ∈ [0, c₀], |Im κ| ≥ b₁, constants certified in Arb.
  Fallback: route (a) L1–L3 with K₀ explicit. Ledger: P4-9 (the
  high-frequency exclusion) at CERTIFIED with the (c₀, b₁) or K₀ actually
  achieved.
- **S4-4 (bridge, ~1 session + hours of compute).** If (b₁, c₀) ⊄ R (or
  route (a)'s K₀ > 15): extend the winding-number run to R′ ⊇ R covering the
  leftover (`modecount` with sonic series order K ≈ 0.91·|κ|_max + 10,
  near-axis tile widths per §5(v)); winding 2 on ∂R′ ⇒ no zeros in R′ \ R. Ledger: P4-6 evidence
  extended (winding on ∂R′), stored under `results/theorem_c/`.
- **S4-5 (assembly + review, ~1 session).** Theorem C = S4-2 + S4-3 + S4-4 +
  Theorem B′; adversarial review of the chain (as for A and B); FORMULATION
  v1.5 records Theorem C and the unconditional linear statement (§1).
  Ledger: P4-10 (Theorem C, CERTIFIED after review) superseding the
  "[open — research item]" of FORMULATION §2; P4-2's remaining scope note
  updated.

Honest failure modes: the 2nd-order reduction has no good weight AND L1
stalls → Theorem C stays a conjecture, but S4-2/S4-4 still yield a much
larger certified zero-free region (Re κ > c₀ plus R′), which is worth
recording either way.
