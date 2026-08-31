# P4 S4-1 — Theorem C: exact 2nd-order reduction, large-|κ| law of E, sign survey

Status: S4-1 done (2026-08-31). Code (new; nothing pre-existing modified):
`src/p4/validated/{tc_reduce,tc_reduced_eigs,tc_leading}.py`; tests
`uv run pytest problems/P4/tests/test_tc_reduction.py -q` (9 tests, ≈ 21 s).
The reduction identities of §1 are **exact** (`fmpq_mpoly`); §§2–4 are float
(nothing here is a certificate). Inputs: theorem-c-plan.md §§2, 4.2, 4.1;
s2-theorem-b.md §§2–3; the P4-3 diagnostics (numerics-report.md §D).

## 1. The exact 2nd-order reduction (route (b) prerequisite)

### 1.1 Why a quotient is forced, and the gauge identity (exact)

At the sonic point the analytic family of the 4D linearised system is 3-dim
(exponents {0,0,0,σ(κ)}), and the constraint functional c(p) of an analytic
solution vanishes identically once c(0) = 0 (s2 §3.8), so analytic ∩ Σ is
2-dimensional — but Σ-solutions on an interval form a 3-dim space, so no scalar
2nd-order ODE can hold for all of them in the fluid variables alone: a
1-dimensional invariant direction must be quotiented out. That direction exists
for every κ: **the gauge vector g = (A′, N′ + κN, W′, V′) solves the 4D
linearised system identically** — `tc_reduce.gauge_residual` builds
(SΔ̃)³[P g′ − (DQ − Ψ − κP_s)g] with u′, g′ by Cramer and gets the *zero
polynomial* in (A, N, W, V, κ), before any constraint elimination; c(g)
eliminates to 0 under A = Anum/S. This closes the "ball identity to order 38"
caveat of Theorem B (s2 §3.7/3.8 rigor caveat (1); ledger P4-6 updated) and is
the engine of the reduction.

### 1.2 The pair (χ, η) and the closure identities (exact)

N-clock gauge-invariant fluid perturbations (N′ + κN = N(F_N + κ) never
vanishes off the window in §1.4; W′, V′, A′ each vanish at an interior point,
so N is the only usable clock):
    χ := Δ̃[(N′ + κN) V_p − V′ N_p] = Δ̃ (N′ + κN) ΔV,   ΔV = V_p − V′ ξ,
    η := Δ̃[(N′ + κN) W_p − W′ N_p] = Δ̃ (N′ + κN) ΔW,   ξ := N_p/(N′ + κN),
polynomial covectors ℓ = (0, −PV, 0, e), m = (0, −PW, e, 0) with e = Δ̃(Q₁+κN),
PV = Δ̃V′, PW = Δ̃W′; both annihilate g exactly. Δv_rel = ∂_V v_rel · ΔV (the
N-clock kills the N-part), so χ IS the v_rel/V-perturbation up to the positive
factor (N²−1)/(N+V)² — the plan's first candidate works. The closure
(`tc_reduce.reduce_pair`, exact after elimination; the Ap/Wp/Vp components
match by construction, the Np component is the identity):
    S²Δ̃² dn e · χ′ = Achi χ + Bchi η + Cchi c(p),   dn := (κ − A)S,
    S²Δ̃² dn e · η′ = Aeta χ + Beta η + Ceta c(p),
so on Σ (c = 0) the pair closes into a 2×2 first-order system with
coefficients a = Achi/D2 etc., D2 = S²Δ̃²·dn·e. Eliminating η gives the scalar
2nd-order equation for χ (§1.3). Flint factorisation: D2 = −(64/3)·N W³ D³ S⁶ ·
(A−κ)·3(F_N+κ); the numerators carry κ-degree 3. Kernel of (ℓ, m) on Σ =
span(g): the pair is the flow on Σ/span(g); at κ = κ̄ the admissible solution
IS g and χ ≡ 0.

### 1.3 Scalar form, the sonic weight, and the κ-structure

χ″ − Tχ′ − Uχ = 0 with T = a + d + b′/b, U = a′ + bc − a(T − a) (b = Bchi/D2
etc. — exact rational functions). Partial fractions in κ (float, residuals
≤ 3e−15, `tc_reduced_eigs`/prototype):
    a + d = −κ(μ₊+μ₋) + q₀(x) + d/dx log[(A−κ)(F_N+κ)],
with the pole residues *exactly* A′ and F_N′ — so in T the (A−κ), (F_N+κ)
log-derivatives **cancel** against b′/b (b ∝ N·fb/((A−κ)(F_N+κ)WD)):
    T = −κ(μ₊+μ₋)(x) + q₀(x) + d/dx log[N·fb/(W·D)],
fb = the κ-cubic irreducible factor of Bchi (39 terms). q₀(x) is not a finite
combination of log-derivatives of the obvious factors (fit residual 4e−3 with
5 factors), so exp(−∫T) is left transcendental and we take the plan's form
with the **algebraic weight**
    (𝔴 χ′)′ + a₁ χ′ + a₀ χ = 0,   𝔴 := −D = −Δ̃/(4SW),
    a₁ = −(𝔴T + 𝔴′) = 𝔴·κ(μ₊+μ₋) − 𝔴[q₀ + (log(N fb/W))′],   a₀ = −𝔴U.
Exact weight identities (tests): Δ̃ = 4SWD and
    **𝔴 = −D = 3(N+V)² (1/3 − v_rel²) = 3(N+V)² · w_GHJS** —
the sonic weight is GHJS's w times the positive factor 3(N+V)². It vanishes
simply at x = 0 with slope −D₁ = −7.5606 (𝔴(−0.05) = 0.3765 ✓) and is positive
on the whole cone. κ-structure of a₀: U = −det C·κ² + u₁κ + u₀ + rational
corrections with the poles (A−κ)⁻¹, (F_N+κ)⁻¹ and the three κ-roots of fb, so
    a₀ = −(3−V²) κ² + a₀₁(x) κ + a₀₀(x) + O(1/κ)-rational,
(𝔴·det C = D det C·(−1)·(−1) = −(3−V²) exactly, since det C = (3−V²)/D): the
plan's "a₀ quadratic in κ" holds for the polynomial part; the extra rational
terms are the intrinsic price of the constraint slaving A_p ∝ 1/(κ−A) and are
uniformly O(1/κ) for |Im κ| ≥ 1 or Re κ > 2. a₁ is κ-affine plus the same
type of O(1) rational correction (fb′/fb).

### 1.4 Apparent singularities (inventory)

(i) x = 0: the genuine regular singular point, 𝔴 ∝ D → 0; scalar exponents
{1, 1+σ(κ)} for χ, i.e. {0, σ(κ)} for ΔV = χ/(Δ̃(N′+κN)) (Fuchs sum checked at
κ = 3: −0.82 = 1 + τ₀ ⟹ σ = −2.820 ✓ = linsonic4's σ(3)). (ii) κ = A(x):
dn = 0; ℓ, m, c become dependent on Σ there (Achi(κ=A) does NOT eliminate to
0), the same apparent singularity as the 3D system's (κ−A) factor (s2 §2.5(c));
real κ ∈ (1, A_max = 1.89] only. (iii) κ = −F_N(x): the N-clock chart
degenerates (ℓ ∥ m); real κ ∈ [κ̄, 1) — note −F_N has an interior dip (min
0.3549 at x ≈ −0.03), so at κ̄ itself the chart is singular at x* ≈ −0.06.
(iv) fb(x; κ) = 0: movable apparent singularities of the *scalar* form (b = 0;
absent from the 2×2 pair); no zeros found on the imaginary axis
(min|b|/τ ≈ 0.50 at τ = 15, 30, 60). All of (ii)–(iii) sit on the real-κ
segment [κ̄, 1.89] ⊂ R where Theorem B already governs; none affect the
exclusion regions of route (b).

## 2. Float validation of the reduction (`tc_reduced_eigs.py`)

Reduced 2×2 integration (coefficients evaluated from the exact polynomials,
DOP853 rtol 1e−12, co-integrated background) vs the full linearised system
(`perturb`), (χ, η)(−8):

| κ | rel. diff χ | rel. diff η |
|---|---|---|
| 2.5 | 2.9e−12 | 1.6e−07 |
| κ₁ = 2.8105525488 | 2.6e−05 (χ ≈ 1.2e2 vs scale 1e10: abs ≈ 3e−13 of scale) | 9.7e−11 |
| 0.34 (complex-x detour) | 1.4e−09 | 2.0e−04 |
| 5 + 7i | 2.7e−12 | 1.4e−08 |

For real κ ∈ (0.35, 1.05) the chart pole (iii) is crossed by a complex-x
detour (sin-bump of height 0.04 on [−0.25, −0.05]; background and solution are
analytic in x, path-independence checked to ≤ 3e−9). Eigenvalues of the
reduced problem, E_red(κ) := χ(x_end)e^{3x_end}, secant zeros (x_end = −11):
    κ₁_red = 2.8105525488268,  |κ₁_red − κ₁_certified| = 3.4e−13  (≤ 1e−9 ✓)
    κ̄_red  = 0.35569920371096, |κ̄_red − κ̄_certified|  = 7e−15    (≤ 1e−9 ✓)
(x_end matters: at −8 the regular part contaminates E_red at relative e^{3x_d}
and shifts the zero by 2.9e−9; the shift scales as e^{3x_end} as predicted.)
At κ̄ the zero is the quotient structure itself: χ(p̃(κ̄)) ≡ 0 since p̃(κ̄) ∝ g —
the reduced problem sees the gauge zero through the kernel, which is why its
reproduction is essentially exact.

## 3. Sign-condition survey for route (b) (float pass, certified tube x-grid)

Grid: the 344 step abscissae of the certified tube ([−8, −0.05]; the tube JSON
is decoded to floats). `tc_reduced_eigs survey`:

| quantity | result on the grid | GHJS analogue / verdict |
|---|---|---|
| 𝔴 = −D | > 0 everywhere, 0.3765 … 1.2e7 | weight positive ✓ |
| 𝔴′/𝔴 = D′/D | **< 0 everywhere**, range [−19.9, −2] | their w′/w < 0 holds *globally*, not just near the sonic point ✓ |
| 𝔴 vs 3(N+V)²w_GHJS | identical (3.7e−9, float) | P4-3 monotonicity transfers verbatim ✓ |
| det C = (3−V²)/D | < 0 everywhere | a₀'s κ²-density −(3−V²) sign-definite ✓ |
| tr C = μ₊+μ₋ | < 0 everywhere (−(μ₊+μ₋) = 1.88 at −0.4 → 0.0009 at −4) | transport term one-signed ✓ |
| Hq := (3−V²) + D(tr C)²/4 | > 0 for x ≤ x_H = −0.354; **< 0 on (−0.354, 0)** (single sign change; interior min +0.034 at −0.364) | naive square completion fails on the sonic shoulder where tr C ~ μ₋ ~ (1−σ₀)/x — exactly the region GHJS cure with the m·w′/w commutator shift |
| U/κ² → −det C at κ = iτ | rel. deviation ≤ 20/(τ·distance-to-0); ≤ 0.11 for x ≤ −0.5 at τ = 60 | the κ²-part dominates uniformly away from the sonic layer ✓ |
| b-zeros on the axis | none: min|b|/τ = 0.50 (τ = 15, 30, 60) | scalar form regular on the axis ✓ |

## 4. Route (a): the large-|κ| law of E (`tc_leading.py`, ledger P4-7)

E on three rays θ = 0, π/4, π/2, |κ| ∈ [30, 300] (13 log-spaced samples each;
sonic start δ = min(0.04, 8/|κ|), K = 44; results/theorem_c/tc_leading.json).
Joint fit of ln|E| = ln C + Φ₊ Re κ − q ln|κ| + Re(d/κ):
    **Φ₊ = 0.861612, q = 0.9971, C = 0.6326, max residual 4.9e−3**,
d_real = −8.50, d_diag = −4.28(1+i), d_imag ≈ 0 — all consistent with a single
real d ≈ −8.5 (the fit only sees Re(d/κ)). C_est = E κ^q e^{−Φ₊κ} at |κ| = 300:
0.614, 0.620+0.013i, 0.631+0.007i on the three rays (drift ≤ 4e−2, the size of
d/κ): one WKB term, one complex constant. **Resolved law and the p ≈ 1 vs 0.7
puzzle** (plan §2.3): q = 1 within 3e−3; dropping the d/κ term drags the fitted
q to 0.911 on [30, 300] and to ≈ 0.7 on the plan's shorter [15, 60] real-ray
window — the "0.7" was the d/κ correction aliasing into the power. So
    E(κ) = C κ^{−1} e^{Φ₊κ} (1 + d/κ + O(1/κ²)),  C ≈ 0.633,  d ≈ −8.5,
uniformly on the three rays — the 1/κ is the constraint-slaving power, as the
plan conjectured. Factorisation (Levinson frame: T-columns = unit scaled fluid
eigenvectors with positive ṽ-component, exponent integrals based at 0):
Φ₊ = 0.861662 by quadrature (matches the fit to 5e−5); centre determinant
D_∞ = lim_{t→0} det[r₁, r₂, e₊]_{(Â,ñ,ṽ)} = **2.038906**, κ-free (= (2/3)·ñ₂(0)
= (2/3)·3.058359; r₁'s ñ-component is κ-dependent but drops out of the det);
zone factor J = ∫_{−∞}^{−0.05}(B̃₄₄+3)dx = 6.2736 (B̃₄₄ → −3.00003 at −16 ✓;
inside |x| < 0.05, B̃₄₄ has a ≈ +0.0285/x tail whose logarithm belongs to the
sonic layer and caps the layer/zone split ambiguity at κ^{±0.03});
    c_s := C e^{J}/D_∞ = 164.6   (in the stated conventions, layer cut 0.05).
The route-(a) shape E = c_s D_∞ e^{−J} κ^{−1} e^{Φ₊κ}(1+O(1/κ)) is thereby
float-validated end to end; what remains analytic is Lemma L1 (the layer
constant from uniform large-order asymptotics) and explicit Levinson constants.

## 5. GATE VERDICT (for S4-2)

**Route (b) looks viable — keep it primary.** The plan's risk (i) (no
sign-definite weight from the 4D reduction) is *retired*: the reduction exists
as an exact identity, and the weight is 𝔴 = −D = 3(N+V)²·w_GHJS — positive on
the cone, simple zero at the sonic point with the D₁ = 7.5606 factor, and with
the commutator sign 𝔴′/𝔴 < 0 holding on the whole certified grid (GHJS only
need it near the sonic point). P4-3's certified-to-be v_rel/ρ̂ monotonicity is
literally the monotonicity of this weight. The κ²-density of a₀ is the
sign-definite −(3−V²). Known costs, none fatal: (a) the Hq < 0 sonic shoulder
(−0.354, 0) means S4-2/S4-3 must use the GHJS commutator mechanism, not a
pointwise inequality — expected; (b) coefficients are rational, not polynomial,
in κ (poles confined to real κ ∈ [κ̄, 1.89] ⊂ R, O(1/κ) uniformly in the
exclusion regions); (c) the scalar form's fb-factor — prefer the 2×2 pair
(χ, η) for the energy identities, which has no fb. Route (a) stays the
fallback, now with the law pinned (q = 1, Φ₊ = 0.8616, C = 0.633, D_∞ = 2.039)
and the compact leftover for route (c) sized by K₀ ≈ 30 being plausible.

## 6. Reproduction

```
PYTHONPATH=problems/P4/src uv run python -m p4.validated.tc_reduce            # exact identities, ~2 s
PYTHONPATH=problems/P4/src uv run python -m p4.validated.tc_reduced_eigs all  # §2 + §3, ~3 min (P4_TUBE_CACHE=<tube.json> for the tube grid)
PYTHONPATH=problems/P4/src uv run python -m p4.validated.tc_leading scan      # ~30 min; then `fit`, `factor`
uv run pytest problems/P4/tests/test_tc_reduction.py -q                       # 9 tests, ~21 s
```
