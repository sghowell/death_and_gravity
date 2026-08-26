# P4 S1 numerics report — Evans–Coleman CSS solution at k = 1/3

Status: S1 (CONJECTURED-level reproduction), 2026-08-25. Code: `problems/P4/src/p4/`
(numpy/scipy only at runtime; sympy was used once, ephemerally via `uv run --with sympy`,
for the derivations quoted below). Tests: `uv run pytest problems/P4/tests -q`.
Plots: `problems/P4/results/`. Everything below is double-precision, non-rigorous.

Modules: `tps.py` (truncated power series + dual numbers), `css.py` (Stage A: KHA
system, sonic data, first-order branches), `taylor.py` (high-order sonic expansions,
background and linearized), `shoot.py` (Stage B), `perturb.py` + `spectrum.py`
(Stage C), `diagnostics.py` (Stage D).

## Stage A — the CSS system (KHA polar-areal variables, γ = 4/3)

**Formulation used.** KHA95 eq. 18 (digest §1.1) with unknowns (A, N, ω, V) of
x = ln(−r/t), sonic point gauge-fixed at x = 0. Before building on the digest's
transcription I re-derived the full time-dependent system from scratch
(sympy: Einstein tensor of ds² = −α²dt² + a²dr² + r²dΩ², perfect fluid p = ρ/3,
∇·T = 0; change of variables t = −e^{−s}, r = e^{x−s}). Result:

- Rows 1, 2 (Hamiltonian constraint, slicing condition) agree with the digest:
  A_x/A = F_A := 1 − A + 2ω(1+V²/3)/(1−V²),  N_x/N = F_N := −2 + A − 2ω/3.
- The momentum constraint is (A_s + A_x)/A = G := −(8/3)NVω/(1−V²); for CSS
  solutions (A_s = 0) this is exactly KHA99 eq. (211), F_A = G.
- The digest's row 4 agrees exactly (reading "N_{,x}" literally as N·F_N).
- **The digest's row 3 has a transcription slip**: its last term must read
  2NV(1 + 4ω/(9(1−V²))), not 2N(1 + …) (a factor V dropped). With that fix both
  fluid rows are exact linear combinations of the derived ∇·T = 0 equations
  (residual 0 symbolically), *with the momentum constraint already used to
  eliminate a_t*, so the rows contain no A_s term. The corrected rows are what
  `css.coeffs` implements; the digest file itself was not modified.

**Sonic point.** det of the principal part is
Δ = −4(3N²V² − N² + 4NV − V² + 3)/(3ω(V²−1)); Δ = 0 ⇔ the fluid speed relative
to the x = const line, v_rel = (1+NV)/(N+V), equals ±1/√3. The digest's closed
forms N₀ = (√3−V₀)/(1−√3V₀), A₀ = (7+2√3V₀−3V₀²)/(4(1−V₀²)),
ω₀ = (3/8)(1−2V₀/√3−V₀²)/(1−V₀²) were verified symbolically to satisfy Δ = 0,
the row-proportionality condition af − ce = 0 (KHA99 215) and F_A = G.

**First-order branches.** With A₁ = A₀F_A(y₀), N₁ = N₀F_N(y₀) fixed, the
first-order coefficient V₁ satisfies the quadratic (derived with sympy, pasted
into `css.first_order_quadratic`)
c₂V₁² + c₁V₁ + c₀ = 0, c₂ = 12√3(1−V₀²)², c₁ = 12V₀(V₀²−1)(3V₀³−2√3V₀²−5V₀+2√3),
c₀ = V₀(9V₀⁶−9√3V₀⁵−21V₀⁴+10√3V₀³+39V₀²+3√3V₀−27), and ω₁ = −(e₀+b₀V₁)/a₀.
Discriminant: 432V₀(1−V₀²)⁴(3V₀³−5√3V₀²+3V₀+3√3), positive exactly for
V₀ ∈ (−1, −1/√3) ∪ (0, 1/√3), negative on (−1/√3, 0) (focus, no real analytic
crossing). The desingularised flow's eigenvalues along the two branches
(`css.sonic_eigenvalues`) have opposite signs on (0, 1/√3) — a **saddle**, as H03
says for EC at k = 1/3 — and equal signs on (−1, −1/√3) — a node. The flat
Friedmann solution sits exactly at the degenerate point V₀ = −1/√3 (A₀ = 1/(1−V₀²)
⇔ 3V₀²−2√3V₀−3 = 0), where both branches coincide.

**Checks (tests/test_css.py).** Zeroth-order data satisfy Δ = 0, af − ce = 0 and
the constraint to 1e−13; both branches satisfy the O(x⁰) equations and ℓ·O(x¹)
(ℓ = left null vector of the principal matrix) to 1e−14; the two branches are
distinct; the constraint identity F_A − G stays ≤ 1e−13 along solutions of the 4D
ODE integrated from sonic data; the order-40 sonic Taylor series agrees with direct
ODE integration to 1e−15 (radius of convergence ≈ 0.99 at the EC value of V₀).

**Regular centre.** In the scaled variables Â = (A−1)e^{−2x}, N̂ = Ne^{x},
ω̂ = ωe^{−2x}, V̂ = Ve^{−x} the regular fixed point is Â = 2ω̂/3, N̂V̂ = −1/2
(KHA99 220–222; both reproduced by the numerics). The irregular direction is a
central point mass, A − 1 ~ m̃ e^{−x} (equivalently NV + 1/2 ~ e^{−3x}); the
mismatch used for shooting is m̃(x_end) = (1−1/A)e^{x}, which is monotone in x.
**Pitfall found:** the constraint surface F_A = G is invariant but *repelling* for
the backward-x flow (the 4D fixed point has a second unstable direction transverse
to it), so naive 4D integration drifts off the constraint and never reaches the
fixed point. All shooting therefore uses the constraint-reduced 3D system
(N, ω, V) with A − 1 = 2ω(1 + V²/3 + (4/3)NV)/(1−V²) (`css.rhs_scaled3`).

## Stage B — shooting on V₀ (`p4/shoot.py`)

**Method.** Start at x = −δ (δ = min(0.1, R/6), R = estimated radius of the
order-40 sonic series; linear data with δ = 1e−5 when R < 0.02), integrate the
reduced scaled system with DOP853 (rtol 1e−13, atol 1e−16) toward x_end, with
terminal events for 2m/r → 1 (A → ∞), |V| → 1 and a second sonic point
(|Δ| relative < 1e−7, from either side). Coarse scan on 58 points of
V₀ ∈ (0, 1/√3) per branch (x_end = −8), brackets from sign changes of m̃,
brentq refinement (x_end = −15, xtol 1e−14).

**Results (saddle range V₀ ∈ (0, 1/√3)).**
- Branch 0 (V₁ < 0): no sign change of m̃ anywhere; every shot ends in A → ∞ with
  m̃ > 0 or at a second sonic point. No regular-centre solution on this branch.
- Branch 1 (V₁ > 0): one genuine root,
  **V₀ = 0.112439401388092** (reproducible to 1e−14 between runs; dm̃/dV₀ ≈ −0.49),
  sonic data (A₀, N₀, ω₀) = (1.8614267226, 2.0113168913, 0.3256888894),
  first-order (A₁, N₁, ω₁, V₁) = (−0.37029024, −0.71542382, −0.48249731, +0.48472202),
  desingularised sonic eigenvalues (−2.8247, +31.348) (saddle; EC lies on the
  branch with Δ increasing through x = 0, V₁ > 0; the other branch is the
  V₁ = −0.29653 direction, the LP/"other" root).
  A second sign change of m̃ near V₀ ≈ 0.1529 is a discontinuity (m̃ jumps from
  −0.02 to +0.38 as the termination reason switches), not a root.
- Regular-centre check for the root: at x = −8, Â − 2ω̂/3 = 1.1e−5, N̂V̂ + 1/2 = 7e−7,
  m̃ = 1.5e−10; limits N̂_∞ = 1.2365999612, ω̂_∞ = 5.82098013 (= 4πρ_c t² in
  polar-areal time), V̂_∞ = −0.40433391, Â_∞ = 3.880663. Beyond x ≈ −9 the
  double-precision limit on V₀ (≈1e−16, amplified by e^{−3x}) makes the solution
  leave the fixed point; this is inherent to one-sided shooting in double
  precision and does not affect Stages C–D (decided by x ≈ −6).
- Zeros of V: exactly **one**, at x = −0.2509 (inside the sound cone; V < 0 near
  the centre, V₀ = +0.112 at the sonic point, V → +0.8968 as x → +∞). Hunter (a)
  type ⇒ this is the Evans–Coleman solution.
- Gauge eigenvalue in the sonic-point gauge: κ̄ = −dN̄_ss/dx(0) = 2 − A₀ + 2ω₀/3
  = **0.3556992037**. The digest quotes KHA99 §V.7.2 as ≃ −0.35699; our value is
  0.355699 — the two agree if the quoted figure lost a digit ("0.35699" vs
  "0.355699"); a genuine 0.35699 would need V₀ smaller by ≈ 9e−4, incompatible
  with the sharp sign change of m̃. Stage C's independent eigenvalue search
  confirms 0.35570 (see below).
- Outward continuation (reduced plain system, DOP853 to x = 30, no events fired):
  a = √A → **1.065340**, 2m/r → 0.118903 (m/r → **0.059451**), Ω = 4πr²ρ = ω/A →
  **9.17397e−3**, V → 0.896838, N ∝ e^{−0.871x} ⇒ α ∝ r^{0.129} on t = 0.
  EC94: a → 1.07, m/r → 0.0596, Ω → 9.56e−3, α ~ r^{0.129} (their PDE-code
  values; ours agree to 0.5 %, 0.3 %, 4 %, and the exponent exactly).
  No zeros of V outside the sonic point; no second sonic point out to x = 30.
- Nodal range V₀ ∈ (−1, −1/√3): the sonic series has tiny radii (near-resonant
  node), the scan is dominated by early terminations, and no clean root was found
  besides the flat Friedmann point at the boundary V₀ = −1/√3; not pursued
  further (GRLP does not exist at k = 1/3; HM01).

## Stage D — monotonicity diagnostics on the sound cone (`p4/diagnostics.py`)

Profile on x ∈ [−8, 0] (regular to ≤ 1e−5 there); 2001-point dense grid, finite
differences; files `results/ec_profile.png`, `results/ec_monotonicity.png`,
`results/ec_profile.npz`, `results/ec_monotonicity.json`. Definitions: V = KHA's
fluid velocity relative to the static (r = const) observer (= HM01's V_R up to
sign); v_rel = (1+NV)/(N+V) = fluid velocity relative to the x = const (z = const)
line (= HM01's |V_z|, = 1/√3 at the sonic point); ρ̂ = 4πρt² = ωe^{−2x}/A;
w = 1/3 − v_rel² (GHJS's signed distance to the sonic point).

| quantity on [−8, 0] | behaviour | GHJS-style hypothesis |
|---|---|---|
| v_rel (HM01 V_z) | **strictly increasing** from 0 to 1/√3 (min d/dx = +1.4e−4 at x = −8, max 0.49) | holds (|V_z| increasing) |
| ρ̂ = 4πρt² = ωe^{−2x}/A | **strictly decreasing**, 5.821 → 0.178 (max d/dx = −9e−6) | holds (ρ̂' < 0) |
| ωe^{−2x} (no 1/A) | also strictly decreasing | holds |
| w = 1/3 − v_rel² | **strictly decreasing** (d/dx ∈ [−0.556, −3.7e−8]), w > 0 inside | holds ⇒ w'/w < 0 on the open cone |
| V (KHA, ≈ V_R) | **not monotone**: 0 → min −0.1095 at x ≈ −0.85 → zero at x = −0.2528 → +0.1124 at x = 0 | fails (one zero inside the cone) |
| \|V\| | not monotone (vanishes at x = −0.2528) | fails |
| A, 2m/r, Ω = 4πr²ρ | increasing except just inside the sonic point (A has its maximum ≈ 1.92 at x ≈ −0.3; A₁ = −0.370 < 0) | — |

Conclusion for Theorem C: the two monotonicity properties that GHJS actually use
(density decreasing outward, and the relative velocity to the similarity lines
increasing so that w'/w < 0 inside the cone) **do hold** for the EC profile at
k = 1/3 on the closed sound cone, in the variables v_rel = V_z and ρ̂. What fails is
monotonicity of the areal-frame velocity V (= V_R): it has its single zero inside
the cone at x = −0.2528 (r/(−t) = e^{x} = 0.777 of the sonic radius). So the
GHJS energy-method template is applicable provided it is formulated with V_z (the
velocity relative to z = const lines), not V_R; the collapsing/expanding sign
change of V_R lives inside the cone but does not enter the w'/w < 0 mechanism.

## Stage C — linear eigenvalue problem (`p4/perturb.py`, `p4/spectrum.py`)

**Linearization.** h = H_ss(x) + ε h_p(x) e^{κs} in the corrected time-dependent
rows of Stage A (∂_s → κ on perturbations; the rows already have a_t eliminated by
the momentum constraint, so no κA_p appears in them). The linearization is done
*exactly* with dual numbers around the numerical background (no sympy at runtime):
the same `css.coeffs` code evaluated on `Dual` objects gives DRes·y_p. The
linearized momentum constraint, κA_p = dC·y_p with C(y) := A(G − F_A) (so that
A_s = C), is solved algebraically for A_p, and the reduced linear system for
(N_p, ω_p, V_p) is co-integrated with the reduced scaled background (DOP853,
rtol 1e−12, atol 1e−14) for a whole batch of κ values at once.

**Sonic point (regular singular point).** Taylor expansion of the linearized
rows with the same ℓ·[Res]_n trick as for the background (`taylor.perturbation_series`):
at order 0 the unknowns (A_p, N_p, ω_p, V_p)(0) are fixed by (i) row
proportionality ℓ·[Res_lin]₀ = 0, (ii) the linearized constraint, (iii) the gauge
N_p(0) = 0, up to the normalisation A_p(0) = 1 — one free parameter, as in KHA99
§V.4. Orders n ≥ 1 follow linearly; the order-n 2×2 system is singular only at the
resonant values κ_n^res = −0.0990 − 1.0990 n (n = 1, 2, …; second Frobenius
exponent σ(κ) = n), all with Re κ < −1.19, so the matching function is analytic
on {Re κ > −1.19} ⊃ the search box. The linearized constraint is propagated by
the recursion to 1e−11 (check). Start at x = −0.1 with K = 36 terms.

**Centre condition.** With the asymptotic background (N̂, ω̂, V̂)_∞ the scaled
perturbation system (N_p e^{x}, ω_p e^{−2x}, V_p e^{−x}) has constant-coefficient
limit with exponents {−3, 0, 0} for every κ tested (0, 0.36, 2.81, 5, 2+3i,
10+10i; at κ = 1 exactly the linearized constraint degenerates at the centre,
κ − ∂_A C → 0, the origin-gauge value of KHA99 §VI, and the exponent set is
{−1, 0, 0}). Hence there is exactly **one** irregular direction — the central
point-mass mode A_p ~ e^{−x} — and the matching function is
E(κ) := A_p(x_end) e^{x_end}, evaluated at x_end = −9 (E changes by < 1e−7
between x_end = −8, −9, −10). E(κ) is real for real κ, E(0) = 1.048.

**Real eigenvalues.** Secant refinement (`spectrum.refine_zero`):
- **κ₁ = 2.8105525487765** (|E| < 1e−15); KHA99 Table 2: 2.81055255 — agreement
  to 1.2e−9, i.e. to all nine quoted digits; **β = 1/κ₁ = 0.355801922**
  (KHA95/99: 0.35580192). Robustness: x_end = −8 → 2.8105525478, −10 →
  2.81055254882; rtol 1e−13 → …87765; K = 30, δ = 0.05 → …87765. Combined
  numerical uncertainty ≈ 1e−9 (dominated by x_end/tolerance), consistent with
  the 1.2e−9 offset from KHA's value.
- **Gauge eigenvalue κ̄ = 0.3556992036** (|E| < 1e−14), equal to
  −dN̄_ss/dx(0) = 2 − A₀ + 2ω₀/3 = 0.3556992037 to 7e−11; the sonic-point data
  of the E-solution at κ̄ coincide with the pure-gauge generator
  (A_ss', N_ss' + κ̄N_ss, ω_ss', V_ss')/A_ss'(0) to 1e−8 (test). So the mode at
  0.35570 is the pure gauge mode of KHA99 eq. 257/285 in the sonic-point gauge;
  the digest's "0.35699" is a dropped digit.

**Mode count in the KHA box (argument principle, `spectrum.winding_number`).**
E sampled on the rectangle boundary (60 points per side, adaptively refined until
every phase increment is < π/2; E is analytic on Re κ > −1.19, so the winding
number equals the number of zeros):

| rectangle (Re κ) × (Im κ) | winding = #zeros | min |E| on contour |
|---|---|---|
| [0, 15] × [−14, 14] (KHA95 search box) | **2** | 4.5e−2 |
| [0, 5] × [−5, 5] | 2 | 0.11 |
| [2, 4] × [−1, 1] | 1 (κ₁) | 0.11 |
| [0, 1] × [−1, 1] | 1 (gauge κ̄) | 0.13 |

Hence the zeros of E in 0 ≤ Re κ ≤ 15, |Im κ| ≤ 14 are exactly κ̄ = 0.35570 (pure
gauge) and κ₁ = 2.81055255: **exactly one non-gauge eigenvalue**, as KHA95/99
state ("no other relevant or marginal eigenmodes" in that box). Relevant mode
count = 1 ⇒ γ_PBH = β = 1/κ₁ = 0.355802.

**Next eigenvalues (Re κ < 0).** E has simple poles at the sonic-point resonances
κ = −0.0990 − 1.0990 n, n = 0, 1, 2, … (the n = 0 one, σ(κ) = 0, comes from the
normalisation A_p(0) = 1: (κ + 0.0990)E(κ) → 0.14 as κ → −0.099, confirmed
numerically), so counts there are Z − P with P known.
- Strip −1.1 ≤ Re κ ≤ 0, |Im κ| ≤ 2: winding −1 = Z − 1 ⇒ Z = 0; on the real
  axis E < 0 throughout (−1.1, −0.1). No eigenvalue between the gauge mode and the
  first resonance.
- Real axis (−2.25, −1.25): one sign change; secant ⇒ **κ₂ = −1.5456213652**
  (|E| < 1e−14). This is KHA99's "eigenmode with Re κ ≲ −1.4" (§V.7.3, from their
  "less complete search" in −1.5 ≤ Re κ < 0, |Im κ| < 2 — our value lies just
  outside that window, which explains their inequality).
- Real axis (−3.35, −2.35): no sign change; E > 0 with a shallow minimum ≈ 0.755 at
  κ ≈ −2.98, suggestive of a nearby complex pair (KHA99's "fourth/fifth mode, a
  complex pair"); see the winding numbers below.

## Known-answer ledger (HEURISTIC → CONJECTURED)

| quantity | this work | literature | status |
|---|---|---|---|
| sonic-point data closed forms (KHA99 §IV.3.2) | verified symbolically at γ = 4/3 | digest §1.1 | ✓ |
| sonic point type at k = 1/3 | saddle on the EC branch range V₀ ∈ (0, 1/√3) | H03: saddle for k ≲ 0.41 | ✓ |
| EC solution parameter | V₀ = 0.1124394013881 (±3e−14), branch V₁ = +0.4847 | not tabulated | new |
| zeros of V (= V_R) | 1, at x = −0.2528 | Hunter (a): exactly one | ✓ |
| −dN̄_ss/dx(0) | 0.3556992037 | KHA99: "0.35699" | ✓ up to a dropped digit |
| t = 0 limits | a → 1.06534, m/r → 0.059451, Ω → 9.174e−3, α ∝ r^{0.129} | EC94: 1.07, 0.0596, 9.56e−3, r^{0.129} | ✓ (0.3–4 %) |
| κ₁ | 2.8105525488 (±1e−9) | KHA99: 2.81055255 | ✓ (1.2e−9) |
| β = 1/κ₁ | 0.355801922 | 0.35580192 | ✓ |
| gauge eigenvalue (sonic-point gauge) | 0.3556992036 = −dN̄_ss/dx(0) | KHA95/99 | ✓ |
| zeros in [0,15]×[−14,14] | 2 (gauge + κ₁) | KHA95: none other | ✓ |
| next real eigenvalue | κ₂ = −1.5456213652 | KHA99: Re κ ≲ −1.4 | ✓ |
| monotonicity on the cone | v_rel ↑, ρ̂ ↓, (1/3 − v_rel²) ↓; V not monotone | needed for Theorem C | new |

## What S1 settles for the formulation's [open] items

1. **Branch choice (Def. of EC, §1 [open]):** the EC solution crosses the sonic
   point along the first-order root with V₁ > 0 (Δ increasing through x = 0,
   desingularised eigenvalue +31.3); the other root (V₁ < 0, eigenvalue −2.82)
   carries no regular-centre solution at all on (0, 1/√3). In GHJ's language the
   EC branch is the root of Q that is *not* LP-type; here both roots are real for
   all V₀ ∈ (0, 1/√3) and only one of them admits a regular centre.
2. **Gauge (Theorem B [open]):** the KHA sonic-point gauge N_p(0) = 0 is
   well-conditioned when the momentum constraint is used to eliminate A (and A_p):
   the linearized constraint is a pointwise algebraic relation, the sonic-point
   expansion has one free amplitude, and the centre condition is a single scalar
   (one irregular direction, exponent −3). The gauge eigenvalue is
   κ̄ = 2 − A₀ + 2ω₀/3, an explicit function of the sonic data.
3. **Rectangle (Theorem B):** R = [0, 15] × [−14, 14] contains exactly κ̄ and κ₁;
   E is analytic on Re κ > −1.19 (nearest resonance/pole at κ = −1.198), so a
   Krawczyk/argument-principle certification on R faces no poles.
4. **Theorem C hypotheses:** the GHJS-type monotonicity holds on the closed sound
   cone for v_rel (= V_z) and ρ̂, and for w = 1/3 − v_rel²; it fails for V_R. The
   energy-method template should be set up in HM01/KHA variables with V_z.
5. **Analytic variable at the centre (A2 [open]):** the regular solution approaches
   the fixed point with corrections in integer powers of e^{x} (the scaled
   variables converge; the sub-leading corrections seen numerically are O(e^{2x})),
   consistent with the digest's remark that the HM01 exponent 2(1+3k)/(3(1+k)) = 1
   at k = 1/3.

## Caveats and follow-ups
- Double precision throughout; V₀ is limited to ~1e−13 by the e^{−3x} amplification
  of the one-sided shooting; the profile beyond x ≈ −9 is not trusted. A two-sided
  (centre ↔ sonic) matching or an mpmath integration would remove this; for S2 the
  validated centre expansion will be used anyway.
- The nodal range V₀ ∈ (−1, −1/√3) was only scanned coarsely (series radii are
  tiny there); the flat Friedmann solution is at its boundary.
- E(κ) is computed with the normalisation A_p(0) = 1; the resonant values
  κ = −0.099 − 1.099n are poles of E, not eigenvalues.

## Reproduction
```
uv run pytest problems/P4/tests -q                       # 28 tests, ~25 s
PYTHONPATH=problems/P4/src uv run python -m p4.diagnostics problems/P4/results
```
Stage C driver (≈ 2 min for the KHA box):
```python
from p4 import perturb, spectrum
prob = perturb.Problem(0.112439401388092, 1, K=36, x_end=-9.0)
spectrum.refine_zero(prob, 2.81, 2.812, real=True)        # kappa_1
spectrum.winding_number(prob, (0.0, 15.0, -14.0, 14.0))   # -> 2
```
