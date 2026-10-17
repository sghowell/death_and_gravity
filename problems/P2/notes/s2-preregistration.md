# P2 S2: pre-registration of the part-(i) feasibility probe

Written 2026-10-17, **before any production run** of the S2 certificates (house
rule: criteria first). Only float design studies preceded it; they are
summarized in §2 so that the route choice is on record. Conventions:
FORMULATION v0.3; the point is P_B (ΛM² = 3/50, Q/Q_max = 999/1000,
δ = 1/20, β := (½+δ)κ₋ = 0.3268γ₀, so c* = 0.33).

## 1. What S2 must decide (s0-gate §7, coordinator brief 2026-10-17)

- Certify the barrier-top window at P_B for l = 20, 50 and 100. Use the
  inequalities of a part-(i) route (pointwise, scalable in h = 1/(l+½)), not
  part-(iii) zero counting.
- Fit the error constant C and project l₀.
- **Kill rule (as given).** If no route certifies l = 100, or the projected l₀
  at P_B is above 100: stop. P2 is then NO-GO-for-now, with the certified
  finite-l and part-(ii) results kept as partial evidence, explicitly not the
  entry point.

## 2. Route choice (recorded before production)

The brief names M2 (Riccati invariant tubes with an M1-quality centre) as the
leading route and M3 (complex deformation plus coercivity) as the alternate.
S2 runs **M3 first**:

- M2 at O(h) quality needs a centre that is uniform through the coalescence of
  two complex turning points: the complex-parameter Liouville–Weber map of
  route M1. That construction is the open analysis itself (s0-gate §4,
  "complex coalescence parameter").
- M3 reduces exclusion to a one-dimensional pointwise inequality whose h → 0
  limit is the harmonic oscillator, so the trapping gap is built in.

Float design study (scratch, 2026-10-17; reproduced by the S2 code). The
contour is x = x_P + e^{iπ/4}s through the photon sphere, with variable phase
φ(s) = φ_∞ s⁴/(s₁⁴+s⁴) and ε₁ = η|φ′|. The quantity is the best ground-state
margin divided by μ_HO := 2(l+½)Ω_Pγ₀:

| l | worst corner of the window | at the l-th PS mode (must fail) |
|---|---|---|
| 20 | 0.58 | −0.020 |
| 50 | 0.65 | −0.004 |
| 100 | 0.66 | −0.001 |

The ideal value at the strip bottom is 1 − β/γ₀ = 0.673. The bound fails at true
QNMs, so the control discriminates.

## 3. Window (frozen here)

For l ∈ {20, 50, 100}, with ω_l := (l+½)Ω_P:

  **𝒲_l = {ω : ω_l − 8γ₀ ≤ Re ω ≤ ω_l + 6γ₀, −β ≤ Im ω ≤ 0}.**

- In the scaled variable ζ = (h²ω² − Ω_P²)/h, this is Re ζ/k ∈ [−8, 6], where
  k = 2Ω_Pγ₀ is the oscillator quantum. So the window runs from 8 quanta below
  the barrier top to 6 above.
- Its upper end must reach the P2-2 threshold a_l: a_l ≤ ω_l + 6γ₀ is checked
  (certified) for each l, so 𝒲_l ∪ {Re ω ≥ a_l} has no gap above the top. For
  l = 100, a_l is computed with the S1 code (`lg_bounds`).
- Below ω_l − 8γ₀ lie the sub-barrier and low-frequency regimes (S4), which S2
  does not address.

## 4. Certificate structure (route M3)

1. **Contour geometry (independent of l and ω).**
   - Validated integration of dr/ds = e^{iπ/4}F(r) from r(0) = r_P gives
     enclosures of r(s) on s-cells for |s| ≤ S.
   - Beyond S: certified bounds |r(s) − r_h| ≤ d(s) from the exact identity
     Re r*(r(s)) = x_P + s cos θ.
   - The sectors between the real x-axis and the contour carry no singularity of
     V(x). This uses a fan of validated rays plus the horizon-funnel inverse
     function bounds.
2. **Coercivity certificate per l** (one ground state per window). Pointwise in
   s:

     E(s) := q + c′G′ − cG′² + cG″ ≥ μ > 0,

   where c = cos φ − η|φ′|, q = Re(e^{−iφ}e^{2iθ}(V_l(r(s)) − ω²)) − |φ′|/(4η),
   and the ground state is g = e^{−G}, with G′ continuous and piecewise linear.
   - The inequality is checked on every cell by the mean-value form
     (Arb, outward rounding), for ω² ranging over the interval hull of
     {ω² : ω ∈ 𝒲_l}. E is affine in ω², so the hull adds no overestimation in ω².
   - The tails |s| > S are checked by the bound in 1.
3. **Implication** (written, notes/s2-m3.md). A QNM in 𝒲_l would give an L²
   solution of −u″ + e^{2iθ}(V − ω²)u = 0 on the contour, since the Jost
   solutions continue through the sectors and decay along the rays. But
   Re⟨e^{−iφ}(−u″ + Qu), u⟩ ≥ μ‖u‖² > 0. Hence W_l(ω) ≠ 0 on 𝒲_l.

## 5. Acceptance criteria (KAT-style; fixed now)

- **A0.** Geometry certificate passes:
  - every ray of the fan stays in |r| ≥ r_min > 0;
  - the ray spacing is below half the certified analyticity radius;
  - the funnels are certified;
  - the tail bounds hold.
- **A1.** For each l ∈ {20, 50, 100}, every cell and both tails satisfy
  E ≥ μ_cert with μ_cert > 0 (certified), on all of 𝒲_l.
- **A2 (negative control).** With the same design, the float margin at the
  l-th PS fundamental mode (from the S1 Frobenius solver) is ≤ 0. The rigorous
  checker, run on a small box around that mode, must not certify.
- **A3 (consistency).** μ_cert(l) ≤ the float best margin λ_min(l) for the same
  design.
- **A4 (replay).** `python -m p2.m3 verify` recomputes every certificate
  verdict.

## 6. Error constant and projection (fixed now)

Let m(l) := μ_cert(l)/μ_HO(l), with μ_HO(l) = 2(l+½)Ω_Pγ₀.
- Least-squares fit of m(l) = m_∞ − C/(l+½) on l = 20, 50, 100.
- Projected l₀ := max(0, ⌈C/m_∞ − ½⌉), the l at which the fitted margin
  vanishes.
- Also reported: the fit m(l) = m_∞ − C₂/(l+½)², since the design study
  suggests O(l⁻²).

**Verdict rule.**
- GO for S3 if A0, A1 (all three l), A2, A3 and A4 hold and the projected
  l₀ ≤ 100.
- NO-GO-for-now if A1 fails at l = 100 or the projected l₀ > 100.
- Anything else is reported as INCONCLUSIVE, with the failing item named.
