# P2 S1: P2-2 — high-frequency exclusion by Liouville–Green with weighted error control

Status: S1, 2026-10-16. This is the written proof behind `src/p2/lg_bounds.py` and the 102
certificates `certificates/lg/lg_{PB,PA}_l{0..50}.json`.
Replay: `python -m p2.lg_bounds verify PB 50` (and `PA`), about 3 s each.

## Statement (certified for X ∈ {P_B, P_A}, 0 ≤ l ≤ 50)

Let β = (½+δ)κ₋. With the exact dyadic W_l stored in `lg_X_l{l}.json`:

  W_l(ω) ≠ 0 for every ω with |Re ω| ≥ W_l and Im ω ≥ −β.

- The part Im ω ≥ 0 is P2-0.
- The part Re ω ≤ −W_l follows from the exact symmetry W_l(−ω̄) = conj W_l(ω). This holds
  because the radial coefficients are real, R₊(−ω̄) = conj R₊(ω), and c₊, c_c are real. It is
  checked in Arb in `tests/test_arbw.py`.
- What remains is a = Re ω ≥ W_l and −β ≤ b = Im ω ≤ 0. That is what is proved below.

## Setting

x = r*, with dx = dr/F and r* = Σ_i A_i ln|r − r_i| (an exact antiderivative). Φ = rR satisfies
Φ″ + p²Φ = 0, where p² = ω² − V(x), V = V_l and ′ = d/dx.

Because a² − b² − V(x) ≥ a² − β² − V_max > 0, we have Re p² > 0 on ℝ. Hence p = √(ω² − V) (principal
branch) is continuous, nonvanishing, and tends to ω at ±∞. As x → ∓∞,
V = O(e^{2κ₊x}) (resp. O(e^{−2κ_c x})) together with its derivatives.

Liouville–Green functions: y^± = p^{−1/2}e^{∓iS}, S′ = p. They solve y″ + (p² − θ)y = 0 with

  θ = p^{1/2}(p^{−1/2})″ = V″/(4p²) + 5V′²/(16p⁴),

and W(y⁺, y⁻) = 2i. With p² = A + iB, A = a² − b² − V > 0 and B = 2ab ≤ 0, we get Im p = B/(2 Re p) ≤ 0
and |Im p| ≤ a|b|/√A.

## Lemma 1 (Volterra equations with horizon references)

Set

  Φ₊ = y⁺(1 + h₊),  h₊(x) = −(1/2i) ∫_{−∞}^{x} [e^{2i(S(x)−S(t))} − 1] (θ/p)(t)(1 + h₊(t)) dt,
  Φ_c = y⁻(1 + h_c),  h_c(x) = −(1/2i) ∫_{x}^{∞} [e^{−2i(S(x)−S(t))} − 1] (θ/p)(t)(1 + h_c(t)) dt.

Here S is normalized by S(x) = ωx + ∫_{−∞}^{x}(p − ω) for Φ₊, and by S(x) = ωx − ∫_x^∞(p − ω)
for Φ_c. Then:

- Φ₊ and Φ_c solve the radial equation.
- Φ₊ ~ ω^{−1/2}e^{−iωx} as x → −∞, and Φ_c ~ ω^{−1/2}e^{iωx} as x → +∞.
- Φ₊ is a nonzero multiple of rR₊, and Φ_c of rR_c.
- W(Φ₊, Φ_c) = (r₊r_c)·(nonzero)·W_l(ω), where the nonzero factor is a product of the
  normalization constants and the entire, nonvanishing Jost exponential.

**Proof.** Variation of parameters for y″ + (p² − θ)y = −θy with the kernel
G(x,t) = [y⁺(t)y⁻(x) − y⁺(x)y⁻(t)]/(2i) gives the first equation; the second is its mirror image.

Kernel bounds. For t ≤ x,

  |e^{2i(S(x)−S(t))}| = exp(2∫_t^x |Im p|),

and symmetrically on the c side. Fix X. For x ≤ X, the kernel ½|e^{…} − 1| is bounded by
K_X(t) = ½(1 + e^{2∫_t^X|Im p|}), which is nondecreasing in X.

Convergence of the kernel integral. As t → −∞, |Im p| → |b| ≤ β < κ₊ and |θ/p| = O(e^{2κ₊t}), so
∫_{−∞}^X K_X|θ/p| < ∞. Picard iteration and Grönwall then give a unique solution with

  |h₊(x)| ≤ exp(Ψ₊(X)) − 1 for x ≤ X,   Ψ₊(X) := ∫_{−∞}^{X} K_X(t)|θ/p|(t) dt.

Applying the same estimate with X = x gives h₊(x) = O(e^{2κ₊x}).

Identification with rR₊. The Frobenius solution satisfies rR₊ = r₊e^{iωc₊}e^{−iωx}(1 + O(e^{2κ₊x})),
because r − r₊ = e^{2κ₊(x−c₊)}(1 + O(r − r₊)) and A₊ = 1/(2κ₊). Hence
D := Φ₊ − (ω^{−1/2}/(r₊e^{iωc₊}))·rR₊ is a solution with |D| = O(e^{−(2κ₊−|b|)|x|}). A nonzero
solution decays at most like the recessive e^{iωx}, i.e. like e^{−|b||x|}. Since 2κ₊ − |b| > |b|,
D ≡ 0. The c side is identical with κ_c (|b| ≤ β < κ_c).

The Wronskian identity is W(rR₊, rR_c) = Δ(R₊R_c′ − R₊′R_c), together with the FORMULATION
normalization of W_l. ∎

## Lemma 2 (the Wronskian cannot vanish)

Let

  Ψ̃₊ = ∫_{−∞}^{X} e^{2∫_t^X|Im p|}|θ/p| dt,

and define Ψ_c, Ψ̃_c on [X, ∞) in the same way. Then:

- |h₊′(X)/p(X)| ≤ e^{Ψ₊}Ψ̃₊. This follows by differentiating the Volterra equation (the integrand
  bracket vanishes at t = x) and using 1 + |h₊| ≤ e^{Ψ₊} on (−∞, X]. The same holds on the c side.
- Direct computation:

  W(Φ₊, Φ_c) = e^{i(S_c − S₊)}·2i·[(1+h₊)(1+h_c) + ((1+h₊)h_c′ − h₊′(1+h_c))/(2ip)],

  and S_c − S₊ = −Θ with Θ = ∫_{−∞}^{∞}(p − ω)dx. Hence

  |W_l(ω)e^{iΘ}/(2iω) − 1| ≤ e^{Ψ₊+Ψ_c}(1 + (Ψ̃₊ + Ψ̃_c)/2) − 1.

- So W_l(ω) ≠ 0 whenever

  e^{Ψ₊+Ψ_c}(1 + (Ψ̃₊ + Ψ̃_c)/2) < 2.        (*)

Float check of the whole chain (normalization, the phase Θ, and the bound):
`tests/test_lg.py::test_bound_dominates_actual_wronskian`. At Re ω = W_l and Im ω = −β, the actual
|W_l e^{iΘ}/(2iω) − 1| is about 0.01, against a certified bound of about 0.98–1.0. The bound is
deliberately crude: it uses absolute values and ignores the oscillation of the kernel.

## Lemma 3 (uniform majorants; what the certificate computes)

For all a ≥ a₀ and b ∈ [−β, 0]:

1. |p|² ≥ Re p² ≥ m(r) := a₀² − β² − V(r).
2. |θ/p| dx = J dr with J ≤ |(F V_r)_r|/(4m^{3/2}) + 5F V_r²/(16m^{5/2}). This uses V′ = F V_r and
   V″ = F(F V_r)_r; the right side is decreasing in a.
3. |Im p| ≤ βG(r), where G = max(1, a₀/√m) if β² + V(r) > 0 and G = 1 otherwise. Reason:
   a/√(a² − c) is decreasing in a for c > 0, increasing in c, and ≤ 1 for c ≤ 0.

So the integrals Ψ, Ψ̃ computed with these majorants at a = a₀ bound the true ones for every ω in
the region, and (*) at a₀ proves the statement for all a ≥ a₀.

Quadrature (`LG.side_bounds`), on r-cells from each horizon to the matching point r_X:

- **Interior cells.** On each cell [r_k, r_{k+1}] (exact dyadic endpoints, about 2.2k cells per
  side: geometric near the horizon, uniform beyond):
  - Arb interval evaluation of the polynomial/rational expressions gives J_k ≥ sup J and
    G_k ≥ sup G, with m > 0 certified on the cell.
  - The cell's r*-length is exact: dx_k = r*(r_{k+1}) − r*(r_k).
  - The weight exponent for the cell is bounded by I_k = Σ_{j ≥ k} βG_j dx_j (cells between the
    cell and r_X).
  - Contributions: ½(1 + e^{2I_k})J_k|cell| to Ψ, and e^{2I_k}J_k|cell| to Ψ̃.
- **End cells** [r₊, r₊ + e] and [r_c − e, r_c], e = 10⁻⁷·span.
  - Write r*(r₁) − r*(r) = A_h ln(|r₁ − r_h|/|r − r_h|) + (regular part), with |A_h| = 1/(2κ_h).
    The oscillation of the regular part on the cell is bounded by B.
  - The weight is ≤ E·e^{2βG B}(|r₁ − r_h|/|r − r_h|)^s, where s = βG/κ_h < 1 is certified
    (s ≈ 0.67 at P_B r₊, 0.25 at P_B r_c, 0.75–0.77 at P_A r₊, 0.52 at P_A r_c).
  - Integration gives Ψ_end ≤ ½J_max e(1 + E e^{2βGB}/(1 − s)) and Ψ̃_end ≤ J_max e E e^{2βGB}/(1 − s).

The threshold search (bisection on a₀ and two candidate matching points, r_P and argmax V_l) is a
design step. The certificate is the final Arb evaluation of (*) at the stored a₀ and r_X, which the
replay recomputes.

## Results

| | W_0/κ₋ | W_1/κ₋ | W_2/κ₋ | W_10/κ₋ | W_50/κ₋ | W_l/((l+½)Ω_P): l=0, 2, 10, 50 |
|---|---|---|---|---|---|---|
| P_B | 5.401 | 10.31 | 15.61 | 56.57 | 255.5 | 2.196, 1.269, 1.095, 1.029 |
| P_A | 2.695 | 4.951 | 7.411 | 25.96 | 115.6 | 2.439, 1.342, 1.119, 1.036 |

- W_l lies below the planning value 1.25(l+½)Ω_P + 2κ₋ of s0-gate §3 for every l except
  P_B, l = 0 (5.40κ₋ against 5.07κ₋).
- For large l it approaches the barrier top √V_max ≈ (l+½)Ω_P.
- The criterion value is in [1.947, 2) in every certificate (bisection to 0.2%).
- Cost: about 1 s CPU per (point, l) for search plus certification, 104 s CPU for all 102.
  Replay takes 0.05 s per l.
