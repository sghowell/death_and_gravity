# P4 S2 — Theorem B (rectangle mode count): validated linearised problem

Status: Stage 1 of 3 done (2026-08-29). Code: `problems/P4/src/p4/validated/{linsys,linsonic}.py`
(new; nothing pre-existing modified). Tests: `uv run pytest problems/P4/tests/test_validated_modes_stage1.py -q`
(7 tests, ≈ 2 s). Every number below is a rigorous enclosure (256-bit balls) unless marked "float".

## 1. Stage 1 — the linearised system and its sonic-point series for complex κ

### 1.1 The linear ODE (`linsys.py`)

Perturbation h = H_ss(x) + ε h_p(x) e^{κs} of the time-dependent KHA system in the cleared
polynomial form of A1 (`systems.sonic_system`: P(u)u_x + P_s(u)u_s = Q(u), u = (A, N, W, V), rows 3–4
multiplied by 3WS resp. WS, so P_s = [[0],[0],[0,0,3S,4VW],[0,0,4VS,4(1+V²)W]]). Linearised exactly:
    P(u) p′ + Ψ(u,u′) p + κ P_s(u) p = DQ(u) p,   Ψ_{rl} = Σ_i u_i′ ∂P_{ri}/∂u_l.
Linearised momentum constraint (C̃ = 0 along the background — proven in A1 — and ∂C̃/∂A = S):
    (κ − A) S A_p = A (C̃_N N_p + C̃_W W_p + C̃_V V_p)          (S1: κA_p = dC·y_p, C = AC̃/S, ∂_A C = A).
Rows 2–4 contain no A_p′ (P's A-column vanishes off row 1), so eliminating A_p and multiplying by
(κ − A)S gives the **linear 3D system for q = (N_p, W_p, V_p)**, polynomial in (u, u′, κ):
    𝒫(x;κ) q′ = 𝒢(x;κ) q,   𝒫 = (κ−A)S P_{rq},   𝒢 = (κ−A)S (DQ − Ψ − κP_s)_{rq} + A DQ_{rA} C̃_q,
    𝒫 = 𝒫₀(u) + κ𝒫₁(u),  𝒢 = Σ_{j≤2} κ^j [Ga_j(u) + Σ_i Gb_{j,i}(u) u_i′]   (exact `fmpq_mpoly` data).
This is exactly the system S1 integrates (`perturb.linearized_rhs`; E(κ) is defined through it).
`LinSystem.coefficient_matrix(u, κ)` returns 𝒜(x;κ) = 𝒫⁻¹𝒢 (3×3 `acb_mat`) from background balls u at
a regular point (u′ by Cramer on the fluid block), `A_p(u, κ, q)` the constraint value; backgrounds:
A1 Taylor model on [−0.05, 0] (`SonicExpansion.eval`), the A3 tube (`background_from_state(tmint.State)`,
interval or point set, via `plain_from_scaled`), the A2 centre family for x ≤ −3.
Checks (test): against S1's dual-number linearisation on S1's float sonic series at x = −0.05, −0.02,
−0.3 (complex κ): |𝒜 − 𝒜_S1| ≤ 1e−12 (entries up to 200; S1's rounding), ball radii ≤ 4e−72; on the
certified backgrounds (A1 model at −0.05, −0.03; A2 family at −3, −4, radii ≤ 6e−40) the balls contain
S1's values evaluated on the midpoints.

### 1.2 Sonic-point series for κ in a box (`linsonic.py`)

Regular singular point: 𝒫(0) = (κ−A₀)S₀ P(u₀) has the background's fluid null vector ℓ, and the level
equations are those of A1, E_n = (F^N_{n−1}, F^{fl}_{n−1}, ℓ·F_n), F = 𝒫θq − 𝒢q, with
    M_n(κ) = n D(κ) + E(κ),  D = (κ−A₀)S₀ D̂ (D̂ the background structure matrix restricted to q),
    E = −(0, 0, ℓ·𝒢₀(κ)),
so D is invertible iff κ ≠ A₀ and M_n(κ) is linear in κ for n ≥ 1 (the resonances of S1). Order 0:
N_p(0) = 0 (gauge), ℓ·𝒢₀ q₀ = 0 (row proportionality) and A₀(C̃_W W_p0 + C̃_V V_p0) = (κ−A₀)S₀
(normalisation A_p(0) = 1): one free amplitude, as in S1/KHA99. Orders 1..K: one ball solve each;
success proves M_n(κ) invertible for **every κ in the box** (non-resonance certificate; a box of
half-width 5e−3 around the first resonance −1.198 is rejected with a singular M₁ — test).
q₀ = (0, 0.09347531, −0.37577332) at κ₁ (S1: same to 1e−8).

**κ-dependence (Taylor model in δ = κ − κ_c, κ_c real).** Naive ball propagation of a box κ
through the recursion would blow up ≈ 8× per order (dependency effect); instead
q_n(κ) = Σ_{k≤m} q_{n,k} δ^k + R_n with the polynomial part from the recursion at the exact centre
(δ a formal variable, level solve triangular in the δ-degree, exact) and, for the remainder, the
(m+1)-th δ-coefficient from the same recursion with the **box** as base point, which encloses
q_n^{(m+1)}(ξ)/(m+1)! for every ξ in the box. For complex δ the Lagrange form is unavailable; the
integral remainder q_n(κ) − Σ_{k≤m} q_{n,k}δ^k = δ^{m+1}/m! ∫₀¹(1−t)^m q_n^{(m+1)}(tδ)dt along the
segment [κ_c, κ] ⊂ box (convex) gives |R_n| ≤ sup_box|q_{n,m+1}| |δ|^{m+1}, |δ| ≤ √2 w. A_p is
treated the same way (series division by (κ−A)S in δ-polynomial arithmetic); plain ball division over
the box wraps (radii 1e7 at n = 40 — the first obstacle met).

**Results (V₀ = 0.112439401388092 point A1 background, K_bg = 41, K = 40, m = 5, w = 1e−6, i.e.
κ ∈ κ_c ± 1e−6 ± 1e−6 i; 256 bits; 0.2 s per box including the certificate).**

| box | ball radii n = 0 / 10 / 20 / 30 / 40 | remainders n = 0 / 20 / 40 | values at x = −0.05 (radii of A_p, N_p, W_p, V_p) |
|---|---|---|---|
| κ₁ = 2.8105525488 | 2.3e−7 / 4.1e−6 / 1.1e−5 / 2.6e−5 / 1.2e−2 | 6.5e−39 / 4.3e−17 / 8.3e−3 | 1.0932642, −0.0995691, 0.1291094, −0.4033548 (3.4e−6, 1.9e−8, 2.6e−7, 9.3e−8) |
| κ̄ = 0.3556992037 | 4.5e−6 / 6.9e−5 / 1.7e−4 / 3.6e−4 / 6.9e−2 | 1.4e−33 / 2.1e−17 / 4.8e−2 | 0.8654719, −0.0046731, 1.3640121, −1.2785223 (2.5e−5, 2.8e−7, 3.8e−6, 4.7e−6) |

The radii up to n ≈ 30 are the true κ-sensitivity (|∂q_n/∂κ|·√2w); at n = 40 the box remainder
dominates (the box run's top coefficient grows ≈ 7.7×/order; m = 7 brings it to 4e−5), which is
irrelevant for values on |x| ≤ 0.05 and for the ν-weighted norms of §1.3. Point runs (w = 0) have
radii ≤ 3e−41. The certified gauge eigenvalue κ̄ = 2 − A₀ + 2W₀/3 = 0.355699203710973 ± 1.3e−16
(A1 balls) lies in the second box; S1's κ₁ = 2.8105525487765 lies in the first.
Test (c): at κ = κ₁ (point) every ball of (A_p, N_p, W_p, V_p)_n, n ≤ 38, contains S1's
`taylor.perturbation_series` float (relative tolerance 1e−10; observed 1.5e−16 at n = 0 to 1.2e−12 at
n = 38, S1's rounding; S1's orders 39–40 suffer its TPS truncation).

**Row 1 of the 4D linearisation** (S A_p′ − 2VA′V_p − DQ₀·p, not used by the 3D system) is checked as a
ball identity on the computed series: all K coefficients contain 0 with |ρ_n| ≤ 1.2e−39 (point runs)
— the Bianchi identity holds order by order, so the constraint-eliminated series is the S1 4D series.

### 1.3 Tail bound on |x| ≤ 0.05 (`certify_linear`) — the argument

Write q = q̄ + z (q̄ = balls n ≤ K over the box, z the tail). For n > K the exact recursion is
M_n z_n = −[E_n(q̄) + Σ_{K<m<n} M_{nm} z_m], M_{nm} = Σ_{(c,r,s)} c[𝒫_{n−m+s+1}·m − 𝒢_{n−m+s}] (row
triples of the level equations), so T(z)_n := −M_n⁻¹[E_n(q̄) + Σ_{K<m<n} M_{nm} z_m] is **affine**
on ℓ¹_ν (‖z‖_ν = Σ_{n>K}|z_n|_∞ν^n). With g := ‖D⁻¹E‖ < K+1 (D, E at the box κ), ‖M_n⁻¹‖ ≤ c/n,
c = ‖D⁻¹‖/(1 − g/(K+1)), and m/n ≤ 1, 1/n ≤ 1/(K+2): the linear part has norm ≤ Z and ‖T(0)‖ ≤ Y with
    Z = c Σ_{k≥1} ‖B̃_k‖_∞ ν^k,  B̃_k[i,l] = Σ|c|(|𝒫^{(r)}_{k+s+1}[l]| + |𝒢^{(r)}_{k+s}[l]|/(K+2)),
    Y = c Σ_{n>K} |E_n(q̄)|_∞ ν^n / n.
The coefficient matrices of the *true* background u = ū + v are split as (polynomial part along the
truncated ū: finite exact sums, Z₁, Y₁) + (increment over the background tail: Z₂, Y₂), using the A1
certificate at ν < ν_u: ‖v‖_ν ≤ ε_u(ν/ν_u)^{K_u+1} =: ε_v and ‖θv‖_ν ≤ ε_θ (`deriv_tail_bound`), and the
Banach-algebra majorants inc(p) = p^{abs}(‖ū‖_ν + ε_v) − p^{abs}(‖ū‖_ν) of the polynomials 𝒫_j, Ga_j and
Σ_i[inc(Gb_{j,i})‖θū_i‖_ν + Gb_{j,i}^{abs}(‖ū‖_ν+ε_v) ε_θ] for the Ψ-part, weighted by |κ|^j over the box:
    Z₂ = c‖[Σ|c|(ν^{−(s+1)} inc𝒫 + ν^{−s} inc𝒢/(K+2))]‖_∞,   Y₂ = (c/(K+1)) Σ_i Σ|c| ν^{−s} Σ_l(inc𝒫_{rl}‖θq̄_l‖_ν + inc𝒢_{rl}‖q̄_l‖_ν).
If Z = Z₁ + Z₂ < 1, T is a contraction of ℓ¹_ν and its unique fixed point is the true tail (the
recursion has a unique solution since every M_n is invertible: n ≤ K by the ball solves, n > K by
g < K+1), so Σ_{n>K}|q_n|ν^n ≤ ε := Y/(1−Z), |q_n| ≤ εν^{−n}, and Σ_{n>K}|q_n||x|^n ≤ ε(|x|/ν)^{K+1}.
Everything is evaluated with the box κ, so the certificate holds for every κ in the box.

| box | ν | c | g | Z₁ | Z₂ | Y | ε | tail on |x| ≤ 0.05 |
|---|---|---|---|---|---|---|---|---|
| κ₁ | 0.09 | 1.280 | 6.82 | 0.902 | 1.4e−38 | 4.3e−39 | 4.4e−38 | 1.5e−48 |
| κ̄ | 0.08 | 0.740 | 3.71 | 0.939 | 1.6e−41 | 1.9e−41 | 3.1e−40 | 1.3e−48 |

(ε_v ≈ 5e−43, ε_θ ≈ 3e−39 from the A1 certificate ν_u = 0.1, ε_u = 3.8e−41, K_u = 41.) Consistency
test: the K = 40 bound |q_n| ≤ εν^{−n} holds for the exactly computed orders 41–50 of a K = 50 run.
`LinSonicExpansion.eval(x)` returns (A_p, N_p, W_p, V_p)(x) over the box with both tails (q and the
A1 background) included — the validated starting data for Stage 2 at x₀ = −0.05.

### 1.4 Timings, what is enclosed, obstacles

Timings (M-series laptop, 256 bits): A1 background K = 41 + certificate 0.08 s; one box expansion
(point run + box run + A_p models + row-1 check) 0.2 s; certificate < 0.05 s; full test file 2.2 s.
Enclosed for every κ in each box: the unique analytic solution of the 3D linearised system with the
gauge/normalisation above (coefficients n ≤ 40 as `acb` balls, a convergent tail bound on |x| ≤ 0.09
resp. 0.08), and its values at x = −0.05; non-resonance of every order for the whole box.
Obstacles / caveats: (1) A_p by plain series division over a box κ wraps (fixed by the δ-model);
(2) the 3D reduced system carries the factor (κ − A(x)): for real κ ∈ (1, A₀) it has an *apparent*
singularity at the x where A(x) = κ (the 4D system has none) — harmless on ∂R (|Im κ| = 14 or
Re κ ∈ {0, 15}) and near κ₁, κ̄, but Stage 2 must keep the tube away from it or switch to the 4D form
for such κ; (3) the certified ν (0.08–0.09) is ≈ 1/10 of the empirical radius (≈ 1), the same loss as
in A1, and needs ν < ν_u = 0.1 for the background's derivative tail; (4) row 1 is verified to order K as
a ball identity, not as an exact polynomial identity (the E(κ) of Theorem B is defined by the 3D
system as in S1, so nothing depends on it); (5) larger κ-boxes: the box run's blow-up (≈ 7.7×/order)
with rem_n ∝ w^{m+1} means w = 1e−4, m = 5 gives rem₄₀ ≈ 1e10 (useless as coefficients, harmless
on |x| ≤ 0.05 where only rem_n·0.05^n enters); for Stage 3 (rectangle contour) many small boxes or a
higher m are needed.

Reproduction:
```
uv run pytest problems/P4/tests/test_validated_modes_stage1.py -q
PYTHONPATH=problems/P4/src uv run python -c "
from p4.validated import sonic, linsonic
bg = sonic.sonic_expansion('0.112439401388092', K=41); bg.certify()
ex = linsonic.linear_sonic_expansion(bg, '2.8105525488', width=1e-6, m=5, K=40); print(ex.certify())
print([z.str(10, radius=True) for z in ex.eval(-0.05)])"
```
