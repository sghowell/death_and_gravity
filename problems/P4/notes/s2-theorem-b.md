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

## 2. Stage 2 — validated propagation, the centre condition and E(κ)

Status: Stage 2 of 3 done (2026-08-29). Code (new, nothing pre-existing modified):
`problems/P4/src/p4/validated/{linscaled,lintube,linstep,linprop,lincentre,linmatch}.py`; tests
`uv run pytest problems/P4/tests/test_validated_modes_stage2.py -q` (6 tests; the three fast ones ≈ 20 s, the
full pipeline ≈ 10 min — or ≈ 5 min with `P4_TUBE_CACHE=<tube.json>` saved by `Tube.save`). Every number is
a rigorous 256-bit ball unless marked float.

### 2.1 The linear systems in the centre-scaled variables (`linscaled.py`)

Propagation in the plain variables is hopeless toward the centre (‖𝒜‖ ~ N/W ~ e^{-3x}), so the Stage-1
system is rewritten exactly for q̃ = Λ⁻¹q, Λ = diag(T⁻¹, T², T), T = e^x (S1's scaled perturbation): with
u = (A, N, W, V) = ((S + 2wT²T̃)/S, n/T, wT², vT), S = 1 − v²T², T̃ = 1 + v²T²/3 + 4nv/3, u′ expressed
through (z, z′) (A′ by differentiating the constraint), every row multiplied by the nonvanishing
S^{dA+2}T^{dN+2} and the common T-power of each row stripped (exact `fmpq_mpoly` division):
    Σ_j κ^j P_j(z) q̃′ = Σ_j κ^j G_j(z, z′) q̃,   z = (n, w, v, T), z′ = (n′, w′, v′),  j ≤ 1 (P), j ≤ 2 (G),
G_j = [Ga_j + Σ_i Gb_{j,i} u_i′]Λ − P_jΛΓ, Γ = diag(−1, 2, 1). Check (test): P⁻¹G = Λ⁻¹𝒜Λ − Γ against the
Stage-1 matrix at a point (|diff| ≤ 5e−73). The **4D system** (A_p kept: P p′ = (DQ − Ψ − κP_s)p, κ-degree
1, no (κ − A) factor) is built the same way with Λ₄ = diag(T², T⁻¹, T², T) (Â_p = A_p e^{-2x}), Γ₄ =
(2, −1, 2, 1); restricted to the linearised constraint surface it reproduces the 3D system (5e−73, test).
The κ-derivative system for (q̃, ∂_κ q̃) is assembled block-wise, P̂ = [[P, 0], [P₁, P]],
Ĝ = [[G, 0], [G₁ + 2κG₂, G]], again κ-graded polynomial matrices (`augmented`), and is propagated with the
weighted norm (derivative block × λ = 2⁻¹², as tmint's y-block) and a blockwise ball inverse of P̂
(interval LU of the 6×6 box matrix leaves a wide zero block that 1/λ amplifies — the first blow-up met).

### 2.2 The certified background along the whole range (`lintube.py`)

Tube region [−3, −0.05]: tmint is run for the V0 box c* ± 1e−16, c* = c + 6.3e−15 = 0.1124394013880983,
which contains the certified EC value (A3: V0* − c ∈ 6.3e−15 ± 7.5e−17), so the interval set contains the
EC background. Per step the exact point Taylor coefficients z_K (K = 28) through the reference midpoint
m_k are recorded together with ε_z = tail_u(h) + rk(sup|Y_K| + bound) (rk the weighted radius of the sets
before the step, `bound` tmint's Grönwall bound of the step Jacobian, both captured from the integrator as
arb; Y_K the 4D fundamental-matrix series — the mean-value argument of `tmint.jacobian_step`'s refinement),
so that every solution of the interval set satisfies |z(s) − z_K(s)| ≤ ε_z on the step: 244 steps, 189 s
(384 bits), ε_z = 3.6e−16 (x = −0.05) → 2.1e−14 (x = −3), interval-set width 2.1e−14 at −3 (207 × w).
Centre region [−8, −3]: the A2 series at the certified (a*, μ*) (balls, radii 9e−16, 5e−14; K = 30,
ν = 0.08, ε = 1.1e−15) gives z(x_k) as balls; the Taylor coefficients through the ball contain the true
ones and ε_z = the A2/A3 Banach tail at h (h = 0.05, 100 steps, ε_z ≈ 1e−35).
κ-independent data per step (`derive`, `system_data`; 0.15 s/step for the 3D system, 0.3–0.4 s for the 4D
and 6D ones, once per process): sub-boxes Z_j of the tube (8 exact pieces of [−h, 0], + ε_z), z′ = f(Z_j),
L_bg = sup‖Df‖, Pinv_bg, the residual sup R of the truncation, ε_dz = L_bg ε_z + Pinv_bg R (bound on
|z′ − z_K′|), coefficient sums a = sup|z_K|, sup|z_K′|; per linear system the s-polynomial matrices
P_j(z_K(s)), G_j(z_K(s), z_K′(s)) exact to order K with the discarded orders bounded through the
coefficientwise majorant p^{abs}(Σ|z_i|s^i) (tails ≈ 1e−70), box values on the Z_j, and the increment
majorants inc P_j, inc G_j = p^{abs}(a + ε) − p^{abs}(a) (per-variable ε = ε_z, ε_dz). The tube is saved
and loaded exactly (dyadic mantissa/exponent of every ball).

### 2.3 The validated linear step and the propagation (`linstep.py`, `linprop.py`)

For fixed κ (point) a step x_k → x_k − h: (1) Y_K(s) = Σ_{i≤K} Y_i s^i, Y₀ = I, from
(i+1)P₀Y_{i+1} = Σ G_iY_{n−i} − Σ (n+1−i)P_iY_{n+1−i}; (2) the defect along the *true* background,
    E = Y_K′ − 𝒜(z)Y_K = P(z)⁻¹[D_K + (P(z) − P(z_K))Y_K′ − (G(z, z′) − G(z_K, z_K′))Y_K],
    D_K = P(z_K)Y_K′ − G(z_K, z_K′)Y_K (exact polynomial in s: orders < K rounding only, ≥ K truncation),
    sup|E| ≤ Pinv[Dsup + (δP + |inc P|)Y′sup + (δG + |inc G|)Ysup],  Pinv, L = sup|P⁻¹|, sup|P⁻¹G| over the Z_j;
(3) Grönwall on Δ = Y − Y_K: ‖Y(−h) − Y_K(−h)‖ ≤ h e^{Lh} sup|E| =: bound (∞-norm, weighted for the 6D
system; every quantity an arb upper bound); (4) Lohner set {c + A r} in C^d: c + Ar → Y_K(−h)(c + Ar) +
box(bound·‖y‖_∞), complex QR re-orthogonalisation. Typical values: L = 216 (x = −0.05) → 50 (−0.9) → 4.3
(x ≤ −2), Pinv = 15 → 0.55, Dsup ≤ 1e−35 on the tube steps (K = 28, h ≤ 0.02) and ≈ 1e−17 on the centre
steps (h = 0.05), bound = 4e−13–6e−13 per tube step (dominated by the increment terms, i.e. by ε_dz = L_bg ε_z ≈
4e−14: the V0 tube, not the Taylor truncation). Widths of the propagated sonic solution: 4e−12 (x = −0.08),
5e−10 (−0.65), 1.4e−8 (−1.25), 2.3e−7 (−3), 1.0 at −8 — relative to |q̃| ∝ e^{−3x} a constant ≈ 5e−9.
Test: on a 15-step tube the propagated Stage-1 solution and Φ(−0.08, −0.05)q(−0.05) both land inside the
Stage-1 series' own enclosure at x = −0.08 (radii 5e−12 vs 1e−13).
κ **boxes** are propagated at the box centre κ_c with the additional per-step bound (Δ = Φ(κ) − Φ(κ_c),
Δ′ = 𝒜(κ)Δ + (𝒜(κ) − 𝒜(κ_c))Φ(κ_c)): bound_κ = h e^{L_box h} sup‖P(κ)⁻¹‖(‖G(κ) − G(κ_c)‖ + ‖P(κ) − P(κ_c)‖‖𝒜(κ_c)‖)(Ysup + bound)
from the cached box data — rigorous, but crude (§2.5); putting the box κ into the Taylor recursion instead
blows up by ≈ 6×/order (the Stage-1 dependency effect with h/ν ≈ 0.75).

### 2.4 The centre condition and E(κ) (`lincentre.py`, `linmatch.py`)

At t = e^x → 0 the scaled **4D** system along the A2 series is Fuchsian, P₀ invertible for every κ, and
B₀ = P₀⁻¹G₀ has spectrum {−3, −3, 0, 0} for every κ (numerically at 0, 0.356, 1, 2.81, 1 ± 0.5i, 5 + 7i;
the reduced 3D form has {−3, 0, 0} and degenerates at κ = 1, P₀ ∝ (κ − 1)). Rows A, W of G₀ vanish
identically (as polynomials), so rank G₀ = 2 and the exponent-0 (regular) family is 2-dimensional and
analytic in κ, including κ = 1; no positive integer is an exponent. Basis: (w̃, ṽ)(0) = (1, 0), (0, 1) with
(Â, ñ)(0) from rows A, W of G₀ (its (Â, ñ) minor is κ-independent, singular values 28.2, 3), then
(nP₀ − G₀)p_n = Σ_{k≥1}[G_k − (n−k)P_k]p_{n−k} (K = 50). Tail: the Euler-form affine contraction of §1.3
with ‖M_n⁻¹‖ ≤ c/n, Z = c Σ_k(‖G_k‖/(K+1) + ‖P_k‖)ν^k + Z₂, background tail via ℓ¹_ν majorants with the A2
certificate (ε_v = tail_bound(ν), ε_θ = ν·deriv_tail_bound(ν)): at ν = 0.06, g = 5.0, c = 1.11, Z₁ = 0.545,
Z₂ = 2e−17, ε = 5e−17, 1.5e−16 (uniform in κ over the tests); the basis values at x = −3 satisfy the
linearised constraint as a ball identity (test). With the propagated vector at x_d = −8 (Â_p from the
constraint when propagated in 3D) and T_d = e^{x_d}:
    E(κ) := e^{3x_d} det[r₁, r₂, p̃(x_d)] restricted to the rows (Â_p, ñ_p, ṽ_p)  ( = ± det[r₁, r₂, p̃, e_W] ),
zero iff p̃(x_d) lies in the regular plane, analytic in κ on the whole rectangle (4D data), same zeros as
S1's E. (The 3D determinant det[r₁ᵠ, r₂ᵠ, q̃] would vanish identically at κ = A(x_d) ≈ 1 + 3e−7, where the
constraint surface projects degenerately to q̃ — a spurious zero inside R; e_W is never in the constraint
surface since ∂C̃/∂W = −2T̃ ≠ 0.) E_fin(κ) := e^{x_d}A_p(x_d) is S1's function at x_end = x_d, enclosed
exactly. dE/dκ: the augmented propagation with ∂_κ q̃(x₀) from the Stage-1 δ-polynomials (box remainder
(m+1) sup|q_{n,m+1}| |δ|^m, tail by Cauchy's estimate over the box enlarged by ρ = 1e−6) and the augmented
regular family (d = 8, same certificate); checked against finite differences (rel. 1e−6).

### 2.5 Results, the apparent singularity, timings, and the plan for Stage 3

Pipeline (`linmatch.Matcher`): Stage-1 series over the V0 box (the A1 box expansion re-packaged with a
ball-valued fluid null vector, `box_background`) → q̃(x₀ = −0.05) → 3D (or 4D) propagation through the 344
steps → Â_p(x_d) from the constraint → 4D regular family at t_d = e^{−8} → E, E_fin. S1 float = `perturb.Problem(V0, x_end=-8).E`.

| κ | E(κ) (ball) | E_fin(κ) (ball) | S1 float E(x_end = −8) | dE/dκ | t (s) |
|---|---|---|---|---|---|
| κ₁ = 2.8105525488 | 0 ± 6.2e−10 | 0 ± 1.5e−9 | 2.0e−10 | 0.021939 ± 6e−7 | 124 (6D, incl. 120 s derivation) |
| κ̄ = 0.3556992037 | 0 ± 1.1e−7 | 0 ± 1.0e−6 | −7.8e−10 | −0.060 ± 7e−4 | 6 |
| 2.5 | −0.00574942 ± 2.6e−9 | −0.05412792 ± 7.4e−9 | −0.054127923 | — | 1.5 (cached) / 53 (first 3D call) |
| 3.0 | 0.004615063 ± 3.6e−10 | 0.043448503 ± 5.3e−10 | 0.0434485029 | — | 1.5 |
| 1 ± 0.5i | −0.0156235 ∓ 0.0016114i (± 6e−8) | −0.147088 ∓ 0.015171i (± 5e−7) | −0.1470878 ∓ 0.0151706i | — | 2.7 |
| 1.2 (4D) | −0.01406685 ± 3.3e−9 | −0.132432397 ± 7e−10 | −0.13243184 (S1 crosses the singularity: off by 5.6e−7) | — | 21 (incl. 4D derivation) |
| 0 | 0.111315859 ± 7.5e−10 | 1.04798363 ± 5.6e−9 | 1.04798363 | — | 1.5 |
| 14i | 0.00050543 − 0.00472644i (± 8e−9) | 0.0047585 − 0.0444981i (± 6e−9) | 0.004758505 − 0.04449806i | — | 2.7 |
| 15 | 1054.879 ± 5.6e−4 | 9930.934 ± 2.0e−4 | 9930.93396 | — | 1.5 |

E_fin/E = 9.4109 (± 1e−6) at every κ tested: at x_d = −8 the regular part of p̃ is O(e^{3x_d}) relative, so
S1's function is the analytic determinant up to a constant factor (the basis normalisation) — the same zeros,
and S1's sign pattern on the real axis (E(2.5) < 0 < E(3.0), E(0) > 0) and the conjugate symmetry are
reproduced by E. **Zeros:** E(κ₁) ∋ 0 with radius 6.2e−10 and |dE/dκ| = 0.0219, so the zero of E lies within
|E|/|dE| ≤ 3e−8 of S1's κ₁ *if* dE varies little over that distance — the rigorous localisation (interval
Newton) needs dE over a κ box, i.e. the Stage-3 Taylor model (below); the same at κ̄ (|E/dE| ≤ 2e−6; the
gauge run is ≈ 200× less tight, its Lohner widths at x_d being 5e6 vs 4e3 for κ₁ — the gauge mode is
not dominated by the point-mass direction early on). 4D vs 3D at κ = 2.5: E = −0.00574942 ± 1.5e−9 (4D)
overlaps the 3D ball; E_fin identical to 5e−10.

**(c) The apparent singularity, real κ ∈ (1, A_max).** Along the background A(x) rises from A₀ = 1.861 to
A_max ≈ 1.89 (x ≈ −0.4) and then decreases to 1.0095 at x = −3 and to 1 + 3e−7 at x = −8; the reduced 3D
system has P ∝ (κ − A(x)), so for real κ ∈ (1, A_max) some sub-box of the tube contains a singular P and the
3D propagation stops with a singular box matrix (κ = 1.2: at x ≈ −1.47, test). The 4D system (A_p kept)
has no such factor; its solution with the sonic data (on the linearised constraint surface, which is
invariant) is analytic through x_κ, and its q̃-part *is* the 3D solution wherever (κ − A) ≠ 0 (the 4D and
3D right-hand sides coincide on the constraint surface: exact polynomial identity checked as balls, §2.1).
So E is defined on the whole rectangle by the 4D propagation, which is what `Matcher.E(..., system="4D")`
does (κ = 1.2: table). Analyticity in κ inside R: q̃(x_d; κ) is entire in κ from the 4D system (polynomial
κ-dependence, regular on [x_d, x₀]) given analytic initial data, and r₁, r₂ are analytic (4D Frobenius,
uniform in κ, §2.4); hence E is analytic on R wherever the *sonic* data are. Two caveats for Stage 3 (both
interior to R, the contour is unaffected): (i) Stage 1's sonic recursion is the 3D form, whose leading
matrix D ∝ (κ − A₀) degenerates at κ = A₀ = 1.861 — its certificate already fails for real κ ≈ 1.3–2.3
(Z₁ = 1.7 at κ = 1.5), so the sonic-side data near A₀ must come from a 4D sonic-point series (regular
singular point of the 4D linearised system) or from a removability argument; (ii) the disk around κ = 1
where the 3D centre form degenerates is harmless (E uses the 4D family), but the 3D-to-4D conversion at
x_d has the factor 1/(κ − A(x_d)) — for |κ − 1| < 1e−6 propagate in 4D from x_c or x₀.
On ∂R the Stage-1 certificate holds with ν ≥ 0.05 at all sampled points; at the corner 15 + 14i ν = 0.05
exactly (Z₁ = 0.991 at 0.05, 1.001 at 0.0505), which does not cover x₀ = −0.05 → build the Stage-3 tube from
x₀ = −0.045 (the A1 model is certified on |x| ≤ 0.1; everything else unchanged).

**κ boxes.** With the per-step Lipschitz bound of §2.3 a box w = 1e−8 around κ = 2.5 gives E_fin =
−0.05 ± 5.9e−3, w = 1e−9: ± 3.0e−4 (κ₁, w = 1e−8: E = ± 3.7e−5), i.e. relative radius ≈ 1e7·w: the
per-step bound is only ≈ 5× pessimistic (‖∂𝒜/∂κ‖ ≈ 3100 near the sonic point), but E is a ≈ 1e−2 component
of the propagated solution, so every relative error in q̃ is amplified ≈ 70×, and the box error
accumulates over 344 steps. Rigorous, but useless beyond w ≈ 1e−9; the box result trivially contains the
point result (test). The κ-box treatment for the contour must therefore be a **Taylor model in κ**: the
order-m augmented system (block-Toeplitz: P̂ with P on the diagonal and P₁ below, Ĝ with G, G₁ + 2κG₂, G₂ on
the three diagonals; `augmented` is m = 1) propagated at the box centre gives E_k = ∂_κ^k E/k! exactly (the
weighted-norm bounds of §2.3 with λ per block), and the remainder from the (m+1)-th derivative over the
box with the crude bound: relative remainder ≈ 1e7·w·w^{m+1}|E^{(m+1)}/E|/(m+1)!. With |E^{(k)}/E| = O(1)
(E varies on the scale 1 in κ; min|E| on ∂R is 4.8e−3 = 0.045/9.41 from S1) m = 3, w = 0.02 gives 1e−3 and
m = 5, w = 0.05 gives 1e−5: 860–2150 boxes for the 86-long contour at ≈ 10–20 s each (one-time derivation
of the 3(m+2)-dimensional system's step data, reusable if assembled from the 3D blocks) — a few hours.
Phase bookkeeping: consecutive boxes overlap, |E| ≥ 4e−3 on each, phase steps ≪ π/2.

**Timings** (M-series laptop). Tube: 189 s (tmint, 384 bits, 244 steps) + 6 s centre extension (100 steps);
κ-independent step data 0.15 s/step (3D), 0.3 (4D), 0.4 (6D) — 53 / 100 / 120 s once per process; per κ
afterwards 1.5 s (3D), 2.7 s (complex κ), 3 s (6D), plus 0.2 s Stage-1 series and 0.1 s centre family.
**What is enclosed.** For the EC background (V0* ∈ c* ± 1e−16, (a*, μ*) as certified balls) and each κ: the
unique solution of the linearised system with the Stage-1 gauge/normalisation, its value at x_d = −8 as a
Lohner set (relative width ≈ 1e−8 for point κ), the 2-dimensional regular family at the centre with
certified tails, and E(κ) = e^{3x_d} det[r₁, r₂, p̃(x_d)]_{(Â,ñ,ṽ)} as a ball (radii 4e−10–6e−8 on the
samples), with dE/dκ when requested; E_fin = S1's function at x_end = −8 as a ball.
**Stage 3 plan.** (1) tube from x₀ = −0.045; (2) order-m κ-Taylor model of the propagation and of the
regular family, remainder from the crude box bound; (3) contour of R = [0, 15] × [−14, 14] in ≈ 10³ boxes,
winding number 2 (S1) ⇒ exactly κ̄ and κ₁ in R once (4) the interior analyticity is closed: 4D sonic-point
series near κ = A₀, 4D propagation on (1, A_max), the κ = 1 disk; (5) interval Newton at κ₁ and κ̄ with the
box derivative from (2) (point data: |E(κ₁)| ≤ 6.2e−10, dE = 0.0219; |E(κ̄)| ≤ 1.1e−7, dE = −0.060).

Reproduction:
```
PYTHONPATH=problems/P4/src uv run python -c "
from flint import arb; from p4.validated import lintube, centre, linmatch; from p4.validated.arbseries import precision
tube = lintube.Tube.build('0.1124394013880983', 1e-16, x_c=-3.0)                       # 190 s
with precision(256):
    a = arb('-0.2123656467659762832750918714540807905889') + arb('7.031e-12') + arb(0, 8.87e-16)
    mu = arb('8.901323275379966931515526907200000000000') + arb('2.3757e-9') + arb(0, 4.61e-14)
    ce = centre.centre_expansion(mu * (2*a).exp(), nhat=(-a).exp(), K=30); ce.certify()
    tube.extend_centre(ce, -8.0, hmax=0.05)
    M = linmatch.Matcher(tube, ce, linmatch.box_background('0.1124394013880983', 1e-16))
    r = M.E('2.8105525488', deriv=True); print(r['E'], r['E_fin'], r['dE'])"
```

## 3. Stage 3 — the mode count in R: winding number, interior analyticity, Krawczyk, gauge mode

Status: Stage 3 of 3 done (2026-08-29). Code (new, nothing pre-existing modified):
`problems/P4/src/p4/validated/{linsonic4,lintaylor,modecount,analyticity,krawczyk_kappa}.py`; tests
`uv run pytest problems/P4/tests/test_validated_modes_stage3.py -q` (8 fast tests ≈ 45 s; 2 slow ones with
`P4_TUBE_CACHE=<tube.json>`; the stored full-R results `results/theorem_b/{winding_R,analyticity_R}.json` are
re-verified by 2 further tests). Revised 2026-08-30 after an adversarial review (§3.8: what changed).
Every number below is a rigorous 256-bit ball unless marked float. The tube is the Stage-2 A3 tube
[−3, −0.05] (244 steps, V0 ∈ c* ± 1e−16), *without* the centre extension: E is evaluated at x_d = −3 (§3.3).

### 3.1 The 4D sonic-point series and its κ-Taylor model (`linsonic4.py`)

Stage 2's caveat (i) is real: the constraint-reduced 3D recursion has D = (κ−A₀)S₀D̂ and a rank-one E, so
det M_n^{3D}(κ) = n²(κ−A₀)²S₀² det D̂ · [n(κ−A₀)S₀ + tr(D̂⁻¹E(κ))] with E quadratic in κ — one root of the
bracket is the resonance κ_n^res, the other is a *ghost* root that tends to A₀ as n → ∞. The ghost roots
ruin every box computation near the real axis: the Stage-1 box remainder at κ₁ is 87 (relative) already
for w = 0.05, m = 5, and the box runs at w = 0.2 meet singular level matrices (measured). The 4D system
P(u)p′ = [DQ − Ψ − κP_s]p (A_p kept, P κ-free, G affine in κ) has none of this: its level equations are
exactly A1's, E_n = (F^A_{n−1}, F^N_{n−1}, F^row_{n−1}, ℓ·F_n), so M_n(κ) = nD₄ + E₄(κ) with
    D₄ κ-independent (block-triangular: P is diag(S, 1, fluid block)), invertible;  E₄(κ) = E⁰ + κE¹ of rank 1,
hence M_n is singular iff n = σ(κ) := −tr(D₄⁻¹E₄(κ)), a single *affine* exponent:
σ(κ₁) = −2.64741095423 ± 7e−14, σ(κ̄) = −0.41375455226, σ(0) = −0.09010597051, σ(15+14i) = −13.7385 − 12.7385i
(S1's resonances κ = −0.0990 − 1.0990n: σ(κ) = −(κ + 0.0990)/1.0990). Order 0: A_p(0) = 1, N_p(0) = 0,
ℓ·G₀p₀ = 0 and the linearised constraint (a 2×2 solve for (W_p, V_p)(0), singular iff σ(κ) = 0, the pole of
E at −0.099); the constraint surface is invariant (exact identity, §3.8) and the series stays on it (`constraint_residual`:
a ball identity, all orders contain 0; proof in §3.8). Tail certificate as in §1.3 with the rank-one sharpening
‖M_n⁻¹‖ ≤ (c/n), c = ‖D₄⁻¹‖(1 + ‖D₄⁻¹E₄‖/(K+1 − max Re σ)) (no g < K+1 needed).
Checks (test): at κ₁, κ̄, 0, 15+14i every coefficient ball (n ≤ 40) overlaps the Stage-1 3D ball and the
values at x₀ = −0.05 overlap (4D radii ≤ 1e−11); the 4D certificate is at least as good as the 3D one:
ν = 0.09 at κ₁, κ̄, 0, **and at κ = 1.5 and κ = A₀ = 1.861** (Z₁ = 0.74, where the 3D form fails), ν = 0.065 at
the corner 15+14i (3D: 0.05) — so x₀ = −0.05 and the existing tube stay valid on the whole contour.
The 4D coefficient balls widen ≈ 5×/order (rounding amplified along the off-constraint formal direction of
the recursion; radii 3e8 at n = 40 for a point run) — harmless for |x| ≤ 0.05 and for the ν-weighted
certificate (5e−32 ≤ ε ≤ 1e−24 on the samples).
**κ-Taylor model** (point run in δ-polynomial arithmetic + (m+1)-th coefficient of the box-based run,
exactly Stage 1's construction, `delta_model(x₀, r)`; remainder includes the box tail). Relative remainder
at x₀ (float):

| κ_c | w = 0.25, m = 8 | w = 0.25, m = 12 | w = 0.5, m = 12 | w = 1e−7, m = 2 |
|---|---|---|---|---|
| 15 + 14i | 9.4e−8 | 3.7e−12 | 1.1e−6 | 1.7e−26 |
| 14i | 3.3e−9 | 1.3e−13 | 6.6e−8 | 3.6e−26 |
| 15 | 1.9e−11 | 1.6e−16 | 2.6e−11 | 2.5e−26 |
| 7.5 + 14i | 1.4e−6 | 1.9e−10 | 2.5e−4 | 2.6e−26 |
| κ₁ | 2.5e−8 | 4.6e−12 | 1.4e−5 | 8.9e−24 |
| κ̄ | (w = 0.05) 1.5e−7 | (w = 0.1) 2.2e−2 | — | 1.2e−20 |
| 0 | (w = 0.02) 3.0e−3 | (w = 2⁻⁷) < 1e−5 | — | 1.2e−18 |

The pole at −0.099 limits the boxes on the side Re κ = 0 (w ≈ dist/10, §3.4); elsewhere w = 0.25 with m = 8
is tight to 1e−6. 0.15 s per box model.

### 3.2 Taylor-model propagation of the fundamental matrix (`lintaylor.py`)

The scaled 4D system P(s)y′ = (G_c(s) + δG₁(s))y is propagated as a Taylor model in δ = κ − κ_c, |δ| ≤ r:
(1) block-Toeplitz Taylor recursion for Y_k(s) = ∂_κ^k Y/k! (k ≤ m; the order-m κ-derivative system is lower
block-triangular Toeplitz, so only its first block column is computed: (m+1)·K² products of 4×4 matrices);
(2) Grönwall bound of the block system in the weighted norm ‖ŷ‖_λ = max_k λ^k|y_k| (defect D_k = PY_k′ − G_cY_k −
G₁Y_{k−1} along the true background, ‖𝒜̂‖_λ ≤ L + λL₁ with L₁ = sup‖P⁻¹G₁‖; λ = 1, the κ-scale of the
solution: every block then carries the same absolute error as the value — λ = r loses 4 % on dE at Krawczyk
scale, measured); (3) **per-step Cauchy remainder**: Y(−h; κ) is entire, sup_{|δ|≤ρ}‖Y(−h;κ)‖ ≤ e^{h(L+ρL₁)}, so
the degree-m truncation error is ≤ e^{h(L+ρL₁)}(r/ρ)^{m+1}/(1 − r/ρ) with ρ = (m+1)/(hL₁). This is what makes
the approach work: in the 4D scaled variables L₁ = 20.4 at x = −0.05 (the 3D form had ‖∂𝒜/∂κ‖ ≈ 3100 — its
1/(κ−A) factor), h = 0.00187 there, so hL₁ ≤ 0.04 everywhere and ε_Y = 3e−23 per step at m = 8, r = 0.25;
(4) composition (Σ_j Y_jδ^j + R_Y)(S₀ + Σ_k S_kδ^k): block-Toeplitz map of the polynomial part (m+1 Lohner sets
S_k in C⁴, the lower blocks entering block k as balls), the degrees m+1..2m and R_Y·y added to S₀.
Three measured facts fixed the design. (a) An *initial* remainder box (the sonic model's 9.4e−8) propagated as
a set wraps ≈ 3000× along the 244 steps; hence the Taylor model of the **fundamental matrix** (exact
unit-vector data) is propagated and composed with the sonic model afterwards by an exact Taylor-model product
(`tm_apply`). (b) At Im κ = 14 the solutions rotate fast and an axis-aligned complex box loses up to √2 per
step: the Lohner radii are **discs** {c + Ar : |r_j| ≤ ρ_j} (rotation-invariant, and the ∞-norm Grönwall
error is a disc anyway) — 140× tighter at x = −1.3 on the tile 2.25 + 14i. (c) Even so a separate
isotropic remainder set wraps faster than the thin sets S_k (e^{8/unit x} vs e^{5.5/unit x} against the
solution's e^{3/unit x}); every δ-independent error is therefore added to the thin set S₀, whose width stays
at the level of the Grönwall boxes: at 2.25 + 14i the total model width went from 2.8e−4 to 7e−9 relative.
At 15 + 14i, r = 0.25, m = 8: |Φ_k|r^k = 1.0e7, 2.0e6, 1.9e5, 1.2e4, 5.8e2, 22, 0.73, 0.03, 0.015 (κ-convergence
scale ≈ 2.5), 38 s per box (244 steps, 0.16 s/step); the Stage-2 point propagation at δ = 0, 0.25, 0.25i and
−0.175 + 0.075i lies inside the model (test, on a short tube; relative width there 1.2e−7, the sonic remainder).

### 3.3 E at x_d = −3, the centre Taylor model, and why x_d does not matter

E(κ) := e^{3x_d} det[r₁, r₂, p̃(x_d)]_{(Â,ñ,ṽ)} = ± e^{3x_d} det[r₁, r₂, p̃, e_W] exactly as in §2.4, but with
x_d = −3 (the tube end; the regular family is certified there with ν = 0.06 > e^{−3}, tails ≤ 1e−20) and the
sonic data from §3.1. **x_d-independence.** For x_d′ < x_d: r_i(x_d′) = Φr_i(x_d), p̃(x_d′) = Φp̃(x_d), Φ = Φ(x_d′, x_d)
entire in κ, det Φ ≠ 0, and Φ⁻¹e_W = α(κ)e_W + σ with σ ∈ Σ(x_d) (the linearised constraint surface, a
κ-analytic hyperplane not containing e_W since ∂C̃/∂W = −2T̃ ≠ 0) and α = ℓ_Σ(Φ⁻¹e_W)/ℓ_Σ(e_W) ≠ 0 because Σ is
invariant (e_W ∉ Σ(x_d′) = ΦΣ(x_d)). Since r₁, r₂, p̃ ∈ Σ (proven for every κ in §3.8: exact invariance identity + exponent
conditions at the sonic point and at the centre), det[r₁, r₂, p̃, σ] = 0 (four vectors in a 3-space) and
    E_{x_d′}(κ) = e^{3(x_d′−x_d)} det Φ(x_d′, x_d; κ) α(κ) · E_{x_d}(κ),
an analytic nonvanishing factor: zeros, multiplicities and the winding number are the same for every
x_d ≤ −3. (Float check: E₃/E_fin^{S1}(x_end = −3) = 0.1109 + 0.0677i, 0.1119 + 0.0689i, 0.1129 + 0.0700i at
15 + 14i − 0.25, 15 + 14i, 15 + 14i + 0.25 — the smooth factor, not a constant, as expected at x_d = −3.)
**Centre model** (`RegularFamilyTM`): the §2.4 recursion in δ-polynomial arithmetic at κ_c (the (Â, ñ)
minor is κ-free, p₀ is affine in δ), the (m+1)-th coefficient from the box-based run (encloses
r_i^{(m+1)}/(m+1)! on the box), tail from `RegularFamily(box κ).certify(0.06)` (box runs converge since
(8·e^{−3})ⁿ → 0). **E** is then the 3×3 determinant in scalar Taylor-model arithmetic (`TM`: products
truncated at degree m, degrees m+1..2m and the cross terms into the remainder). At 15+14i, w = 0.25, m = 8:
E(15+14i) = 639.2 − 853.8i (± 0.02), |c_k|w^k = 1.1e3, 2.2e2, 23, 1.6, 0.084, 3.6e−3, 2.3e−4, 1.1e−4, 1.1e−4,
remainder 9.3e−3 (relative 9e−6) before (b)–(c), 2.4e−7 at the tile 2.25 + 14i after; 38 s.

### 3.4 The winding number of E around ∂R (`modecount.py`)

Contour: ∂R counterclockwise (bottom Im κ = −14, right Re κ = 15, top, left Re κ = 0), each side parametrised by
s ∈ [0, ℓ] and tiled by segments [s, s + 2w] with dyadic half-widths w ≤ 0.25 (w ≤ dist(κ, −0.099)/10 next to the
pole, i.e. 2⁻⁷ = 0.0078 around κ = 0), consecutive tiles sharing their endpoints *exactly* (all positions are
dyadic floats). On a tile: the E model of §3.3 on |δ| ≤ w (m = 8); the tile is split into n_sub = 8 pieces; on
each piece the enclosure B_j = E(piece) must satisfy Re(B_j e^{−iθ_j}) > 0 with θ_j = arg mid B_j (a ball
comparison: 0 excluded and the continuous argument on the piece is the principal argument relative to θ_j,
whose variation is < π); the tile's argument change is telescoped, Σ_j(b_j − a_j) = b_{n−1} − a_0 +
Σ_j[arg rot_j − arg rot_{j+1} + 2πℓ_j] with the integers ℓ_j fixed by the ball b_j − a_{j+1} (radius < π
required), so that only the two endpoint arguments contribute to the radius (the first run summed the piece
increments and lost 6.5 rad on the tiles near Re κ = 2 where |E| ≈ 0.005). A tile whose increment radius is not
< 2e−3 or whose sonic model is not tight (relative remainder > 1e−6) is halved. The increments are stored as
exact dyadic arb balls and summed in arb; the winding number is the integer N with |Σ/2π − N| < 1/2.
Run (8 workers, spawn; 172 chunks of length 0.5; the κ-free step data of the 244 steps derived once per worker):

| side | tiles | half-widths | min |E| on the side | max increment radius | max propagation remainder / |p̃| | mean s per tile |
|---|---|---|---|---|---|---|
| bottom (Im κ = −14) | 42 | 0.25 ×18, 0.125 ×24 | 4.82e−3 (at Re κ ≈ 1.3) | 4.3e−4 | 1.0e−6 | 32 |
| right (Re κ = 15) | 56 | 0.25 ×56 | 1039 | 1.4e−5 | 1.6e−7 | 29 |
| top (Im κ = 14) | 42 | 0.25 ×18, 0.125 ×24 | 4.82e−3 | 4.3e−4 | 1.0e−6 | 30 |
| left (Re κ = 0) | 150 | 0.125 ×102, 2⁻⁴ ×10, 2⁻⁵ ×10, 2⁻⁶ ×12, 2⁻⁷ ×16 | 4.81e−3 (at Im κ ≈ ±13.6) | 9.1e−5 | 4.3e−6 | 29 |

**Total: Σ Δarg / 2π = 2.00000 ± 8.4e−4 ⟹ winding number 2**, certified; 290 tiles, min |E| on ∂R = 4.81e−3
(re-run 2026-08-30 with the fixes of §3.8 from `python -m p4.validated.modecount`, stored as
`results/theorem_b/winding_R.json` with every tile increment as an exact dyadic ball: 1.99999996 ± 8.39e−4, the same
290 tiles, the same per-side table to the digits shown, 1231 s wall on 8 workers; re-summed by
`test_full_rectangle_result_if_available`)
(S1's float minimum 4.5e−2 was for E_fin; E₃ = E_fin × the factor of §3.3), 1255 s wall (8665 CPU-s), the
conjugate symmetry E(κ̄) = E(κ)̄ visible tile by tile (the bottom and top rows are conjugates to all digits).
Corner/side values (floats of the endpoint balls): E₃(0) = 0.113077, E₃(15) = 1125.08, E₃(±14i) = 7.31e−4 ∓
4.77e−3 i, E₃(15 ± 14i) = 639.19 ∓ 853.80 i (radii ≤ 2e−2 at the corners, 6e−7 where |E| is small).
Cross-checks: the small rectangle [2.5, 3.1] × [−0.3, 0.3] gives 1.0000 ± 1.9e−4 (κ₁ only; 8 tiles, 66 s on 4
workers, test); the fallback rectangle [0, 5] × [−5, 5] gives 2.00000 ± 1.6e−4 (158 tiles, min |E| = 1.2e−2, 36 min on 2 workers) — consistent with the S1 float table of `numerics-report.md`.
The first full run (isotropic remainder set, piece-wise summation) gave 1.99 ± 5.2 — winding "2" but *not*
certified — which is what forced the three fixes of §3.2 (b), (c) and the telescoping above; the second
run's tiles are the same 290 (no extra halving was needed once the models were tight).

### 3.5 E is analytic on a neighbourhood of R (`analyticity.py`) — the hypotheses of the argument principle

E(κ) = e^{3x_d} det[r₁(κ), r₂(κ), Φ(x_d, x₀; κ)p(x₀; κ)] is analytic wherever its three ingredients are:
(S) **sonic data.** On each box of a cover of R the 4D box run of §3.1 certifies: M_n(κ) invertible for
1 ≤ n ≤ K = 40 (ball solves) and for n > K (rank-one bound with max Re σ < K + 1) — non-resonance on the
box — and Z < 1 with ν > |x₀|, uniformly on the box. The p_n(κ) are rational without poles on the box, the
partial sums converge uniformly, so p(x₀; ·) is analytic on the open box (Weierstrass). No transport is
needed (ν ≥ 0.065 > 0.05 everywhere).
(P) **propagation.** P is κ-free and ball-invertible on all 8 sub-boxes of all 244 steps (`tube_regular`,
max ‖P⁻¹‖ = 11.7), G is affine in κ: Φ(x_d, x₀; ·) is entire.
(C) **regular family.** On each box `RegularFamily(box κ).certify(0.06)` proves P₀ invertible, (nP₀ − G₀(κ))
invertible for 1 ≤ n ≤ 50 (ball solves) and for n > 50 (g = ‖P₀⁻¹G₀‖ ≤ 6.12 < 51: no positive-integer
exponent), rank G₀ = 2 (κ-free minor) and Z < 1: the exponent-0 family is 2-dimensional for every κ and its
basis r₁, r₂ (normalised by (w̃, ṽ)(0) = (1,0), (0,1)) is analytic in κ on the box.
Cover: R tiled by 420 squares of half-width 0.5, bisected where a certificate fails (only next to the pole
at −0.099: 432 boxes in the end, min half-width 0.125); every tile is certified on the box enlarged by 5 %,
so the open certified boxes cover the closed rectangle and E is analytic on a neighbourhood of R.
Results: min ν = 0.065, max Re σ = 0.105 (< 1: no resonance in R; the resonances are at σ ∈ {1, 2, …}, i.e.
Re κ ≤ −1.198, and the pole σ = 0 at −0.099 ∉ R), max Z_sonic = 0.99974, max Z_centre = 0.572, max g_centre
= 6.12; 76 s with 2 workers (re-run 2026-08-30 with the certify4 fix of §3.8, `results/theorem_b/analyticity_R.json`:
the same 432 boxes and the same numbers to the digits shown, max ‖P⁻¹‖ = 11.74, 79 s; the JSON also records the
constraint-exponent certificate of §3.8 and is re-checked by `test_stored_analyticity_cover`). Consequences used below: (i) for every κ ∈ R the admissible sonic solution with
N_p(0) = 0 is unique up to scale and has A_p(0) ≠ 0 (σ(κ) ≠ 0), so the normalisation A_p(0) = 1 loses
nothing; (ii) E(κ) = 0 ⟺ p̃(x_d; κ) ∈ span(r₁, r₂) ⟺ κ is an eigenvalue in the class of §3.8 (using e_W ∉ Σ(x_d)
and r₁, r₂, p̃ ∈ Σ, both proven for every κ ∈ R in §3.8, and r₁(x_d), r₂(x_d) independent: a dependence
would give a solution vanishing at x_d, hence ≡ 0 near x_d (P(t_d) invertible, `tube_regular`) and on the
whole disc |t| < 0.06 by analyticity, against (w̃, ṽ)(0) = (1, 0), (0, 1)); (iii) E is analytic on a neighbourhood of R, so the winding number of
E|∂R equals the number of zeros in R counted with multiplicity.

### 3.6 Krawczyk enclosures of the two zeros (`krawczyk_kappa.py`)

With the Taylor model of §3.3 on |δ| ≤ w around κ_c (m = 3): E(κ_c) = c₀, E′(B) = Σ_k k c_k B^{k−1} + R′ on the
box B = κ_c ± w/2 with |R′| ≤ rem/(w − √2·w/2) (Cauchy for the remainder's derivative; |δ| ≤ √2·w/2 on the square box), Y = 1/mid(c₁),
K(B) = −Y E(κ_c) + (1 − Y E′(B))·(B − κ_c), E(κ_c) = c₀ + R(0) ∈ c₀ ± rem. K(B) ⊂ int B (acb boxes over-approximate the real 2×2 interval
products, so the test is sound) ⟹ exactly one zero of E in B, contained in κ_c + K(B). κ_c is first
Newton-refined from the float S1 values (one extra model, 20–40 s).

| zero | κ_c (after Newton) | w/2 | E(κ_c) | E′(κ_c) | E′(B) | K(B) radius | enclosure |
|---|---|---|---|---|---|---|---|
| κ₁ | 2.8105525488271472 | 1e−7 | ± 1.0e−10 (incl. R(0) ≤ 1.6e−23) | 0.022490638 ± 4.3e−10 | 0.02249064 ± 4.9e−9 | 4.4e−9 | 2.81055254883 ± 5.61e−9 (Im ± 4.4e−9) |
| κ̄ | 0.3556992037109642 | 1e−5 | ± 2.2e−10 (incl. R(0) ≤ 7.4e−18) | −0.061293328 ± 8.8e−10 | −0.06129 ± 6.4e−6 | 4.6e−9 | 0.35569920371 ± 8.35e−9 (Im ± 4.6e−9) |

Outward-rounded: **κ₁ ∈ [2.81055254439, 2.81055255326]**, **κ̄ ∈ [0.35569919907, 0.35569920835]**, both with
|Im| ≤ 4.7e−9, both simple (E′(B) ∌ 0), and **γ = 1/κ₁ ∈ [0.3558019218, 0.3558019231]** (ball 0.35580192 ±
3.0e−9). The enclosure radii are set by |E(κ_c)|/|E′| with rad E(κ_c) = 1e−10 (the Lohner width of the
propagation, 5e−10 relative), not by the Taylor model (remainders 1.6e−23, 7.4e−18). Re-run 2026-08-30 with the corrected operator
(§3.8 (iii)); the outward-rounded enclosures are unchanged; `results/theorem_b/krawczyk.json`. S1's κ₁ =
2.8105525487765 (float, ± 1e−9) and KHA99's 2.81055255 lie inside; γ agrees with KHA's 0.35580192 to all
digits given. The two boxes are disjoint, so with the winding number 2 of §3.4 they are *all* the zeros.

### 3.7 The zero at κ̄ is the pure-gauge mode

(i) κ̄ = 2 − A₀ + 2W₀/3 = 0.355699203710964 ± 1.6e−16 from the certified A1 data (V0 box) lies in the Krawczyk
box of §3.6, which contains exactly one zero. (ii) The pure-gauge perturbation g = (A′, N′ + κ̄N, W′, V′)/A′(0)
of the background (the residual coordinate freedom x → x + εe^{κs} of the time-dependent KHA system,
whose lapse picks up κεe^{κs}N; it is admissible in the sonic-point gauge iff N_p(0) = N′(0) + κN(0) = 0, i.e.
iff κ = −N′(0)/N(0) = κ̄) is checked as a ball identity on the certified series (`gauge_checks`, 0.1 s): the
4D linearised residual [Pg′ − G(κ̄)g]_j contains 0 for all j ≤ 38 (widths ≤ 3.2e−9, set by the V0-box
widths of the background coefficients at high order), the linearised constraint residual contains 0 (≤
1.8e−10), and every coefficient ball of the 4D sonic series over a box containing κ̄ overlaps g_n (n ≤ 40;
distance ≤ the ball widths for n ≤ 10). g is analytic at the sonic point and regular at the centre (the
background is analytic in t = e^x there, and Â_p = A′e^{−2x}, ñ_p = (N′ + κ̄N)e^x, w̃_p, ṽ_p are O(1)), so
g ∈ span(r₁, r₂) at x_d and E(κ̄) = 0 exactly, given the gauge identity; by uniqueness (§3.5 (i)) g *is* the
admissible sonic solution at κ̄. Hence the zero in the box is κ̄ and it is the gauge mode, of multiplicity
one (E′ ≠ 0 on the box). What is verified computationally is the ball identity to order 38 plus the
Krawczyk uniqueness; the exact statement "g solves the linearised system for every κ" is the symmetry
argument (S1 verified it numerically to 1e−8; it is not re-derived symbolically here).

### 3.8 What is proven (Theorem B, computer-assisted part), caveats, timings

**Theorem B′ (as certified).** Let R = [0, 15] × [−14, 14] ⊂ ℂ and let E: R → ℂ be the matching function
E(κ) = e^{3x_d} det[r₁(κ), r₂(κ), p̃(x_d; κ)]_{(Â, ñ, ṽ)}, x_d = −3, of the linearised KHA system about the
certified EC background (V0* ∈ c* ± 1e−16, (a*, μ*) the A3 balls), with p̃ the solution analytic at the
sonic point with N_p(0) = 0, A_p(0) = 1 on the linearised constraint surface and r₁, r₂ the regular
(exponent-0) family at the centre (§2.4, §3.3). Then E is analytic on a neighbourhood of R and has
**exactly two zeros in R, both simple and real: the gauge zero κ̄ ∈ [0.35569919907, 0.35569920835], which is
the pure-gauge mode (κ̄ = 2 − A₀ + 2ω₀/3 = 0.355699203710964 ± 1.6e−16), and κ₁ ∈ [2.81055254439,
2.81055255326]** (E′(κ₁) = 0.022490638 ± 4e−10, E′(κ̄) = −0.061293328 ± 1e−9). Hence there is exactly one
non-gauge eigenvalue with Re κ ≥ 0 in R, it is real, and the mass-scaling exponent is
**γ = 1/κ₁ ∈ [0.3558019218, 0.3558019231]** (γ = 0.35580192 ± 3.0e−9; KHA: 0.35580192).
Reality: E(κ̄) = E(κ)̄ (real data), so the unique zero in each box, symmetric about the real axis, is real.
Proof structure: §3.5 (analyticity on a neighbourhood of R and E(κ) = 0 ⟺ eigenvalue) + §3.4 (winding
number 2 = number of zeros with multiplicity) + §3.6 (two disjoint Krawczyk boxes, one simple zero each)
+ §3.7 (the zero near κ̄ is the gauge mode). FORMULATION's targets κ₁ ∈ [2.81055, 2.81056] and
γ ∈ [0.355801, 0.355802] are met with ≈ 3 more digits.

**Function-space caveat (what "eigenvalue" means here).** A zero of E is a κ for which the *4D* linearised
KHA system (spherical symmetry, sonic-point gauge N_p(0) = 0, time dependence e^{κs}) has a solution that is
(i) analytic in x at the sonic point (radius ≥ 0.065), (ii) on the linearised momentum-constraint surface,
(iii) regular at the centre, i.e. analytic in t = e^x with the scaled variables (Â_p, ñ_p, w̃_p, ṽ_p) bounded
(the point-mass direction ∝ t^{−3} excluded), (iv) normalised by A_p(0) = 1 — no loss in R since the
admissible sonic solution is unique up to scale and has A_p(0) ≠ 0 for σ(κ) ≠ 0, §3.5 (i). Kink modes
(non-analytic at the sonic point), non-spherical modes and Re κ < 0 are outside this class, and Theorem C
(no zeros with Re κ ≥ 0 outside R) remains the open analytic item; the count refers to the KHA sonic-point
gauge (the HM01 gauge would move κ̄ but not κ₁).
**Rigor caveats (as revised 2026-08-30).** (1) The gauge identity "g = (A′, N′ + κN, W′, V′)/A′(0) solves the linearised system for
every κ" is used as the symmetry argument and verified as a ball identity to order 38 (§3.7), not as an exact
polynomial identity; the Krawczyk statement (exactly one simple zero in the box containing κ̄) is
independent of it. (2) *Closed.* r₁, r₂, p̃ ∈ Σ is now proven for every κ ∈ R by the exact identity and the two exponent
conditions below (previously: exponent count {−3, 0, 0} of the reduced form, ball-checked). (3) The A1/A2/A3
certificates of Theorem A are inherited (V0 box, centre balls, tube); an independent review of the whole rigor
chain (Stage 1–3) is still required before quoting a theorem. (4) All winding-number arithmetic is exact
dyadic/arb; the per-tile models are Taylor models with rigorous remainders; the sub-segment criterion
Re(E e^{−iθ}) > 0 is a ball comparison.

**Invariance of the constraint hyperplane and the equivalence E = 0 ⟺ eigenvalue (exact, all κ).** Let
c(x) := (κ − A)S A_p − A(C̃_N N_p + C̃_W W_p + C̃_V V_p) along the background (Σ(x) = ker c(x, ·); c = AS × the
linearised E_tr, given row 1). `linsys.linear_constraint_propagation` establishes, in `fmpq_mpoly` exactly as
Lemma 2 of `modelling-scope.md` (u′ and p′ by Cramer, A eliminated by the background constraint A = (S + 2WT̃)/S,
exact division with remainder 0, gcd cancellation, denominator checked to be a constant multiple of S·D):
    16 S D · dc/dx = Λ_lin · c,   D := Δ̃/(4SW) = 3N²V² − N² + 4NV − V² + 3,   Λ_lin ∈ ℚ[N, W, V] (19 terms, κ-free),
for every solution p of the 4D linearised system P p′ = (DQ − Ψ − κP_s)p over a background solving the 4D system
with C̃ = 0 (both facts hold for the A1 series, the A3 tube and the A2 series, `modelling-scope.md` Lemma 4).
Hence Σ is invariant wherever Δ̃ ≠ 0, in particular along the whole tube (`tube_regular`: P invertible on every
sub-box ⟺ SΔ̃ ≠ 0), and the ODE has the two singular ends handled by exponents (`analyticity.constraint_exponents`,
ball certificate, recorded in `analyticity_R.json`): (a) *sonic point.* D(u₀) = 0 exactly (A1's closed forms,
`test_sonic_closed_forms_exact_in_Q_V0_sqrt3`, with S₀W₀ ≠ 0), D₁ = 7.5606 ± 7e−15 ≠ 0, so c′ = λc has a simple
pole with residue γ_lin = Λ_lin(u₀)/(16 S₀ D₁) = 0 ± 1.3e−15 ∉ ℕ⁺; an analytic c with c(0) = 0 (the order-0
constraint of the 4D recursion) then has (nD₁·16S₀ − Λ_lin(u₀))c_n = (lower orders) and c ≡ 0 on |x| < ν: p̃ ∈ Σ,
transported to x_d by invariance. (b) *centre.* With N = n/t, W = wt², V = vt (t = e^x) one has c = t²c̃ with c̃
analytic for a solution analytic in t, and θc̃ = (λ(t) − 2)c̃, λ = Λ_lin/(16SD): t²D = −n² + O(t²) and Λ_lin has
no N³ monomial (exact), so λ is analytic at t = 0 with λ(0) = 16n₀²/(16·(−n₀²)) = −1 — the constraint exponent
is ρ = −3 (ball −3 ± 4e−15; the off-Σ point-mass direction) ∉ ℕ₀, hence (k − ρ)c̃_k = (lower orders) forces c̃ ≡ 0
for every solution analytic in t: r₁, r₂ ∈ Σ on |t| < 0.06 ∋ e^{−3}. (c) ℓ_Σ·e_W = 2AT̃ = 0.6778 ± 1e−14 ≠ 0 at
x_d: e_W ∉ Σ(x_d). All three are κ-independent. With §3.5 (i) (uniqueness of the admissible sonic solution: the
order-0 2×2 solve M₀ is a ball solve on every box of the cover) and (C) (the exponent-0 family is exactly
span(r₁, r₂)), E(κ) = 0 ⟺ κ is an eigenvalue in the class above, for every κ ∈ R — in particular at κ₁ and κ̄.

**What changed in the review revision (2026-08-30).** (i) The exact identity and exponent certificate above
(new: `linsys.linear_constraint_propagation`, `analyticity.constraint_exponents`, test). (ii) `linsonic4.certify4`:
residual orders up to j = K + max(l_P, l_G) + 1 and Z₁ summed to k = max(l_P, l_G) (previously l_P + K + 1 and
k < max; with l_P = l_G = 42 one term of size ~1e−9 was omitted; all certificates re-run, numbers unchanged).
(iii) Krawczyk: E(κ_c) = c₀ + R(0) includes the model remainder and the Cauchy denominator is w − √2·w_B for
the square box (§3.6 re-run: enclosures unchanged to the digits shown). (iv) `modecount.Context` refuses a cached
tube whose `info` (V0, box width, x₀) differs from the certified constants. (v) The full-rectangle winding and
the analyticity cover are re-runnable from `python -m p4.validated.{modecount,analyticity} TUBE` and their
results (all tile increments as exact dyadics, box parameters, tube info) are stored in `results/theorem_b/` and
re-checked by tests.

**Timings** (M-series laptop, 256 bits). 4D sonic series + certificate 0.05–0.1 s; δ-model (m = 8, box run)
0.15 s; centre model 0.4 s; fundamental-matrix Taylor model over the 244 tube steps 38 s at m = 8 (0.16
s/step, of which the block recursion ≈ 40 %, the 10 Lohner sets ≈ 30 %), 9–20 s at m = 3; one contour tile
≈ 40 s (+ 15 s per worker for the κ-free step data); the whole contour (§3.4) in ≈ 20 min on 8 workers (1231 s); the
interior cover 79 s (2 workers); Krawczyk 39 s (κ₁), 19 s (κ̄) including the Newton re-centring; fast
tests 45 s, slow tests ≈ 5 min with the cached tube; the exact identity of §3.8 0.01 s, its exponent certificate 1 s.

Reproduction (tube from §2.5's `Tube.build` saved with `tube.save(path)`, ≈ 190 s; results under
`problems/P4/results/theorem_b/`, re-checked by the tests):
```
PYTHONPATH=problems/P4/src uv run python -c "
from p4.validated import lintube, modecount; lintube.Tube.build(modecount.V0_EC, modecount.W_V0, x_c=-3.0).save('tube_full.json')"
PYTHONPATH=problems/P4/src uv run python -m p4.validated.modecount tube_full.json --workers 8      # winding_R.json, ~20 min
PYTHONPATH=problems/P4/src uv run python -m p4.validated.analyticity tube_full.json --workers 2    # analyticity_R.json, ~80 s
PYTHONPATH=problems/P4/src uv run python -c "
from p4.validated import modecount, krawczyk_kappa, lintube; from flint import arb
ctx = modecount.Context(lintube.Tube.load('tube_full.json'))                                        # krawczyk.json numbers
for kc, w in ((modecount.KAPPA1, 2e-7), (modecount.KGAUGE, 2e-5)):
    k = krawczyk_kappa.krawczyk(ctx, arb(kc), w); print(k['ok'], k['zero'], k['dE'], k['dE_box'])
print(krawczyk_kappa.gauge_checks(ctx.bg))"
P4_TUBE_CACHE=tube_full.json uv run pytest problems/P4/tests -q -k stage                            # Stage 1-3, slow tests included
```
