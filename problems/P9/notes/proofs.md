# P9(a) — Proofs backing the certified H₀ bound

Status: 2026-08-29. This note states and proves the lemmas on which the
certificate chain of FORMULATION.md §3.6–3.7 rests, and matches each lemma to
the code that implements it (`file:function`, paths relative to
`problems/P9/src/p9/`). It is a companion to FORMULATION.md (the statement) and
`results/summary.md` (the numbers); nothing here changes either.

Conventions. Every number stored as an IEEE-754 double (grid nodes, data values,
whitening factors, L, box endpoints, brackets, T) is read as the exact dyadic
rational it represents. Expressions such as e^{x_k}, log₁₀(e^{x_k} − 1) are then
real numbers, enclosed in Arb balls on the verification path. fl↑(·), fl↓(·)
denote a double that is ≥, resp. ≤, the real argument (`verify._endpoint`).

## 1. Setting and notation

**Grid.** G = G₃₇₃ is the midpoint-refined (r = 2) geometric grid
0 = x₀ < x₁ < … < x_N = X = ln 3.5, N = 372, h_k := x_{k+1} − x_k ≤ 0.00495
(`geometry.geometric_grid`, `model.ClassSpec.x`). Base grid: x₁ = ln 1.005,
x_{k+1} = 1.1 x_k while the spacing is < 0.02, then spacing 0.02 (94 nodes);
r = 2 inserts midpoints twice. All nodes are doubles.

**Box.** ũ_lo := fl↓(c/(H_max r_hi)), ũ_hi := fl↑(c/(H_min r_lo)) with
c = 299792.458 km s⁻¹, H_min = 20, H_max = 2000 km s⁻¹ Mpc⁻¹, r_lo = fl(146.57),
r_hi = fl(147.61) Mpc (`model.ClassSpec.u_box`; outward rounding so that the
float box contains {c/(r_d H) : H ∈ [H_min, H_max], r_d ∈ [r_lo, r_hi]}).
Note fl(146.57) = 146.57 − 3/439804651110400 < 146.57.

**Class.** C(G, L) is the set of continuous ũ: [0, X] → (0, ∞), linear on each
segment [x_k, x_{k+1}], with node values ũ_k := ũ(x_k) satisfying

    ũ_lo ≤ ũ_k ≤ ũ_hi                                            (box)
    |ũ_{k+1} − ũ_k| ≤ L h_k min(ũ_k, ũ_{k+1})                     (min-form slope bound)

for all k. A member is identified with its node vector ũ ∈ ℝ^{N+1}. Physically
ũ(x) = c/(r_d H(z)), x = ln(1+z).

**Observables.** For x ∈ [0, X]:

    D̃_H(x) := ũ(x),      D̃_M(x) := ∫₀ˣ ũ(x′) e^{x′} dx′,      D̃_i := D̃_M(x_i)  (D̃₀ = 0).

**Data.** BAO rows r = 1..n_b (n_b = 12): redshift z_r, value d_r, kind ∈
{D_M, D_H}; x_r := ln(1 + z_r); P_r(ũ) := D̃_M(x_r) for kind D_M, D̃_H(x_r) for
kind D_H. SN rows j = 1..n (n = 1580): m_j, zHD_j, zHEL_j; x_j := ln(1 + zHD_j).
The redshifts are doubles; the x's are transcendental (§13.2).

**Covariance convention (C̃).** W_b (12×12) and W_s (n×n) are the float64
lower-triangular inverse Cholesky factors of the released C_BAO and C_SN
(`model.whitener`); they are the recorded objects and define
C̃ := (WᵀW)⁻¹ for each block. Recorded deviations max|W C Wᵀ − I|:
6.2×10⁻¹⁵ (SN), 2.2×10⁻¹⁶ (BAO) (`data/MANIFEST.json`, "whitening").

**Statistic.** For ũ ∈ C(G, L), M′ ∈ ℝ:

    r_j(ũ, M′) := m_j − 5 log₁₀ D̃_M(x_j) − 5 log₁₀(1 + zHEL_j) − M′,
    χ²(ũ, M′)  := ‖W_b (d − P(ũ))‖² + ‖W_s r(ũ, M′)‖²           (`model.Frozen.chi2`, exact: `verify.rigorous_chi2`).

**Feasible set and target.** For a recorded double T,

    F = F(G, L, T) := { (ũ, M′) : ũ ∈ C(G, L), M′ ∈ [0, 40], χ²(ũ, M′) ≤ T },
    λ₀^min := inf { log₁₀ ũ₀ : (ũ, M′) ∈ F }.

Since H₀ = c/(r_d ũ₀) and F does not involve r_d (M′ absorbs 5 log₁₀ r_d; BAO is
in units of r_d), for every (ũ, M′) ∈ F and r_d ≥ r_lo: H₀ ≤ c/(r_lo · 10^{λ₀^min}).
T is fixed from a reference point by Lemma 8: T = fl↑(χ²(ũ*, M′*) + Δ).

**Relaxation variables (FORMULATION §3.6, `lkr_rows.layout_for`).**
c_i := log₁₀(e^{x_i} − 1) for i ≥ 1, c₀ := 0. For a class member define

    λ_i := log₁₀ ũ_i,                          κ₀ := λ₀,   κ_i := log₁₀[D̃_i/(e^{x_i} − 1)]  (i ≥ 1),
    ℓ_j := log₁₀ D̃_M(x_j),                     ρ_i := λ_i − κ_i − c_i = log₁₀(ũ_i/D̃_i)     (i ≥ 1),
    δ_i := κ_{i+1} − κ_i + c_{i+1} − c_i,     s_i := λ_{i+1} − κ_i − c_i  (i ≥ 1),  s₀ := λ₁ − λ₀,
    E_{d,i} := 10^{δ_i},  E_{r,i} := 10^{ρ_i},  E_{s,i} := 10^{s_i},
    y_b := log₁₀ D̃_M(x_r)  (D_M rows),   P_r := P_r(ũ),   U_k := ũ_k (k ∈ E, the D_H segment endpoints),
    w_b := W_b (d − P),   w_s := W_s r.

So δ_i = log₁₀(D̃_{i+1}/D̃_i) for i ≥ 1 and δ₀ = log₁₀(D̃₁/ũ₀); s_i = ρ_i + (λ_{i+1} − λ_i).
The conic program has 5065 variables, 1971 equalities, 12571 inequalities and one
second-order cone of dimension 1593 (for G₃₇₃, n = 1580).

**Brackets.** A bracket vector B = (ρ_lo, ρ_hi, y_lo, y_hi, λ_lo, λ_hi)
(`lkr.Brackets3`) is *valid* if every (ũ, M′) ∈ F satisfies
ρ_i ∈ [ρ_lo,i, ρ_hi,i], y_b ∈ [y_lo,b, y_hi,b], λ_i ∈ [λ_lo,i, λ_hi,i]. F_rel(B) denotes
the feasible set of the relaxed conic program built from B and T by
`lkr_rows.build` (rows listed in §§3–7).

## 2. Theorem and corollary

**Theorem (certified rows; CLAIMS P9-2, P9-3).** Fix the frozen inputs, G = G₃₇₃,
M′ ∈ [0, 40], and a row (L, Δ) with its recorded reference point (ũ*, M′*) and
T = fl↑(χ²(ũ*, M′*) + Δ) (`results/certificates/lkr_<row>/state.json`). Then for
every (ũ, M′) ∈ F(G, L, T): log₁₀ ũ₀ ≥ λ₀^lb, hence for every r_d ≥ 146.57 Mpc,
H₀ = c/(r_d ũ₀) ≤ H_up, with

| L | Δ | T | λ₀^lb | H_up (km s⁻¹ Mpc⁻¹, rounded up) | passes |
|---|---|---|---|---|---|
| 1   | 4 | 1600.6025273225905 | 1.466652998279891  | 69.8430 | 6 |
| 1.5 | 1 | 1380.7251890839852 | 1.4723753834749773 | 68.9288 | 5 |
| 1.5 | 4 | 1383.7251890839852 | 1.4666604946265958 | 69.8418 | 5 |
| 1.5 | 9 | 1388.7251890839852 | 1.4607606987404465 | 70.7971 | 5 |
| 2   | 4 | 1381.149522938679  | 1.4652332794512115 | 70.0717 | 5 |
| 3   | 4 | 1376.8754806108354 | 1.461857760931861  | 70.6184 | 6 |
| 5   | 4 | 1370.248688985442  | 1.4536056769402457 | 71.9731 | 6 |
| 10  | 4 | 1358.456480933817  | 1.4283891947449474 | 76.2758 | 7 |

H_up = ceil₄(c/(fl(146.57) · 10^{λ₀^lb})) evaluated in ball arithmetic (e.g.
69.8417512308 ± 4×10⁻¹¹ for L = 1.5, Δ = 4); using fl(146.57) < 146.57 only
enlarges the bound, so the statement holds for r_d ≥ 146.57 exactly.

*Proof.* Lemma 8 gives (ũ*, M′*) ∈ F and F ⊇ (Δ-sublevel set of χ² on the class).
Lemma 7 (induction over the tightening passes, base = class box) shows that
every bracket vector used in the chain is valid, so by the Proposition of §8 the
lift of F lies in F_rel(B_p) for every pass p, and by Lemma 6 each ball-verified
dual certificate is a lower bound on the corresponding objective over F_rel(B_p),
hence over F. The final certificate (`p<last>_final_lambda0.npz`) has objective
λ₀, giving λ₀ ≥ λ₀^lb on F. Since ũ₀ = 10^{λ₀} and r_d ≥ r_lo, H₀ ≤ c/(r_lo 10^{λ₀^lb}). ∎

**Corollary (two-sided interval; CLAIMS P9-1r, P9-3).** For each row, the point
(ũ_f, M′_f) of `results/certificates/feasible_<row>_r2.json` lies in F (§12), so

    H_lo ≤ sup_{(ũ,M′)∈F} c/(r_lo ũ₀) ≤ H_up,

with H_lo = floor₄(c/(r_lo ũ_{f,0})): 69.6497 (L=1), 68.3373 (L=1.5, Δ=1),
69.5227 (L=1.5, Δ=4), 70.5562 (L=1.5, Δ=9), 69.5876 (L=2), 69.7561 (L=3),
70.2001 (L=5), 71.6439 (L=10). (Two of the lower bounds in `results/summary.md`
and CLAIMS are printed rounded to nearest — 69.5228 and 69.5877 — rather than
down; the exact values are 69.52275852… and 69.58766072…, so the rounded-down
figures above are the ones that follow from the certificates. The upper bounds
there are rounded up, as required.)

**What is claimed / not claimed** (mirrors FORMULATION §2). Claimed: a theorem
about the finite-dimensional set F defined above from the frozen inputs, the
recorded W (i.e. C̃), the recorded T, and the float grid/box/L. Not claimed:
frequentist coverage of Δ; anything about r_d beyond the box (r_d enters only
through H₀ = c/(r_d ũ₀)); anything about the SN absolute calibration (M′ free,
calibrators excluded); anything about the continuous class Lip_L(ln H) beyond
C(G, L) ⊂ Lip_L (Lemma 1(b); §13.5); optimality of the reference point (§11).
Assumed: flatness; D_L = (1 + zHEL) D_M(zHD); the released covariances (through W)
as the definition of the likelihood.

## 3. Lemma 1 — the class in λ = log₁₀ ũ

Notation: lh_k := log₁₀(1 + L h_k); for x ∈ (0, X],

    q(x) := e^x/(e^x − 1),
    I_±(x) := ∫₀ˣ e^{±L(x−x′)} e^{x′} dx′ = e^{±Lx} (e^{(1∓L)x} − 1)/(1 ∓ L)     (= x e^{Lx} for I_+ when L = 1),
    θ(x) := ũ(x)(e^x − 1)/D̃_M(x),    θ_lo(x) := (e^x − 1)/I_+(x),    θ_hi(x) := (e^x − 1)/I_−(x),
    J_±(x_i, x_{i+1}) := ∫_{x_i}^{x_{i+1}} e^{±L(x′−x_i)} e^{x′} dx′ = e^{x_i}(e^{(1±L)h_i} − 1)/(1 ± L)   (= e^{x_i} h_i if 1 ± L = 0).

Code: `geometry._I_pm`, `geometry._J_pm`, `geometry.theta_bounds` (floats) and
the same three functions in `lkr_rows` (parameterized arithmetic). The
degenerate-denominator branch (`abs(mid(a)) < 1e-9`) is exact for the L values
used: 1 ∓ L is either exactly 0 (L = 1) or ≥ 0.5 in modulus.

**Lemma 1.** Let ũ be continuous, linear on each segment of G, with positive node
values, and λ_i := log₁₀ ũ_i.

(a) *(Equivalence.)* For each k, |ũ_{k+1} − ũ_k| ≤ L h_k min(ũ_k, ũ_{k+1}) iff
|λ_{k+1} − λ_k| ≤ lh_k. Hence C(G, L) = {ũ : box, |λ_{k+1} − λ_k| ≤ lh_k ∀k}.

(b) *(Lipschitz.)* If ũ ∈ C(G, L), then |(ln ũ)′| ≤ L on every open segment and
ln ũ is L-Lipschitz on [0, X]; equivalently |d ln H/d ln(1+z)| ≤ L.

(c) *(θ-range.)* If ũ ∈ C(G, L) and x ∈ (0, X], then D̃_M(x) ∈ ũ(x)·[I_−(x), I_+(x)]
and θ(x) ∈ [θ_lo(x), θ_hi(x)]. Moreover θ_lo is nonincreasing and θ_hi nondecreasing
on (0, X], both tend to 1 as x → 0⁺, so θ_lo ≤ 1 ≤ θ_hi. In the variable t = e^x:
θ_hi = (1+L) t^L (t−1)/(t^{1+L} − 1), θ_lo = (1−L)(t−1)/(t − t^L) (L ≠ 1).

(d) *(Ratio bounds.)* If ũ ∈ C(G, L): D̃_{i+1} − D̃_i ∈ ũ_i·[J_−, J_+]; for i ≥ 1,
D̃_{i+1}/D̃_i ∈ [1 + J_−/I_+(x_i), 1 + J_+/I_−(x_i)] =: [r_lo,i, r_hi,i]; for i = 0,
D̃₁/ũ₀ ∈ [J_−, J_+] =: [r_lo,0, r_hi,0]. Hence δ_i ∈ [log₁₀ r_lo,i, log₁₀ r_hi,i].

*Proof.* (a) Put r := ũ_{k+1}/ũ_k > 0 and h := h_k. If r ≥ 1 the minimum is ũ_k and
the condition reads r − 1 ≤ Lh, i.e. r ≤ 1 + Lh. If r < 1 the minimum is ũ_{k+1}
and it reads 1 − r ≤ Lhr, i.e. r ≥ (1 + Lh)⁻¹. Together: (1+Lh)⁻¹ ≤ r ≤ 1 + Lh,
i.e. |log₁₀ r| ≤ log₁₀(1 + Lh) = lh_k.

(b) On (x_k, x_{k+1}), ũ′ ≡ (ũ_{k+1} − ũ_k)/h_k and, ũ being linear,
ũ(x) ≥ min(ũ_k, ũ_{k+1}); so |ũ′/ũ| ≤ |ũ_{k+1} − ũ_k|/(h_k min(ũ_k, ũ_{k+1})) ≤ L.
ln ũ is continuous and piecewise C¹ with |derivative| ≤ L on finitely many
segments, hence L-Lipschitz on [0, X]. ln H = const − ln ũ.

(c) By (b), for 0 ≤ x′ ≤ x: ũ(x) e^{−L(x−x′)} ≤ ũ(x′) ≤ ũ(x) e^{L(x−x′)}. Multiply by
e^{x′} > 0 and integrate over [0, x]: ũ(x) I_−(x) ≤ D̃_M(x) ≤ ũ(x) I_+(x). Dividing
ũ(x)(e^x − 1) by these gives θ ∈ [θ_lo, θ_hi]. The closed forms follow from
I_±(x) = e^{±Lx} ∫₀ˣ e^{(1∓L)x′} dx′.
Monotonicity: substituting s = x − x′, I_−(x) = e^x ∫₀ˣ e^{−(1+L)s} ds,
I_+(x) = e^x ∫₀ˣ e^{(L−1)s} ds, and e^x − 1 = e^x ∫₀ˣ e^{−s} ds, so

    θ_hi(x) = ∫₀ˣ e^{−s} ds / ∫₀ˣ e^{−(1+L)s} ds,     θ_lo(x) = ∫₀ˣ e^{−s} ds / ∫₀ˣ e^{(L−1)s} ds.

Claim: if f, g > 0 are continuous on [0, X] and f/g is nondecreasing
(nonincreasing), then R(x) := ∫₀ˣ f / ∫₀ˣ g is nondecreasing (nonincreasing) on
(0, X] and R(x) → f(0)/g(0) as x → 0⁺. Indeed, with F := ∫₀ˣ f, G := ∫₀ˣ g,
R′ = (fG − Fg)/G² and

    f(x)G(x) − F(x)g(x) = ∫₀ˣ [f(x)g(s) − f(s)g(x)] ds = ∫₀ˣ g(s) g(x) [f(x)/g(x) − f(s)/g(s)] ds,

whose sign is that of the monotonicity of f/g; the limit is l'Hôpital. Apply
with f = e^{−s}: for θ_hi, g = e^{−(1+L)s} and f/g = e^{Ls} is increasing; for
θ_lo, g = e^{(L−1)s} and f/g = e^{−Ls} is decreasing. Both limits equal 1. The
t-forms are the closed forms with t = e^x.

(d) By (b) from the node x_i forward: for x′ ∈ [x_i, x_{i+1}],
ũ_i e^{−L(x′−x_i)} ≤ ũ(x′) ≤ ũ_i e^{L(x′−x_i)}; multiply by e^{x′} and integrate:
D̃_{i+1} − D̃_i = ∫_{x_i}^{x_{i+1}} ũ e^{x′} dx′ ∈ ũ_i [J_−, J_+]. For i ≥ 1, (c) at
x = x_i (where ũ(x_i) = ũ_i) gives D̃_i ∈ ũ_i [I_−(x_i), I_+(x_i)], and
D̃_{i+1}/D̃_i = 1 + (D̃_{i+1} − D̃_i)/D̃_i lies in [1 + J_−/I_+, 1 + J_+/I_−] (the scale
ũ_i cancels). For i = 0, D̃₀ = 0. Finally δ_i = log₁₀(D̃_{i+1}/D̃_i) for i ≥ 1 and
δ₀ = log₁₀(D̃₁/ũ₀) by the definitions of §1 (κ₀ = λ₀, c₀ = 0). ∎

Rows (`lkr_rows.build`): class rows ±(λ_{k+1} − λ_k) ≤ lh_k with
lh_k = log₁₀(1 + L h_k) (ball); λ box rows λ_i ∈ [λ_lo,i, λ_hi,i] where the base
bracket is [fl↓(log₁₀ ũ_lo), fl↑(log₁₀ ũ_hi)]; δ bracket rows
−(κ_{i+1} − κ_i) ≤ −(lo_d − dc_i), κ_{i+1} − κ_i ≤ hi_d − dc_i with
lo_d = fl↓(log₁₀ r_lo,i), hi_d = fl↑(log₁₀ r_hi,i), dc_i = c_{i+1} − c_i (ball).
The exact class check of the reference and feasible points uses the min-form in
rational arithmetic (`certify_feasible.in_class_exact`); by (a) this is the same
class as the λ-form used by the relaxation.

## 4. Lemma 2 — interpolation of κ

Let ũ ∈ C(G, L), κ as in §1 and F := (ln 10) κ, i.e. F(x) = ln D̃_M(x) − ln(e^x − 1)
for x ∈ (0, X], F(0) := ln ũ₀. Put a := ũ′/ũ on open segments, and

    B₀ := L/6 + L² (1/2 + x₁/6)²,
    B_k := q(x_k) [ L θ_hi(x_{k+1}) + dθ_{k+1} ( q(x_k)(1 + θ_hi(x_{k+1})) − 1 ) ]    (k ≥ 1),
    dθ_{k+1} := max( 1 − θ_lo(x_{k+1}), θ_hi(x_{k+1}) − 1 ),
    e_k := h_k² B_k / (8 ln 10).

Code: `geometry.kappa_second_derivative_bound` (B_k), `geometry.kappa_interp_slack`
(e_k), `lkr_rows.curvature_bound` (B_k in the verifier's arithmetic); `build`
stores e_seg[k] = fl↑(e_k).

**Lemma 2.** (a) F ∈ C¹((0, X]) and on each open segment,

    F′ = q(θ − 1),    F″ = q [ θa − (θ − 1)(q(1 + θ) − 1) ],    with q′ = −q(q − 1),  θ′/θ = a + q − qθ.

(b) For k ≥ 1, |F″| ≤ B_k on (x_k, x_{k+1}).

(c) On [0, x₁]: F(x) = ln(ũ₀ + s φ(x)) with s := (ũ₁ − ũ₀)/h₀ and
φ(x) := ((x − 1)e^x + 1)/(e^x − 1) = x q(x) − 1, φ(0) := 0. For all x ≥ 0,
0 ≤ φ″(x) ≤ 1/6 and 1/2 ≤ φ′(x) ≤ 1/2 + x/6; and |s|/(ũ₀ + sφ(x)) ≤ L on [0, x₁].
Consequently |F″| ≤ B₀ on (0, x₁).

(d) *(C^{1,1} interpolation inequality.)* If f ∈ C¹[a, b] and f′ is Lipschitz with
constant B, and p is the linear interpolant of f at a, b, then
|f(x) − p(x)| ≤ B (x − a)(b − x)/2 ≤ B (b − a)²/8 on [a, b].

(e) For every k and x ∈ [x_k, x_{k+1}], with t := (x − x_k)/h_k and κ_k := κ(x_k):
|κ(x) − (1 − t) κ_k − t κ_{k+1}| ≤ e_k. In particular, for SN j on segment k = k(j),
|ℓ_j − (1 − t_j) κ_k − t_j κ_{k+1} − log₁₀ zHD_j| ≤ e_k, and for a D_M BAO row on
segment k, |y_b − (1 − t) κ_k − t κ_{k+1} − log₁₀ z_r| ≤ e_k.

*Proof.* (a) D̃_M ∈ C¹[0, X] with D̃_M′ = ũ e^x (ũ is continuous) and D̃_M > 0 on
(0, X], so F ∈ C¹((0, X]) with F′ = ũe^x/D̃_M − e^x/(e^x − 1) = qθ − q, using
ũe^x/D̃_M = θ e^x/(e^x − 1) = qθ. On an open segment ũ is linear, so θ, F are
smooth there. From q = 1 + (e^x − 1)⁻¹, q′ = −e^x/(e^x − 1)² = −q(q − 1). From
ln θ = ln ũ + ln(e^x − 1) − ln D̃_M, θ′/θ = a + q − qθ. Therefore
F″ = q′(θ − 1) + qθ′ = −q(q−1)(θ−1) + qθ(a + q − qθ)
   = q[θa + qθ(1 − θ) − (q − 1)(θ − 1)] = q[θa − (θ − 1)(qθ + q − 1)] = q[θa − (θ−1)(q(1+θ) − 1)].

(b) On (x_k, x_{k+1}), k ≥ 1: |a| ≤ L by Lemma 1(b); by Lemma 1(c) and its
monotonicity, θ(x) ∈ [θ_lo(x_{k+1}), θ_hi(x_{k+1})], so 0 < θ ≤ θ_hi(x_{k+1}) and
|θ − 1| ≤ dθ_{k+1}; q is decreasing (q′ < 0), so 1 < q(x) ≤ q(x_k), and
q(1+θ) − 1 > 0. Hence |F″| ≤ q[Lθ + |θ − 1|(q(1+θ) − 1)] ≤ B_k, every factor being
nonnegative and nondecreasing in q and in θ.

(c) On [0, x₁], ũ(x′) = ũ₀ + s x′, and ∫₀ˣ x′ e^{x′} dx′ = (x − 1)e^x + 1, so
D̃_M(x) = ũ₀(e^x − 1) + s[(x − 1)e^x + 1], i.e. D̃_M/(e^x − 1) = ũ₀ + sφ(x) and
F = ln(ũ₀ + sφ). Also (x − 1)e^x + 1 = x e^x − (e^x − 1) gives φ = xq − 1.
Bounds on φ: with y := x/2, (e^x − 1)⁻¹ = (coth y − 1)/2, hence
φ = x + x/(e^x − 1) − 1 = y + y coth y − 1. Differentiating in x (d/dx = ½ d/dy):
φ′ = ½[1 + coth y − y/sinh²y], φ″ = ¼[−2/sinh²y + 2y cosh y/sinh³y] = (y cosh y − sinh y)/(2 sinh³y).
The numerator vanishes at 0 and has derivative y sinh y ≥ 0, so φ″ ≥ 0. Put
h(y) := sinh³y − 3(y cosh y − sinh y): h(0) = 0 and h′ = 3 sinh y (sinh y cosh y − y) ≥ 0
because sinh y cosh y = ½ sinh 2y ≥ y; so φ″ ≤ 1/6. As sinh y cosh y − y = O(y³),
φ′(0⁺) = ½, and φ′(x) = ½ + ∫₀ˣ φ″ ∈ [½, ½ + x/6]. φ(0) = 0 since xq → 1.
Ratio bound: g(x) := ũ₀ + sφ(x) = D̃_M(x)/(e^x − 1) = ∫₀ˣ ũ e^{x′}dx′ / ∫₀ˣ e^{x′}dx′ is a
weighted mean of ũ over [0, x], so g ≥ min_{[0,x]} ũ ≥ min(ũ₀, ũ₁) (ũ linear on
[0, x₁]); and |s| = |ũ₁ − ũ₀|/h₀ ≤ L min(ũ₀, ũ₁) by the class condition. Thus
|s|/g ≤ L. Now F′ = sφ′/g and F″ = sφ″/g − (sφ′/g)², so
|F″| ≤ Lφ″ + L²φ′² ≤ L/6 + L²(½ + x₁/6)² = B₀ on (0, x₁).

(d) Let e := f − p. Then e(a) = e(b) = 0 and e′ = f′ − const is Lipschitz with
constant B, hence absolutely continuous with |e″| ≤ B a.e. Two integrations by
parts give the Peano-kernel identity e(x) = −∫_a^b G(x, s) e″(s) ds with
G(x, s) = (s − a)(b − x)/(b − a) for s ≤ x and (x − a)(b − s)/(b − a) for s ≥ x.
Since G ≥ 0 and ∫_a^b G(x, s) ds = (x − a)(b − x)/2, |e(x)| ≤ B(x − a)(b − x)/2 ≤ B(b−a)²/8.

(e) Apply (d) to f = F on [x_k, x_{k+1}]. F is C¹ on the closed segment (for k = 0
including the endpoint 0, by the closed form of (c): φ extends analytically to
0), and F′ is Lipschitz there with constant B_k because |F″| ≤ B_k on the open
segment and F′ is continuous on the closed one. So |F(x) − interp| ≤ h_k²B_k/8;
divide by ln 10. For the SN/BAO forms use ℓ_j = log₁₀ D̃_M(x_j) = κ(x_j) + log₁₀(e^{x_j} − 1)
with e^{x_j} − 1 = zHD_j exactly (and e^{x_r} − 1 = z_r). ∎

Rows (`lkr_rows.build`, "SN interpolation" and "BAO D_M rows"):
ℓ_j − (1−t_j)κ_k − t_jκ_{k+1} ≤ c_j + e_k and −ℓ_j + (1−t_j)κ_k + t_jκ_{k+1} ≤ −c_j + e_k,
with c_j = log₁₀ zHD_j and t_j = (x_j − x_k)/h_k as balls, e_k = fl↑(h_k²B_k/(8 ln 10));
the same pair for y_b with c = log₁₀ z_r. The segment index k(j) is chosen from
the float log1p(zHD_j) (`geometry.segment_index`); (e) requires the *real* x_j to
lie in [x_k, x_{k+1}]. This was verified in 200-bit ball arithmetic for all 1580 SN
and 12 BAO redshifts on G₃₇₃ (smallest distance of any x_j to a node: 6.8×10⁻⁷,
far above the float error of x_j); see §13.2. `geometry.kappa_second_derivative_bound`
refuses x₁ > 0.1; by (c) the bound B₀ needs no such restriction (the guard is
merely conservative). Monte-Carlo test of (e) over adversarial class members:
`tests/test_kappa.py::test_interpolation_slack_is_rigorous`.

## 5. Lemma 3 — segment identities

For segment i put h := h_i and

    a_i := (e^{x_{i+1}} − (h + 1) e^{x_i}) / h,      b_i := ((h − 1) e^{x_{i+1}} + e^{x_i}) / h.

Code: `lkr_rows.build` (lists `a`, `b`), `socp2.full_segment_coeffs`, and the
partial-segment version in `geometry.dm_matrix` / `verify.rigorous_chi2` (`dm`).

**Lemma 3.** Let ũ be continuous and linear on each segment of G.
(a) D̃_{i+1} − D̃_i = a_i ũ_i + b_i ũ_{i+1} for every i.
(b) a_i > 0, b_i > 0 and a_i + b_i = e^{x_{i+1}} − e^{x_i}.
(c) *(Scale-free form.)* For i ≥ 1, 10^{δ_i} − 1 = a_i 10^{ρ_i} + b_i 10^{s_i}, i.e.
E_{d,i} − a_i E_{r,i} − b_i E_{s,i} = 1. For i = 0, 10^{δ₀} = a₀ + b₀ 10^{s₀}, i.e.
E_{d,0} − b₀ E_{s,0} = a₀.
(d) For x ∈ [x_k, x_{k+1}],
D̃_M(x) = Σ_{i<k} (a_i ũ_i + b_i ũ_{i+1}) + ũ_k [(x_{k+1} − x + 1)e^x − (h_k + 1)e^{x_k}]/h_k + ũ_{k+1} [(x − x_k − 1)e^x + e^{x_k}]/h_k,
with nonnegative coefficients.

*Proof.* (a) On [x_i, x_{i+1}], ũ(x′) = ũ_i (x_{i+1} − x′)/h + ũ_{i+1} (x′ − x_i)/h.
Integrating by parts,
∫_{x_i}^{x_{i+1}} (x_{i+1} − x′) e^{x′} dx′ = [(x_{i+1} − x′)e^{x′}] + ∫ e^{x′} = −h e^{x_i} + e^{x_{i+1}} − e^{x_i} = h a_i,
∫_{x_i}^{x_{i+1}} (x′ − x_i) e^{x′} dx′ = [(x′ − x_i)e^{x′}] − ∫ e^{x′} = h e^{x_{i+1}} − e^{x_{i+1}} + e^{x_i} = h b_i.
(b) a_i and b_i are integrals of the positive functions (x_{i+1} − x′)e^{x′}/h,
(x′ − x_i)e^{x′}/h over an interval of positive length; their sum is
∫ e^{x′} dx′ = e^{x_{i+1}} − e^{x_i}. (Equivalently h a_i = e^{x_i}(e^h − 1 − h) and
h b_i = e^{x_i} ψ(h) with ψ(h) = (h − 1)e^h + 1, ψ(0) = 0, ψ′ = h e^h.)
(c) For i ≥ 1 divide (a) by D̃_i > 0: D̃_{i+1}/D̃_i − 1 = a_i ũ_i/D̃_i + b_i ũ_{i+1}/D̃_i.
Here D̃_{i+1}/D̃_i = 10^{δ_i}, ũ_i/D̃_i = 10^{ρ_i}, and ũ_{i+1}/D̃_i = 10^{λ_{i+1} − κ_i − c_i} = 10^{s_i}
(since log₁₀ D̃_i = κ_i + c_i). For i = 0, D̃₀ = 0 and D̃₁/ũ₀ = a₀ + b₀ ũ₁/ũ₀ =
a₀ + b₀ 10^{λ₁ − λ₀}, while D̃₁/ũ₀ = 10^{δ₀} by definition.
(d) Same computation with upper limit x instead of x_{k+1}; the two partial
coefficients are integrals of nonnegative functions. ∎

Rows (`lkr_rows.build`, "segment identities"): E_{d,0} − b₀ E_{s,0} = a₀ and
E_{d,i} − a_i E_{r,i} − b_i E_{s,i} = 1 (i = 1..N−1), coefficients as balls.
(Segment N−1 is the last one; there is no identity for i = N.)

## 6. Lemma 4 — sandwiches for convex 10^y

Fix lo < hi, K = 4 (`lkr_rows.N_TANGENTS`), tangent points
y_t^{(k)} := lo + k (hi − lo)/(K + 1), p_k := 10^{y_t^{(k)}}, and the chord slope
σ := (10^{hi} − 10^{lo})/(hi − lo).

**Lemma 4.** (a) For all real y and each k: 10^y ≥ p_k (1 + (y − y_t^{(k)}) ln 10).
(b) For y ∈ [lo, hi]: 10^y ≤ 10^{lo} + σ (y − lo).
(c) *(Range.)* If a pair (E, y) satisfies the K inequalities E ≥ p_k(1 + (y − y_t^{(k)}) ln 10),
the inequality E ≤ 10^{lo} + σ(y − lo), and y ∈ [lo, hi], then
E ∈ [max_k p_k (1 + (lo − y_t^{(k)}) ln 10), 10^{hi}].

*Proof.* (a) e^{u} ≥ e^{v}(1 + u − v) for all u, v (convexity of exp, or
e^{u−v} ≥ 1 + (u − v)); take u = y ln 10, v = y_t ln 10. (b) The chord of a convex
function lies above its graph between the abscissae. (c) Each tangent value is
increasing in y (slope p_k ln 10 > 0), so at y ≥ lo it is ≥ its value at lo; the
chord is increasing and equals 10^{hi} at y = hi. ∎

Remark. The docstring of `lkr_rows.sandwich_range` calls the lower end
positive; that holds iff (y_t^{(1)} − lo) ln 10 < 1 and fails for wide brackets
(e.g. the base λ box, width 2.003), where the lower end is negative. Only
|E| ≤ max(|E_lo|, |E_hi|) is used (Lemma 6), so the sign is immaterial.

Rows (`lkr_rows.sandwich_rows(iE, lin, lo, hi)` produces, for y = Σ lin·vars,
the K tangent rows −E + p_k ln 10 · y ≤ p_k (y_t^{(k)} ln 10 − 1) and the chord row
E − σ y ≤ 10^{lo} − σ lo; the caller adds the bracket rows y ∈ [lo, hi] and
rescales the E coefficient). The three scaled families in `build`:

- E_{d,i} = 10^{δ_i} = 10^{dc_i} · 10^{y′}, y′ := κ_{i+1} − κ_i, dc_i := c_{i+1} − c_i,
  y′ ∈ [lo_d − dc_i, hi_d − dc_i] (Lemma 1(d) bracket rows). The E coefficient is
  divided by 10^{dc_i}, i.e. the rows sandwich E_{d,i}/10^{dc_i} = 10^{y′}.
  Range: E_{d,i} ∈ 10^{dc_i}·[E_lo, E_hi].
- E_{r,i} = 10^{ρ_i} (i ≥ 1): y := λ_i − κ_i = ρ_i + c_i ∈ [ρ_lo,i + c_i, ρ_hi,i + c_i]
  (ρ bracket rows −(λ_i − κ_i) ≤ −(ρ_lo,i + c_i), λ_i − κ_i ≤ ρ_hi,i + c_i). Since
  10^y = 10^{c_i} E_{r,i}, the E coefficient is multiplied by 10^{c_i}.
  Range: E_{r,i} ∈ 10^{−c_i}·[E_lo, E_hi].
- E_{s,i} = 10^{s_i}: for i ≥ 1, y := λ_{i+1} − κ_i = s_i + c_i with
  s_i = ρ_i + (λ_{i+1} − λ_i) ∈ [ρ_lo,i − lh_i, ρ_hi,i + lh_i] =: [lo_s, hi_s] by
  Lemma 1(a); for i = 0, y := λ₁ − λ₀ = s₀ ∈ [−lh₀, lh₀], shift 0. Bracket rows on y
  are added; the E coefficient is multiplied by 10^{c_i} (10^0 for i = 0).
  Range: E_{s,i} ∈ 10^{−c_i}·[E_lo, E_hi].

All tangent/chord rows hold for the lifted variables of any class member whose
δ, ρ, s lie in the stated brackets (Lemma 4(a),(b)); the brackets themselves are
class-only for δ (Lemma 1(d)), valid-by-induction for ρ (Lemma 7), and derived
for s.

## 7. Lemma 5 — BAO rows

Let ũ ∈ C(G, L), and for BAO row r let k = k(r) be the segment with
x_r ∈ [x_k, x_{k+1}] and t_r := (x_r − x_k)/h_k ∈ [0, 1] (`geometry.segment_index`
on the float log1p(z_r); membership of the real x_r verified, §13.2).

**Lemma 5.** (a) *(D_M rows.)* With y_b := log₁₀ D̃_M(x_r):
|y_b − (1 − t_r) κ_k − t_r κ_{k+1} − log₁₀ z_r| ≤ e_k, and P_r = 10^{y_b}. Given a valid
bracket y_b ∈ [y_lo, y_hi], the pair (P_r, y_b) satisfies the sandwich rows of
Lemma 4 on [y_lo, y_hi] and P_r ∈ [E_lo, E_hi] of Lemma 4(c).
(b) *(D_H rows.)* D̃_H(x_r) = ũ(x_r) = (1 − t_r) ũ_k + t_r ũ_{k+1} = (1 − t_r) U_k + t_r U_{k+1}
exactly, with U_k := 10^{λ_k}. Given λ_k ∈ [λ_lo,k, λ_hi,k], (U_k, λ_k) satisfies the
sandwich rows on that bracket, U_k ∈ [E_lo, E_hi](λ_lo,k, λ_hi,k), and
P_r ∈ [min(U_lo,k, U_lo,k+1), max(U_hi,k, U_hi,k+1)].
(c) *(Whitened residuals.)* With off_j := 5 log₁₀(1 + zHEL_j):
w_b + W_b P = W_b d ⇔ w_b = W_b (d − P), and
w_s + 5 W_s ℓ + (W_s 𝟙) M′ = W_s (m − off) ⇔ w_s = W_s r with r_j = m_j − 5ℓ_j − off_j − M′;
then χ²(ũ, M′) = ‖w_b‖² + ‖w_s‖², so χ² ≤ T ⇔ (√T, w_b, w_s) ∈ SOC.

*Proof.* (a) Lemma 2(e) at x = x_r, with e^{x_r} − 1 = z_r; P_r = D̃_M(x_r) = 10^{y_b};
Lemma 4. (b) ũ is linear on [x_k, x_{k+1}]; Lemma 4 for (U_k, λ_k); a convex
combination of U_k, U_{k+1} lies between the smallest lower and largest upper
bound. (c) Definitions of §1: 5 log₁₀ D̃_L = 5ℓ_j + off_j. ∎

Rows (`lkr_rows.build`): "BAO: w_b + Wb P = Wb d" (n_b equalities), the SN dense
block (assembled in `lkr2.LKRModel2.__init__` for the solver and in
`verify3.Verifier3.certify` for the verifier, from `B.sn_rhs` = m − off), the
"D_H rows" P_r − (1 − t_r)U_k − t_r U_{k+1} = 0, the "BAO D_M rows" (interpolation
pair, y_b bracket pair, sandwich of P_r on [y_lo, y_hi]), and the "U nodes"
sandwiches on [λ_lo,k, λ_hi,k] for k ∈ E = {201, 202, 226, 227, 251, 252, 288, 289,
302, 303, 361, 362} (the D_H segment endpoints on G₃₇₃, `layout_for`). The
second-order cone is (√T, w_b, w_s) with √T computed as a ball from the double T.

## 8. Proposition — validity of the relaxation (the lift)

**Proposition.** Let B be a valid bracket vector (§1) and T the recorded
threshold. For (ũ, M′) ∈ F(G, L, T) define the lift y = lift(ũ, M′) ∈ ℝ^{5065} by

    λ_i = log₁₀ ũ_i;  κ₀ = λ₀, κ_i = log₁₀[D̃_i/(e^{x_i} − 1)] (i ≥ 1);  M′;  ℓ_j = log₁₀ D̃_M(x_j);
    E_{d,i} = 10^{δ_i};  E_{r,i} = 10^{ρ_i} (i ≥ 1), E_{r,0} := 0;  E_{s,i} = 10^{s_i};
    y_b = log₁₀ D̃_M(x_r) (D_M rows);  P_r = P_r(ũ);  U_k = ũ_k (k ∈ E);
    w_b = W_b (d − P(ũ));  w_s = W_s r(ũ, M′).

Then y satisfies every row of the relaxed conic program `lkr_rows.build(fr, B, T)`
(with the real, unrounded coefficients), i.e. lift(F) ⊆ F_rel(B). The objective
and the tightened quantities are linear in y: λ₀ = e_{λ₀}·y, ρ_i + c_i = λ_i − κ_i,
y_b, λ_i.

*Proof.* Row groups in build order.
Equalities: (E1) w_b + W_b P = W_b d and (E2) w_s + 5W_sℓ + (W_s𝟙)M′ = W_s(m − off):
Lemma 5(c). (E3) κ₀ − λ₀ = 0: definition. (E4) segment identities: Lemma 3(c).
(E5) D_H rows: Lemma 5(b).
Inequalities: (I1) λ_i ∈ [λ_lo,i, λ_hi,i]: `build` uses the larger of fl↓(log₁₀ ũ_lo)
and B.λ_lo,i (resp. the smaller of fl↑(log₁₀ ũ_hi), B.λ_hi,i); both are bounds on
λ_i, by the class box and the validity of B. (I2) class rows: Lemma 1(a).
(I3) SN interpolation pair: Lemma 2(e) (e_k rounded up only weakens the row).
(I4) δ bracket rows: Lemma 1(d) (endpoints rounded outward). (I5) E_d sandwich:
Lemma 4(a),(b) applied to y′ = κ_{i+1} − κ_i ∈ [lo_d − dc_i, hi_d − dc_i], with
E_{d,i}/10^{dc_i} = 10^{y′}. (I6) ρ bracket rows (i ≥ 1): validity of B. (I7) E_r
sandwich: Lemma 4 with 10^{λ_i − κ_i} = 10^{c_i} E_{r,i}. (I8) s bracket rows:
s_i = ρ_i + (λ_{i+1} − λ_i) with (I6) and (I2); for i = 0, s₀ = λ₁ − λ₀ with (I2).
(I9) E_s sandwich: Lemma 4 with 10^{λ_{i+1} − κ_i} = 10^{c_i} E_{s,i}. (I10) BAO D_M
rows: interpolation pair by Lemma 5(a), bracket by validity of B, sandwich by
Lemma 4. (I11) U sandwiches: Lemma 4 with y = λ_k ∈ [λ_lo,k, λ_hi,k] by (I1).
(I12) M′ ∈ [0, 40]: definition of F.
Cone: (√T, w_b, w_s) ∈ Q since ‖w_b‖² + ‖w_s‖² = χ²(ũ, M′) ≤ T (Lemma 5(c)).
Linearity of the objective and of λ_i − κ_i, y_b, λ_i is immediate. ∎

## 9. Lemma 6 — weak duality with residual absorption

Consider the conic program

    minimize qᵀy   subject to   A y + s = b,   s ∈ K := {0}^{m₀} × ℝ₊^{m₁} × Q_{m₂},

Q_m := {(s₀, s̄) ∈ ℝ × ℝ^{m−1} : ‖s̄‖₂ ≤ s₀} the second-order cone, and
F_rel := {y : b − Ay ∈ K}. In our instance m₀ = 1971, m₁ = 12571, m₂ = 1593, the
SOC block has rows s = (√T, w_b, w_s) (A has −I on the w columns and b = (√T, 0)).

**Lemma 6.** Suppose Y_k ≥ sup_{y∈F_rel} |y_k| for every coordinate k. Let
z = (z₀, z₁, z₂) with z₁ ≥ 0 componentwise and z₂ ∈ Q_{m₂} (z₀ arbitrary), and
ρ := Aᵀz + q. Then for every y ∈ F_rel,

    qᵀy ≥ −bᵀz − Σ_k |ρ_k| Y_k.

*Proof.* Q is self-dual: for s, z ∈ Q, zᵀs = z₀s₀ + z̄ᵀs̄ ≥ z₀s₀ − ‖z̄‖‖s̄‖ ≥ 0. Hence
z ∈ K* := ℝ^{m₀} × ℝ₊^{m₁} × Q_{m₂} and zᵀs ≥ 0 for all s ∈ K. For y ∈ F_rel,
s = b − Ay ∈ K gives zᵀb ≥ zᵀAy. Therefore
qᵀy = ρᵀy − (Aᵀz)ᵀy = ρᵀy − zᵀAy ≥ ρᵀy − bᵀz ≥ −bᵀz − Σ_k |ρ_k||y_k| ≥ −bᵀz − Σ_k |ρ_k| Y_k. ∎

No optimality, feasibility or accuracy of z is required: any z ∈ K* yields a
valid lower bound, and a poor z only makes the bound weak (through −bᵀz being
far from the optimum and through the loss term).

**Implementation (`verify3.Verifier3.certify`).**
(i) *Cone membership.* z₁ = z_le ≥ 0 is checked exactly on the floats. For the
SOC multiplier, ‖z̄₂‖ is computed as a ball; if z₂,₀ ≥ ‖z̄₂‖ is not provable, z₂,₀
is replaced by fl↑(‖z̄₂‖) (`_endpoint(·, +1)`), which puts the modified z in Q; the
bound is then evaluated for the modified z, so it is valid by the lemma.
(ii) *Ball evaluation.* The rows are produced by `lkr_rows.build(fr, B, T, ArbArith())`:
the same code path as the solver's, but with every real coefficient
(e^{x_k}, c_i, a_i, b_i, lh_k, e_k, tangent/chord data, √T) enclosed in a 160-bit
ball; doubles (grid, data, W, brackets, L, T, z) are exact. ρ_k = Σ_rows A_{rk} z_r + q_k
and val = −bᵀz are accumulated as balls enclosing the exact real quantities of
the real program; loss = Σ_k |ρ_k| Y_k with Y_k = max(|lo_k|, |hi_k|) from
`build.var_lo/var_hi`; finally lb = fl↓(val − loss). The SN dense block is
handled as matrices: ρ_ℓ += 5 W_sᵀ z_s, ρ_{w_s} += z_s, ρ_{M′} += (W_s𝟙)ᵀ z_s,
val −= (W_s(m − off))ᵀ z_s, matching the row w_s + 5W_sℓ + (W_s𝟙)M′ = W_s(m − off).
For the SOC block: val −= √T · z₂,₀ and ρ_w −= z̄₂. Length assertions on z_eq and
z_le tie the multipliers to the rows in build order (eq = [BAO | SN block | κ₀ = λ₀
| segment identities | D_H rows], as stacked in `lkr2.LKRModel2.__init__`).
Conclusion: qᵀy ≥ lb for every y ∈ F_rel(B), where F_rel(B) is the feasible set of
the real relaxed program; the solver's floating-point copy of the program is
never trusted.
(iii) *Solver inaccuracy.* The multipliers come from Clarabel
(`lkr2.LKRModel2.solve_dual`, tolerance 1e-10, relaxed to 1e-6 on retry) and are
merely proposals; inexact optimality/feasibility shows up as ρ ≠ 0 and is paid
for by the loss term. Recorded gaps |solver value − lb| per pass are ≤ 2.4×10⁻⁵
(`lkr_cert_L1.5_D4_r2.log`; e.g. 1.46666049669 vs 1.46666049463 for the final λ₀).
(iv) *Outward rounding (`verify._endpoint`).* For a ball x = (m, r) the endpoint
m ± r is formed in ball arithmetic, printed with 25 significant decimal digits,
converted to a correctly rounded double and moved two ulps outward. The
conversion error is < 1 ulp (½ ulp from decimal→double, ~10⁻²⁵ relative from
printing, ~10⁻⁴⁸ from the 160-bit arithmetic), so the result is ≥ (resp. ≤) the
real endpoint.

**Variable boxes valid on F_rel(B)** (`lkr_rows.build`, "variable bounds"). Each
bound below is implied by rows of the program, hence holds on all of F_rel(B),
as Lemma 6 requires (not merely on lift(F)):
λ_i by (I1); κ_i ∈ [λ_lo,i − ρ_hi,i − c_i, λ_hi,i − ρ_lo,i − c_i] by (I1) and (I6)
(κ₀ = λ₀ with ρ₀ := 0); M′ by (I12); ℓ_j: |ℓ_j| ≤ max_{i∈{k,k+1}} max(|κ_lo,i|, |κ_hi,i|)
+ |c_j| + e_k by (I3) and 0 ≤ t_j ≤ 1; E_{d,i}, E_{r,i}, E_{s,i} by Lemma 4(c) with the
scalings of §6; y_b by the (I10) bracket rows; P_r (D_M) by Lemma 4(c) on
[y_lo, y_hi]; P_r (D_H) by (E5) as a convex combination of the U bounds; U_k by
Lemma 4(c) on the λ bracket; w_b, w_s by |w_·| ≤ ‖(w_b, w_s)‖ ≤ √T from the SOC row.
E_{r,0} appears in no row and not in the objective, so ρ_{E_r,0} = 0 and its
placeholder box (±10⁶) contributes nothing.

## 10. Lemma 7 — bound tightening (the certificate chain)

The chain (`run_lkr_certified.main`, `run_lkr_certified.certified_pass`):
B₀ := `lkr.initial_brackets3(fr)`. In pass p (p = 0, 1, …), for each objective in
the job list J — (ρ, i, lo/hi) for i = 1..N, (λ, i, lo/hi) for i ∈ E ∪ {0}, and
(y_b, r, lo/hi) for the six D_M rows; |J| = 2·372 + 2·13 + 2·6 = 782 on G₃₇₃ — with

    q = e_{λ_i} − e_{κ_i} (ρ lo),  −(e_{λ_i} − e_{κ_i}) (ρ hi),  ±e_{λ_i},  ±e_{y_b},

the solver proposes a dual vector for min qᵀy over its float copy of F_rel(B_p),
and Lemma 6 certifies lb_q over the real F_rel(B_p). All jobs of a pass use the
same B_p. Then B_{p+1} := B_p ∩ (new bounds), and one more certified solve gives
λ₀^{lb,p} := certified lower bound of λ₀ over F_rel(B_{p+1})
(`p<p>_final_lambda0.npz`). The loop stops when |λ₀^{lb,p} − λ₀^{lb,p−1}| < 2×10⁻⁵ or
the pass budget is exhausted; the reported λ₀^lb is the last λ₀^{lb,p}.

**Lemma 7.** (a) B₀ is valid.
(b) If B_p is valid and lb is a certified lower bound of qᵀy over F_rel(B_p), then
the updates
   ρ_lo,i ← max(ρ_lo,i, fl↓(lb − c_i)),   ρ_hi,i ← min(ρ_hi,i, fl↑(−lb − c_i)),
   λ_lo,i ← max(λ_lo,i, lb),   λ_hi,i ← min(λ_hi,i, −lb),   y_lo ← max(y_lo, lb),   y_hi ← min(y_hi, −lb)
(for the respective q) produce a valid bracket vector; hence B_{p+1} is valid.
(c) Every B_p is valid, and λ₀^{lb} ≤ log₁₀ ũ₀ for all (ũ, M′) ∈ F.
(d) A job whose solve fails leaves its bracket entry unchanged, which is valid.
(e) Validity does not depend on convergence or on the number of passes; more
passes only shrink F_rel(B_p) ⊇ lift(F).

*Proof.* (a) For i ≥ 1, ρ_i = log₁₀(ũ_i/D̃_i) = log₁₀(θ(x_i)/(e^{x_i} − 1)) and Lemma 1(c)
gives θ(x_i) ∈ [θ_lo(x_i), θ_hi(x_i)]; `initial_brackets3` evaluates
log₁₀(θ_lo/(e^{x_i} − 1)), log₁₀(θ_hi/(e^{x_i} − 1)) as balls and rounds outward.
For y_b: D̃_M(x_r) = ∫₀^{x_r} ũ e^{x′} dx′ ∈ [ũ_lo, ũ_hi]·(e^{x_r} − 1) = [ũ_lo z_r, ũ_hi z_r] by the
box, so y_b ∈ [fl↓ log₁₀(ũ_lo z_r), fl↑ log₁₀(ũ_hi z_r)]. For λ: the box,
[fl↓ log₁₀ ũ_lo, fl↑ log₁₀ ũ_hi].
(b) By the Proposition lift(F) ⊆ F_rel(B_p); by Lemma 6, qᵀy ≥ lb on F_rel(B_p),
hence on lift(F). For q = e_{λ_i} − e_{κ_i}, qᵀy = λ_i − κ_i = ρ_i + c_i, so ρ_i ≥ lb − c_i
on F; `certified_pass` forms fl↓(lb − c_i) with c_i a ball
(`_endpoint(arb(lb) − c_ball[i], −1)`), a double ≤ the real lb − c_i. For the
negated objective, ρ_i + c_i ≤ −lb. For λ_i and y_b the objective is the variable
itself and lb is already a double. The max/min with the previous entry is the
intersection of two valid brackets, which is valid.
(c) Induction on p using (a), (b); the final job is (b) with q = e_{λ₀}, and
λ₀ = log₁₀ ũ₀ on lift(F).
(d) `_job` returns None on a solver failure and `certified_pass` skips it.
(e) Each pass is a separate application of Lemma 6 to a separate (valid)
relaxation; nothing is inherited except the brackets, which are valid by (c). ∎

Row-level remark. Between passes the row *structure* is fixed and only the
bracket numbers change (bracket rows, tangent points y_t^{(k)}, chord slopes);
`build` regenerates all rows from B_p in both arithmetics, so the solver's
program and the verifier's program are the same program up to enclosure.

## 11. Lemma 8 — the reference point and T

**Lemma 8.** Let ũ* ∈ ℝ^{N+1} (doubles) pass `certify_feasible.in_class_exact`
(exact rational check of the box and of the min-form slope bound with the double
grid and L), let M′* ∈ [0, 40], let [χ_lo, χ_hi] be a ball enclosure of
χ²(ũ*, M′*) (`verify.rigorous_chi2`), and T := fl↑(χ_hi + Δ)
(`_endpoint(ball + arb(Δ), +1)`). Then
(a) (ũ*, M′*) ∈ F(G, L, T);
(b) F(G, L, T) ⊇ { (ũ, M′) ∈ C(G, L) × [0, 40] : χ²(ũ, M′) ≤ χ²_min + Δ },
    χ²_min := inf_{C(G,L)×[0,40]} χ²;
(c) T is an upper endpoint: the theorem is a statement about F(T) exactly as
    defined, and a non-optimal reference only enlarges F(T), so the certified
    upper bound on H₀ is conservative for the Δ-sublevel set in (b).

*Proof.* (a) ũ* ∈ C(G, L) by the exact check (the min-form is the class, Lemma 1(a));
χ²(ũ*, M′*) ≤ χ_hi ≤ T − Δ ≤ T. (b) χ²_min ≤ χ²(ũ*, M′*) ≤ χ_hi, so χ² ≤ χ²_min + Δ
implies χ² ≤ χ_hi + Δ ≤ T. (c) F(T) is monotone in T. ∎

`verify.rigorous_chi2` evaluates D̃_M at the SN and BAO abscissae by Lemma 3(d)
in 160-bit balls (x_j = log(1 + zHD_j) as a ball; e^{x_k} as balls; segment chosen
from a 20-digit midpoint, correct by §13.2), the exact log₁₀, the offsets
5 log₁₀(1 + zHEL_j), W_b, W_s as exact doubles, and sums the squares of the whitened
residuals. Recorded for L = 1.5: χ²(ũ*, M′*) ∈ 1379.7251890839848043 ± 9.65×10⁻¹⁸,
M′* = 16.403897…, T = 1383.7251890839852 (Δ = 4); the per-row T values are in the
table of §2 and are identical in `state.json` and `feasible_<row>_r2.json`.
The reference point itself is a numerically optimized class minimizer
(`classmin.minimize_chi2_over_class`, best of two starts, in
`run_lkr_certified.main`); its optimality is neither claimed nor needed.

## 12. Feasible-point certificates (the lower bounds)

**Lemma 9.** Let ũ_f ∈ ℝ^{N+1} (doubles) pass `in_class_exact`, M′_f ∈ [0, 40], and
let [χ_lo,f, χ_hi,f] be a ball enclosure of χ²(ũ_f, M′_f) with fl↑(χ_hi,f) ≤ T. Then
(ũ_f, M′_f) ∈ F(G, L, T), and therefore

    sup_{(ũ,M′)∈F} c/(r_lo ũ₀) ≥ c/(r_lo ũ_{f,0}),

i.e. at r_d = r_lo the class contains a member with χ² ≤ T and H₀ = c/(r_lo ũ_{f,0}).

*Proof.* Membership by Lemma 1(a) and χ²(ũ_f, M′_f) ≤ χ_hi,f ≤ T. ∎

The certificate is the point itself (`certify_feasible.main` with `--refine 2
--T <chain T>`; output `results/certificates/feasible_<row>_r2.json`, fields
`u_feasible`, `Mp_feasible`, `chi2_feasible_enclosure`, `in_class_exact`,
`certified`). How the point was found (`feasible.max_H0_point`, a numerical
maximization of ũ₀⁻¹ along the χ² = T boundary) is irrelevant to validity. For
L = 1.5, Δ = 4: χ²(ũ_f, M′_f) ∈ 1383.7241890930381910 ± 3.5×10⁻¹⁷ ≤ T = 1383.7251890839852,
M′_f = 16.4300727…, c/(r_lo ũ_{f,0}) = 69.5227585253 ± 4×10⁻¹¹ (ball arithmetic), so
H_lo = 69.5227 after rounding down. Together with the Theorem this gives the
Corollary of §2. Both bounds refer to the same class G₃₇₃, the same M′ box and
the same T; this consistency was the substantive fix of the adversarial review
(CLAIMS P9-1r/P9-2).

## 13. Remarks

**13.1 Covariance convention.** The statistic is defined through the recorded
W (C̃ = (WᵀW)⁻¹), not through the released C; max|W C Wᵀ − I| is 6.2×10⁻¹⁵ (SN) and
2.2×10⁻¹⁶ (BAO). Transferring the theorem to C would require a perturbation bound
on χ²_C − χ²_C̃ over F: with W C Wᵀ = I + E, χ²_C − χ²_C̃ ≈ −w̃ᵀE w̃ for the whitened
residual w̃ (‖w̃‖² = χ² ≈ 1400), so crudely |χ²_C − χ²_C̃| ≲ ‖E‖₂ χ² ≤ n·max|E|·χ² ≈ 2×10⁻⁸,
which is *larger* than the ulp slack in T; a C-statement would therefore need
T enlarged by such a margin. That is not part of the claim; the claim is about
C̃ exactly.
`data/MANIFEST.json` pins the raw files and the derived covariance cache by sha256
(`data.verify_manifest`, `data.load_pantheon_cov`).

**13.2 Floats as exact rationals; transcendental abscissae.** All doubles are
exact. The grid is a double vector, so e^{x_k}, c_k = log₁₀(e^{x_k} − 1), a_k, b_k are
transcendental and appear only as balls; the class constraints, the box and the
exact class check use only rational arithmetic on doubles. The abscissae
x_j = ln(1 + zHD_j), x_r = ln(1 + z_r) are irrational (Lindemann), so no x_j equals a
node; the segment index is chosen from a double approximation
(`geometry.segment_index`, `verify.Verifier._segment`, `rigorous_chi2.dm`). For
Lemma 2(e), Lemma 5(b) and Lemma 3(d) this choice must agree with the real x_j.
Check (200-bit balls, all 1580 SN and 12 BAO redshifts, G₃₇₃): the real x_j lies in
the assigned closed segment for every j, and the smallest distance from any x_j to
a node is 6.8×10⁻⁷, which is 10⁹ times the error of any double approximation, so
every float-based assignment (log1p or a 20-digit midpoint) agrees. Outward
roundings on the verification path: base brackets and λ box (`initial_brackets3`,
`build`), e_k (`build`), δ brackets (`build`), lb − c_i (`certified_pass`),
T (`main`), the final lb (`Verifier3.certify`), the class box (`ClassSpec.u_box`).
The last conversion H₀ = c/(r_lo · 10^{λ₀^lb}) is done in doubles in
`run_lkr_certified.main`; the values in §2 were recomputed in ball arithmetic
(error ~10⁻¹¹) and rounded up at 10⁻⁴, so the double evaluation is immaterial.
Presentation: `results/summary.md` (`report_certified.py`) prints the feasible
lower bounds with `:.4f`, i.e. rounded to nearest; for L = 1.5, Δ = 4 and L = 2 this
shows 69.5228 and 69.5877 where the rounded-down values are 69.5227 and 69.5876.
The certificates themselves are unaffected; §2 uses the rounded-down values.

**13.3 ΛCDM reference.** FORMULATION §1 foresees a second statement with the
flat-ΛCDM best fit as reference. Nothing in §§3–10 changes: one needs ũ_ΛCDM (the
node vector of the fitted H(z)) to pass `in_class_exact` — plausible for L ≥ 1.5
since d ln H/d ln(1+z) = (3/2)Ω_m(z) ≤ 1.43 at z = 2.5 for Ω_m ≈ 0.31, but it must be
checked exactly on G₃₇₃ — and T_ΛCDM := fl↑(χ²(ũ_ΛCDM, M′) + Δ). Since χ²_ΛCDM exceeds
the class minimum by ≈ 20 (1399.7 vs 1380.2 on G₉₄), F(T_ΛCDM) is strictly larger,
and every certificate must be regenerated (T enters the SOC row of every solve);
the resulting statement is "every class member fitting at least as well as ΛCDM
up to Δ has H₀ ≤ …". No such row is certified at present.

**13.4 Loose rows (L ≥ 5).** The certified intervals widen from 0.32 (L = 1.5) to
1.77 (L = 5) and 4.63 km s⁻¹ Mpc⁻¹ (L = 10); the final ρ-bracket widths are
0.48 and 1.01 dex (vs 0.089 at L = 1.5). The gap is relaxation slack, not a
property of F: B_k grows with L (θ_hi − θ_lo ∝ L for small x, and the L² term in B₀),
so e_k ∝ L² h_k² and the per-SN slack is used coherently by the relaxed minimizer;
the chord-minus-tangent gap of each sandwich grows like (bracket width)², and the
wide class lets the ρ brackets stall. Remedies inside the present framework: a
further midpoint refinement (e_k ∝ h²), more tangents (K > 4), tightening λ at all
nodes, or a Jensen-type correction. For L = 1 the class minimum is 217 χ² units
above ΛCDM's: the class is strongly disfavoured by the data, but the row is still
a theorem about its F.

**13.5 Discretization (C(G, L) vs Lip_L).** By Lemma 1(b), C(G, L) ⊂ Lip_L(ln H),
and the theorem concerns C(G, L) only. Of the lemmas, only Lemma 1(c),(d) and
Lemma 2(a),(b) use nothing beyond the Lipschitz bound; Lemma 2(c) (first segment),
Lemma 3 (segment identities) and Lemma 5(b) (D_H by interpolation) use linearity
on segments. An extension to Lip_L would need either (i) an approximation lemma —
for every ũ ∈ Lip_L a member of C(G, L′) (possibly L′ > L) whose D̃_M, D̃_H differ by
explicit O(h²), O(h) amounts, turning χ² ≤ T into χ² ≤ T + ε(h) for the projected
member, so that max over Lip_L(Δ) ≤ max over C(G, L′, Δ + ε) — or (ii) a relaxation
in which the exact identities of Lemma 3 and Lemma 5(b) become inclusions with
class-only widths. Both are deferred (FORMULATION §6, v2 items).

**13.6 Code state.** Function references are to the committed tree at
`2cd323e` (the chain that produced every row of §2). At the time of writing the
working tree carries uncommitted, additive extensions (a D_V/BGS row option
`ClassSpec.use_dv`, DES-SN loaders, `lkr_rows.obbt_objectives`,
`Brackets3.apply_bound`); with `use_dv = False`, which is the case for every
certified row, the variable layout and row set are unchanged by construction
(empty blocks), and the arguments above apply verbatim.

## Appendix — lemma-to-code map

| Statement | Code |
|---|---|
| Class, box, grid (§1) | `model.ClassSpec` (`x`, `u_box`), `geometry.geometric_grid`; `certify_feasible.in_class_exact` |
| χ² (exact) | `model.Frozen.chi2` (float), `verify.rigorous_chi2` (balls) |
| Lemma 1 (a) class rows | `lkr_rows.build` "class" rows; `lh = log10(1 + L h)` |
| Lemma 1 (c) θ-range | `geometry._I_pm`, `geometry.theta_bounds`; `lkr_rows._I_pm`, `lkr_rows.theta_bounds`; base ρ bracket in `lkr.initial_brackets3` |
| Lemma 1 (d) ratio bounds | `geometry._J_pm`, `geometry.kappa_difference_bounds`; `lkr_rows.build` δ bracket rows (`_J_pm`, `_I_pm`, `ar.dn/ar.up`) |
| Lemma 2 (b),(c) B_k | `geometry.kappa_second_derivative_bound`; `lkr_rows.curvature_bound` |
| Lemma 2 (e) e_k, rows | `geometry.kappa_interp_slack`; `lkr_rows.build` "SN interpolation", "BAO D_M rows"; `tests/test_kappa.py` |
| Lemma 3 a_i, b_i, identities | `lkr_rows.build` (`a`, `b`, "segment identities"); `socp2.full_segment_coeffs`; `geometry.dm_matrix` |
| Lemma 4 sandwiches, ranges | `lkr_rows.sandwich_rows`, `lkr_rows.sandwich_range`; scalings 10^{dc}, 10^{c_i} in `build` |
| Lemma 5 BAO rows | `lkr_rows.build` "BAO", "D_H rows", "BAO D_M rows", "U nodes"; `layout_for` (enodes) |
| Proposition (lift) | variable layout `lkr_rows.layout_for`; rows `lkr_rows.build` |
| Lemma 6 weak duality | `verify3.Verifier3.certify`; `verify._endpoint`; boxes `build` "variable bounds"; solver `lkr2.LKRModel2.solve_dual` |
| Lemma 7 chain | `run_lkr_certified.certified_pass`, `run_lkr_certified.main`; `lkr.initial_brackets3`; `state.json`, `p<p>_*.npz` |
| Lemma 8 reference, T | `run_lkr_certified.main` (reference, `rigorous_chi2`, `_endpoint(ball + Δ, +1)`); `classmin.minimize_chi2_over_class` |
| Lemma 9 feasible point | `certify_feasible.main`; `feasible.max_H0_point`; `feasible_<row>_r2.json` |
| Results table | `results/summary.md` (`report_certified.py`); CLAIMS P9-1r, P9-2, P9-3 |
