# P4 S2 — validated local expansions for Theorem A (items A1, A2)

Status: S2 items A1–A2 done (CERTIFIED-level ingredients), 2026-08-25. Code:
`problems/P4/src/p4/validated/` (python-flint / Arb ball arithmetic; imports the S1 modules
`p4.css`, `p4.taylor`, `p4.shoot` only for cross-checks). Tests:
`uv run pytest problems/P4/tests/test_validated_*.py -q` (18 tests, ≈ 3 s).
Every number below is a rigorous enclosure unless marked "float". Nothing committed.

Modules: `arbseries` (truncated series over `arb`), `biseries` (bivariate series in
(x, δ), Kronecker-packed), `polysys` (exact `fmpq_mpoly` systems, evaluation, majorant
evaluation), `systems` (the two systems), `recursion` (certified level-by-level solve),
`tailbound` (Banach fixed-point tail bound), `sonic` (A1), `centre` (A2).

## 1. The systems in cleared-denominator polynomial form

Both expansions are computed for a system  **P(t,u)·θu = Q(t,u)** with P a d×d matrix and
Q a d-vector of *exact* rational polynomials (verified against S1's float right-hand sides
`css.rhs_plain`/`css.rhs_scaled3` at random points to 4e−15, `test_validated_core.py`).

**Sonic point** (t = x, θ = d/dx, u = (A, N, W, V), S = 1 − V², Φ_A = S·F_A =
(1−A)S + 2W(1+V²/3), F_N = −2 + A − 2W/3): rows 1–2 are S·A' = AΦ_A and N' = N F_N;
rows 3–4 are 3WS × (S1 row 3) and WS × (S1 row 4):
  3(1+NV)S W' + 4W(N+V) V' = NVWΦ_A − 4NVWS F_N − 6NVWS − (8/3)NVW²,
  (4V+N+3NV²)S W' + 4W(1+V²+2NV) V' = −WN[SΦ_A + 4(1+V²)S F_N + 2(1+3V²)S].
Degrees ≤ 7 in u. The fluid principal block has determinant
Δ̃ = −4(V−1)(V+1)W·(3N²V² − N² + 4NV − V² + 3) (exact factorisation), whose last
factor is KHA's sonic condition.

**Regular centre** (t = ε = e^x, θ = t d/dt = d/dx, u = (n, w, v) = (N e^x, W e^{−2x},
V e^{−x}), constraint-reduced with A − 1 = 2WT/S, T = 1 + V²/3 + (4/3)NV, S = 1 − v²t²).
Substituting N = n/t, W = w t², V = v t and F_A = G = −(8/3)NVW/S into S1's rows,
multiplying row 3 by 3wS and row 4 by t·w·S (this clears all negative powers of t):
  S n' = (4/3) n w t² (1 + 2nv + v²t²),
  3S(1+nv) w' + 4w(n + v t²) v' = −[6wS(1+nv) + 4wv(n+vt²) + nvw(2S + (8/3)w t²(3 + 5nv + 2v²t²))],
  (n + 4vt² + 3nv²t²)S w' + 4w t²(1 + v²t² + 2nv) v' = −[2w(n+4vt²+3nv²t²)S + 4wvt²(1+v²t²+2nv)
        + nw(S(−2 + 2v²t² − (8/3)nvw t²) + (16/3)w t²(1+v²t²)(1+2nv+v²t²))].
At t = 0: row 4 gives w' = 0, row 3 gives v' = −3(1+2nv)/(2n), row N gives n' = 0, so the
fixed-point set is {nv = −1/2} (S1's N̂V̂ = −1/2; Â = 2w/3 follows from the constraint) and
the Jacobian P₀⁻¹DQ₀ of the t = 0 system has spectrum {−3, 0, 0} (float check; the −3 is
the point-mass mode e^{−3x}, the two zeros are the parameters n_∞, w_∞).

## 2. The coefficient recursion and its certificate

For a series u = Σ u_n tⁿ the residual coefficients F_m = [P(t,u)θu − Q(t,u)]_m are used
through *level equations* E_n: at the sonic point E_n = (F¹_{n−1}, F²_{n−1}, F^{fl}_{n−1},
ℓ·F_n) with ℓ = (P⁴_W, −P³_W)(u₀) the left null vector of the (rank-one) fluid block and
"fl" the better-scaled fluid row; at the centre E_k = F_k. E_n is affine in u_n with the
exact matrix (`recursion.py`, derived from (θu)_j = s_j u_{j+σ}, s_j = j + σ, σ = 1 for
d/dt, 0 for t d/dt)
    M_n = n·D + E,   D_row = Σ c P^{(r)}_{s+σ}(u),   E_row = Σ c (Ψ^{(r)}_s − DQ^{(r)}_s)(u),
where P_k, DQ_k are Taylor coefficients of P(t,u(t)), ∂Q/∂u(t,u(t)) and
Ψ = Σ_i (θu_i)·∂P_{·i}/∂u; all depend on u₀, u₁ only. The formula is checked against a
finite-difference extraction at several orders (test). Each order is then one ball solve
`arb_mat.solve(M_n, −E_n(u_{<n}))`, which succeeds only if M_n is *proven* invertible —
this is the non-resonance certificate at every order ≤ K — and the balls enclose the exact
coefficients by induction. Self-check: every level residual of the truncated series
contains 0.

**Zeroth order (sonic).** KHA99's closed forms A₀ = (7+2√3V₀−3V₀²)/(4(1−V₀²)),
N₀ = (√3−V₀)/(1−√3V₀), W₀ = (3−2√3V₀−3V₀²)/(8(1−V₀²)) satisfy Δ̃(u₀) = 0, ℓ·Q₀ = 0
(KHA99 214–215) and the constraint C̃ := (A−1)S − 2W(1+V²/3+(4/3)NV) = 0 (KHA99 211)
**exactly**: verified as polynomial identities in Q[V₀, √3] with `fmpq_mpoly`
(homogenised numerators reduced modulo s² − 3; `test_sonic_closed_forms_exact_in_Q_V0_sqrt3`).

**First order (sonic).** A₁ = A₀F_A(u₀), N₁ = N₀F_N(u₀); the fluid rows at order 0 give
W₁ = α + βV₁; ℓ·F₁ = 0 is an exact quadratic in V₁ whose coefficients are extracted from
three evaluations (exact since the map is quadratic). Discriminant 8.520048 > 0 (ball);
roots V₁ ∈ {−0.2965262, +0.4847220}; EC branch = the root with V₁ > 0. The extracted
quadratic is proportional to S1's closed form `css.first_order_quadratic` to 1e−12.

**Momentum constraint along the 4D series.** Exact identity (fmpq_mpoly division,
remainder 0): along the 4D CSS flow  S·Δ̃·dC̃/dx = Λ·C̃  with Λ a degree-10 polynomial
(`systems.sonic_constraint_propagation`). Δ̃(u(x)) has a simple zero at x = 0
(Δ̃₁ = 9.72505345987352 ± 6e−16 on the EC branch) and C̃(u₀) = 0 exactly, so the Taylor
coefficients c_n of C̃(u(x)) obey (n S₀Δ̃₁ − Λ₀)c_n = (combination of c_j, j < n); since
γ := Λ₀/(S₀Δ̃₁) = 0 (ball [±6e−76]) is not a positive integer, c_n = 0 for all n:
**the analytic sonic solution satisfies the momentum constraint identically**, and the
constraint-reduced 3D solution used at the centre and for shooting is the same object.

## 3. Parameter dependence: interval V₀ handled as a Taylor model

Naive ball propagation of an interval V₀ of width 1e−10 through the recursion is useless:
the radii grow ≈ 5× per order (dependency/wrapping effect; the rate is the inverse of the
contraction radius of §4): radius 6 at n = 10, 1.5e8 at n = 20, 3e21 at n = 40 (float
diagnostic). Instead the recursion is run on **bivariate series in (x, δ)**, δ = V₀ − c,
truncated at δ^{m+1} (`biseries.py`; exact for the retained δ-coefficients because the
recursion is triangular in the δ-degree; the closed forms and the quadratic root are
computed as δ-series by series division/sqrt). Two runs:
(i) base point c (exact centre): the polynomial part Σ_{k≤m} u_{n,k} δ^k, radii ≤ 1e−45
    (K = 40) / 5e−33 (K = 60) — pure rounding at 256 bits;
(ii) base point the whole interval X = [c−w, c+w]: its top coefficient encloses
    u_n^{(m+1)}(ξ)/(m+1)! for all ξ ∈ X, so by Taylor's theorem with Lagrange remainder
    |R_n| ≤ rem_n := sup_X |u_{n,m+1}| · w^{m+1}; the interval blow-up of run (ii) is
    harmless because it multiplies w^{m+1}.
The enclosure over X is u_n(X) ⊂ Σ_{k≤m} u_{n,k}[−w, w]^k + [−rem_n, rem_n].
Results for c = 0.112439401388092, w = 5e−11 (width 1e−10):
  K = 40, m = 3 (1.4 s): ball radii 1.5e−10 (n=0) → 1.5e−8 (n=40) — the true sensitivity
    ≈ n|u_n|·w — remainders 2e−40 (n=0) → 1.6e−10 (n=40);
  K = 60, m = 5 (5.2 s): ball radii 1.5e−10 (n=0) → 5.1e−8 (n=60), remainders
    2e−60 (n=0) → 3.3e−13 (n=60). The δ-Taylor coefficients grow ≈ 20× per δ-degree
  (δ-analyticity radius ≈ 0.05), so m = 5 keeps rem_n ≪ ball radius up to K = 60.
General widths are supported (`width=`, `m=`); rem_n ∝ w^{m+1}·5ⁿ·20^{m+1} gives the
usable (w, m, K) range, e.g. w = 1e−6, m = 4, K = 12 has radii < 1e−3.
S1's float coefficients (`taylor.background_series`, K = 40) lie inside every ball (test (i));
the point run agrees with them to 4e−14 (S1's own rounding).

## 4. Tail bound: Banach fixed point on the tail in ℓ¹_ν  (rigorous geometric majorant)

**Lemma.** Let ū = (u_n)_{n≤K} be the (exact, enclosed) coefficients, D, E as in §2,
ν > 0, |·| the max norm on R^d, ‖v‖_ν = Σ_{n>K}|v_n|νⁿ. Suppose
 (a) g := ‖D⁻¹E‖_∞ < K+1, and put c := ‖D⁻¹‖_∞/(1 − g/(K+1)) (then ‖M_n⁻¹‖ ≤ c/n, n > K,
     and M_n is invertible for every n > K, so the formal series is unique);
 (b) Y := Σ_{n>K} (c/n)|E_n(ū)|νⁿ  (a finite sum: the residual of a polynomial);
 (c) Z₁ := c Σ_{k≥1} ‖B̃_k‖_∞ νᵏ with B̃_k the entrywise bounds of the off-diagonal blocks
     ∂E_n/∂u_m at ū (k = n − m), using s_{m−σ}/n ≤ 1 and 1/n ≤ 1/(K+2):
     B̃_k[i,l] = Σ_(c,r,s)|c|(|P^{(r)}_l[k+s+σ]| + (|Ψ^{(r)}_l[k+s]| + |DQ^{(r)}_l[k+s]|)/(K+2));
 (d) Z₂(ε) := c‖[ν^{−(s+σ)}inc(P^{(r)}_l) + ν^{−s}(K+2)⁻¹(Σ_i‖θū_i‖_ν inc(∂P^{(r)}_{il}) + inc(DQ^{(r)}_l))
     + ν^{−(s+σ)} ε Σ_i |∂P^{(r)}_{il}|^{abs}(‖ū‖_ν + ε)]‖_∞, inc(p) = p^{abs}(‖ū‖_ν+ε) − p^{abs}(‖ū‖_ν)
     (Banach-algebra Lipschitz bounds; the last term is the θv contribution, whose factor
     s_j ≤ n cancels the 1/n).
If Z₁ + Z₂(ε) < 1 and Y + (Z₁+Z₂(ε))ε ≤ ε, then the map T(v)_n = v_n − M_n⁻¹E_n(ū+v) is a
contraction of the closed ball B_ε ⊂ ℓ¹_ν into itself; its fixed point is a tail solving
the recursion, hence (uniqueness) the true tail. **Conclusion:** Σ_{n>K}|u_n|νⁿ ≤ ε, so
|u_n| ≤ ε ν^{−n} (n > K), the series converges on |x| < ν, its sum is an analytic solution,
and for |x| ≤ r ≤ ν:  Σ_{n>K}|u_n||x|ⁿ ≤ ε (r/ν)^{K+1},
Σ_{n>K} n|u_n||x|^{n−1} ≤ (ε/ν) q^K (K+1−Kq)/(1−q)², q = r/ν.
All quantities are evaluated in ball arithmetic over the interval enclosures of ū, so the
certificate holds for every V₀ in the interval. (`tailbound.py`; `Certificate.tail_bound`.)

**Sonic results (EC branch, V₀ = 0.112439401388092 point; identical ν for the width-1e−10
interval, see tests).** c ≈ 1.03–1.04, ‖D⁻¹E‖ ≈ 1.9.
| K | ν | Z₁ | Y | ε | tail on [−0.05, 0] | tail on [−0.02, 0] |
|---|---|---|---|---|---|---|
| 40 | 0.1 | 0.799 | 7.8e−41 | 5.9e−40 | 2.7e−52 | 1.3e−68 |
| 60 | 0.1 | 0.780 | 7.2e−61 | 4.9e−60 | 2.1e−78 | 1.1e−102 |
The interval version (width 1e−10, m = 3 or 5) gives the same certificate, so for **every**
V₀ in the interval, with the tail included,
  (A, N, W, V)(−0.05) ∈ ([1.8787434870 ± 6.8e−11], [2.047337137 ± 2.5e−10],
                          [0.3503962287 ± 6.7e−11], [0.0884703313 ± 6.7e−11])
— the validated starting data for the sonic-side integration of A3 (test).
ν = 0.13 fails (Z₁ = 1.07); the certified radius 0.1 is ≈ 1/10 of the empirical radius
0.99 — the usual loss of Cauchy-majorant arguments at the derivative coupling m·P_{n−m},
not a property of the solution. Test (iii): with K = 40 the bound |u_n| ≤ εν^{−n} and
Σ_{41≤n≤60}|u_n|νⁿ ≤ ε hold for the exactly computed orders 41–60.

## 5. Regular-centre expansion (A2)

**Recursion.** With the system of §1, u = (n, w, v)(t), t = e^x, the level-k equation is
(kP₀ − DQ₀)u_k = (known from u_{<k}); P₀ = P(0, u₀) is invertible (det = −4w₀n₀²) and
P₀⁻¹DQ₀ has spectrum {−3, 0, 0} (float; the certificate does not use it: each order is a
proven ball solve and orders > K are covered by ‖D⁻¹E‖ < K+1). Hence: **for every point
(n_∞, w_∞, −1/(2n_∞)) of the fixed-point set there is exactly one formal series solution,
and it converges** (Briot–Bouquet, made explicit by the §4 lemma with θ = t d/dt).
Only even powers of t occur (odd coefficients are exactly 0), i.e. the corrections are in
e^{2x}, consistent with S1's observation.

**Parameters and the translation symmetry.** The family has two parameters (n_∞, w_∞).
Since the CSS system is autonomous in x, x ↦ x + a maps solutions to solutions and acts as
n_k ↦ n_k e^{(k−1)a}, w_k ↦ w_k e^{(k+2)a}, v_k ↦ v_k e^{(k+1)a}; the invariant is
μ = n_∞² w_∞. `centre.py` computes the normalised family n_∞ = 1, w_∞ = μ as a Taylor model in
δ = μ − c (same bivariate machinery as §3) and `rescale(n_∞)` gives any member; the identity
"direct run at (n_∞, w_∞) = rescaled normalised run at μ" is verified as balls (test, 5e−55).
For the matching (A3) both parameters are needed because the sonic point is pinned at x = 0.

**Results at S1's parameters** (n_∞, w_∞) = (1.2365999612, 5.82098013), v_∞ = −1/(2n_∞) =
−0.4043344782 (S1's quoted V̂_∞ = −0.40433391 differs by 5.7e−7: it is the x = −8 value,
which carries S1's e^{−3x} drift; the fixed-point relation is exact), K = 30 (0.05 s): |u_k| ≈ 2.4e2, 2.2e4, 8.8e5 at k = 10, 20, 30 (empirical radius
in t ≈ 0.63, i.e. x ≈ −0.46); point radii ≤ 2e−46. The series reproduces S1's integrated
profile (`shoot`, DOP853 from the sonic point) at x = −4 and −6 to 1e−11 (n), 1.6e−9 (w),
1.3e−9 (v) — the size of S1's parameter/profile error — while at x = −8 S1's profile has
drifted by 5e−7 in v (its report's e^{−3x} amplification), the series has not (test (iv)).
Tail certificate (§4 lemma, Euler form): ν = 0.08 (Z₁ = 0.994, Y = 7.7e−29, ε = 2.0e−26,
c = 1.15), so the series converges for e^x < 0.08, i.e. **valid for x ≤ x_c := −3** with
Σ_{k>30}|u_k|e^{kx} ≤ 8.4e−33 at x = −3, 2.9e−46 at x = −4, 3.4e−73 at x = −6.
(Normalised family at μ = n_∞²w_∞: ν = 0.05, Z₁ = 0.77, ε = 2.0e−31.) Parametrised run
(μ interval of width 2e−8, m = 3, K = 30, 0.6 s): radii 1e−8 (k=0) → 0.63 (k=30) — again
the true sensitivity ≈ k|u_k|·w·2.4 — remainders ≤ 2.3e−6.
Â = (A−1)e^{−2x} is provided as the exact truncated series 2wT/S (`Ahat_series`); its tail
follows from those of (n, w, v) by the constraint (not separately certified).

## 6. What is established, what remains for A3, obstacles

**Established (CAP-grade, every step a ball computation or an exact polynomial identity):**
- (A1) For every V₀ ∈ [c − 5e−11, c + 5e−11], c = 0.112439401388092: the KHA 4D CSS system
  has a unique formal power-series solution at the sonic point x = 0 with the closed-form
  data u₀(V₀) and first-order data on the branch V₁ > 0 (V₁ ∈ 0.4847220 ± 4e−8); its
  coefficients u_n, n ≤ 60, are enclosed as degree-5 Taylor models in V₀ − c (radii ≤ 5.1e−8
  over the interval, ≤ 5e−33 at the point); the series converges on |x| < 0.1 and
  Σ_{n>60}|u_n||x|ⁿ ≤ 4.9e−60·(10|x|)^{61}; its sum is real-analytic, solves all four KHA
  equations and the momentum constraint on (−0.1, 0.1), and the sonic point is a saddle
  crossing with Δ̃₁ = 9.725… ≠ 0 (simple zero of the sonic determinant).
- (A2) For the fixed point (n_∞, w_∞, −1/(2n_∞)) with S1's values, the constraint-reduced
  scaled system has a unique analytic solution in t = e^x, coefficients enclosed to k ≤ 30,
  convergent for e^x < 0.08, with the tail bounds above (x ≤ −3); the two-parameter family is
  the translation orbit of a one-parameter Taylor-model family in μ = n_∞²w_∞.

**Remaining for A3 (validated shooting/matching).** Unknowns (V₀, n_∞, w_∞), matching
(N, W, V) at some x_m ∈ [−1, −0.3] (A follows from the constraint on both sides — proven
above). Recommended tool: a rigorous Taylor-method ODE integrator in Arb (Picard/Taylor
models with the same `Series`/`fmpq_mpoly` machinery: the reduced 3D system is polynomial
after clearing S = 1 − V² and W, so the Taylor-step remainder can be bounded by the same
Banach-algebra majorants as §4), rather than VNODE-LP: everything else in this pipeline is
already in python-flint, the ODE is smooth on [x_c, −δ] ∪ [−δ, 0)… and 256-bit precision is
free, which matters because the matching Jacobian must be enclosed on a parameter box.
Krawczyk/interval Newton on the 3D box: the series side is ready (Taylor models in V₀ and
μ give the parameter derivatives too); the ODE side must propagate the same Taylor models
(or, cheaper, use a box small enough — width ≲ 1e−10 — that plain ball integration with
the mean-value form is contractive).
**Obstacles found.** (1) Naive interval propagation through either recursion blows up ≈ 5×
per order (the inverse contraction radius) — solved by the Taylor-model-in-parameter runs;
the same will apply to the ODE integration and forces either Taylor models in the
parameters or very small boxes. (2) The certified radii (0.1 at the sonic point, 0.08 at
the centre) are ≈ 1/10 of the empirical ones because of the derivative coupling m·P_{n−m}
in the ℓ¹_ν operator norm; this is enough here (start the integration at x = −0.05 with tail
2e−78, and at x = −3…−6 from the centre) but a sharper argument (Toeplitz-structured
tail inverse or a weighted norm) would extend the sonic domain to |x| ≈ 0.3.
(3) The centre series has empirical radius e^x ≈ 0.63 (x ≈ −0.46), so the centre side must
be integrated numerically from x_c ≈ −3 to x_m; the sonic side from −0.05 to x_m.

Reproduction: `uv run pytest problems/P4/tests/test_validated_*.py -q` (18 tests, ≈ 3 s);
`sonic.sonic_expansion("0.112439401388092", K=60, width=5e-11, m=5).certify()`;
`centre.centre_expansion("5.82098013", nhat="1.2365999612", K=30).certify()`.
