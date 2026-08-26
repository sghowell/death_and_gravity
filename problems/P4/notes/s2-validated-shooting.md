# P4 S2 — validated shooting and matching for Theorem A (item A3)

Status: S2 item A3, 2026-08-26 (this note is appended to incrementally; the last section
states exactly what is proven).  Code: `problems/P4/src/p4/validated/` — new modules
`shootsys`, `variational`, `lintail`, `tmint`, `matching`, `a3driver` (the A1/A2 modules
are imported, not modified).  Tests: `problems/P4/tests/test_validated_shoot*.py`.
Every number below is a rigorous ball enclosure unless marked "float"/"empirical".
Nothing committed.

## 1. Set-up: what is matched

Unknowns (V0, a, μ): the sonic velocity V0 (A1), and the regular-centre family written as
the translation orbit of A2's normalised family (n∞ = 1, w∞ = μ):
    C(a, μ)(x) = (e^{−a} ñ(t′), e^{2a} w̃(t′), e^{a} ṽ(t′)),   t′ = e^{a+x},
so that a ≈ −ln n∞ = −0.21236 and μ = n∞² w∞ ≈ 8.90132 at S1's values.  The sonic side
Φ(x; V0) = (n, w, v)(x) = (N e^x, W e^{−2x}, V e^{−x}) is the A1 sonic series evaluated at
x0 = −0.05 (with tails) and then integrated with a validated integrator to x_c.  The
matching map and its Jacobian are
    F(V0, a, μ) = Φ(x_c; V0) − C(a, μ)(x_c) ∈ R³,     F′ = [ ∂Φ/∂V0 | −∂C/∂a | −∂C/∂μ ],
    ∂C/∂a = (e^{−a}(θñ − ñ), e^{2a}(θw̃ + 2w̃), e^{a}(θṽ + ṽ)),  θ = t′ d/dt′  (the x-derivative,
    from the series with the certified derivative tail), ∂C/∂μ = (e^{−a}, e^{2a}, e^{a}) ⊙ ∂Ỹ/∂μ.
A zero of F is a solution of the constraint-reduced CSS system on (−∞, x_c] ∪ [x_c, 0.1)
that is analytic through the sonic point on the EC branch (A1) and regular at the centre
(A2): an Evans–Coleman solution.  Krawczyk on a box X ∋ m: K(X) = m − Y F(m) +
(I − Y F′(X))(X − m) ⊂ int X proves existence and uniqueness of the zero in X.

**System integrated.**  The constraint-reduced scaled system of A2 (`systems.centre_system`,
whose t d/dt is d/dx) with T = e^x adjoined (T′ = T): a 4D autonomous polynomial system
P(T,u) u′ = Q(T,u) (`shootsys.shoot_system`, exact fmpq_mpoly re-mapping).  Chosen over the
plain (N, W, V) system because (i) S1 showed the 4D unreduced system is off-constraint
unstable toward the centre, (ii) the scaled variables stay O(1) on the whole range
(N ~ e^{−x}, W ~ e^{2x} in the plain ones), (iii) the regular-centre condition and A2's
series are native in these variables.  Checked against S1's `css.rhs_scaled3` (test).

**Derivative with respect to V0.**  F′ needs ∂Φ/∂V0 = y, which solves the variational
equation.  The augmented system z = (u, y) is again exact-polynomial
(`variational.augment`: P̃ = [[P,0],[Σ_l y_l ∂P/∂u_l, P]], Q̃ = (Q, DQ·y); dT/dV0 ≡ 0 is
dropped, 7 unknowns), so all A1/A2 machinery (recursion, tail bound) applies to it.

## 2. Tails of the parameter derivatives (new certificate, `lintail`)

A1/A2 certify |u_n| for n > K but say nothing about |∂u_n/∂V0|, which F′ needs.  For
n ≤ K the derivative Taylor model is the shifted δ-polynomial plus the Lagrange term
(m+1)·sup_X|u_{n,m+1}|·w^m from A1's interval-base-point run (`variational.derivative_balls`).
For n > K: the y-rows of the augmented level equations are *affine in y with the same level
matrices M_n = nD + E as the u-block* (P̃_yy = P, Ψ̃_yy = Ψ, DQ̃_yy = DQ), so for the fixed
(unknown) u-tail v with ‖v‖_ν ≤ ε_u the map T(z) = z − M^{-1}E^y(ū+v, ȳ+z) on the y-tail is an
affine contraction with the *u-certificate's* constant Z = Z1 + Z2(ε_u); hence
    ‖y-tail‖_ν ≤ (Y_y + Γ ε_u)/(1 − Z),
Y_y the y-residual of the truncation beyond K, Γ the ℓ¹_ν operator-norm majorant of
∂E^y/∂(u-tail) over the ball (same Banach-algebra bounds as Z2, including the θv term).
The direct 8D certificate is useless here (ν ≤ 0.015 even after rescaling y, because the
max-norm couples the blocks); the affine argument recovers the u-block's ν.  Checked against
explicitly computed orders 41–60 at the sonic point and 29–44 at regular points (tests).

Sonic level equations for the y-block: ℓ·(rows 7,8)_n + ℓ′·(rows 3,4)_n with ℓ′ = dℓ/dV0
(over the interval, from the δ-series of ℓ): the ℓ′ term cancels the u_{n+1}-contribution
Ψ̂_0 (n+1) u_{n+1} exactly (from d/dV0[ℓ P_0] = 0); without it the order-K residual is
nonzero, with it every order's residual contains 0 (test).

Results (V0 = c ± 1e−9, K = 40): u-tail ν = 0.1, ε_u = 5.9e−40; y-tail ε_y = 2.2e−36
(Γ = 725), tails at x0 = −0.05: u ≤ 2.7e−52, y ≤ 9.8e−49.  Centre (μ ± 1e−8, K = 30):
ν = 0.05, ε_u = 2.0e−31, ε_y = 4.9e−30; at t′ = e^{−3.21}: u-tail 1.9e−33, y-tail 4.4e−32.

## 3. The validated integrator (`tmint`) — what was tried and what works

**Obstacle (measured, decisive).**  Interval Taylor coefficients of this system blow up
by ≈ 12× per order at x = −0.5 for a box of width 1e−6 (radius 1e21 at order 25), by
≈ 250× per order for A1's interval-base-point trick applied to the step recursion at
x ≈ −0.55 (d²u/dV0² "enclosure" 4.6e72 at order 28), and the crude Picard box over one step
gives ‖Df̃‖ ≈ 2500 (e^{Lh} ≈ 1e6).  The majorant radius of the system (≈ 0.03–0.2, from
`tailbound`) is 10–30× below its empirical radius of convergence (0.83 at x = −0.5, 3.8 at
x = −4).  Consequences: the textbook Lagrange-remainder-over-a-Picard-box Taylor method and
a Taylor model in V0 propagated with A1's Lagrange trick are both useless here.

**Design that works (all balls tight, no interval Taylor coefficients ever formed).**
One reference trajectory m(x) (V0 = c) is integrated by *point* Taylor data; two Lohner
sets ride on it: the point set (V0 = c: tails + rounding) and the interval set (all V0 in
[c ± w]).  Per step x → x − h:
1. Taylor coefficients z_i of the 7D augmented solution through m (K = 28, 384 bits, block
   forward substitution on the block-triangular P̃; the point recursion amplifies rounding by
   up to 500×/order near the sonic point, hence 384 bits);
2. x-tails: the u-block by the A1/A2 Banach certificate on the *point* balls (ν from 0.005
   at x = −0.05 to 0.2 at x = −3; ε ~ 1e−41), the y-block by the affine-contraction bound of
   §2; step h = min(0.02, 0.75 ν) with (h/ν)^{K+1} ε ≤ 1e−45 per step;
3. the Jacobian J of the step map over the hull of the interval set (radius r_k): the point
   fundamental-matrix series Y_K(s) (tight) plus a Grönwall bound of its defect along the
   certified tube Z(s) = Σ z_i s^i + tails + ρ_R, ρ_R = r_k e^{Lh} (fixed point in ρ_R):
       ‖Y(−h; z) − Y_K(−h)‖ ≤ h e^{Lh} ( ‖P̃^{-1}‖ sup|D̂| + ‖D²f̃‖ ρ_R sup|Y_K| ),
   D̂(s) = P̃(z_K(s)) Y_K′(s) − G(z_K(s)) Y_K(s) an *exact polynomial* in s (its coefficient
   sum bounds sup|D̂| without any s-interval dependency blow-up: ~1e−41–1e−60), L, ‖P̃^{-1}‖,
   ‖D²f̃‖ (exact second derivatives of P̃, Q̃) from one-shot box evaluations on 8 sub-intervals
   of the step, all in a weighted norm (y scaled by 1/max|y|).  Typical values: L = 40–130,
   ‖D²f̃‖ = 600–6000, bound 1e−5–1e−3 (relative to |J| ~ 1);
4. mean-value form Φ̃(m + d) ∈ Φ̃(m) + Y_K(−h) d + e with ‖e‖_w ≤ bound·‖d‖_w for both sets:
   the bound is applied as the *operator-norm* bound it is (a box of radius bound·‖d‖_w/Sc_i),
   not entrywise — the entrywise interval matrix costs a dimension factor ≈ 50 and made the
   set width obey a Riccati-type recursion (blow-up at width ≈ 1e−5 around x = −1.3 in the
   first runs); the exact block structure is used (the u-rows of e do not see d_y, the
   T-row of e is 0), tails and rounding are added, and each set is QR re-orthogonalised
   (Lohner).
The step halves h whenever e^{Lh} > 4 or a box solve is not provably regular (needed only
near the sonic point, where ‖P̃^{-1}‖ ∝ 1/|x|).  The tube enlargement ρ_R is refined once
with the fundamental matrix (ρ_R = r_k(sup|Y_K| + bound) instead of r_k e^{Lh}).

**Shape of the initial set.**  The V0-tube at x0 is a *curve*: the interval set is
initialised in Lohner form with the tangent (a1, 2 a2) = ∂(u, y)/∂V0 as its first basis
direction (r1 = |τ| w) and only the O(w²) remainders of the degree-1 Taylor models in the
transverse directions (r_j ≈ 1e−17 for w = 1e−9).  A box-shaped initial set (all directions
1e−9 wide) was amplified ≈ 1000× by the flow's unstable direction and by ∂y/∂u0 = D²Φ·y, and
the run blew up at x ≈ −1.9.

Validation (tests): the point set contains S1's DOP853 trajectory (to S1's 1e−10) and its
central-difference V0-derivative at x = −0.3 and −1; the interval set contains the point
set; the point set has width 1e−47 at x = −0.3 and 1e−26 at x = −1 (starting from S1 data
at −0.5); the interval set for w = 1e−9 has width 3e−8 at x = −0.3.  Cost: ≈ 0.5–0.9 s per
step (63 steps to reach x = −0.3 from −0.05, ≈ 135 more to x = −3).

## 4. The certified trajectory from the sonic point to x_c = −3 (w = 1e−9)

Run: A1 Taylor model at x0 = −0.05 for V0 ∈ [c ± 1e−9], c = 0.112439401388092 (initial
point set 3e−52 (u), 1e−48 (y); interval set a segment with r1 = |τ| w = 1.25e−8 and
transverse radii ≤ 1.5e−17), integrated to x_c = −3 in 297 steps (229 s; h = 0.002–0.007
for x > −0.2, 0.01 for −1.7 < x < −0.3, 0.02 beyond; ν = 0.005 → 0.2).

| x | point-set width | interval set: u | interval set: y | max\|∂u/∂V0\| (float) |
|---|---|---|---|---|
| −0.33 | 1.5e−47 | 1.7e−9 | 5.5e−9 | 0.9 |
| −0.71 | 2.4e−47 | 9.3e−9 | 3.8e−8 | 1.6 |
| −1.59 | 1.0e−45 | 1.7e−8 | 4.6e−8 | 3.5 |
| −2.07 | 1.7e−45 | 7.4e−8 | 6.6e−7 | 13 |
| −2.67 | 2.8e−45 | 1.9e−7 | 3.1e−6 | 77 |
| −3.00 | 9.0e−45 | 5.9e−7 | 6.7e−6 | 207 |

The true width of the V0-tube is |∂u/∂V0|·w (2.1e−7 at x = −3 for u; ≈ 3e−6 for y from
S1's second differences), so the enclosure over-estimates by ≈ 3× (u) and ≈ 2× (y) after
297 steps; the Jacobian bound stays at 1e−6–6e−5 per step.  At x_c = −3:
    (n, w, v)(−3; c) ∈ ([1.23662180055708 ± 3.6e−15], [5.77448322178866 ± 1.2e−15],
                        [−0.402996853997845 ± 3.3e−16])          (S1: 1.2366218, 5.77448322, −0.40299685)
    ∂(n, w, v)/∂V0(−3; c) ∈ ([13.8621536540 ± 4.4e−12], [−149.553641573 ± 2.2e−10],
                             [−207.190706750 ± 1.5e−10])         (S1 central differences: 13.86, −149.55, −207.19)
    ∂(n, w, v)/∂V0(−3; V0) ∈ ([13.862154 ± 6.4e−7], [−149.55364 ± 7.4e−6], [−207.1907 ± 1.4e−5])
                                                                  for every V0 ∈ [c ± 1e−9],
i.e. the Jacobian column of the matching map is enclosed with relative radius 6.5e−8 on the
whole V0 box — the Krawczyk contraction factor is ~1e−7, far from the limit.

## 5. The Krawczyk certificate (Theorem A, item A3) — what is proven

Box X = [c − 1e−9, c + 1e−9] × [a_c − 1e−8, a_c + 1e−8] × [μ_c − 1e−8, μ_c + 1e−8] with
    c = 0.112439401388092,  a_c = −0.21236564676597628…(= −ln 1.2365999612),  μ_c = 8.90132327537996693…
    (= 1.2365999612² × 5.82098013),  x_c = −3,  x0 = −0.05.
Pipeline (`a3driver.run`, 232 s, 297 integration steps, 384 bits): A1 sonic Taylor model
(K = 40, m = 5) with the y-tail certificate → `tmint` to x_c → A2 centre Taylor model
(K = 30, m = 5, normalised family in μ) with the y-tail certificate → F(m), F′(X) → Krawczyk
with Y = mid(F′(m))^{-1} (floats, converted exactly).

    F(m) = ( [−8.7699e−12 ± 4.4e−19], [1.6078e−9 ± 2.7e−14], [−1.0808e−12 ± 4.4e−17] )
    K(X) − m = ( [6.3e−15 ± 6.9e−17], [7.031e−12 ± 7.9e−16], [2.3757e−9 ± 4.4e−14] )
              ⊂ int( [±1e−9] × [±1e−8] × [±1e−8] ).

**Theorem (certified).**  F has exactly one zero (V0*, a*, μ*) in X, and it lies in K(X):
    V0* = 0.11243940138809834 ± 7e−17,   a* = −0.21236564675894539 ± 8e−16,
    μ*  = 8.9013232777556436 ± 4.4e−14,
    i.e. n∞ = e^{−a*} = 1.2365999611913 ± 6e−15, w∞ = μ* e^{2a*} = 5.8209801316354 ± 6e−14,
    v∞ = −1/(2n∞) = −0.40433447815922 ± 2e−15  (S1: 1.2365999612, 5.82098013, −0.40433391 at x = −8).
Consequently there is a solution of the constraint-reduced CSS system (S1's `rhs_scaled3`,
equivalently the KHA system with the momentum constraint) on (−∞, 0.1) which
(i) coincides for x ≥ −0.05 with the analytic sonic-point solution of A1 on the EC branch
    (V1 > 0) with sonic velocity V0*, and hence crosses the sonic point analytically;
(ii) coincides for x ≤ −3 with the regular-centre solution of A2 with parameters (a*, μ*):
    (n, w, v) → (n∞, w∞, −1/(2n∞)) as x → −∞ with corrections in integer powers of e^{2x},
    so the centre is regular (A − 1 = O(e^{2x}), no point mass);
(iii) on [−3, −0.05] lies in the certified tube of §4 (point-set width ≤ 9e−45).
It is the Evans–Coleman solution; uniqueness holds in the sense that no other (V0, a, μ) ∈ X
gives such a solution (Krawczyk).  With S1's sign data this is Theorem A's existence claim
with x_sp ↔ V0* pinned to 16 digits.

**A4 (zero count of V_R = V).**  Along the certified solution, V = v e^x has exactly one
zero on (−∞, 0] (`matching.sign_certificate_v`, test): the A1 series gives V > 0 on
[−0.05, 0] (10 sub-intervals); the integration log gives exactly one step whose certified
tube contains v = 0, x ∈ [−0.25813, −0.25063], with v′ > 0 there and v of fixed sign on every
other step (S1: zero at −0.2509); the A2 series gives ṽ ≤ −1/2 + 1.7e−3 < 0 for x ≤ −3.
(The zero at the centre is the limit V → 0 as x → −∞, not a zero at a finite point.)
Hunter type (a) on the closed sound cone is therefore certified as well.

## 6. What remains, and recommendations for Theorem B

**Remaining for Theorem A.**
- Local uniqueness is uniqueness of the zero of F in the box X (V0, a, μ).  "Unique V0
  in [c ± w]" additionally needs that every regular-centre solution meeting the sonic-side
  tube at x_c has centre parameters in the (a, μ) box, i.e. a global injectivity/graph
  statement for the regular-centre manifold at t′ ≤ e^{x_c + a}; A2's series gives
  v = −1/(2n) + O(t′²) uniformly, so this is plausible but not certified here.
- Extension of the enclosure past the sonic point (x > 0) and to t = 0 (x → +∞): the same
  integrator applies to the plain 3D system outward from x = +0.05 (the A1 series covers
  |x| < 0.1); not run.
- Larger V0 boxes: the interval set's width is ∝ w and the Jacobian bound is ∝ (set width).
  A run with w = 1e−7 (w_a = w_μ = 1e−6) failed at x ≈ −0.97 ("Jacobian enclosure failed"):
  the tube enlargement ρ_R = r_k e^{Lh} uses the *weighted* radius of the whole 7D set, so the
  (100× larger) y-width is fed into the u-tube, the box evaluations of L and ‖D²f̃‖ widen,
  and e^{Lh} > 4 persists down to h ≈ 1e−5.  Remedy (not implemented): block-wise enlargement
  (the u-deviation depends on d_u only, with L_u = ‖Df‖ ≈ 40) and, for w ≳ 1e−6, a
  second-order (Taylor-model) representation of the set or a log-norm in an adapted norm.  The certified radius of the A1 sonic series (0.1, vs empirical 0.99) forces
  the start at x0 = −0.05, where ‖P̃^{-1}‖ ≈ 19 makes the first ~60 steps tiny (h ≈ 0.002–0.007);
  the sharper sonic tail bound suggested in the A1/A2 note would remove this.
- Performance: ≈ 0.5–0.9 s per step, ≈ 200 steps to x_c = −3; Newton refinement + certificate
  ≈ 10 min.  The bivariate machinery is not used by the integrator at all.

**Recommendations for Theorem B (validated integration of the linearised system, complex κ).**
1. The perturbation system is linear in (N_p, ω_p, V_p) with coefficients along the certified
   background: exactly the situation of `tmint.jacobian_step` (a linear ODE with certified
   coefficients handled by a point series plus a Grönwall bound of its defect along the
   background tube).  Reuse it column-wise: the fundamental matrix of the linearised system
   is what E(κ) needs (Wronskian/mismatch at x_c), and the Grönwall bound with
   L = sup‖A(x; κ)‖ over the background tube is uniform in κ on a box.
2. Complex κ: `arbseries.Series`/`arb_mat` have exact `acb` twins in python-flint; the
   recursion/tail-bound modules are ring-generic except for `abs_upper`/`norm_inf` (use
   `acb.abs()` upper bounds).  The affine-contraction tail argument (`lintail`) applies
   verbatim to the linearised sonic/centre expansions: they are affine in the perturbation
   with the background's level matrices shifted by κ, M_n(κ) = nD + E + κ(·); its
   invertibility for all κ in a box is the resonance check (S1: poles at κ = −0.099 − 1.099n
   only), so on Re κ > −1.19 a single certificate per box suffices.
3. E(κ) enclosure on a box: propagate a *point* fundamental matrix at the box centre and
   bound the κ-variation by the same Grönwall estimate with ‖∂A/∂κ‖·(box radius); then
   the argument principle on ∂R with interval-enclosed E and a Krawczyk enclosure of κ1,
   as in S1's `spectrum.winding_number`, but with balls.  The gauge mode κ̄ = 2 − A0 + 2ω0/3
   is available exactly from the certified sonic data (A1 balls).
4. Precision: 384 bits were necessary near the sonic point for the point recursion (rounding
   amplification 500×/order at x = −0.05); the same will hold for the linearised recursion.

## 7. Files, tests, reproduction

New modules (`problems/P4/src/p4/validated/`, nothing pre-existing modified):
- `shootsys.py` — the 4D polynomial shooting system, Taylor/variational coefficients with
  block forward substitution, box-safe polynomial evaluation (`arb ** int` is NaN on
  zero-containing balls in python-flint 0.9), Hessian norm;
- `variational.py` — augmented (variational) polynomial systems, derivative Taylor models,
  augmented level equations (with the ℓ′ terms);
- `lintail.py` — the affine-contraction tail certificate for parameter derivatives;
- `tmint.py` — the Lohner-type integrator (two sets, Grönwall–defect Jacobian, QR);
- `matching.py` — sonic-side initial sets, centre side with derivatives, Krawczyk, the A4
  sign certificate;  `a3driver.py` — Newton refinement + certificate driver.
Tests: `test_validated_shoot.py` (4 tests, ≈ 60 s: system/Jacobian vs S1, y-tail
certificate vs explicit orders, ℓ′ level equations, short validated integration vs S1) and
`test_validated_shoot_matching.py` (1 test, 232 s: the full certificate incl. A4 — passed).
Result file: `problems/P4/results/a3_midpoint.json` (the certified midpoint and K(X) − m).

Reproduction:
```
uv run pytest problems/P4/tests/test_validated_shoot.py -q                 # 4 tests, ~60 s
uv run pytest problems/P4/tests/test_validated_shoot_matching.py -q        # Krawczyk + A4, ~5 min
PYTHONPATH=problems/P4/src uv run python -m p4.validated.a3driver -3.0 1e-9 1e-8 1e-8 2
```
