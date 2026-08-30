# Modelling scope of Theorem A: the KHA CSS system is the self-similar reduction of Einstein–Euler

Status: 2026-08-29. Closes the gap named in `notes/s2-validated-shooting.md` §5 ("Not part of
the certificate: the derivation of the KHA ODE system from the Einstein–Euler equations").
Code: `src/p4/derivation/` (new; nothing pre-existing modified). Tests:
`uv run --with sympy pytest problems/P4/tests/test_derivation.py -q`.

## 0. Claim

Let (M, g, ρ, u) be a spherically symmetric solution of the Einstein–Euler system
G_{μν} = 8πT_{μν}, ∇_μT^{μν} = 0, T^{μν} = (ρ+p)u^μu^ν + p g^{μν}, p = ρ/3, written in
polar-areal coordinates ds² = −α²(t,r)dt² + a²(t,r)dr² + r²dΩ² (areal radius r; this
gauge exists wherever ∇r is spacelike, in particular on the region covered by the
certificate, where A = a² ≥ 1 is finite), with u^t = 1/(α√(1−V²)), u^r = V/(a√(1−V²)).
Suppose it is continuously self-similar: in s = −ln(−t), x = ln(−r/t) the KHA unknowns
N = α/(a e^x), A = a², W = 4πr²a²ρ, V depend on x only. Then (A, N, W, V)(x) solves the
KHA CSS system (rows 1–4 of `css.coeffs` with ∂_s = 0) together with the momentum
constraint F_A = G (KHA99 eq. 211); conversely any solution of the constraint-reduced 3D
system `css.rhs_scaled3` on which the constraint holds at one point defines, through the
formulas above, a CSS solution of Einstein–Euler. The precise statements are Lemmas 1–4.

## 1. What is computed and how it is checked (deliverable 1: `derivation/einstein_euler.py`)

`derive()` returns every intermediate object as a sympy expression, computed by explicit
index loops (no `sympy.diffgeom`, no tensor package): g, g⁻¹, Γ^λ_{μν}, R_{μν}, R, G_{μν};
T^{μν}, T_{μν} (only u^μu^ν products occur, so no square roots), u·u (= −1, checked);
U^ν = ∇_μT^{μν}; E_{μν} = G_{μν} − 8πT_{μν}. Sanity facts verified exactly: U^θ = U^φ = 0;
E_{tθ} = E_{tφ} = E_{rθ} = E_{rφ} = E_{θφ} = 0; α_t is absent from U^t, U^r (the Euler
equations never see the lapse's time derivative); E_tt, E_tr, E_rr, U^t, U^r contain
only first derivatives of (α, a, ρ, V) (asserted in `to_jet`).

The five first-order equations are then rewritten as rational functions of the
"(s,x)-jet" symbols (a, N, W, V, a_s, a_x, N_s, N_x, W_s, W_x, V_s, V_x, e^s, e^x) via
α = N a e^x, ρ = W e^{2s}/(4π e^{2x} a²), ∂_t = e^s(∂_s + ∂_x), ∂_r = e^{s−x}∂_x
(`to_similarity`; the chain rule is done by sympy on Function objects, then the
derivatives are replaced by symbols and e^s, e^x by positive symbols es, ex, so that
every subsequent identity is a decision `cancel(·) == 0` on rational functions — an
exact procedure). With A = a² and the rows of `css.coeffs` (`kha_rows_symbolic`; its two
float literals 8.0/3.0, 4.0/3.0 are restored to 8/3, 4/3 and asserted to be the only
floats):

    row1 := A_x/A − F_A,   row2 := N_x/N − F_N,   rowM := (A_s+A_x)/A − G,
    row3 := sa W_s + sb V_s + a W_x + b V_x + e,   row4 := sc W_s + sd V_s + c W_x + d V_x + f,

`identities()` establishes, each with residual exactly 0:

    I1   E_tt = N² e^{2s} · row1                                (Hamiltonian constraint)
    I2   E_tr = e^{2s−x} · rowM                                  (momentum constraint)
    I3   E_rr = e^{2s−2x} · (2 row2 + row1)                       (slicing condition)
    I4   after a_s → aG/2 − a_x (rowM = 0), a_x → aF_A/2 (row1 = 0), N_x → N F_N (row2 = 0):
         (U^t, U^r)ᵀ = M · (row3, row4)ᵀ   identically (no use of F_A = G is needed), with
         M = [[W e^{3s}/(4πN²a⁴e^{4x}),  V W e^{3s}/(12πN²a⁴e^{4x}(1−V²))],
              [0,                       W e^{3s}/(12πN a⁴e^{3x}(1−V²))]],
         det M = W² e^{6s}/(48π²N³a⁸e^{7x}(1−V²)) ≠ 0 for W ≠ 0, |V| ≠ 1.

So on any open x-interval where W ≠ 0 and |V| < 1: {E_tt = E_tr = E_rr = 0, U^t = U^r = 0}
⟺ {row1 = row2 = rowM = 0, row3 = row4 = 0}, for time-dependent fields as well (the s-derivative
terms are included; this is the full KHA95 eq. 18 with its momentum constraint, not only the CSS
specialisation). The CSS system is the case W_s = V_s = A_s = 0.

**Published equations reproduced.** KHA95 eq. 18 (γ = 4/3): rows 1, 2 and 4 exactly as
transcribed in `notes/literature-digest.md` §1.1; row 3 exactly as transcribed *except* its
last term, which the derivation gives as 2NV(1 + 4W/(9(1−V²))) — the digest's display has
2N(1 + …), a dropped factor V. This confirms S1's finding (`notes/numerics-report.md`,
Stage A) and is what `css.coeffs` implements; the digest file is unchanged (its §1.1 should be
read with this correction). KHA99 eq. 211, F_A = G, is I2 with A_s = 0 combined with I1.

**Link to the certified object.** `identities_with_certified_system()` converts the exact
`fmpq_mpoly` data of `validated.systems.sonic_system` (the 4D polynomial system P·u′ = Q whose
Taylor models A1 certifies) to sympy and shows P·u′ − Q = (A S row1, N row2, 3WS row3, WS row4)
exactly at W_s = V_s = 0, S = 1 − V². The existing exact test
`test_centre_system_is_the_reduced_sonic_system_exactly` (`tests/test_validated_core.py`) links
`sonic_system` to the constraint-reduced `centre_system` (= `css.rhs_scaled3`), the object of the
Krawczyk certificate. The chain Einstein–Euler → KHA rows → `sonic_system` → `centre_system` is
therefore exact end to end.

**The angular Einstein equation** (`derivation/bianchi.py`). E_θθ (and E_φφ = sin²θ E_θθ,
checked) is second order and is not among the five equations above. Lemma 1 below shows it is
implied; the two facts used are verified exactly: ∇_μG^{μν} ≡ 0 for the general metric
(`contracted_bianchi`, third derivatives of α, a, all four components cancel to 0), and for a
spherically symmetric symmetric tensor with vanishing (t,r)-block, ∇_μE^{μr} = −(2r/a²) E^{θθ}
and ∇_μE^{μt} = 0 (`angular_equation_coefficient`).

## 2. Lemmas

Notation: u = (A, N, W, V); S = 1 − V²; the 4D polynomial CSS system P(u)u′ = Q(u) of
`validated.systems.sonic_system` (rows: S A′ = AΦ_A, N′ = N F_N, and the fluid block
B(u)(W′, V′)ᵀ = (Q₃, Q₄)ᵀ), with Δ̃ = det B = −4(V−1)(V+1)W(3N²V² − N² + 4NV − V² + 3);
constraint C̃(u) = (A−1)S − 2WT, T = 1 + V²/3 + (4/3)NV (KHA99 eq. 211 cleared of
denominators; ∂C̃/∂A = S). The reduced 3D system is rows 2–4 with A := 1 + 2WT/S substituted
(`css.rhs_plain3`/`rhs_scaled3`; in polynomial form `validated.systems.centre_system`).

**Lemma 1 (reduction of the PDE system).** On the region r > 0, W ≠ 0, |V| < 1, N ≠ 0 the
Einstein–Euler system for the metric and fluid of §0 is equivalent to
{E_tt = E_tr = E_rr = 0, U^t = U^r = 0}, and for CSS fields this set is equivalent to
{row1 = row2 = rowM = 0, row3 = row4 = 0} with ∂_s = 0, i.e. to the KHA rows 1–4 together with
F_A = G.
*Proof.* E_{tθ}, E_{tφ}, E_{rθ}, E_{rφ}, E_{θφ}, U^θ, U^φ vanish identically (spherical symmetry;
verified). E_θθ = 0 (hence E_φφ = 0) follows from the other five: ∇_μG^{μν} ≡ 0 and ∇_μT^{μν} = 0
give ∇_μE^{μν} = 0; if E^{tt} = E^{tr} = E^{rr} = 0 on an open set, the r-component of this
divergence reduces to −(2r/a²)E^{θθ} (verified), so E^{θθ} = 0. The equivalence with the rows is
I1–I4 of §1 (the factors N²e^{2s}, e^{2s−x}, e^{2s−2x} and det M are nonzero on the region;
I4's eliminations are the equations rowM = row1 = row2 = 0 themselves). ∎

**Lemma 2 (invariance of the constraint surface).** Let u(x) be a C¹ solution of the 4D system
on an interval I on which S Δ̃ ≠ 0. Then C̃(u(x)) = C̃(u(x₀)) exp∫_{x₀}^{x} Λ/(SΔ̃) du; in
particular C̃(u(x)) ≡ 0 on I iff C̃(u(x₀)) = 0 for one x₀ ∈ I.
*Proof.* `sonic_constraint_propagation` verifies the exact polynomial identity
S Δ̃ · (∇C̃ · f) = Λ · C̃, where f is the vector field obtained from P u′ = Q by Cramer's rule
(A′ = AΦ_A/S, N′ = N F_N, W′ = P_W/Δ̃, V′ = P_V/Δ̃) and Λ is an explicit polynomial (degree 10;
exact division, remainder 0, `test_constraint_is_invariant_of_4d_flow`). Hence c(x) := C̃(u(x))
satisfies the scalar linear ODE c′ = (Λ/(SΔ̃))(u(x)) c with continuous coefficient, whose solution
through c(x₀) is unique. ∎

**Lemma 3 (reduced ⇒ full).** Let y(x) = (N, W, V)(x) be a C¹ solution of the reduced 3D system
on an interval I with N ≠ 0, W ≠ 0, |V| < 1, and let A(x) := 1 + 2W T/S. If the zeros of
Δ̃(u(x)) on I are isolated (e.g. u analytic with Δ̃(u) ≢ 0, or Δ̃ ≠ 0 throughout), then
u = (A, N, W, V) solves the full 4D CSS system on I, the constraint holds identically, and
(by Lemma 1) the fields α = N a e^x, a = √A, ρ = W e^{2s}/(4πe^{2x}A), V solve Einstein–Euler on
{(t, r): t < 0, r > 0, ln(−r/t) ∈ I}. Conversely (Lemma 2), a 4D solution on I with C̃ = 0 at
one point has C̃ ≡ 0 and its (N, W, V) solve the reduced system.
*Proof.* Rows 2–4 hold for u: the reduced rows are exactly rows 2–4 with A = 1 + 2WT/S
substituted (`test_centre_system_is_the_reduced_sonic_system_exactly`, exact fmpq_mpoly
quotients 1, t³S, t²S after clearing denominators). C̃(u(x)) ≡ 0 by construction, so
0 = S·Δ̃·(d/dx)C̃(u) = S·Δ̃·(S A′ + C_N N′ + C_W W′ + C_V V′). By the adjugate form of Cramer's
rule, rows 3–4 give Δ̃ W′ = P_W and Δ̃ V′ = P_V with no division; with N′ = Q₂ and the identity of
Lemma 2 written out (C_A Q₁ Δ̃ + C_N Q₂ S Δ̃ + S(C_W P_W + C_V P_V) = Λ C̃, C_A = S) this becomes
S Δ̃ (S A′ − Q₁) = Λ C̃ = 0. Hence row 1, S A′ = Q₁, holds wherever Δ̃ ≠ 0 and, by continuity
of S A′ − Q₁, at the isolated zeros. The rest is Lemma 1. ∎

**Lemma 4 (through the sonic point; the certified object).** The Theorem A object consists of
(i) the A1 sonic series u(x) = Σ u_n xⁿ, convergent on |x| < 0.1, satisfying P(u)u′ = Q(u) as an
identity of convergent series — this is the meaning of "analytic crossing", Δ̃(u(0)) = 0 being a
simple zero — and satisfying C̃(u(x)) ≡ 0 there (exponent argument, `s2-validated-expansions.md`
§2: (nS₀Δ̃₁ − Λ₀)c_n = combination of c_j, j < n, with Λ₀/(S₀Δ̃₁) ∉ ℕ⁺); (ii) the validated
integration of the reduced system on [−3, −0.05] started from the series data at −0.05; (iii) the
A2 regular-centre series for x ≤ −3, matched at −3 (Krawczyk). Claim: the glued function is
real-analytic on (−∞, 0.1), satisfies all five equations of Lemma 1 pointwise on (−∞, 0.1) (rows
3–4 in the polynomial form B(W′,V′)ᵀ = (Q₃,Q₄)ᵀ, which at x = 0 is the analytic-crossing
condition), and the associated fields solve Einstein–Euler on {t < 0, r > 0, −r/t < e^{0.1}}.
*Proof.* On (−0.1, 0.1) the series is a 4D solution with C̃ ≡ 0; the equations E = 0, U = 0 are
equivalent to the rows there by I1–I4, whose factors are regular at the sonic point (M involves
only W, N, a, V, e^s, e^x; the polynomial rows are 3WS·row3, WS·row4 with 3WS ≠ 0), so the
singularity at Δ̃ = 0 is a property of the *solved-for* vector field, not of the equations. On
[−3, −0.05] the reduced solution is analytic and Δ̃ ≠ 0 (every step of `tmint` bounds ‖P̃⁻¹‖ by a
ball inverse of the principal block along its tube, `validated/tmint.py` l. 255, so its success
certifies invertibility); on (−∞, −3] the A2 series is analytic in t = e^x with reduced block
determinant −4wn² + O(t²) ≠ 0 near t = 0, so Δ̃(u) ≢ 0 and its zeros are isolated. Lemma 3
therefore gives the 4D solution and the constraint on both pieces. Gluing: at −0.05 the integrated solution and the series
are reduced solutions with the same data at a regular point, hence coincide on a neighbourhood
(Picard) and, both being analytic, on their whole common interval; the same at −3 (the Krawczyk
zero is exactly the statement that the data agree). Hence the glued function is analytic and
each of the five equations, an identity of analytic functions on each piece, holds on the
union. ∎

**Domain and what "solution" means at r = 0.** The spacetime statement is for r > 0 (x finite).
The regular-centre assertion (A − 1 = O(e^{2x}), N V → −1/2, corrections in integer powers of
e^{2x} = r²/t²) is the A2 certificate and is what makes the fields even in r; C^∞/analytic
extension to r = 0 in Cartesian coordinates is not claimed here (FORMULATION §3).

## 3. Independent check (deliverable 2: `derivation/independent_check.py`, `derivation/qjet.py`)

Route chosen: exact rational pointwise identity testing on the jet space with forward-mode AD
over ℚ (`QJet`: truncated second-order Taylor polynomials in (t, r, θ, φ) with `Fraction`
coefficients; θ₀ is a point with rational sin/cos, (3/5, 4/5)). At each random rational point
(t₀ < 0, r₀ > 0, α₀, a₀, ρ̄₀ > 0, V₀ ∈ (−0.9, 0.9), random first and second partials) the module
recomputes Γ, R_{μν}, G_{μν} from the metric *on jets*, T̄^{μν} = 4πT^{μν} (so π never appears)
and U^ν = ∇_μT̄^{μν}, and evaluates the KHA rows from a *second* transcription — KHA95 eq. 18 as
displayed in the digest, with the row-3 factor V, written directly in Fraction arithmetic,
without `css.coeffs` and without the a…f decomposition. It verifies exactly (`run`, 10 points in
the test, 0.5 s):
(a) E_tt = N²e^{2s}·row1, E_tr = e^{2s−x}·rowM, E_rr = e^{2s−2x}(2·row2 + row1) at unconstrained
    points — the hand-derived factors; since the second-order jet coefficients are random this
    also shows they cancel from the three constraint/slicing equations;
(b) the Euler equations are independent of α_t;
(c) at points where a_r, a_t, α_r are solved from E_tt = E_tr = E_rr = 0 (each affine in its
    variable; the residuals are checked to be exactly 0 afterwards, and row1 = row2 = rowM = 0
    there), (U^t, U^r) is affine in (ρ̄_t, ρ̄_r, V_t, V_r) (checked), the 2×2 matrix M is
    reconstructed from the (ρ̄_t, V_t) columns alone, and the (ρ̄_r, V_r) columns and the source
    terms then match exactly, with det M ≠ 0 — so U = 0 ⟺ rows 3–4 = 0 at the point;
(d) Misner–Sharp: m = r(1 − a⁻²)/2 obeys m_t = r²T̄_t{}^r, m_r = −r²T̄^t{}_t at such points
    (Hayward's covariant form of HM01 eqs. 3–4), anchoring E_tt, E_tr to the comoving
    formulation of digest §1.2 without any coordinate change.
`test_pointwise_M_and_rows_match_symbolic_exactly` closes the loop: the symbolic rows and
4πM of (1), evaluated at the rational points of (2), equal the Fractions of (2) exactly.

Why this is independent of (1): no sympy, no symbolic differentiation or substitution, no
`cancel`/`simplify`, no shared source; the chain rule to (s, x) is done by explicit rational
formulas (f_x = r f_r, f_s = −t f_t − r f_r) rather than by sympy on Function objects; the rows
come from a different transcription; and the pass/fail decision is equality of rationals. What it
shares with (1) by necessity are the definitions (metric, T^{μν}, u^μ, p = ρ/3, the KHA
variables) and the textbook formulas for Γ and R_{μν}; a common error there would not reproduce
KHA's published rows nor the Misner–Sharp identity, both of which are matched. The pointwise
evidence is probabilistic in the Schwartz–Zippel sense (an identity of rational functions of
14 variables false as a polynomial identity would be caught at a random point with probability
≥ 1 − deg/|range|); the deterministic proof is (1). A float cross-check (`cross_check_css`)
confirms the Fraction transcription agrees with `css.fluid_residuals` to 3e−15.

## 4. Status of the modelling gap

Rigorous now (exact, machine-checked, `tests/test_derivation.py`, ≈ 5 s):
- Einstein–Euler (polar-areal, p = ρ/3) ⟺ KHA95 eq. 18 (with the row-3 factor V) + KHA99 eq. 211,
  including the time-dependent terms (I1–I4), and the implied angular equation (Bianchi);
- KHA rows ⟺ `css.coeffs` (up to the two float literals, restored) ⟺ `validated.systems.sonic_system`
  (exact fmpq_mpoly) ⟺ `centre_system` (existing exact test) — the certified object;
- Lemmas 1–4: a reduced solution with the constraint at one point is a full CSS Einstein–Euler
  solution, through the sonic point by the series, on x ∈ (−∞, 0.1).

Still relying on reading rather than on an exact identity:
- that the metric/fluid/velocity conventions of §0 are KHA95's (they reproduce KHA95's rows as
  digested, up to the row-3 slip, which is strong evidence; the paper itself was not consulted
  here, only `notes/literature-digest.md`, so whether the slip is the digest's or the paper's is
  not settled — irrelevant to the certificate, which now rests on the derivation);
- on [−3, −0.05], invertibility of the principal block is read off from the success of the
  validated integrator (its ball inverse of P̃ at every step), not re-verified here; on x ≤ −3
  Lemma 3 needs only isolated zeros of Δ̃, which analyticity of the A2 series supplies;
- `css.coeffs`'s float literals are the correctly rounded 8/3, 4/3 (asserted via `nsimplify`);
  the certificate itself never uses these floats.
Not touched: the remaining items of `s2-validated-shooting.md` §6 (global uniqueness in V₀,
continuation past x = 0.1 and to t = 0).
