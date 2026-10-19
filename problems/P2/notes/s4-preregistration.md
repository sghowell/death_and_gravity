# P2 S4 pre-registration: part (i) at P_B, uniformly in l

Status: written 2026-10-18, before any production run, on branch `p2/s4-part-i` (from main at
`dc5bfd27`). Only float design studies (§9) were run before this document. Any later change to
anything fixed here goes in the dated **Amendments** section at the end, never by a silent edit:
test points, thresholds, l ranges, design constants and grids all count.

Notation. h = 1/(l+½). The strip depth is written **d_δ := (½+δ)κ₋** (= 0.55κ₋ at P_B). β is
kept for D3's gap α/κ₋ (S3 review, W3).

## 1. Target and scope

- **Target (FORMULATION P2-1 at P_B).** An explicit l₀ such that for every l ≥ l₀, W_l has no
  zero in S_δ.
- **Reductions.**
  - Im ω > 0 is P2-0, cited at P2-0's level as in P2-1w and P2-3B.
  - Re ω < 0 follows by conjugation symmetry, W_l(−ω̄) = conj W_l(ω) (FORMULATION §1.2; S3,
    zeq_identity Z4).
  - What remains is the closed half-strip H⁺ = {Re ω ≥ 0, −d_δ ≤ Im ω ≤ 0}, for all l ≥ l₀.
- **Target l₀ = 20** (h₀ = 1/20.5). S3 certified l ≤ 19 (P2-3B), so l₀ = 20 closes P2-3 at P_B.
- **Fallback (§7, K4).** Suppose the certificates close only for h ≤ h₁ < h₀. Then
  l₀ := ⌈1/h₁ − ½⌉. If l₀ ≤ 30, the S3 pipeline (unchanged, `p2.part3`) runs for
  20 ≤ l < l₀: about 7.9 core-h at l₀ = 30 (S3 cost model × 1.2). If l₀ > 30, stop and report
  (K4); those l must then wait for an approved budget.
- **Out of scope.**
  - P_A (§10 says only whether the design transfers).
  - P2-6.
  - The assembly (P2-4, P2-5, the ledger rows), which is S5.

## 2. The regime cover of {l ≥ l₀} × H⁺

Write ω = Re ω + iν, with ν ∈ [−d_δ, 0], and the scaled frequency ζ := h·Re ω. The design
constants (§5) are:
- Ω_L = 8κ₋ and Ω_S = 6κ₋ (absolute);
- ζ_SU = 0.9Ω_P, ζ_TL = 0.8Ω_P, ζ_TU = 2.0Ω_P and ζ_A = 1.5Ω_P (scaled; Ω_P/κ₋ = 4.9188).

| regime | Re ω range | method |
|---|---|---|
| **L** (low frequency, near the horizon throats) | 0 ≤ Re ω ≤ Ω_L | each throat by a validated scaled Frobenius series in τ = (r − r_h)/h² (h² an interval, h = 0 included), then Riccati invariance through the forbidden bulk |
| **S** (sub-barrier, including the large-ω throat) | Ω_S ≤ Re ω ≤ ζ_SU/h | each throat by a half-line M3 certificate on a curved contour passing below the turning point (h → 0 and ζ → 0 closed by parameter boxes), then the same bulk Riccati invariance |
| **T** (barrier top) | ζ_TL/h ≤ Re ω ≤ ζ_TU/h | global M3 on the S2 contour x = x_P + e^{iπ/4}s, with an h-independent design and the exact semiclassical splitting h²E = E₀ + hE₁ + h²E₂ (§3.3) |
| **A** (above the barrier, to ∞) | Re ω ≥ ζ_A/h | global M3 on the S2 contour, with phase φ ≤ φ₀ < 0; positive at leading order at ζ = ζ_A, then monotone in Re ω (§3.4) |

**Overlaps.**
- L ∩ S = [6, 8]κ₋.
- S ∩ T = [0.8, 0.9]Ω_P/h.
- T ∩ A = [1.5, 2.0]Ω_P/h.

For every l ≥ 20, S is non-empty (ζ_SU/h ≥ 0.9·Ω_P·20.5 = 90.8κ₋ > Ω_L), so L ∪ S ∪ T ∪ A = H⁺.
The cover is checked as exact rational inequalities on the constants in `certificates/part1/cover_PB.json`.

**Closing the limits.**
- **h → 0.**
  - L: h² ∈ [0, h₀²] is an interval; h = 0 is the Bessel limit.
  - S: the parameters are ζ ∈ [0, ζ_SU] and ε := 1/Re ω ∈ [0, min(1/Ω_S, h₀/ζ)]. h enters only
    through ζ and ε (§3.2), and the boxes include ζ = 0 and ε = 0.
  - T: a per-cell bound valid on all of h ∈ (0, h₀] (§3.3).
  - A: positive at leading order (h⁰), with explicit O(h) terms.
- **ω → ∞.** Only in A, by monotonicity in Re ω along horizontal lines (§3.4).
- **ζ → 0** (Re ω fixed, l → ∞) is the L regime. **ζ → 0 with Re ω → ∞** is S's closing box at
  ζ = 0: the exponential-throat model, through a near-horizon Taylor model in
  ρ = e^{2κ_h(x − c_h)}.

The pole sets E₊ ∪ E_c lie at Im ω ≤ −κ₊ < −d_δ, so they are outside every domain.

## 3. Methods per regime

### 3.1 L: validated scaled Frobenius at each horizon, and the bulk

- **Scaling.** Near r₊ put t = r − r₊ = h²τ. Then h²·(radial equation) is analytic in (τ, h²),
  with a regular singular point at τ = 0 (exponents ±iω/(2κ₊)). Its h = 0 limit is the modified
  Bessel equation.
- **The incoming solution.** It is the Frobenius series at τ = 0 with exponent −iω/(2κ₊), and
  its coefficients are polynomials in h².
- **Evaluation.** It is evaluated in Arb at the exit τ_a, with ω in boxes and h² ∈ [0, h₀²], with
  majorant tail bounds (the S1/S3 recurrence style). The series radius in τ is (r₊ − r₋)/h² ≥ 40
  at h₀, against τ_a ≤ 8.
- **Continuation.** The checker also continues by Taylor steps in τ beyond the Frobenius disc.
  Production does not need this, but the control N3 (small l) does.
- **Exit rule.** h⁻²U(r_a) = max(4|ω|², C_W), with U = F/r² and
  C_W := 4·max(κ₊², max_{[r₊,r_c]}(−W̃)). Here W̃ = FF′/r − U/4, and V_l = h⁻²U + W̃ exactly;
  float max(−W̃) = 0.0079, and C_W is enclosed rigorously.
- **Claims certified.**
  - Re y₊(x_a) > 0, where y = Φ_x/Φ and Φ = rR₊.
  - At r_c, symmetrically, Re y_c(x_a^c) < 0.
- **Bulk.**
  - U is strictly increasing on (r₊, r_P) and strictly decreasing on (r_P, r_c) (certified once).
  - W̃ ≥ −max(−W̃).
  - Hence Re(V − ω²) > 0 on [x_a, x_a^c].
- **Riccati invariance** (`s0-gate.md` §4). It carries Re y₊ > 0 to x_P and Re y_c < 0 back to
  x_P, and Φ₊ has no zero there. So y₊(x_P) ≠ y_c(x_P), and W_l(ω) ≠ 0.

### 3.2 S: half-line M3 throat certificates

- **Normalized form.** For ζ = h·Re ω and ε = 1/Re ω:
  V − ω² = (Re ω)²[Ũ − 1 − 2iεν + ε²(ν² + W̃)], with Ũ = U/ζ². This is exact, so h appears only
  through ζ and ε.
- **Local variable.** Y = k_T(x − x_T), with U(x_T) = ζ² and k_T = |U_x(x_T)|/(2ζ²) (= κ_h as
  ζ → 0).
- **Contour (r₊ side).** x = x_T + Y(t)/k_T, with Y(t) = t + iη(t).
  - It is C², a ray of angle π/4 for t ≤ t₁ (tail toward the horizon, Im x < 0).
  - It passes below the turning point at distance ≥ d_min > 0.
  - It lies on the real axis for t ≥ t_e, and ends at the exit y_a.
  - The r_c side is the point reflection x = x_T^c − Y/k_T (tail in Im x > 0, where the outgoing
    solution decays).
  - Both lie in S2's singularity-free sectors Σ₋ and Σ₊ (Lemma W1).
- **Phase.** φ(t) is C², goes from φ_∞ ∈ [−π/3, −π/4] in the tail to 0, and is 0 for t ≥ t_e.
- **Design selection.** For each ζ-box the parameters (bend width w, exit y_a, φ_∞, phase centre
  and width) are chosen by a deterministic float maximin over a fixed grid; the grid is committed
  in code before production. They are then certified.
- **Certified.** On every t-cell, for all (ζ, ε, ν) in the box:
  - the leading-order term Re(e^{−iφ}x′²(Ũ − 1)) ≥ m_box > 0;
  - its corrections: the −2iεν term, the ε² terms, and the first-order term (x″/x′)u′ and the
    phase costs, both absorbed by Young's inequality with explicit cost;
  - the kinetic coefficient c ≥ 0;
  - the tail beyond the last cell, by disc bounds as in S2 Lemma 2;
  - at the exit, Ũ(x_a) ≥ 1 + δ_a.
- **Conclusion (Lemma W3).** The half-line coercivity gives
  Re(e^{−iφ(t_a)}u′ū)(t_a) ≥ μ‖u‖² > 0. With φ(t_a) = 0 and x′(t_a) = 1, this is
  Re y₊(x_a) > 0 (and Φ₊(x_a) ≠ 0).
- **Bulk.** As in L. Here Ũ ≥ 1 + δ_a between the exits, by the monotonicity of U, and
  ε²·max(−W̃) < δ_a is certified.
- **ζ → 0.** For ζ ≤ ζ₁ (small; fixed in code), Ũ is a Taylor model in ρ_T = e^{2κ_h(x_T − c_h)}
  around the exponential throat, with an explicit remainder.

### 3.3 T: the barrier top, uniformly in h

On the S2 contour, with Q = i(V − ω²), phase φ, ε₁ = η|φ′| and the semiclassical ground state
g = e^{−G} with G′ = Γ/h (Γ independent of h), the identity is exact:

    h²E = E₀(s; ζ) + h E₁(s; ζ, ν) + h² E₂(s; ν),
    E₀ = Re(e^{−iφ} i (U − ζ²)) − cΓ²,
    E₁ = 2ζν cos φ + (cΓ)′,
    E₂ = Re(i e^{−iφ}(W̃ + ν²)) − |φ′|/(4η),

where E is S2's E = q + c′G′ − cG′² + cG″ and c = cos φ − η|φ′|.

- **Design (pre-registered).**
  - φ = φ_∞ s⁴/(s₁⁴ + s⁴), with φ_∞ = −1.4375, s₁ = 7 and η = ½.
  - Γ = (1 − ε_g)·θ(s)·sign(s)·√(q₀(s; Ω_P)/c), with ε_g = 0.05, where
    q₀(s; Ω_P) = Re(e^{−iφ}i(U − U_P)).
  - θ = 1 for |s| ≤ 2, a C¹ cosine taper on 2 ≤ |s| ≤ 4, and θ = 0 for |s| ≥ 4.
  - Near s = 0, Γ ≈ a s with a = (1 − ε_g)√(|U″|/2). In σ = s/h^{1/2} this is the oscillator
    ground state e^{−aσ²/2}.
- **Certified, per s-cell and per (ζ, ν)-box.**
  - E₀ ≥ 0, and E₀/h₀ + E₁ + min(0, h₀E₂) ≥ μ_T > 0.
  - Since E₀/h ≥ E₀/h₀ for h ≤ h₀, this gives h⁻¹·(that bound) ≤ E for every h ∈ (0, h₀]. It is
    uniform in h with no σ boxes.
  - Near s = 0, E₀ = s²·Ẽ₀(s) is used, with Ẽ₀ from a Taylor model of U at x_P, because
    E₀(0) = 0 exactly.
  - The tails, c ≥ 0 and the contour enclosure reuse S2's geometry certificate
    (`certificates/m3/m3_PB_geometry.json`, replayed).
- **Departure from the coordinator's S2-scaling recipe** (s2-report §3), stated plainly:
  - **No spline nodes.** Γ is an explicit smooth formula evaluated in Arb, so there are no nodes,
    and so no sawtooth, rather than node spacing ∝ (l+½)^{−1/2}.
  - **A wider window.** The window is a fixed scaled interval [0.8, 2.0]Ω_P, not a fixed number of
    quanta. It contains every quanta-sized window, and it is needed for the overlaps with S and A.
  - **The h-closing.** The exact splitting above replaces closing boxes in σ.
  - **No hull loss.** E₀ is affine in ζ² and E₁ bilinear in (ζ, ν), so endpoint evaluation is
    exact (S2 review §12).

### 3.4 A: above the barrier, up to ∞

- **Design.** The same contour and tail as T, with φ = φ₀ + (φ_∞ − φ₀)s⁴/(s₁⁴ + s⁴),
  φ₀ = −0.2, and no ground state.
- **Certified.** At ζ = ζ_A, for h ∈ (0, h₀] and ν ∈ [−d_δ, 0]:
  - q₀(s; ζ_A) ≥ m_A > 0 at leading order;
  - E ≥ h⁻²(m_A − h·2ζ_A d_δ − h²·max|E₂|) > 0.
- **Monotonicity (Lemma W7).** ∂q/∂(Re ω) = 2|ω| sin(arg ω − φ(s)), and the ground-state terms do
  not depend on ω. Since φ(s) ≤ φ₀ < −arctan(d_δ h₀/ζ_A) ≤ arg ω, q is nondecreasing along
  horizontal lines. So E ≥ μ at Re ω = ζ_A/h implies E ≥ μ for every larger Re ω.
- Lemma 3's decay condition holds throughout A (Re ω > d_δ).

### 3.5 Written lemmas, and which are non-routine

| lemma | content | status |
|---|---|---|
| W1 | sector analyticity of r(x) on Σ₋ ∪ Σ₊ (S2 Lemma 1, unchanged); every S4 contour lies in it | **non-routine** (as in P2-1w) |
| W2 | continuation of the Jost solutions to the contours, and their decay along the tails (S2 Lemma 3), extended to the curved throat contours with π/4 tails | **non-routine** |
| W3 | half-line coercivity on a curved contour, with the first-order term (x″/x′)u′ absorbed by Young, giving Re y(x_a) > 0 at the exit | routine (integration by parts with one boundary term) |
| W4 | global coercivity plus the ground-state bound with C¹ compactly supported G′ (S2 §3–4, generalised from piecewise-linear G′) | routine |
| W5 | Riccati invariance in the forbidden bulk, and no zeros of Φ | routine |
| W6 | identification of the scaled Frobenius series with the Jost solution, uniformly in h | routine |
| W7 | monotonicity in Re ω for A | routine |
| W8 | the exact splitting V_l = h⁻²U + W̃ and the per-cell h-bound | routine algebra (sympy-checked) |
| W9 | the cover, the conjugation symmetry and P2-0 for Im ω > 0 | routine bookkeeping and citation |

**Expected level: WRITTEN_PROOF.**
- W1 and W2 are non-routine prose, exactly as in P2-1w, and the Im ω > 0 part is at P2-0's level.
- CERTIFIED would need W1 and W2 machine-checked (2D validated continuation of r(x) over the
  sectors). That is not planned in S4.
- The proof note will be `notes/s4-proof.md`.

## 4. Acceptance criteria (KAT-style) and controls

- **K1 (geometry replay).** S2's contour certificate replays (`p2.m3 verify` geometry part).
- **K2 (two Wronskian codes).**
  - L's scaled-Frobenius Jost values, converted to W_l, overlap zsol and arbw (S3's C0 wrapper) at
    12 points: l ∈ {20, 30, 50} × ω ∈ {0, 2κ₋ − 0.2iκ₋, 5κ₋ − 0.55iκ₋, 8κ₋ − 0.3iκ₋}.
  - Both balls must have relative radius below 1e−10.
- **K3 (overlaps).** At 3 sample (ζ, ν, h) points in each overlap (L∩S, S∩T, T∩A), both adjacent
  certificates hold. This is required by the cover anyway, and reported explicitly.
- **K4 (float consistency).** At 200 deterministic random points of the cover (l ∈ [20, 400]), the
  float W_l is nonzero, and the float S0/S1 solver finds no QNM in H⁺.
- **N1 (near-mode, PS).** The T checker is run with ν-boxes containing the l-th PS fundamental
  mode (from the S1 float solver) at l = 20, 50 and 100. The boxes are ±γ₀/50 about the mode, and
  STRAD, Im ω ∈ [−1.05, −0.95]γ₀. Each **must not certify**.
- **N2 (an h-box holding true QNMs).** The T checker on h ∈ [1/50.5, 1/20.5], with the ζ-box
  [Ω_P − 0.001κ₋, Ω_P] and the ν-box [−1.69, −1.68]κ₋. These contain the PS fundamental modes of
  every l ∈ [20, 50] (KAT-2: E_l ≈ −0.112, D_l ≈ 0.0685). It **must not certify**.
- **N3 (near-mode, NE).** The L checker at l = 1 and l = 2 (h = 2/3, 2/5), with ν-boxes containing
  the NE modes −1.911925iκ₋ and −2.891324iκ₋. These lie between the poles −nκ₊; the box is checked
  to avoid E₊. Each **must not certify**.
- **N4 (code path, S).** The S checker with the exit placed before the turning point, in the
  allowed region, **must not certify**.
- **Replay.** `python -m p2.part1 verify` recomputes every certificate's verdict from its stored
  design (S: the design per ζ-box). It recomputes a deterministic 5% sample of cells exactly, and
  re-checks the cover file.
- **Tests.** `tests/test_part1.py`: fast replays of samples, and P2_SLOW full replays of one box per
  regime, plus N1 and N3.

## 5. Design constants (fixed here)

- **The h range.** h₀ = 1/20.5.
- **Regime boundaries (§2).**
  - Ω_L = 8κ₋, Ω_S = 6κ₋;
  - ζ_SU = 0.9Ω_P, ζ_TL = 0.8Ω_P, ζ_TU = 2.0Ω_P, ζ_A = 1.5Ω_P.
- **T.** φ_∞ = −1.4375, s₁ = 7, p = 4, η = ½, ε_g = 0.05, taper [2, 4]. The s-range is S2's
  [−64, 32] plus tails. Cells are bisected to depth ≤ 12.
- **A.** φ₀ = −0.2, the same tail, Γ = 0.
- **S.** Tail angle π/4. The grid is: w ∈ {0.15, 0.3, 0.6, 1.0}, y_a ∈ {0.05, 0.1, 0.2, 0.4, 0.8},
  φ_∞ ∈ {−π/4, −π/3}, phase centre ∈ {−1, −0.5, 0}, phase width ∈ {0.25, 0.5, 1} (§9).
  - ζ-boxes: [0, ζ₁] with ζ₁ = 0.02Ω_P, then width ≤ 0.02Ω_P. They are halved on failure, down to
    0.0025Ω_P.
  - ε-boxes: ε ∈ [0, min(1/Ω_S, h₀/ζ)], halved on failure, down to 4 levels.
- **L.** ω-boxes of side ≤ 0.05κ₋ on [0, 8κ₋] × [−d_δ, 0] (halved on failure, down to 0.0125κ₋),
  h² ∈ [0, h₀²], and the exit rule of §3.1.
- **Precision.** 128 bits; 256 bits for K2.

## 6. Certificates (paths)

Directory `problems/P2/certificates/part1/`:
- `T_PB.json`, `A_PB.json`;
- `S_PB_plus.json`, `S_PB_c.json`, `L_PB_plus.json`, `L_PB_c.json`;
- `bulk_PB.json` (U monotone; the W̃ bound; C_W);
- `cover_PB.json`, `controls_PB.json`, `replay_PB.json`.

Every file records the code commit, the precision, the design and the timing.

Code: `src/p2/part1.py`, with `part1_T.py`, `throat.py` and `lowfreq.py` as needed. Runs are
resumable, with JSONL working files and launch records under `results/part1/`.

## 7. Kill criteria

- **K-unsound.** Any of N1–N3 certifies. Stop: the checker is not trusted.
- **K-QNM.** A checker fails, and the float solver finds a QNM of some l ≥ l₀ in H⁺ at the failing
  box. Stop and report; this would refute part (i) at the target l₀.
- **K-gap.** A regime cannot be certified over its domain within the §5 refinement limits and
  design grids. Its boundaries may then move only by a dated amendment that keeps the §2 cover
  (for example ζ_SU, ζ_TL). If no admissible cover exists, report the uncovered set; P2-1 is not
  established.
- **K4 (l₀ too large).** The achieved uniform h₁ < 1/30.5. Stop and report: the S3 fallback for
  20 ≤ l < l₀ would exceed about 8 core-h.
- **K-compute.** The projection exceeds 40 core-h, or the total reaches the 48 core-h cap.

## 8. Compute budget and projection

- **Projection for the S4 certificates.**
  - T: about 0.3 core-h (≈ 10⁴ s-cells × about 24 ζ-boxes, with bisection).
  - A: < 0.05 core-h.
  - S: about 3 core-h for both sides (about 100 ζ-boxes × ≈ 5·10³ contour cells × ε-boxes, plus
    complex-contour enclosures).
  - L: about 0.2 core-h.
  - Controls and K-checks: about 1 core-h. Replays and tests: about 1 core-h.
  - Development runs: about 3 core-h.
  - **Total about 9 core-h.**
- **Fallback.** If 20 < l₀ ≤ 30: S3 for 20 ≤ l < l₀, up to about 7.9 core-h.
- **Cap.** 48 core-h, at most 4 workers at nice 10. Runs are detached and resumable, with launch
  records, and outputs are committed as they land.

## 9. Pilot (float design studies, before this document)

Scratch scripts, not committed; about 3 CPU-min in all.

- **T design (§5).** Float min over s of h²E/(h·a*) (a* = √(|U″|/2) = 2Ω_Pγ₀), at ν = −d_δ:

  | ζ/Ω_P | 0.8 | 0.9 | 1.0 | 1.5 | 2.0 | 2.5 |
  |---|---|---|---|---|---|---|
  | l = 20 | 0.175 | 0.243 | 0.320 | 0.460 | 0.296 | 0.133 |
  | l ≥ 50 | 0.689 | 0.656 | 0.623 | 0.460 | 0.296 | 0.133 |

  The h → 0 limit at ζ = Ω_P is 0.623 = 0.95 − 0.327, the S2 ideal with 5% slack. The S2 phase
  with the S2 support (s₁ = 5, taper [3, 6]) fails below ζ ≈ 0.95Ω_P; the narrower support fixes
  that.
- **A design.** At ζ_A = 1.5Ω_P, min_s q₀ = 0.0109 > 0, attained at s = 0. The margin is ≫ μ_HO.
- **S** (both sides, best grid design per ζ). The leading-order margin min_t q₀/|Y′|² is:

  | ζ/Ω_P | 0.02 | 0.3 | 0.6 | 0.8 | 0.85 | 0.9 | 0.95 |
  |---|---|---|---|---|---|---|---|
  | margin | 0.56 | 0.51 | 0.34 | 0.17 | 0.17 | 0.087 | 0.042 |

  The r_c side agrees to 0.01. ζ_SU = 0.9 leaves margin 0.087; S ∩ T starts at 0.8.
- **L** (mpmath, the true equation). min over Re ω ≤ 8κ₋ and Im ω ∈ [−d_δ, 0]:
  - Re y₊(x_a)/max(|ω|, κ₊) = 1.44;
  - −Re y_c(x_a^c)/max(|ω|, κ_c) = 0.85;
  - for l = 20, 50, 200 and 1000 alike (the Bessel limit).
- **Fixed functions.** min W̃ = −0.0079 (at r ≈ 3.95). U is monotone on both sides of r_P (float).

## 10. P_A (not in S4)

The design transfers structurally. The binding constant is T's h → 0 margin, which drops from
0.623 to about 1 − 0.667 − 0.05 = 0.283 μ_HO. L's horizon condition Im ω/κ₊ ≥ −0.75 > −1 still
holds. l₀ at P_A would be larger; no run is planned.

## 11. Outcome statements

- **Success.** "At P_B, for every l ≥ l₀ (target 20) and every ω with Re ω ≥ 0 and
  −d_δ ≤ Im ω ≤ 0, W_l(ω) ≠ 0." With conjugation symmetry and P2-0 this gives P2-1 at P_B. Level:
  **WRITTEN_PROOF** (W1, W2 non-routine; Im ω > 0 at P2-0's level).
- **With S3** (l ≤ 19, or ≤ l₀ − 1 after the fallback): P2-3 at P_B closes, and with P2-2 and P2-0
  the strip S_δ at P_B has no QNM other than ω = 0 (l = 0). The assembly is S5.
- **Partial outcome** (K-gap): the certified regimes are reported with their exact domains, and
  the uncovered set is named.

## Amendments

### Amendment 1 (2026-10-18, before any production run; requested by the coordinator)

**W1 for the new S contours, certified contour by contour.**
- S2's Lemma 1 reasons as follows. Along a vertical continuation from the real axis,
  Re r*(r(x)) = Re x, and Σ A_i = 0. So the continuation can fail only where r → 0, on the line
  Re x = Re r*(0) = −2.29664, or where r → ∞, on the line Re x = 0.
  - Off those two lines, the vertical continuation exists at every height. The ODE argument of
    Lemma 1's last bullet makes it analytic jointly in x.
  - On the two lines, S2 checked the segments up to heights 0.1413 and 2.438, respectively. They
    hold for Im x < 0 as well, by conjugation (r(x̄) = conj r(x)).
- **The vertical-strip form of W1, used in S4.** r(x) is analytic on each of three open vertical
  strips, at every height: {Re x < −2.29664}, {−2.29664 < Re x < 0} and {Re x > 0}. The strips are
  joined across the two lines only through the checked segments.
- **For every S design actually used, the certificate checks one of two things:**
  - the contour lies in Σ₋ ∪ Σ₊ (S2's sectors); or
  - the contour crosses the lines Re x = −2.29664 and Re x = 0 only within the checked segments
    (|Im x| ≤ 0.1413, resp. |Im x| ≤ 2.438), from enclosures of the crossing heights.
- The r₊-side contours have Re x ≤ x_a⁺ < x_P < −2.29664, so they lie in the left strip; this is
  checked too.
- The result is recorded as a field `w1_check` in `S_PB_plus.json` and `S_PB_c.json`.
- W2 (Lemma 3's representation and decay) carries over in the same strip form: its disc argument
  uses only Re x ≤ X₋, at any height.
- This extension is stated in W1 (§3.5) and in `notes/s4-proof.md`.

### Amendment 2 (2026-10-18, before any production run; requested by the coordinator)

**New negative control N5, on L's production path.** The L checker is run with h² ∈ [0, h₀²],
ω-boxes in L's domain, and no Taylor continuation, with the exit placed in the allowed region
(h⁻²U(r_a) < (Re ω)²). It must not certify. (N3 runs at l = 1 and 2 through the continuation
path, which production does not use.)

### Amendment 3 (2026-10-18, before any production run; requested by the coordinator)

**Corrected wording of §2, "Closing the limits".** The sentence "ζ → 0 (Re ω fixed, l → ∞) is
the L regime" holds only for Re ω ≤ Ω_L. The correct statements are:
- ζ → 0 with Re ω ≤ Ω_L is the L regime;
- ζ → 0 with Re ω fixed above Ω_L, or with Re ω → ∞, is S's closing box ζ ∈ [0, ζ₁], with
  ε = 1/Re ω ∈ [0, 1/Ω_S].

### Amendment 4 (2026-10-18, before any production run; T design)

**T's ground state is linear in s, not a square root.** §3.3 and §5 fixed
Γ = (1 − ε_g)·θ(s)·sign(s)·√(q₀(s; Ω_P)/c). It is replaced by

    Γ(s) = (1 − ε_g)·a_d·s·θ(s),   a_d = 32178552/2³⁰ = 0.02996861190…,

an exact dyadic just below a* = √(|U″|/2) = 2Ω_Pγ₀ = 0.02996861191… (a* − a_d ≈ 3.7·10⁻¹¹).
ε_g = 0.05 and the taper θ (C¹ cosine, 1 on |s| ≤ 2, 0 on |s| ≥ 4) are unchanged, and so is
everything else in T: the contour, the phase (φ_∞ = −1.4375, s₁ = 7, p = 4), η = ½, the
ζ-range [0.8, 2.0]Ω_P, h₀ and the per-cell bound.
- **Reason.** The √ form vanishes like |s| at s = 0, where q₀/c vanishes quadratically.
  Enclosing it and its first two derivatives on s-cells near 0 would need a Taylor model of
  q₀/c at s = 0 for each derivative (the square root of an interval containing 0 has an
  unbounded derivative). The linear form is a polynomial times the taper: Γ, Γ′ and Γ″ are
  evaluated directly in Arb, and E₀ = s²Ẽ₀ near 0 needs only the second derivative of U, as §3.3
  already plans.
- **Validity is unaffected.** Any C¹ compactly supported Γ is admissible in the ground-state
  bound (W4); the choice changes only the margin. §3.3's remark that Γ ≈ a s near 0 now holds
  exactly on |s| ≤ 2.
- **Float support** (the §9 scripts, same grid). On 0 < |s| ≤ 4,
  0.993 ≤ q₀(s; Ω_P)/(c·a*²s²) ≤ 1.057, so the two forms differ by under 6% on the support. The
  h-uniform margin min_s F/a* (F = E₀/h₀ + E₁ + h₀·min(0, E₂), ν = −d_δ):

  | ζ/Ω_P | 0.8 | 0.9 | 1.0 | 1.5 | 2.0 |
  |---|---|---|---|---|---|
  | linear Γ (this amendment) | 0.156 | 0.224 | 0.301 | 0.460 | 0.296 |
  | √ form (§3.3 as registered) | 0.175 | 0.243 | 0.320 | 0.460 | 0.296 |

  E₀ ≥ 0 on the whole contour for both, with E₀/s² → a*² − (1 − ε_g)²a_d² ≈ 0.0975a*² at s = 0.

### Amendment 5 (2026-10-18, before any production run of L, S and the controls; methods fixed in development)

Development runs (sample boxes only, about 0.5 core-h in all) showed where §3.1, §3.2 and §5 must be made
precise or changed. No production certificate had been written. T and A are as in Amendment 4, with
cell bisection target a*/10.

**S (§3.2, §5).**
1. **Contour family.** For the grid parameters (w, y_a): Y(t) = t + iη(t), with η = t for t ≤ −w,
   η = −w·m₀(−t/w) on |t| ≤ w, and η = 0 for t ≥ w. Here m₀ is the C³ smoothed max(0, u) on [−1, 1]
   (m₀″ = (15/16)(1 − u²)², m₀(0) = 5/32). The exit is t_a = w + y_a.
   - So the contour is the π/4 ray for t ≤ −w and is exactly real for t ≥ w. It passes below the
     turning point at depth (5/32)w.
   - The §9 pilot used a softplus bend, which is never exactly real. W3 needs the exit on the real axis
     (φ(t_a) = 0 and x′(t_a) real).
2. **Grid.** w ∈ {0.15, 0.3, 0.6, 1.0} is extended by {2, 3, 4}.
   - A compact half-width w plays the role of about three softplus widths.
   - At small ζ the best designs need half-widths 2 to 4: the float margin at ζ = 0.02Ω_P is 0.37 (+)
     and 0.63 (c) with them, against 0.06 and −0.03 without.
   - The other grid values are unchanged.
   - The phase is φ_∞ for t ≤ p_c − p_w; φ_∞(1 − S₅(v)) with S₅ = 10v³ − 15v⁴ + 6v⁵ and
     v = (t − p_c + p_w)/(2p_w); and 0 after.
   - Designs with p_c + p_w > t_a are inadmissible.
3. **One contour per ζ-box.** x(t) = X_a + s(Y(t) − t_a)/k (s = ±1). k = k_T and X_a = x_T + s·t_a/k
   are taken at the geometric mean of the box's ζ endpoints and rounded to multiples of 2⁻³⁰. 1/ζ²
   enters q̂ affinely, so the box is handled by endpoint evaluation.
4. **Young weight.** η = cos φ/κ_tot pointwise, so the kinetic coefficient is 0 and
   q̂ = Re(e^{−iφ}Y′²[U/ζ² − 1 + 2iε|ν| + ε²(ν² + W̃)]) − (kε)²κ_tot²/(4 cos φ), with
   κ_tot = |Y″/Y′ − iφ′|.
5. **Cells.**
   - The mean-value form in t, with t-derivatives by acb_series, and separable interval bounds in 1/ζ²,
     ε, ε|ν| and ε²ν². Family polynomials are evaluated in centred form.
   - r(t) by validated Taylor integration (order 20) of dr/dt = x′(t)F(r), backward from the real exit
     r_a = r(X_a) (interval Newton).
   - Steps are halved until the a-priori box exists and the remainder is below 10⁻³⁰. The Lohner
     factor is F(B)/F(r₀), which also holds for the non-autonomous r′ = x′(t)F(r).
   - t_L = −8. The tail t ≤ t_L is bounded by S2's horizon disc (Lemma 2), with Re x linear in t on
     the ray.
6. **Exit and order.** Certified per box: U(r_a)/ζ_hi² − 1 > ε_max²·W_bound, and X_a < x_P (+) resp.
   X_a > x_P (c). Between the two exits, U ≥ min(U(x_a⁺), U(x_a^c)) by monotonicity.
7. **Refinement.**
   - If the selected design's float margin is ≤ 0, the ζ-box is halved at once.
   - Otherwise the ε-box is halved first (≤ 4 levels), then the ζ-box (down to 0.0025Ω_P, as
     registered).
   - Cells start at the contour steps and are bisected to depth ≤ 8, with target
     max(margin_float/4, 10⁻³) rounded to a dyadic and stored for replay.
8. **Closing box ζ ∈ (0, 0.02Ω_P].**
   - The contour is x_T(ζ) + s·Y(t)/κ_h (k = κ_h exactly). Along it,
     ρ := e^{±2κ_h(x − c_h)} = ρ_T e^{2Y}.
   - G(ρ) = U(r(ρ)) is taken as its Taylor polynomial of order 14 (by series reversion of
     ρ(r) = s(r − r_h)f(r)), plus a Cauchy remainder on |ρ| < R = m_f·d_r. By Rouché, r(ρ) stays in
     |r − r_h| ≤ d_r there (d_r = 3/50 at r₊ and 1 at r_c).
   - Ũ = e^{2Y}P(ρ_T e^{2Y})/P(ρ_T), with P = G/ρ, ρ_T ∈ [0, ρ₁], and ρ₁ certified from G(ρ₁) ≥ ζ₁².
   - P on the radial segment is bounded by g₁ + Σ gₙzⁿ⁻¹·[0, 1].
   - |W̃| ≤ M_W(|ρ|/R)/(1 − |ρ|/R).
   - Closing-box designs need t_a ≤ 2.1, so that ρ₁e^{2t_a} < R.
9. **W1 (Amendment 1)** is checked per box as specified. In the closing box the contour stays in the
   near-horizon disc, whose r* range lies in the left strip (+) resp. the right strip (c); this is
   certified through r*(r_h ± d_r).

**L (§3.1, §5).**
10. **Arithmetic.** The scaled Frobenius recurrence runs in first-order affine forms (`affine.py`) in
    the deviations of ω (a disc about the box centre) and of h² (a disc about a sub-interval's centre).
    Second-order terms go into the remainder.
    - Plain Arb rectangles wrap by up to √2 per complex multiplication.
    - Over a box of side κ₋/20 at Re ω = 3κ₋, the series lost all accuracy.
11. **h².** [0, h₀²] is covered by 16 equal sub-intervals; the certificate is their conjunction.
12. **Exit, per box.**
    - τ_a is the smallest multiple of 1/64 at which the box minimum of h⁻²U = D/r⁴ (over the 16
      sub-intervals) is at least the box maximum of max(4|ω|², C_W). This is the box form of §3.1's
      rule.
    - C_W = 4·max(κ₊², W_bound), with W_bound = 99/12500 = 0.00792 ≥ max(−W̃) = 0.0079101…
      (`bulk_PB.json`).
13. **N3.**
    - It runs the continuation path in τ (affine Taylor steps, as §3.1 plans) on boxes
      Re ω ∈ [0, 1/20000]κ₋, Im ω = mode ± 1/40000 κ₋. The continuation runs against the recessive
      direction and amplifies a box's deviation about 30-fold.
    - A mode-free box at half the mode's depth is reported alongside.

**Controls and replay.**
14. **N4.** The design selected for ζ ∈ [0.3, 0.32]Ω_P, with X_a moved to x_T ∓ 1/(2k) (allowed
    region). It runs with ε-levels 0 and cell depth 2: q̂ < 0 on an interval near that exit, so
    refinement cannot change the verdict.
15. **K4.** 200 points with seed 20261018:
    - l uniform in [20, 400], the regime uniform in {L, S, T, A}, Re ω uniform in its range
      (A: [1.5, 3]Ω_P/h), and ν uniform in [−d_δ, 0];
    - mpmath at 20 digits;
    - a secant search (≤ 12 steps) seeded at each point.
16. **Replay.**
    - T and A: 5% of cells from the stored partitions, plus each box's worst cell, plus the tails.
    - S: every 20th ζ-leaf, recomputed from its stored design and target.
    - L: every 20th ω-leaf.
    - Verdicts and the stored bounds must agree exactly.

### Amendment 6 (2026-10-18, after the first K4 results; post hoc, disclosed as such)

**What "the float solver finds a QNM" means in K4.** The registered wording (§4; Amendment 5,
item 15) was a secant search seeded at each point, with no criterion for a found zero.
- **What happened.** After 14 of the 200 points, three secant runs had stopped inside H⁺, each within
  10⁻⁶ (relative) of its seed:
  - l = 122 (T), Re ω ≈ 1047κ₋;
  - l = 136 (A), Re ω ≈ 1502κ₋;
  - l = 253 (A), Re ω ≈ 2222κ₋.
- **They are not zeros.** At these points |W| equals its value at the seed to six digits (83.2271 and
  123.14, at both 20 and 40 digits). The rigorous zsol code gives |W| = 83.227148 ± 3.6·10⁻⁷ at the
  l = 122 point.
- **Why the stopping rule fires.** The secant stops when the relative step falls below
  10^{−dps+12}. At large |ω| the phase of W turns fast (d log W/dω is large), so the step f/W′ is
  tiny although W is not small.
- **The K4 evaluation from now on.** A secant stop in H⁺ counts as a QNM only if
  |W(root)| < 10⁻⁸·|W(seed)| (float, the run's settings) and the zsol ball at the root contains 0.
  - Every flagged stop is checked both ways.
  - The raw number of stops in H⁺ is reported alongside.
  - The K4 run itself (the points, the solver and its settings) is unchanged.

### Correction to Amendment 6 (2026-10-19, after the independent review, finding S2)

Amendment 6's text stands as written. This corrects it on two points.
- **A rule did exist at launch.** The registered *text* gave no criterion for a found zero, but the
  code did.
  - At `35201758`, and at `ca66270b`, from which K4 was launched, `assemble_K4` set
    `pass = … and not any(r["root_in_H_plus"])`. So every secant stop in H⁺ counted as a zero.
  - **Under that launch-time rule K4 fails**, with 61 of the 200 points stopped in H⁺. Amendment 6
    replaced the rule after the first results.
- **The mechanism is overshoot and return, not a fast-turning phase.** The review traced the secant
  at three of the stops.
  - d arg W/d(Re ω) at the seeds is only −3.35, −1.87 and −3.68.
  - The first secant step overshoots to Im ω ≈ −7.5, −13.0 and −6.9 κ₋, where |W| ≈ 10¹⁵ to 10⁴⁸.
  - The next step returns to the previous iterate to 16 digits.
  - The step after that is at most 3·10⁻¹⁶ relative, and the step-only test fires at the perturbed
    seed x₁ = x₀(1 + 10⁻⁶) + 10⁻⁹i. This is why every stop lies within 10⁻⁶ of its seed.
- **The new criterion stands.** A secant stop is a zero only if W is small there.
  - At all 61 stops the residual ratio |W(root)|/|W(seed)| is 1.0000008–1.0000013, and the
    rigorous zsol balls exclude 0. The reviewer's independent code gives D = 1 at 13 of them.
  - K4 is a float sanity check; the proof does not use it.

### Errata (2026-10-19, after the independent review, items W1 and W5)

- **W1 (Amendment 4).** "a* = 0.02996861191… (a* − a_d ≈ 3.7·10⁻¹¹)" is wrong in its last digits.
  - In Arb, a* = 2Ω_Pγ₀ = 0.0299686149113, and a_d = 32178552/2³⁰ = 0.0299686118960, so
    a* − a_d = 3.0·10⁻⁹. The 3.7·10⁻¹¹ came from the float pilot's finite-difference a*.
  - This is harmless: any a_d < a* is admissible. The best 30-bit dyadic below a* would have been
    32178555/2³⁰.
- **W5.** "K4" names two different things in this document: the float acceptance check (§4) and the
  kill criterion h₁ < 1/30.5 (§7). The S4 report calls the latter "kill criterion K4 (l₀ too large)".
