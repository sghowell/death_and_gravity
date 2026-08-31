# P9(b) — Minimal EFT-of-DE operator content for a stable phantom crossing

Status: formulation v1 (scoping), 2026-08-31. Nothing below is certified yet;
every claim id is a ledger *target*. Changes after the first derivation lands
go in the Revision log, never by silent edit (discipline of `../FORMULATION.md`).
Pre-registered statement: `docs/problems/open-problems-theoretical-cosmology-2026.tex`
§P9(b). Literature basis: `notes/literature-digest.md` (D-numbers cited below).

## 1. Frozen conventions

### 1.1 Units, signature, background
c = 1; reduced Planck mass M_Pl; signature (−,+,+,+); flat FLRW
ds² = −dt² + a²(t)δ_ij dxⁱdxʲ; H = ȧ/a; overdot = d/dt; x = ln a.
Matter sector: pressureless dust (ρ_m > 0, p_m = 0), minimally and
universally coupled. Radiation may be added in S1 if the fiducial needs it;
if so, by revision log. No other matter.

### 1.2 The EFT action (the frozen object)
Unitary gauge, single scalar degree of freedom, universal coupling. The
frozen quadratic action is the ADM α-basis action of
Gleyzes–Langlois–Vernizzi, arXiv:1411.3712 (GLV14), their eq. (86) — the
five-function form in {δN, δK_ij, δR^(3)} with coefficients built from
{M*², α_K, α_B, α_T, α_M, α_H} — together with their α definitions
(α_H: eq. (82)); background equations re-derived in S1 from the same action.
Reporting convention for the braiding is Bellini–Sawicki, arXiv:1404.3713
(BS14): α_B^BS = −2 α_B^GLV (map verified as a derivation identity in S1;
any mismatch is a formulation bug fixed by revision log). α_M ≡ d ln M*²/dx
(BS14 eq. (3.4)). Everything downstream (D, c_s², all no-go identities) is
*derived* from the frozen action by exact sympy in `derivation/` (house style
of `problems/P4/derivation`-equivalent: every step an exact identity test);
published formulas BS14 (3.12)–(3.13) and GLV14 (83), (85) are known-answer
tests, not inputs.

GPV operator dictionary (Gubitosi–Piazza–Vernizzi, arXiv:1210.0201): the
unitary-gauge operators map as M₂⁴ ↔ α_K, m₃³ ↔ α_B, m̄₄²/m̃₄² ↔ α_T/α_H,
f(t) ↔ α_M. The exact mapping table (with signs and H-factors) is an S1
deliverable, verified by sympy; it is documentation, not an input.

### 1.3 Frozen background family and fiducial
Family B: flat FLRW with dust as in §1.1, H ∈ C³(J), H > 0 on a compact
interval J = [t₁, t₂] (equivalently in x). The effective DE fluid is defined
by subtraction at fixed gravitational coupling:

    ρ_DE ≡ 3 M_Pl² H² − ρ_m,      p_DE ≡ −M_Pl²(2Ḣ + 3H²) − p_m,
    w_DE ≡ p_DE/ρ_DE,             f ≡ ρ_DE + p_DE = −2M_Pl²Ḣ − (ρ_m + p_m),

with M_Pl² ≡ M*²(t₀) when α_M is switched on. This is the standard
(DESI-compatible) convention: all modification-of-gravity terms are lumped
into the DE fluid. No-go claims quantify over all of B; existence claims use
the frozen fiducial

    F*: CPL w_DE(a) = w₀ + w_a(1−a), (w₀, w_a) = (−4/5, −4/5), Ω_m0 = 3/10,

round rationals inside the DESI-preferred quadrant (w₀ > −1, w_a < 0), with
crossing at a_c = 3/4 (z_c = 1/3) exactly. DESI DR2 posterior values are
motivation only (§4); the fiducial is deliberately NOT a fitted number.

### 1.4 Definitions (frozen)
- **D1 (crossing).** w_DE crosses −1 at t_c ∈ int J if f ∈ C¹(J), f(t_c) = 0,
  f changes sign at t_c (∃ε > 0 with f strictly negative on one side interval
  and strictly positive on the other), and ρ_DE(t_c) > 0. Transversal
  crossing: additionally ḟ(t_c) ≠ 0. Theorems are stated for sign change;
  where a proof needs transversality this is stated in the claim.
- **D2 (stable on J).** (i) No ghost: D ≡ α_K + (3/2)(α_B^BS)² > 0 on J and
  α_B^BS ≠ 2 on J (the α_B = 2 degeneracy locus is excluded from the ladder);
  (ii) no gradient instability: c_s² ≥ 0 on J (weak inequality — this makes
  the no-gos strongest), where c_s² is the scalar sound speed derived from
  the frozen action; existence witnesses must satisfy the strict version
  c_s² ≥ c_min > 0 with an explicit rational c_min; (iii) tensor sector:
  c_T² = 1 + α_T ≥ 0 and M*² > 0 on J (active only when α_T or α_M is on).
  No tachyon condition is part of "stable": the effective-mass term is
  gauge- and scale-convention dependent; a tachyon diagnostic with threshold
  |μ²| ≤ 10 H² may be *reported* for witnesses but enters no claim.
- **D3 (operator content).** Ladder subsets S ⊆ A = {α_K, α_B, α_T, α_M, α_H}
  with α_K ∈ S always. "S switched on" means: the α's in S are free functions
  in C¹(J) with frozen bounds 0 < α_K ≤ 10³, |α_i| ≤ 10 (i ≠ K) on J; the
  α's not in S vanish identically; M*² ≡ M_Pl² unless α_M ∈ S. The
  boundedness bounds are part of every statement: they close the known
  degenerate escape α_K → ∞ ⟺ c_s² → 0 at crossing (D6, D8 in the digest);
  the boundary case is discussed, not smuggled in. S is *minimal admitting
  stable crossing* if STABLE-CROSSING(S) holds and fails for every proper
  subset. The ladder is the 16 subsets of A containing α_K.

## 2. Theorem ladder (claim templates for CLAIMS.md)

### (b0) k-essence no-go, exact-algebra form  [target: CERTIFIED]
**P9b-0.** In the frozen EFT with S = {α_K} the exact background identity

    I1:   α_K c_s² = 3 Ω_DE (1 + w_DE),        Ω_DE ≡ ρ_DE/(3M_Pl²H²),

holds (equivalently α_K c_s² H² M_Pl² = f). Consequently there is NO pair
(α_K, background) with α_K ∈ C⁰(J), 0 < α_K ≤ 10³, background ∈ B, w_DE
crossing −1 at t_c ∈ int J (D1), and D > 0, c_s² ≥ 0 on J (D2): on the
phantom side f < 0 forces α_K c_s² < 0, contradicting α_K > 0, c_s² ≥ 0.
This is the α-basis analogue of Vikman 2005 (digest D1) and of the quintom
no-go (D3): the machine part is I1 as an exact sympy identity from the
frozen action plus the k-essence dictionary check (c_s² = P_X/(P_X+2XP_XX)
maps to I1); the sign argument is three lines of stated analysis.
Evidence path: `derivation/` chain + `tests/`. Remark recorded with the
claim: the only escapes are α_K → ∞ (degenerate, excluded by the frozen
bound; the c_s² → 0 crossing of Creminelli et al. 2008, digest D6) or a
ghost phase (quintom-A/B two-field route, out of scope §4).

### (b1) Classification over the operator ladder  [target: CERTIFIED per row]
For each of the 16 ladder rows S, exactly one of:
- **(E) Existence certificate.** Explicit witness α-functions (rational in a)
  on the fiducial F*, with a validated-numerics certificate (python-flint
  Arb, as in `problems/P4/src/p4/validated/`) that on J: D ≥ D_min > 0,
  c_s² ≥ c_min > 0, c_T² ≥ 0, M*² > 0 (interval enclosures along the
  trajectory, outward-rounded), and that f changes sign (certified strict
  signs on stated subintervals left and right of a_c; crossing by IVT).
- **(N) No-go identity.** An exact sympy identity in the frozen conventions
  showing D2 ∧ D1 is contradictory for every background in B and all bounded
  α-functions in S, in the style of I1 (auxiliary hypotheses, e.g.
  transversality or Ω_DE(t_c) > 0, stated per row).

**P9b-1.S** (one ledger row per subset S, ids P9b-1.K, P9b-1.KB,
P9b-1.KT, …): the (E) or (N) statement for that row.
**P9b-1** (the classification): "Over the frozen ladder, bounds, fiducial
F* (existence) and family B (no-gos), the minimal sets admitting stable
crossing are exactly {…}", assembled from the 16 rows; its level is the
minimum level over rows. Expected but not presupposed: {α_K, α_B} is the
unique minimal set (KGB, digest D2, D8); the 8 rows without α_B are the
no-go candidates. Honest flags: (i) rows with α_M or α_H but no α_B are
genuinely open — a stable crossing there would be a discovery, and the
no-go may need extra stated hypotheses; (ii) if for some row neither an (E)
witness nor an (N) identity is found, the row ships as CONJECTURED with the
numerical evidence and the classification P9b-1 is stated over the decided
rows only, saying so.

### (b2) Positivity intersection  [target: CERTIFIED conditional on POS]
There is no theorem transporting Adams et al. flat-space positivity to
FLRW/dS; the digest (§4) documents the gap (boost breaking, massless scalar,
no S-matrix). (b2) is therefore *conditional by construction*:
- **Assumption set POS (frozen).** The Melville–Noller shift-symmetric
  Horndeski bounds (arXiv:1904.05874, eqs. (16)–(17); example-model α-form
  their eq. (10): α_B ≤ 2α_T/(1+α_T), sign convention re-checked in S1),
  applied at each t ∈ J under their stated flat-space-transfer assumption.
  POS is an assumption set, named in every statement; it is not claimed to
  follow from UV physics on this background.
- **P9b-2 (weak, trajectory-level).** For each (b1) existence witness: a
  certified evaluation (Arb interval sweep; the witness is rational in a, so
  each POS inequality is a polynomial sign condition on J) of whether the
  witness satisfies POS at every t ∈ J. Statement shipped either way
  (satisfies / violates, with the violating subinterval certified).
- **P9b-2s (strong, region-level; stretch).** Freeze a finite-dimensional
  parametrization (α's polynomial in a of degree ≤ d with rational
  coefficient boxes); decide whether {stable crossing} ∩ POS = ∅ by exact
  semialgebraic certificate (SOS/Positivstellensatz with rational rounding,
  or an explicit certified point). If this is not honestly reachable at the
  frozen degree, P9b-2s is dropped and the ledger says so; P9b-2 stands.

## 3. Deliverables, ledger plan, DESI touchpoint

### 3.1 Deliverables
- `derivation/`: sympy chain (frozen action → background eqs → quadratic
  action → D, c_s² → I1 and per-row identities), each step an exact test.
- `src/p9b/`: witness search (float proposer) + Arb validators (verifier).
- `certificates/`, `results/`: per-row witness/identity artifacts (JSON) and
  the classification table; `tests/`: known-answer + identity tests.

### 3.2 Known-answer tests (gate: nothing above HEURISTIC until these pass)
1. Re-derivation reproduces BS14 (3.12)–(3.13) and GLV14 (83), (85) exactly
   (up to the verified convention map).
2. k-essence dictionary: I1 ⟺ Vikman's P(φ,X) statement under the map.
3. The DPSV 2010 KGB example (digest D2) reproduces a stable crossing
   numerically with our derived D, c_s².
4. Convention check: f(R) limit gives α_B^BS = −α_M, α_K = 0 (digest D10).

### 3.3 Ledger levels (mirroring part (a))
- CONJECTURED: float pipeline finds the row outcome (witness or numeric
  no-go evidence).
- VERIFIED_N: independent re-derivation/implementation path agrees (second
  CAS route for identities; independent ODE/interval code for witnesses).
- CERTIFIED: exact sympy identity chain (no-gos) / Arb certificate with
  outward rounding (witnesses); certificate files committed.
b0: CERTIFIED target. b1: CERTIFIED per decided row. b2: CERTIFIED
conditional on POS (the condition named in the ledger row itself).

### 3.4 DESI DR2 touchpoint (explicit)
The DESI DR2 w₀–w_a posterior (arXiv:2503.14738; 3.1σ vs CMB alone,
2.8–4.2σ with SN compilations) enters (b) ONLY as motivation for the
fiducial quadrant and crossing redshift. No likelihood is computed, no data
file is frozen, no posterior contour is used in any claim. Part (a) owns the
frozen-data stack; part (b) is theory-side algebra + certificates.

### 3.5 Stages (sized like part (a)'s S0–S3)
- S0 (this doc): formulation freeze; digest committed.
- S1: derivation chain + known-answer tests + convention/dictionary pins
  (the equation-number pins of §1.2 re-verified against the papers).
  → gate for everything.
- S2: P9b-0 and the (N) rows (exact identities). → first CERTIFIED rows.
- S3: (E) rows — witness search + Arb validation on F*. → classification
  table complete or honestly partial.
- S4: P9b-2 (and P9b-2s if reachable); assembly of P9b-1; summary.
- S5: writeup + adversarial review (same protocol as (a): review before the
  ledger row is called final; fixes by re-run, supersession logged).

## 4. What is NOT claimed (scope fences)
- No nonlinear completion, no UV model, no claim any witness admits a
  standard UV completion (that is exactly what POS conditionalizes).
- No screening (Vainshtein or otherwise), no nonlinear/quasi-static
  phenomenology, no G_eff/ISW/growth claims, no observational fitting.
- No multi-field: quintom two-field crossings (ghost pair) are out of scope;
  the classification is over the single-scalar unitary-gauge EFT only.
- No DHOST beyond the α_H toggle; no explicit breaking of spatial diffs; no
  matter beyond §1.1; no claim about backgrounds outside B or bounds outside
  D3; no tachyon/nonperturbative stability claims.
- No claim that flat-space positivity holds on FLRW (b2 is conditional); no
  claim about which model nature realizes; no statement about the DESI
  data dispute (part (b) is well posed regardless, per the pre-registration).

## Revision log
- v1 (2026-08-31): initial scoping formulation. Frozen: GLV14 eq. (86)
  action with BS14 α_B reporting convention; family B; fiducial F* =
  CPL(−4/5, −4/5), Ω_m0 = 3/10; ladder = 16 subsets ⊇ {α_K}; bounds
  α_K ≤ 10³, |α_i| ≤ 10; POS = Melville–Noller (16)–(17)/(10). Equation
  pins recorded from web fetches on 2026-08-31 and re-verified in S1.
- v1.1 (2026-08-31, S1 close-out; `notes/s1-derivation.md`, `derivation/`,
  `tests/`). (a) All §1.2 equation pins re-verified verbatim against ar5iv
  fetches (GLV14 (68), (73)–(89); BS14 (3.4), (3.12)–(3.13); DPSV (40)–(41),
  (59)); the map α_B^BS = −2 α_B^GLV verified as a derivation identity.
  (b) Naming clarification: the frozen operator names (M₂⁴, m₃³, m̄₄², m̃₄²)
  are those of GLV14 eq. (87)/GLPV 1304.4840; arXiv:1210.0201's own names are
  m̄₁³ = m₃³, M̄₂² = −M̄₃² = 2m̄₄², μ₁² = m̃₄²/2 (dictionary content
  unchanged, sympy-verified). (c) Pin sharpening for D2(ii), by exact
  derivation, no silent adaptation: with dust and α_H ≠ 0 the GLV14 (85)
  formula is not the UV eigen-speed of the coupled (ζ, δσ_m) system; the
  frozen meaning of "c_s² derived from the frozen action" is henceforth the
  exact DE eigen-speed, which equals (85) + α_H² ρ_m/(M² H² α) for dust
  (kinetic matter mixing, cf. arXiv:1609.01272) and coincides with
  BS14 (3.13) = GLV14 (85) whenever α_H = 0 (all ladder rows without α_H).
  The §3.2(1) known-answer gate passes in that form: (3.13) exact at
  α_H = 0, and the exact product rule c₁²c₂² = c_m²·c_s²(85) for all five
  α's. (d) §3.2 gate status: all four known-answer tests pass; S2 is GO.
