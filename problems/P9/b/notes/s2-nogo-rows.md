# P9(b) S2 — P9b-0 and the 8 no-alpha_B ladder rows

Status: S2 complete, 2026-08-31.  Conventions: FORMULATION.md v1.1 (frozen
GLV14 (86) action, family B, fiducial F*, D1-D3, and the rev-v1.1 pin: for
alpha_H rows "c_s^2" means the exact dust DE eigen-speed).  Everything below
builds on the S1 chain (`derivation/`, `notes/s1-derivation.md` §6); new code:
`src/p9b/fiducial.py`, `tests/test_s2_rows.py` (25 new tests; suite: 49).
Run: `uv run pytest problems/P9/b/tests -q`.

Notation: fhat ≡ 3 Omega_DE (1 + w_DE) = -2 Hdot/H^2 - rho_m/(M_Pl^2 H^2);
on F*, fhat = 3 Omega_DE (4a-3)/5, so sign(fhat) = sign(a - 3/4).  All row
targets are alpha_K c_s^2 with D = alpha_K (alpha_B = 0 => D = alpha_K, and
alpha_B^BS = 0 != 2).  aM = alpha_M, aH = alpha_H, aH' = d(alpha_H)/dx.

## 1. Verdict table (the 8 rows without alpha_B)

| row | exact target: alpha_K c_s^2 = ... | verdict | evidence |
|---|---|---|---|
| {K} | fhat | **(N) proven** (P9b-0) | I1 + sign lemma §2; 3 machine checks |
| {K,T} | fhat - 2 aT | **E-candidate** | aT = -1/10, J=[0.6,1]: min c_s^2 = 0.076 |
| {K,M} | 2 aM - 2 Hdot/H^2 - rho_m/(M^2 H^2) | **E-candidate** | aM = 1/2, J=[0.6,1]: min c_s^2 = 0.304 |
| {K,H} | 2 aH + 2 aH' - 2(1+aH) Hdot/H^2 - (1+2 aH) rho_m/(M_Pl^2 H^2) | **E-candidate** | aH = 1 const, J=[0.7,1]: min c_s^2 = 0.280 |
| {K,T,M} | {K,M} target - 2 aT | **E-candidate** | joint witness: min c_s^2 = 0.504 |
| {K,T,H} | {K,H} target - 2 aT | **E-candidate** | joint witness: min c_s^2 = 0.480 |
| {K,M,H} | 2 aH + 2 aH' + 2(1+aH)(aM - Hdot/H^2) - (1+2 aH) rho_m/(M^2 H^2) | **E-candidate** | joint witness: min c_s^2 = 1.338 |
| {K,T,M,H} | {K,M,H} target - 2 aT | **E-candidate** | joint witness: min c_s^2 = 1.538 |

All witnesses: alpha_K = 1 (D = 1), crossing a_c = 3/4 interior to J,
fhat sign change transversal, rho_DE > 0, c_T^2 >= 0.9, M*^2 > 0, all
alphas far inside the frozen D3 bounds.  Every target expression is an
exact sympy test (`test_master_target_matches_derived_rows`, all 8 rows,
against `rows.cs2_times_aK_row` from the derived exact dust eigen-speed);
every witness a numeric test (`test_near_witness`).

**Classification impact (honest flag).**  If S3 certifies the {K,T}, {K,M},
{K,H} witnesses, each is a MINIMAL set admitting stable crossing (their only
proper subset {K} is (N)): the expected "{K,B} unique minimal" FAILS over the
frozen ladder.  This is the discovery branch FORMULATION §2(b1) flag (i)
anticipated for aM/aH; it extends to aT.  Physics reading: each rescue
sources zeta-gradient energy without braiding — aT via m4^2 (extrinsic
curvature squared), aM via running M*^2, aH via the delta g^00 delta R^(3)
operator — exactly the "higher/extrinsic-curvature stabilisation near
crossing" mechanism of digest D6 (Creminelli et al. 2008).

## 2. P9b-0 (row {K}): the no-go, complete

**Identity I1 (S1, exact).**  For S = {alpha_K} on family B with dust,
derived from the frozen action (tadpole eqs (88)-(89) -> reduction -> speed):

    alpha_K c_s^2 = 3 Omega_DE (1 + w_DE) = f / (M_Pl^2 H^2),   f = rho_DE + p_DE.

Test: `test_reduction.py::test_I1_exact` (plus the k-essence dictionary
test: I1 <=> Vikman's c_s^2 = P_X/(P_X + 2X P_XX)).

**Lemma b0 (sign).**  Let J be compact, H > 0 on J (family B), alpha_K > 0
on J (D2(i)/D3), c_s^2 >= 0 on J (D2(ii)), and let I1 hold.  Then w_DE does
not cross -1 at any t_c in int J.  *Proof.*  D1 sign change gives a phantom-
side point t- in J with f(t-) < 0.  At t-, I1 gives alpha_K c_s^2 =
f/(M_Pl^2 H^2) < 0 since M_Pl^2 H^2 > 0.  But alpha_K > 0 and c_s^2 >= 0
give alpha_K c_s^2 >= 0 — contradiction.  QED (three lines).

Auxiliary hypotheses: NONE beyond D1 + family B + D2(i)-(ii).  Neither
transversality, nor rho_DE(t_c) > 0, nor the frozen bound alpha_K <= 10^3 is
used: the contradiction holds for every alpha_K > 0 (unbounded alpha_K only
sends c_s^2 -> 0^- on the phantom side, never to 0).  The frozen bound's
role is purely to keep the D6 degenerate boundary (alpha_K -> inf,
c_s^2 -> 0) outside the ladder, as FORMULATION §1.4 D3 states.

**Machine checks of the hypothesis structure** (`test_s2_rows.py`):
1. `test_p9b0_hypotheses_inconsistent_assumption_engine`: with sympy
   assumptions alpha_K > 0, c_s^2 >= 0, fhat < 0, the equation
   Eq(alpha_K*c_s^2, fhat) auto-evaluates to `False` (sympy proves
   (alpha_K c_s^2 - fhat).is_positive is True).
2. `test_p9b0_sign_decomposition_certificate`: the obstruction written as a
   sum of manifestly signed terms, alpha_K c_s^2 - fhat =
   [alpha_K c_s^2 (>= 0)] + [-fhat (> 0)] > 0, each sign machine-verified,
   the decomposition exact — but I1 forces the sum to be 0.
3. `test_p9b0_phantom_side_forces_negative_cs2_fractions`: exact Fraction
   sweep incl. the boundary alpha_K = 10^3: c_s^2 = fhat/alpha_K < 0 on the
   phantom side, i.e. f < 0 forces alpha_K c_s^2 < 0.

Escapes (recorded with the claim, as in FORMULATION §2(b0)): ghost phase
(alpha_K < 0; quintom-A/B territory, out of scope) or the degenerate
c_s^2 -> 0 boundary (excluded by the frozen bounds).  Verdict: **(N)**,
CERTIFIED-grade evidence complete (exact identity chain + stated lemma +
machine checks).  Draft ledger rows in §6.

## 3. Per-row analysis (adversarial, both directions)

Search discipline used before each verdict: for (N) — look for bounded-alpha
escapes, the c_s^2 = 0 boundary, kinetic-factor sign flips; for (E) — check
the FULL frozen D2 (incl. c_T^2 >= 0, M*^2 > 0, the v1.1 eigen-speed pin)
and the D3 bounds on the witness.  D = alpha_K = 1 > 0 in every witness, so
the kinetic factor never flips sign; c_s^2 = target exactly.

**{K,T}: E-candidate.**  Target fhat - 2 aT; the S1 "-2 aT rescue budget" is
real: c_T^2 = 1 + aT >= 0 allows aT down to -1, a budget of +2 against
fhat's dip, and fhat -> 0 at the crossing, so any J on which fhat > -2
(automatic near t_c by continuity) is rescuable.  Witness: aT = -1/10,
alpha_K = 1 on J = [0.6, 1] (min fhat = -0.124 at a = 0.6): c_s^2 =
fhat + 1/5 in [0.076, 0.62], c_T^2 = 0.9, all bounds slack.  Why (N) is
impossible here: an (N) identity would need fhat <= -2 forced somewhere on
every crossing background — false on F* itself.  No extra hypothesis of the
D1/D2 kind restores it.  S3 target: Arb sweep of this exact witness
(already rational; background exp(12/5(1-a)) is Arb-elementary).

**{K,M}: E-candidate** (was "genuinely open"; resolved E).  Target
2 aM - 2 Hdot/H^2 - rho_m/(M^2 H^2) = fhat + 2 aM + (rho_m/H^2)(1/M_Pl^2 -
1/M^2).  Rescue: constant aM > 0.  Adversarial cost checked: normalising
M^2(a=1) = M_Pl^2 (t0 = today, FORMULATION 1.3) makes M^2 = M_Pl^2 a^{aM} <
M_Pl^2 in the past, so the matter term PENALISES the rescue by
3 Omega_m (1 - a^{-aM}) < 0, which grows toward matter domination and kills
naive witnesses at a <= 0.5 (e.g. aM = 1/2 fails at a = 0.5: margin -0.10).
On J = [0.6, 1] with aM = 1/2 the balance is positive: min c_s^2 = 0.304 at
a = 0.6, M^2/M_Pl^2 in [0.775, 1].  Why no (N): the offset 2 aM (budget
+-20) dominates fhat near the crossing; no bounded-alpha obstruction exists.

**{K,H}: E-candidate** (was "genuinely open"; resolved E).  Target (with the
v1.1 EXACT eigen-speed, i.e. including the kinetic-matter-mixing term
+ aH^2 rho_m/(M^2 H^2 alpha) — this is what `cs2_times_aK_row(aH_on=True)`
is built from) equals fhat + 2 aH B + 2 aH', with the background bracket
B ≡ 1 - Hdot/H^2 - rho_m/(M_Pl^2 H^2).  On F*, B > 0 for a >~ 0.65 (B =
0.275 at a_c, 0.76 at a = 1; B(0.6) = -0.046) but B < 0 toward matter
domination (B -> -1/2):
a constant-aH witness needs J inside the B > 0 region.  Two witnesses:
primary aH = 1 constant on J = [0.7, 1] (min c_s^2 = 0.280; uses NO
derivative freedom); variant aH = ln(a/0.6) on J = [0.6, 1] (min c_s^2 =
1.876), showing the 2 aH' = 2 d(aH)/dx term (C^1 freedom, no frozen
derivative bound) is an independent, effectively unbounded rescue.  Why no
(N): both escapes are inside the frozen D3 bounds; an (N) would need
sup_J (2 aH B + 2 aH') < -fhat forced, false on F*.  Note the witness is
checked against the pinned eigen-speed target, not GLV (85).

**Combined rows {K,T,M}, {K,T,H}, {K,M,H}, {K,T,M,H}: E-candidates.**
Two independent grounds: (i) *monotonicity/inheritance*: STABLE-CROSSING(S)
is monotone under S ⊆ S' among no-alpha_B rows, because alpha_i ≡ 0 is an
admissible bounded free function for every i != K, and the pinned eigen-speed
of S' at the extra alphas ≡ 0 equals the pinned eigen-speed of S exactly
(`test_row_monotonicity_substitution`, incl. the aH -> 0 and M^2 -> M_Pl^2
degenerations); so any {K,T}/{K,M}/{K,H} witness embeds.  (ii) *joint
witnesses* with every mechanism of the row switched on (aT = -1/10,
aM = 1/2, aH = 1 on the rows' J; table §1), verifying no destructive
interference — margins ADD, exactly as the offset identities predict
(`test_combined_row_offsets_exact`).

## 4. What blocks the other direction (per-row, for the record)

- (N) for any row except {K} is blocked by explicit in-bounds escapes on F*
  itself (§3): no identity quantified over family B can hold when a family-B
  member + admissible alphas violate it numerically by O(0.1..2) margins.
- (E) at CERTIFIED level is blocked only by validated numerics not yet run:
  the S2 witnesses are float-level (4001-point grids, margins >= 0.076,
  no interval arithmetic).  That is precisely S3's job (FORMULATION §3.5),
  not an open mathematical question.  No row is UNDECIDED.

## 5. What S3 needs

1. Arb (python-flint) interval sweeps along F* for the three minimal-row
   witnesses {K,T}, {K,M}, {K,H} (witnesses already rational; background
   needs exp/log enclosures), certifying D >= 1, c_s^2 >= c_min (suggest
   c_min: 1/20, 1/4, 1/4), c_T^2 >= 9/10, M*^2 >= 3/4 on J, and certified
   strict fhat signs on stated subintervals left/right of a_c = 3/4 (IVT).
2. The 8 alpha_B rows (S3 also owns {K,B} etc. existence per the plan).
3. Combined no-alpha_B rows: cite inheritance (no new certificates needed),
   or optionally re-run the sweep on the joint witnesses.
4. P9b-1 assembly: with S3 certificates, minimal sets = {K,T}, {K,M}, {K,H}
   (+ {K,B} if certified); the "unique minimal" expectation is refuted.

## 6. Draft ledger rows (do NOT edit CLAIMS.md from S2; main session owns it)

| P9b-0 | P9(b) | Frozen EFT (GLV14 (86)), S = {alpha_K}, family B (dust): the exact identity I1: alpha_K c_s^2 = 3 Omega_DE (1+w_DE) holds, and consequently no pair (alpha_K, background) with alpha_K > 0 and c_s^2 >= 0 on J admits a w_DE = -1 crossing (D1) — on the phantom side f < 0 forces alpha_K c_s^2 < 0. Holds for ALL alpha_K > 0 (the frozen bound alpha_K <= 10^3 is not needed); auxiliary hypotheses: none beyond D1, family B, D2(i)-(ii). Escapes recorded: ghost phase or the degenerate c_s^2 -> 0 boundary (outside D2/D3) | CERTIFIED (exact sympy identity chain action -> (88)-(89) -> reduction -> I1, `test_I1_exact`; sign lemma stated in notes/s2-nogo-rows.md §2 and machine-checked three ways in `tests/test_s2_rows.py::test_p9b0_*`) | problems/P9/b/derivation/rows.py; problems/P9/b/tests/test_s2_rows.py; problems/P9/b/notes/s2-nogo-rows.md | 2026-08-31 |
| P9b-1.K | P9(b) | Ladder row S = {alpha_K}: (N) — no stable crossing on any family-B background with the frozen bounds; by P9b-0 (which proves the stronger unbounded-alpha_K statement) | CERTIFIED (= P9b-0) | as P9b-0 | 2026-08-31 |

Suggested CONJECTURED rows for the 7 E-candidates (S3 upgrades to CERTIFIED
(E) or downgrades honestly): statement pattern "Ladder row S: (E)-candidate —
float near-witness on F* with alpha_K = 1, [row alphas], J = [a1, 1],
min c_s^2 = [value], D = 1, c_T^2 >= 9/10, M*^2 > 0, all D3 bounds
satisfied; exact target identity tested; Arb certification pending (S3)",
evidence `src/p9b/fiducial.py` + `tests/test_s2_rows.py::test_near_witness`.

## 7. Test inventory added in S2 (all in tests/test_s2_rows.py)

- `test_master_target_matches_derived_rows[K..KTMH]` (8): the master
  expression driving all numerics == derived alpha_K c_s^2, exact, per row.
- `test_combined_row_offsets_exact`: combined rows = pure-row offsets, exact.
- `test_row_monotonicity_substitution` (5): inheritance-lemma algebra, exact.
- `test_p9b0_*` (3): the P9b-0 machine checks (§2).
- `test_near_witness[KT..KTMH]` (8, incl. KH slope variant): D1/D2/D3
  near-witness checks + recorded-margin regression pins.
Suite: 49 passed (24 S1 + 25 S2), ~17 s.
