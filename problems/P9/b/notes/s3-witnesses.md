# P9(b) S3 — Arb-certified existence witnesses for the 15 rows S != {alpha_K}

Status: S3 complete, 2026-08-31.  Conventions: FORMULATION.md v1.1 (frozen
GLV14 (86) action, family B, fiducial F*, D1-D3; the rev-v1.1 pin applies:
c_s^2 = the exact dust DE eigen-speed, incl. kinetic matter mixing for the
alpha_H rows).  Builds on S2 (`notes/s2-nogo-rows.md`).  New code:
`src/p9b/witnesses.py` (frozen witnesses + braiding master target),
`src/p9b/certify.py` (Arb certifier), `certificates/P9b-1.<ROW>.json` (15),
`tests/test_s3_witnesses.py` (62 new tests; suite: 111, ~17 s).
Run: `uv run pytest problems/P9/b/tests -q`; re-generate certificates with
`PYTHONPATH=problems/P9/b/src:problems/P9/b uv run python -m p9b.certify`.

## 1. Method (the (b1)(E) template, executed)

Single master expression `MASTER_TARGET_B` (witnesses.py) = D c_s^2 for every
row, in scalars (BS braiding b = alpha_B^BS, ' = d/dx, g = 1 - b/2):

    D c_s^2 = -2 g^2 (1+aT) + 2 g (1+aH)(1+aM-hd) + 2 g aH' + (1+aH) b'
              - (1+2 aH) rm,      hd = Hdot/H^2,  rm = rho_m/(M^2 H^2),

matched EXACTLY (sympy) to the derived exact dust eigen-speed
`rows.cs2_times_aK_row` for all 8 braiding rows (`test_master_b_matches_
derived_rows`; the 8 alpha_B = 0 instances were matched in S2), and reducing
exactly to S2's master at b = 0.  It is a dyadic polynomial, so the SAME
lambdified function evaluates on floats (proposer) and on arb balls
(verifier); all transcendentals sit in the background enclosures.

Background F* (documented in every certificate; H0 = 1, M_Pl^2 = 1; exact
rationals -> balls, outward rounding throughout, prec = 128 bits):
1+w = (4a-3)/5; g_de = exp((9/5) ln a + (12/5)(1-a)); E2 = (3/10)a^-3 +
(7/10)g_de; hd = -(3/2)[Om + (1+w)Ode]; fhat = 3 Ode (1+w); m2 = M^2/M_Pl^2
= exp(aM ln a); rm = 3 Om/m2.  D = 1 + (3/2)b^2 (alpha_K = 1); c_s^2 =
(D c_s^2)/D by ball division.  Sweep: adaptive bisection of J with rational
endpoints until EVERY subinterval's ball verdict is strict; subdivision
counts and rigorous endpoints stored (house style: P4 validated/, P9
verify.py::_endpoint outward float conversion).

Per row, the certified checks (all strict ball verdicts):
(i) c_s^2 > c_min on J; (ii) D > D_min (exact D = 1 for b = 0 rows;
sweep for braiding rows) and alpha_B^BS < 2 on J (degeneracy locus);
(iii) c_T^2 = 1 + aT >= 0 exact; M*^2 = m2 > 0 on J (aM rows);
(iv) rho_DE > 0 (Ode > 0) and H^2 > 0 (E2 > 0) on J; (v) IVT crossing:
fhat < 0 on [a1, 3/4 - 1/50], fhat > 0 on [3/4 + 1/50, 1]; with H^2 > 0,
sign f = sign fhat, so f changes sign in the gap (and f(3/4) = 0 exactly on
F*: fhat = (12/5) Ode (a - 3/4); transversality holds analytically but is
not needed — D1 is the sign-change form).

## 2. Frozen witnesses and certified margins (all CERTIFIED (E))

alpha_K = 1 everywhere; constants aT = -1/10, aM = 1/2, aH = 1 where on;
braiding alpha_B^BS(a) = a/2 (rational in a; magnitude from the S1 KGB check,
aB_BS ~ 0.5; slope term b' = a/2).  M^2 = M_Pl^2 a^(1/2) for aM rows.
J = [3/5, 1] (no alpha_H) or [7/10, 1] (alpha_H rows); a_c = 3/4 interior.
"lo(c_s^2)" = min over accepted subintervals of the rigorous lower ball
endpoint (a certified lower bound for c_s^2 on J); n = subintervals/depth of
the c_s^2 sweep; lo(D) likewise (braiding rows; D_min = 11/10).

| row     | J        | c_min | lo(c_s^2)  | n     | lo(D)   |
|---------|----------|-------|------------|-------|---------|
| {K,B}   | [3/5,1]  | 1/10  | +0.101969  | 41/8  | 1.1050  |
| {K,T}   | [3/5,1]  | 1/20  | +0.050112  | 78/9  | 1 exact |
| {K,M}   | [3/5,1]  | 1/4   | +0.255325  | 27/8  | 1 exact |
| {K,H}   | [7/10,1] | 1/4   | +0.251444  | 33/9  | 1 exact |
| {K,B,T} | [3/5,1]  | 1/5   | +0.205044  | 32/7  | 1.1050  |
| {K,B,M} | [3/5,1]  | 7/20  | +0.354089  | 28/9  | 1.1050  |
| {K,B,H} | [7/10,1] | 3/10  | +0.300733  | 34/9  | 1.1669  |
| {K,T,M} | [3/5,1]  | 2/5   | +0.410424  | 21/7  | 1 exact |
| {K,T,H} | [7/10,1] | 2/5   | +0.400706  | 26/8  | 1 exact |
| {K,M,H} | [7/10,1] | 5/4   | +1.272565  | 19/8  | 1 exact |
| {K,B,T,M} | [3/5,1]  | 2/5 | +0.404579  | 18/7  | 1.1050  |
| {K,B,T,H} | [7/10,1] | 2/5 | +0.401453  | 30/9  | 1.1669  |
| {K,B,M,H} | [7/10,1] | 4/5 | +0.814317  | 16/7  | 1.1669  |
| {K,T,M,H} | [7/10,1] | 7/5 | +1.416927  | 16/7  | 1 exact |
| {K,B,T,M,H} | [7/10,1] | 1 | +1.001148  | 22/9  | 1.1669  |

fhat sign sweeps (shared background): max upper endpoint on the left
subinterval -0.01389 (J from 3/5) / -0.01504 (from 7/10); min lower endpoint
on the right subinterval +0.000583.  c_T^2 = 9/10 exact (aT rows), 1
otherwise; m2 in [(3/5)^(1/2), 1], certified > 0.  D3: alpha_K = 1 in
(0, 10^3]; max |alpha_i| = 1 <= 10; |b| <= 1/2 < 2 certified (degeneracy
locus avoided).  All 15 witnesses are genuinely JOINT (every alpha of the
row switched on and nonzero on J) — no certificate leans on inheritance.

## 3. Inheritance lemma (second route, not load-bearing)

Every row was certified directly, so the classification does not use the S2
combined-row inheritance lemma.  It is nevertheless now fully machine-checked
as a second route: setting an extra alpha == 0 (alpha_i ≡ 0 is admissible
under D3 — |0| <= 10, C^1, and alpha_B = 0 is off the alpha_B = 2 locus; for
alpha_M, M^2 ≡ M_Pl^2) maps the superset row's derived target exactly onto
the subset row's (`test_inheritance_lemma_b_rows`, 8 braiding cases, chaining
with S2's `test_row_monotonicity_substitution` to reach every row from a
two-operator row).  So the four two-operator certificates alone would already
decide the other 11 rows; the direct joint certificates are strictly stronger
(they also witness no destructive interference).

## 4. P9b-1 classification (all 16 rows decided)

Over the frozen ladder (16 subsets containing alpha_K), frozen bounds (D3),
family B for the no-go and fiducial F* for existence: row {K} is (N)
(P9b-0, S2, CERTIFIED); the other 15 rows are (E), each with an Arb-certified
witness (§2).  Since every proper ladder subset of a two-operator row is {K},

    **the minimal sets admitting stable crossing are exactly
    {K,B}, {K,T}, {K,M}, {K,H}.**

The pre-registered expectation "{K,B} is the unique minimal set" is REFUTED
over the frozen ladder: braiding is sufficient but not necessary; alpha_T,
alpha_M, alpha_H each rescue the crossing without braiding (the discovery
branch FORMULATION §2(b1) flag (i) anticipated, extended to alpha_T; physics
reading in s2-nogo-rows.md §1).

## 5. Draft CLAIMS rows (main session owns CLAIMS.md; do not edit from S3)

Per-row statement pattern (15 rows, ids P9b-1.KB ... P9b-1.KBTMH):

| P9b-1.<ROW> | P9(b) | Ladder row S = <SET>: (E) — explicit witness on F* (alpha_K = 1, <ALPHAS>, J = <J>) admits a stable phantom crossing at a_c = 3/4: Arb-certified (128-bit, outward, adaptive bisection) c_s^2 (exact dust eigen-speed, rev-v1.1 pin) >= <C_MIN> and D >= <D_MIN> on J, alpha_B^BS != 2, c_T^2 >= 0, M*^2 > 0, rho_DE > 0, and f = rho_DE + p_DE certified strictly negative/positive left/right of a_c (IVT); all D3 bounds hold | CERTIFIED (certificates/P9b-1.<ROW>.json; tests/test_s3_witnesses.py) | problems/P9/b/src/p9b/{witnesses,certify}.py; problems/P9/b/certificates/P9b-1.<ROW>.json | 2026-08-31 |

with <SET>, <ALPHAS>, <J>, <C_MIN>, <D_MIN> from §2 (e.g. P9b-1.KB: S =
{alpha_K, alpha_B}, alpha_B^BS = a/2, J = [3/5, 1], c_min = 1/10,
D_min = 11/10).  Assembled classification row:

| P9b-1 | P9(b) | Classification over the frozen ladder: all 16 rows decided — {alpha_K} is (N) (P9b-0); every other row is (E) with an Arb-certified witness on F*; the minimal sets admitting stable crossing are exactly {alpha_K,alpha_B}, {alpha_K,alpha_T}, {alpha_K,alpha_M}, {alpha_K,alpha_H}; the "unique minimal {alpha_K,alpha_B}" expectation is refuted | CERTIFIED (minimum level over the 16 rows; per-row evidence above) | problems/P9/b/notes/s3-witnesses.md §4; certificates/ | 2026-08-31 |

## 6. S4 handoff (P9b-2: POS evaluation)

All 15 (E) witnesses feed P9b-2: each is rational-in-a alphas on J with the
documented exp/power background, so every Melville–Noller inequality
(1904.05874 (16)-(17); example alpha-form (10): alpha_B <= 2 alpha_T/(1+alpha_T),
sign convention re-checked per S1) is an Arb-sweepable sign condition on J.
Expected shape (S4 must certify, statement shipped either way): braiding
witnesses have alpha_B^BS = a/2 > 0 with alpha_T in {0, -1/10}, so the
example-form bound (<= 0 resp. -2/9) looks violated; {K,M}/{K,H}-type
witnesses sit at alpha_B = alpha_T = 0 (boundary).  The witness JSONs carry
everything needed (exact witness, J, background formulas).  P9b-2s (region
level) unchanged from FORMULATION §2(b2).

## 7. Test inventory added in S3 (tests/test_s3_witnesses.py; 62 tests)

- `test_master_b_matches_derived_rows[8]`: MASTER_TARGET_B == derived
  D c_s^2, exact, all braiding rows (rev-v1.1 pin inherited).
- `test_master_b_reduces_to_s2_master`: b = 0 reduction == S2 master, exact.
- `test_inheritance_lemma_b_rows[8]`: inheritance algebra, braiding rows.
- `test_recertify_coarse[15]`: full certification path re-run at 64 bits.
- `test_certificate_json_valid[15]`: stored certificates match the frozen
  registry, all checks pass, rigorous endpoints on the right side of bounds.
- `test_witness_float_margins[15]`: float scans vs certified bounds + D1/D2/
  D3 regression; rigorous lower bounds never exceed float minima.
Suite: 111 passed (24 S1 + 25 S2 + 62 S3), ~17 s.
