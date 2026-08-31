# P9(b) — certified results summary (2026-08-31)

Frozen conventions: FORMULATION.md v1.2 (GLV14 eq. (86) five-α action, BS14
braiding convention, family B = flat FLRW + dust, fiducial F* = CPL(−4/5, −4/5)
with Ω_m0 = 3/10, crossing at a_c = 3/4; stability D2 = no ghost + exact-eigen-
speed c_s² ≥ 0 (+ c_T² ≥ 0, M*² > 0); bounds D3: 0 < α_K ≤ 10³, |α_i| ≤ 10).
Adversarial review 2026-08-31: chain survives, no flaw, no claim-level gap.

## The classification (P9b-1, CERTIFIED)

| row S | verdict | witness / identity | certified margin |
|---|---|---|---|
| {K} | (N) no-go | I1: α_K c_s² = 3Ω_DE(1+w_DE) (all α_K > 0) | exact |
| {K,B} | (E) | α_B^BS = a/2, J = [3/5,1] | c_s² ≥ 1/10, D ≥ 11/10 |
| {K,T} | (E) | α_T = −1/10 | c_s² ≥ 1/20 |
| {K,M} | (E) | α_M = 1/2, J = [3/5,1] | c_s² ≥ 1/4 |
| {K,H} | (E) | α_H = 1 | c_s² ≥ 1/4 |
| 11 larger rows | (E) | genuinely joint witnesses | c_s² ≥ 1/5 … 7/5 |

**Minimal operator sets admitting a stable phantom crossing: exactly
{α_K,α_B}, {α_K,α_T}, {α_K,α_M}, {α_K,α_H}.** The folklore "braiding is the
unique single-field rescue" is refuted within the frozen conventions: tensor-
speed, Planck-mass running, and the beyond-Horndeski operator each rescue the
crossing on their own. (GW170817's |c_T² − 1| bound is observational and
outside every claim — FORMULATION §4.)

## Positivity (P9b-2, CERTIFIED conditional on POS)

POS = Melville–Noller 1904.05874 (16)–(17) + α-form (10) (quartic example-
model prior, frozen as an assumption), applicable to the 7 Horndeski rows
only (α_H rows outside its domain). Under POS: row {K,T} is excluded exactly
(POS forces α_T ≥ 0; the row identity α_K c_s² = f̂ − 2α_T then forbids any
crossing); {K,M} satisfies by saturation; POS-compatible witnesses re-
certified for {K,B}, {K,B,T}, {K,B,M}, {K,T,M}, {K,B,T,M}. Minimal sets
within POS's domain: {α_K,α_B} and {α_K,α_M}.

## Evidence

Exact sympy derivation chain from the frozen action (`derivation/`), 22 Arb
certificates (`certificates/P9b-{1,2}.*.json`, regenerable via
`python -m p9b.certify` / `python -m p9b.positivity`), 162 tests
(`uv run pytest problems/P9/b/tests -q`). Ledger: CLAIMS.md rows P9b-0,
P9b-1.K, P9b-1.rows, P9b-1, P9b-2, P9b-2s. Notes: s1–s4 under `notes/`.
Eigen-speed citation of record: GLPV 1408.1952 dispersion relation
(1609.01272 eq. (2.18) discrepancy recorded, FORMULATION rev. v1.2).
