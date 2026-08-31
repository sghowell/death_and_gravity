# P9(b) S4 — P9b-2: the positivity intersection, conditional on POS

Status: S4 complete, 2026-08-31.  Conventions: FORMULATION.md v1.1.  Builds
on S3 (`notes/s3-witnesses.md`; the 15 (E) witnesses).  New code:
`src/p9b/positivity.py` (POS pin, applicability map, Arb verdicts,
alternative witnesses), `certificates/P9b-2.<ROW>.json` (7),
`tests/test_s4_positivity.py` (51 tests; suite: 162, ~18 s).
Run: `uv run pytest problems/P9/b/tests -q`; re-generate certificates with
`PYTHONPATH=problems/P9/b/src:problems/P9/b uv run python -m p9b.positivity`.

## 1. The pinned POS inequalities (re-verified against the paper)

Source re-fetched 2026-08-31 (ar5iv render of arXiv:1904.05874, Melville–
Noller, PRD 101, 021502(R); equations transcribed from the LaTeX alttext).
MN conventions: signature (−,+++); Λ₁ = M_Pl, Λ₂² = M_Pl H₀, Λ₃³ = M_Pl H₀²;
X ≡ −∇^μφ∇_μφ/(2Λ₂⁴) dimensionless; overbar = flat-background value.

- MN (4)–(5): 2→2 tree amplitude A = c_ss s²/Λ₂⁴ + c_sst s²t/Λ₃⁶ + …;
  bounds c_ss ≥ 0, c_sst ≥ −(3Λ₃⁴/2Λ₂⁴)c_ss.
- **MN (16)** (general Horndeski, φφ→φφ, their (14)–(15); Λ₂ ≫ Λ₃ assumed):
  c_ss ≥ 0 AND c_sst ≥ 0.
- **MN (17)** (shift-symmetric Horndeski):
  2Ḡ₂,X Ḡ₄,X ≥ −Ḡ₂,XX Ḡ₄  and  2Ḡ₄,XX + 2Ḡ₄,X²/Ḡ₄ ≤ Ḡ₃,X².
- **MN (10)** (example model, their (3): shift-symmetric QUARTIC Horndeski
  L = Λ₂⁴G₂(X) + M_Pl²G₄(X)R + Λ₂⁴G₄,X([Φ]²−[Φ²]), plus small mass):
      pos. prior:   α_B ≤ 2α_T/(1+α_T).
  Derivation chain re-checked by hand from the fetched equations: the c_sst
  part of (16) specialises to Ḡ₄,X²/Ḡ₄ ≤ −Ḡ₄,XX (their (7)); with their (8)
  M²α_T = 4XG₄,X, M²α_B = 8X(G₄,X + 2XG₄,XX), M² = 2(G₄ − 2XG₄,X)
  = 2G₄/(1+α_T), multiplying by 16X²/(M²) > 0 gives exactly (10).  Side
  conditions used (checked per row): X > 0, M² > 0, G₄ > 0 ⟺ 1 + α_T > 0.
- MN (11), the (sub)luminality prior α_T ≤ 0, is a SEPARATE prior and is
  NOT part of POS (FORMULATION froze (16)–(17)+(10) only).  Kept in mind in
  §4 as a remark; it enters no claim.
- MN's transfer caveat (their appendix, "Positivity caveats"): bounds are
  derived for massive particles on flat space; on FLRW they *assume* the
  bounds "continue to hold for the G_i evaluated on the cosmological ⟨φ⟩".
  That assumption, applied pointwise at each a ∈ J, IS our frozen POS.

Sign conventions (the S1-checked map applies): MN's α_B is α_B^BS — their
(8) quartic term M²α_B = 8X(G₄,X + 2XG₄,XX) coincides with BS14 (A.7), and
MN cite BS14 for the α definitions and use the BS14 parametrization
α_i = c_iΩ_DE (their (12)).  The apparent G₃-sign difference in MN's
appendix α_B (−2(φ̇/H)XG₃,X vs BS14's +2φ̇XG₃,X/H) is exactly the L₃ =
±G₃□φ convention (G₃^MN = −G₃^BS, cf. the S1 KGB pin G₃ = −G_ours) — no
discrepancy; irrelevant to (10) anyway (G₃ ≡ 0 in the example model).
So (10) applies verbatim to our reported α_B^BS with NO sign flip.

## 2. Applicability map (part of the theorem statement)

Which rows can POS be evaluated on at all?

1. Only (10) is an α-form.  (16)–(17) constrain flat-background values
   Ḡ₂,X, Ḡ₂,XX, Ḡ₃,X, Ḡ₄, Ḡ₄,X, Ḡ₄,XX that an α-trajectory does not
   determine; MN state the c_ss part has no α-form (only α_K depends on
   G₂, and they drop it).  So POS on an α-witness = pointwise (10).
2. **Horndeski gate.** MN's amplitudes (14)–(15) are computed in Horndeski
   (2).  α_H ≠ 0 requires beyond-Horndeski (GLPV) operators outside that
   computation: POS forces NO inequality on the 8 α_H rows.  Verdict class
   (iii): KH, KBH, KTH, KMH, KBTH, KBMH, KTMH, KBTMH — POS not applicable
   (this is a statement about POS's domain, not about those rows' physics).
3. **Trajectory-level conditionality (recorded, not smuggled).**  Within the
   example model the α's obey ties (their (9): α_M = −(Ẋ/4HX)α_B; α_B
   const if Ẋ = 0).  NO frozen witness — and generically no c_iΩ_DE
   trajectory of MN's own MCMC — is an exact example-model trajectory.
   POS as frozen is (10) applied pointwise to the trajectory, exactly MN's
   own usage (they scan c_B, c_M, c_T independently under prior (10)).
   Rows with α_M on (KM, KBM, KTM, KBTM) rely on that usage one step more
   (α_M free rather than tied); named in the claim's conditionality.
4. Implicit X > 0 (rolling scalar) — automatic in the unitary-gauge EFT.

Applicable: the 7 Horndeski rows KB, KT, KM, KBT, KBM, KTM, KBTM.
({K} has no witness — P9b-0 — and is outside P9b-2's quantifier.)

## 3. P9b-2 (weak) verdicts — Arb, 128-bit, outward rounding

q(a) ≡ 2α_T/(1+α_T) − α_B^BS(a); POS ⟺ q ≥ 0 on J.  All sweeps strict
ball verdicts (house style of certify.py); J = [3/5, 1] for all 7 rows.
"margin" = certified rigorous endpoint of q over J (max upper endpoint for
violations: q below it everywhere; min lower endpoint for satisfaction).

| row  | RHS 2α_T/(1+α_T) | verdict                     | certified margin  |
|------|------------------|-----------------------------|-------------------|
| KB   | 0                | violates at EVERY a ∈ J     | q ≤ −0.299999…    |
| KT   | −2/9             | violates at EVERY a ∈ J     | q ≤ −0.222222…    |
| KM   | 0                | satisfies (exact: 0 ≤ 0)    | q ≡ 0 (saturated) |
| KBT  | −2/9             | violates at EVERY a ∈ J     | q ≤ −0.522222…    |
| KBM  | 0                | violates at EVERY a ∈ J     | q ≤ −0.299999…    |
| KTM  | −2/9             | violates at EVERY a ∈ J     | q ≤ −0.222222…    |
| KBTM | −2/9             | violates at EVERY a ∈ J     | q ≤ −0.522222…    |

The certified violating subinterval is all of J in each violating case
(q = −b for α_T-off rows, q = −2/9 − b for α_T = −1/10; b = a/2 ≥ 3/10 or
b = 0).  Side condition 1 + α_T ∈ {1, 9/10} > 0 holds exactly on all rows.
KM saturates (α_B = α_T ≡ 0): the weak inequality holds with exact
equality at every a; a strict version is unattainable IN-ROW (both α's are
identically zero by the row definition).  Certificates: P9b-2.<ROW>.json.

## 4. Alternative witnesses (one bounded attempt per violating family)

A violating WITNESS does not decide the row.  One bounded search per family
(float proposer + full Arb re-certification of stability AND strict POS,
same background F*, same J, all D3 bounds):

- **B-family** (KB, KBT, KBM, KBTM; violation cause α_B > 0 with α_T ≤ 0):
  make the braiding negative.  α_T-off rows: b(a) = a/2 − 21/40 (< 0 on J,
  POS margin ≥ 1/40); α_T = −1/10 rows: b ≡ −1/4 < −2/9 (margin 1/36).
- **Braidingless-T family** (KT, KTM; violation cause: with α_B ≡ 0, (10) ⟺
  α_T ≥ 0 while the witnesses need α_T < 0): flip to α_T = +1/10.

| row  | alternative (α_K = 1)          | lo(c_s²)  | n      | POS margin ≥ |
|------|--------------------------------|-----------|--------|--------------|
| KB   | b = a/2 − 21/40                | +0.100116 | 69/9   | 1/40         |
| KBT  | b ≡ −1/4, α_T = −1/10          | +0.050242 | 126/9  | 1/36         |
| KBM  | b = a/2 − 21/40, α_M = 1/2     | +0.601410 | 34/9   | 1/40         |
| KBTM | b ≡ −1/4, α_T = −1/10, α_M=1/2 | +0.550905 | 37/9   | 1/36         |
| KTM  | α_T = +1/10, α_M = 1/2         | +0.100187 | 53/12  | 2/11         |

(lo = rigorous lower ball endpoint of c_s² over J; n = subintervals/depth;
D ≥ 1 certified — braiding rows D > 1 strict, lo(D) 1.00047 / 1.09375;
fhat sign sweeps, Ω_DE > 0, E² > 0, M² > 0, c_T² ≥ 0, |α_B| ≠ 2 all
re-certified per row; every alternative is joint — row α's nonzero on J.)

- **KT: no alternative exists — the whole row is excluded (exact).**
  POS with α_B ≡ 0 forces α_T ≥ 0 (side condition 1+α_T > 0).  Exact S2
  row identity (re-verified in test_kt_row_identity_exact):
  α_K c_s² = f̂ − 2α_T with f̂ = f/(M_Pl²H²).  A D1 crossing has f̂ < 0 on
  the phantom side, so α_K c_s² < 0 there; α_K > 0 (D3) gives c_s² < 0,
  violating D2(ii).  Holds for every background in family B and every
  bounded α_T: {stable crossing} ∩ POS = ∅ on row {K,T} — region-level,
  same evidence level as P9b-0 (exact identity + stated sign analysis).
- Remark (no claim): on braidingless-T rows POS forces α_T ≥ 0, hence
  c_T² ≥ 1 — witnesses there are superluminal-tensor whenever α_T ≠ 0
  (KTM-alt: c_T² = 11/10).  Adding MN's separate luminality prior (11)
  would re-close the braidingless-T family entirely (α_T ≡ 0 collapses
  KTM's α_T off) — exactly MN's Fig. 1 pos∧lum squeeze.  Recorded for S5.

Net, conditional on POS: stable crossing survives on 6 of the 7 applicable
rows ({K,M} frozen witness; 5 replacement witnesses) and is excluded on
{K,T}.  Restricted to its applicability domain, POS shifts the minimal-set
classification to: {K,B} and {K,M} minimal ({K,T} excluded, {K,H} outside
POS's domain).  The full P9b-1 classification (unconditional) is unchanged.

## 5. P9b-2s (strong, region-level) — feasibility verdict

Representative rows examined: {K,T} and {K,B}, degree ≤ 2 polynomial-in-a
α's with rational coefficient boxes (FORMULATION §2(b2) freeze).
(a) Nonemptiness of {stable crossing} ∩ POS is ALREADY decided by exact
certificates on all 6 surviving applicable rows: the certified points are
the §4 witnesses (degree ≤ 1 in a) — this is precisely the "explicit
certified point" branch of the FORMULATION's P9b-2s.
(b) Emptiness on {K,T} is ALREADY decided exactly (§4): the row identity is
linear in α_T, and the sign argument needs no SOS — a Positivstellensatz
certificate would be trivial (linear combination).
(c) The only content a full SOS/Positivstellensatz machine would add is
emptiness proofs over entire coefficient boxes on rows where we instead
exhibit points — not needed for any shipped statement.  It would require
certified polynomial envelopes for the transcendental background (f̂, Ḣ/H²,
ρ_m/M²H² involve exp) before semialgebraic methods apply: real work, no
payoff.  **Verdict: P9b-2s ships only in the degenerate forms (a)+(b)
above, which are already certified; the general SOS machinery is DROPPED,
and the ledger row will say so** (as FORMULATION §2(b2) provides).

## 6. Draft CLAIMS rows (main session owns CLAIMS.md; do not edit from S4)

| P9b-2 | P9(b) | Positivity intersection, conditional on POS (frozen: MN 1904.05874 (16)–(17); α-form (10) α_B^BS ≤ 2α_T/(1+α_T), side condition 1+α_T > 0, applied pointwise on J under MN's flat-space-transfer assumption and their trajectory-level (12)-usage; BS14 α_B convention, S1 map). Applicability is part of the statement: POS constrains exactly the 7 Horndeski rows; the 8 α_H rows are outside its domain (GLPV operators absent from MN's amplitudes — no verdict). Arb-certified (128-bit, outward): {K,M}'s witness satisfies POS (exact saturation 0 ≤ 0); the other 6 frozen witnesses violate POS at every a ∈ J (certified margins ≥ 2/9); POS-compatible replacement witnesses are certified (full D1–D3 stability + strict POS) for {K,B}, {K,B,T}, {K,B,M}, {K,T,M}, {K,B,T,M}; row {K,T} is excluded entirely — POS forces α_T ≥ 0 and the exact identity α_K c_s² = f̂ − 2α_T then contradicts D2 at any crossing (whole row, family B). Conditional on POS, stable crossing survives on 6 of 7 applicable rows; minimal sets within POS's domain: {α_K,α_B}, {α_K,α_M} | CERTIFIED conditional on POS (certificates/P9b-2.*.json; tests/test_s4_positivity.py) | problems/P9/b/src/p9b/positivity.py; problems/P9/b/notes/s4-positivity.md | 2026-08-31 |

| P9b-2s | P9(b) | Region-level positivity intersection: shipped only in the certified degenerate forms — exact emptiness of {stable crossing} ∩ POS on row {K,T} (identity-level) and certified-point nonemptiness on the 6 surviving applicable rows (the P9b-2 witnesses, degree ≤ 1 in a); the general degree-≤2 SOS/Positivstellensatz emptiness machinery is dropped per FORMULATION §2(b2) (s4-positivity.md §5) | CERTIFIED conditional on POS (in the stated degenerate forms; no general-box claim) | problems/P9/b/notes/s4-positivity.md §5 | 2026-08-31 |

## 7. S5 handoff (writeup + adversarial review must cover)

- The POS pin: equations re-verified against the paper (fetch of
  2026-08-31); the G₃ sign reconciliation (§1) and the no-flip conclusion
  for α_B^BS; reviewer should re-derive (7)⇒(10) independently.
- The applicability map is load-bearing: defend (i) α_H rows NA (Horndeski
  gate), (ii) trajectory-level usage for α_M rows (MN's own (12) practice),
  (iii) "only (10) has an α-form".  An attack via OTHER positivity results
  (boost-breaking, dS — digest D14/D15) is out of scope: POS is frozen.
- Scope hygiene: "violates POS" is a statement about WITNESSES except for
  {K,T}, where it is a row-level theorem; KM's satisfaction is saturation
  of a weak inequality (α_B = α_T ≡ 0), not an interior point.
- The luminality remark (§4) must stay a remark (no claim); if S5 wants it
  as a corollary it needs its own assumption row (POS+lum).
- P9b-2s: verify the drop decision reads as FORMULATION §2(b2) intended.
- Re-run path: python -m p9b.positivity regenerates all 7 certificates.

## 8. Test inventory (tests/test_s4_positivity.py; 51 tests)

- `test_kt_row_identity_exact`, `test_kt_certificate_records_impossibility`:
  exact sympy backing + certificate content for the {K,T} exclusion.
- `test_applicability_partition`, `test_na_rows_have_no_certificate[8]`:
  the frozen applicability map; NA rows refuse certification.
- `test_pos_rhs_exact`: exact rational RHS values; side-condition guard.
- `test_certificate_json_valid[7]`, `test_certificate_alternative_valid[5]`,
  `test_km_certificate_saturation`: stored-certificate validation.
- `test_recertify_coarse[7]`: full S4 path re-run at 64 bits.
- `test_frozen_verdict_float[7]`, `test_alternative_float_margins[5]`:
  independent numpy re-checks; rigorous endpoints on the right side.
- `test_ext_eval_matches_certify_row_eval[7]`: S4 evaluator vs S3
  certify.row_eval ball overlap on frozen witnesses.
Suite: 162 passed (24 S1 + 25 S2 + 62 S3 + 51 S4), ~18 s.
