# P2 S3 pre-registration: part (iii), validated zero counts on R_l

Status: written 2026-10-17, before any production run, on branch `p2/s3-part-iii` (from main at
`7004375e`, which contains S2). Code: `src/p2/zsol.py`, `src/p2/part3.py`, `src/p2/zeq_identity.py`,
`tests/test_part3.py`. Scope from `s0-gate.md` §7 (S3) and FORMULATION P2-3.

## 1. What S3 does, and for which l and points

- **Target (FORMULATION P2-3).** For each l in the range below:
  - W_l is analytic on a neighbourhood of R_l;
  - W_l ≠ 0 on ∂R_l;
  - the winding number of W_l around ∂R_l is 1 for l = 0 (the simple zero at 0) and 0 otherwise.
- **Point.** P_B only. P_A is not in S3 production.
  - P_A's part (i) must reach 0.667γ₀, and S4 has not yet shown which l₀ it gets there. Counting
    zeros at P_A before that would spend budget on an l-range that may change.
  - The code takes the point as a parameter. P_A enters S3 only as a positive control (C2(d)).
- **l-range: l = 0, 1, …, 19.**
  - Part (i) is not certified at any l yet. The only part-(i)-route certificates are S2's
    barrier-top windows (P2-1w), at l = 20, 50 and 100.
  - S3 therefore covers every l below the lowest l with a part-(i)-route certificate.
  - The S2 projection l₀ = 2 is not used; it is the output of a mis-specified fit (s2-report §3).
  - If S4 ends with l₀ > 20, the same pipeline extends to l = 20, …, l₀ − 1. The cost per l grows
    about like (l+½)² (§5).
- **The region.** R'_l = [−a_l, a_l] × [−b, e], with exact dyadic corners:
  - a_l is the certified P2-2 threshold (`certificates/lg/lg_PB_l{l}.json`, field `W_l_exact`);
  - b is β = (½+δ)κ₋ rounded up to 2⁻³⁰ (b/κ₋ = 0.55000001), so R'_l contains the strip part
    −β ≤ Im ω;
  - e is η = κ₋/10 rounded down to 2⁻³⁰ (0 < e ≤ η). Im ω > 0 is covered by P2-0, and 0 ∈ R'_l.
- **The l = 0 enclosure.** A Krawczyk enclosure of the l = 0 near-extremal zero
  ω_* ≈ −0.909381iκ₋, computed with the S3 Wronskian (independent of the S1 code behind KAT-6).
  Output: `certificates/part3/krawczyk_PB_l0.json`.
  - It gives one side of β: β ≤ 0.909381.
  - The two-sided β enclosure (P2-6) is not attempted in S3. Its strip edge −0.92κ₋ passes E₊'s
    first pole (FORMULATION v0.2 note), which needs new foundations.

## 2. Method

**W_l: phase-factored Jost solutions** (`zsol.py`; proof sketch in its docstring).
- Put Z₊ = e^{iω(r* − c₊)}R₊ and Z_c = e^{−iω(r* − c_c)}R_c. These are the solutions, analytic at
  r₊ resp. r_c with value 1 there, of
  ΔZ″ + (Δ′ − 2isωr²)Z′ − (2isωr + L)Z = 0, with s = +1 resp. −1.
- The S1 Jost-normalized Wronskian is then exactly
  W_l(ω) = [Δ(Z₊Z_c′ − Z₊′Z_c) + 2iωr²Z₊Z_c]/(r₊r_c), at any r ∈ (r₊, r_c).
- sympy checks these identities (`zeq_identity.py`: Z1–Z4, including the indicial factor and the
  conjugation symmetry).
- Series: power series at r₊ and r_c (exponent 0; majorant tail bounds), then Taylor steps to the
  matching point r_P (majorant tails). The algorithms are those of S1's `arbw.py`, applied to the
  factored equation.
- **Why the change.** Pilot runs (§5) show that the S1 box evaluation of W over an ω-disc
  overestimates by about 10⁴⁷ at the corner of R_5. The free phase e^{±iωr*} wraps in ball
  arithmetic. The factored form never evaluates that phase.
- **Enclosure over a disc.** The centre trajectory is propagated with point step matrices. The
  deviation over the disc is kept in a Lohner frame with disc radii: the parallelepiped of the
  centre propagation, and an orthonormal frame when that one degenerates. The ω-variation of each
  step enters additively (step deviation applied to the centre vector), plus a second-order term.

**Counting: Cauchy-interpolation discs** (`part3.py`; the scheme of S1 `winding.py` that the S1
report recommended, §7.3 there).
- By the symmetry W_l(−ω̄) = conj W_l(ω), only the right half H of ∂R'_l is computed: the bottom
  edge, the right edge, the top edge.
- The winding number is Δ_H arg W_l/π.
- H is covered by segments with exact dyadic endpoints, each within ρ/2 of its centre c.
  - ρ ≤ κ₋, and the box [c ± ρ]² must avoid E₊ ∪ E_c: certified, as √2ρ < dist(c, E).
  - The planner uses ρ ≤ dist/(1.5√2).
- On each disc:
  - M ≥ sup|W_l| from the zsol box evaluation (step phase θ = 1/4);
  - K point values on |ω − c| = ρ/2, at 128 bits with θ = 3;
  - K is adaptive: K = max(24, 8⌈log₂(4M/|W_c|/10⁻⁶)/8⌉), capped at 160;
  - the coefficients carry the aliasing bound Mx^K/((1 − x^K)ρ^k), with x = ½;
  - the truncation bound is M(|d|/ρ)^K/(1 − |d|/ρ).
- **Pieces.** Each segment starts as pieces of length ≤ κ₋/16. A piece is bisected, up to depth
  14, until its enclosure excludes 0. Its increment is the principal Arg of the ratio of its end
  enclosures, certified below π in absolute value.
- **Failure handling.** A segment that fails is halved, with its radius, down to ρ = κ₋/256.
- **Certificate.** `certificates/part3/winding_PB_l{l}.json` holds every segment, in order:
  - exact endpoints, centre and radius;
  - the increment as a ball;
  - K, M and timing;
  - the commit.

  The winding number is certified when Σ/π lies within ¼ of an integer.
- **Analyticity.** It is checked disc by disc, against the explicit pole set E₊ ∪ E_c (S1, D1′).
  R'_l lies in Im ω ≥ −b > −κ₊, so W_l is analytic on a neighbourhood of R'_l.

**Deviations from FORMULATION §4**, declared here:
- Cauchy discs instead of ω-Taylor-model tiles.
- The corners are dyadic roundings, b ≥ β and e ≤ η.
- The pole checks are recorded per segment inside the winding certificates. There is no separate
  `analyticity_*.json`.
- The files live under `certificates/part3/`.

## 3. Acceptance criteria (KAT-style) and controls

The gate before production is C0–C3 and C6, run by `python -m p2.part3 controls`, which writes
`certificates/part3/controls_PB.json` and `krawczyk_PB_l0.json`.

- **C0 (two independent W codes).** zsol agrees with S1's `arbw.wronskian` (the balls overlap) at
  24 points:
  - l ∈ {0, 1, 2, 5, 10, 19};
  - for each l: (0.3a_l, −0.55), (0.7a_l, 0.1), (a_l, −0.55) and (1, 0.1), in κ₋ units.
- **C1 (D4).** W₀(0) contains 0, and the Cauchy coefficient W₁ at 0 contains
  i(r₊² + r_c²)/(r₊r_c), S1's exact value.
- **C2 (known zeros are counted).** Full-boundary windings, no symmetry used:
  - (a) P_B, l = 1, [1, 12] × [−2, 0.1]κ₋ → 1 (PS 7.2149 − 1.7533i);
  - (b) P_B, l = 2, [1, 17] × [−2, 0.1]κ₋ → 1 (PS 12.1965 − 1.7079i);
  - (c) P_B, l = 0, [−0.05, 0.05] × [−0.95, −0.87]κ₋ → 1 (NE −0.9094i, between the poles);
  - (d) P_A, l = 1, [3.19, 3.29] × [−0.85, −0.75]κ₋ → 1 (CCDHJ18 Table I).
- **C3 (negative control: discrimination at a zero).** P_B, l = 1:
  - [6.9, 7.5] × [−1.80, −1.743]κ₋ → 1;
  - [6.9, 7.5] × [−1.743, −1.70]κ₋ → 0.

  The shared edge passes 0.0103κ₋ above the PS mode. A wrong certified count kills; a failure to
  certify is reported but does not kill.
- **C4 (symmetry).** At l = 2, the left half of ∂R'_2 is computed explicitly
  (`python -m p2.part3 symmetry`). Its sum must overlap the certificate's right-half sum.
- **C5 (float consistency).** The certified counts at l = 0, 1, 2 and 5 equal the S0 float values
  1, 0, 0 and 0.
- **C6 (the l = 0 enclosure).** The zsol Krawczyk test certifies a unique zero, and its enclosure
  overlaps KAT-6's (−0.90938097899480545552493376262i ± 7e−30, in κ₋ units).
- **C7 (replay).** `python -m p2.part3 verify --l 0-19`:
  - re-sums every certificate exactly from its stored increments;
  - checks that the segments tile H in order, with exact shared endpoints;
  - recomputes a deterministic 5% sample of segments per l (at least 3), whose increments must
    overlap the stored ones.

  `tests/test_part3.py` does the same at a smaller scale.

## 4. Expected outcome and kill criteria

**Expected.** Winding 1 at l = 0, from the simple zero at 0 (D4), and 0 for l = 1, …, 19 (S0 floats
for l = 0, 1, 2, 5).

**K1 (kill).** A certified winding number different from the expected value, at any l.
- Stop production.
- Localize the zero by bisecting R'_l, then enclose it with Krawczyk.
- Compare with the float spectrum, and report.
- A certified extra zero in R'_l would refute β > ½ + δ at P_B.

**K2 (kill).** Any of C0–C4 or C6 fails, or C3 certifies a wrong count. The pipeline is then not
trusted: stop and report.

**K3 (compute).** Stop and report in either case:
- the projection from the measured cost of the completed l exceeds 20 core-h;
- the S3 total reaches the 24 core-h cap.

**Not a kill.** A segment that cannot be certified after halving to κ₋/256 makes its l
INCONCLUSIVE. That l is reported, and the other l continue.

## 5. Compute budget and execution

**Pilot.** Development runs, before this document; no production use.
- Sup-bound surveys at l = 0, 1, 5, 10 and 19, which led to the factored form. The S1 box
  evaluation gave M/|W| ≈ 10⁴⁷ at the corner of R_5; zsol gives 10⁴–10¹⁸ for ρ ≤ κ₋ up to l = 19.
- Full windings at l = 0, 1, 2 (1, 0, 0; 93 segments, 230 s CPU).
- 21 timed segments at l = 5, 10 and 19.
- Smoke tests of C1, C2(c) and C6.
- Total: about 0.5 core-h.

**Cost model.** About 1.3 + 0.12·(Re ω/κ₋) CPU-s per segment.
- The planner gives 2351 segments for l = 0–19: 21 at l = 0, 123 at l = 10, 213 at l = 19.
- The projection is 3.4 core-h, or about 5 core-h with an allowance of 1.5× for splits.
- The controls add about 0.3 core-h, and the replay sample about 0.2 core-h.
- Planned total: about 6 core-h. The cap is 24 core-h.

**Execution.**
- `nohup nice -n 10 python -m p2.part3 run --point PB --l 0-19 --workers 4`, at most 4 workers.
- The run is resumable: one JSONL line per finished segment, keyed by its exact endpoints, under
  `results/part3/`.
- A launch record (`results/part3/launch_*.json`) holds the command, commit, host, pid, start, end
  and failures.
- `python -m p2.part3 assemble --l …` writes the certificates.
- Outputs are committed as the l complete.

## 6. Level

If everything passes, each count is CERTIFIED: Arb certificates, re-checked by tests. The written
steps are routine:
- the argument principle on a rectangle;
- the conjugation symmetry (sympy Z4);
- analyticity off E₊ ∪ E_c (S1, D1′);
- the factored identities (sympy Z1–Z3).

The claim is P2-3 at P_B for l ≤ 19 only. P2-3 itself needs every l < l₀, and l₀ comes from S4.
