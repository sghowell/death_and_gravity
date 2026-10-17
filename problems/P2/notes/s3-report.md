# P2 S3 report: part (iii) at P_B, l = 0–19

Status: S3 exit, 2026-10-17. Branch `p2/s3-part-iii`, from main at `7004375e`.
Pre-registration: `notes/s3-preregistration.md` (commit `89a2336b`), written before any
production run.

**Verdict: PASS.** Every pre-registered criterion holds (C0–C7), and no kill criterion fired.
- The count is 1 at l = 0 (the simple zero at ω = 0) and 0 for l = 1, …, 19.
- With P2-0 and P2-2, the strip S_δ at P_B therefore contains no QNM of degree l ≤ 19 other than
  ω = 0.
- The l = 0 near-extremal zero has an independent Krawczyk enclosure.
- Compute: 4.7 core-h of the 24 core-h cap.

## 1. What was certified

For each l = 0, …, 19 at P_B, the zeros of W_l in R'_l = [−a_l, a_l] × [−b, e] were counted:
- a_l is the P2-2 threshold;
- b ≥ β = (½+δ)κ₋, with b/κ₋ = 0.55000001;
- 0 < e ≤ η = κ₋/10, with e/κ₋ = 0.09999998.

All corners are exact dyadics. The count is the winding number of W_l around ∂R'_l, computed on
the right half H by the conjugation symmetry.

| l | a_l/κ₋ | segments | K max | max M/abs(W_c) | min abs(W) on H ≥ | winding ball | count | CPU s |
|---|---|---|---|---|---|---|---|---|
| 0 | 5.40 | 21 | 64 | 9.0e11 | 0.0218 | 1.0000 ± 6.5e−5 | **1** | 33 |
| 1 | 10.31 | 31 | 80 | 2.4e16 | 0.525 | 0 ± 5.1e−5 | **0** | 64 |
| 2 | 15.61 | 41 | 80 | 2.5e17 | 0.748 | 0 ± 1.1e−4 | **0** | 96 |
| 3 | 20.86 | 51 | 88 | 1.0e18 | 1.00 | 0 ± 1.1e−4 | **0** | 132 |
| 4 | 26.02 | 63 | 88 | 2.2e18 | 1.19 | 0 ± 1.7e−4 | **0** | 179 |
| 5 | 31.19 | 73 | 88 | 6.4e18 | 1.43 | 0 ± 2.3e−4 | **0** | 230 |
| 6 | 36.28 | 83 | 88 | 4.2e19 | 1.60 | 0 ± 2.3e−4 | **0** | 288 |
| 7 | 41.37 | 93 | 88 | 1.9e19 | 1.83 | 0 ± 2.6e−4 | **0** | 353 |
| 8 | 46.43 | 102 | 96 | 1.5e20 | 2.00 | 0 ± 3.9e−4 | **0** | 423 |
| 9 | 51.52 | 113 | 88 | 7.2e19 | 2.22 | 0 ± 3.6e−4 | **0** | 503 |
| 10 | 56.57 | 123 | 96 | 1.8e20 | 2.37 | 0 ± 3.3e−4 | **0** | 592 |
| 11 | 61.60 | 133 | 96 | 1.7e20 | 2.58 | 0 ± 4.1e−4 | **0** | 691 |
| 12 | 66.60 | 143 | 96 | 5.4e20 | 2.72 | 0 ± 4.1e−4 | **0** | 799 |
| 13 | 71.66 | 153 | 96 | 8.9e20 | 2.94 | 0 ± 5.9e−4 | **0** | 907 |
| 14 | 76.70 | 163 | 96 | 2.4e21 | 3.00 | 0 ± 4.9e−4 | **0** | 1037 |
| 15 | 81.62 | 173 | 96 | 7.5e21 | 3.27 | 0 ± 5.3e−4 | **0** | 1169 |
| 16 | 86.72 | 183 | 96 | 4.0e21 | 3.33 | 0 ± 5.9e−4 | **0** | 1309 |
| 17 | 91.70 | 193 | 96 | 1.1e22 | 3.60 | 0 ± 5.9e−4 | **0** | 1467 |
| 18 | 96.66 | 203 | 96 | 1.4e22 | 3.70 | 0 ± 6.5e−4 | **0** | 1635 |
| 19 | 101.71 | 213 | 96 | 1.1e22 | 3.95 | 0 ± 6.2e−4 | **0** | 1838 |

- All 2351 planned segments certified at the first attempt: no segment needed halving.
- All 20 certificates were computed at `c7d1bfbc` (the launch record), on 4 workers at nice 15,
  in 57 min wall.
- The per-segment JSONL working files in `results/part3/` are not committed. The certificates
  carry the same records.
- "min abs(W) on H" is the smallest lower bound of |W_l| over the pieces of H. At l = 0 it is
  attained next to ω = 0 on the top edge, where |W₀| ≈ |W₀′(0)|·e.

- **The l = 0 count is 1.** That zero is the simple zero at ω = 0 (D4, P2-0), so R'_0 contains no
  other zero.
- **Counts 0 for l ≥ 1.** These agree with the S0 float counts at l = 1, 2 and 5 (C5).
- **With P2-0 and P2-2.** For every l ≤ 19, the strip S_δ = {Im ω ≥ −(½+δ)κ₋} contains no QNM of
  degree l except ω = 0 at l = 0:
  - P2-2 covers |Re ω| ≥ a_l on the lower strip;
  - P2-0 covers Im ω ≥ 0;
  - the counts cover the rest.

## 2. Method (as pre-registered)

**W_l from phase-factored Jost solutions** (`src/p2/zsol.py`).
- Z₊ = e^{iω(r*−c₊)}R₊ and Z_c = e^{−iω(r*−c_c)}R_c solve
  ΔZ″ + (Δ′ ∓ 2iωr²)Z′ − (±2iωr + L)Z = 0, with value 1 at their horizons.
- The S1 Jost-normalized Wronskian is exactly
  [Δ(Z₊Z_c′ − Z₊′Z_c) + 2iωr²Z₊Z_c]/(r₊r_c).
- sympy checks the factored equation, this identity, the indicial factor and the conjugation
  symmetry (`certificates/part3/zeq_identity.json`).
- **Why the change.** The S1 code evaluates the free phase e^{±iωr*} in ball arithmetic, so its
  box evaluation over an ω-disc wraps: M/|W| ≈ 10⁴⁷ at the corner of R_5.
- zsol never evaluates that phase. Its disc bounds use Lohner frames with disc radii, and give
  M/|W| ≤ 10¹⁸ for ρ ≤ κ₋ up to l = 19.
- Point values are accurate to about 10⁻³⁰ relative at 128 bits.

**Cauchy-interpolation discs** (`src/p2/part3.py`; the scheme of S1's `winding.py`).
- H is covered by segments of length ≤ κ₋, each within ρ/2 of its disc centre.
- Each disc carries:
  - a box bound M;
  - K adaptive point values on the circle of radius ρ/2;
  - the aliasing and truncation bounds;
  - pieces bisected until each enclosure excludes 0.
- The pole set E₊ ∪ E_c is checked against every disc's box.

## 3. Acceptance criteria and controls

| criterion | outcome |
|---|---|
| C0: zsol against S1's arbw at 24 points (l = 0–19, both edges, up to the corners) | pass: relative difference ≤ 4e−36. arbw needed 512–1024 bits at the high-frequency points |
| C1: D4 (W₀(0) ∋ 0; Cauchy W₁ ∋ i(r₊² + r_c²)/(r₊r_c)) | pass |
| C2: known zeros counted (P_B l = 1, 2 PS; P_B l = 0 NE; P_A l = 1 PS) | pass: 1, 1, 1, 1 |
| C3: negative control, rectangles sharing an edge 0.0103κ₋ above the l = 1 PS mode | pass: 1 and 0 |
| C4: the left half of ∂R'_2 computed explicitly against the right half | pass: 0 ± 5.3e−4 against 0 ± 3.3e−4. An extra, unregistered check at l = 0 gives π(1.000 ± 2.9e−4) on the left half, matching the right half |
| C5: float consistency at l = 0, 1, 2, 5 | pass: 1, 0, 0, 0, as in S0 |
| C6: the l = 0 Krawczyk enclosure (zsol) overlaps KAT-6 | pass |
| C7: replay (exact re-sum, coverage, 5% recomputed sample per l) | pass: `part3 replay: OK` for l = 0–19; 132 segments recomputed; every increment overlaps |

**Deviation (C0 implementation, commit `3e7fca05`).**
- As first run, C0 failed in the reference code.
- S1's `arbw` stops its Frobenius series before its own tail condition N > 2|ρ| + 1 holds at large
  |ω|, and then raises.
- At 256 bits its ball is also too wide for an overlap to mean anything.
- C0 now calls arbw unchanged, through a wrapper that raises the minimum term count, at the lowest
  precision giving a relative radius below 1e−10. It now requires both balls below 1e−10 as well
  as overlapping, which is stricter than the pre-registered "the balls overlap".
- The fix touched only the control. The pipeline code is unchanged since `89a2336b`, apart from a
  diagnostic (`c7d1bfbc`, the lower bound of |W| per piece).

## 4. The l = 0 enclosure

The zsol Krawczyk test (`certificates/part3/krawczyk_PB_l0.json`) gives:
- a unique zero of W₀ in the square of half-width 10⁻²⁰|ω_*| about
  ω_* = −0.909380978994805455524933762622i κ₋ (± 3.2e−31);
- contraction 5.3e−15.

This independently re-derives KAT-6 (P2-K6), with a different Wronskian code. It gives the upper
side of β: β ≤ 0.9093809789948054555249337627. The two-sided enclosure (P2-6) was not attempted;
it remains the stretch goal.

## 5. Compute

| item | CPU |
|---|---|
| pilot (surveys, prototypes, l = 0–2 trial windings, timing samples), before the pre-registration | ≈ 0.3 core-h |
| control gate C0–C3, C6 (first attempt stopped in C0, see §3) | ≈ 0.1 core-h |
| production, 2351 segments | 3.82 core-h (57 min wall, 4 workers) |
| C4, C7 and the extra l = 0 symmetry check | ≈ 0.3 core-h |
| tests (fast 18 s, 76 pass; `P2_SLOW=1` 418 s, 89 pass) | ≈ 0.15 core-h |
| **S3 total** | **≈ 4.7 core-h** (cap 24) |

The cost model of the pre-registration (3.4 core-h) was 12% low, at 3.82 core-h; contention on
the loaded machine accounts for part of the difference.

## 6. What S3 does not establish

- **P2-3 itself.** It needs every l < l₀, and l₀ comes from S4. If S4 ends with l₀ > 20, the
  pipeline runs unchanged for l = 20, …, l₀ − 1, at a cost of about 2.6a_l + 0.12a_l² CPU-s per l.
- **P_A.** It was not counted, apart from the C2(d) control.
- **The two-sided β (P2-6).**

## 7. Suggested CLAIMS text (for the coordinator; not applied)

**New row:**

| P2-3B | P2 | At P_B (ΛM² = 3/50, Q = (999/1000)Q_max, δ = 1/20) and for each l = 0, 1, …, 19, let R'_l = [−a_l, a_l] × [−b, e], with a_l the P2-2 threshold and the dyadics b ≥ (½+δ)κ₋ (b/κ₋ = 0.55000001) and 0 < e ≤ κ₋/10. Then W_l is analytic on a neighbourhood of R'_l and has no zero on ∂R'_l. Its winding number around ∂R'_l is 1 for l = 0 (the simple zero at ω = 0) and 0 for 1 ≤ l ≤ 19. With P2-0 and P2-2: for every l ≤ 19, S_δ contains no QNM of degree l other than ω = 0 (l = 0). This is P2-3 for l ≤ 19 only: P2-3 itself needs every l < l₀, and l₀ comes from S4 | CERTIFIED (Arb, 128-bit; not yet independently reviewed). Cauchy-interpolation discs on the right half of ∂R'_l (conjugation symmetry), 2351 segments; every winding ball lies within 7e−4 of its integer. W comes from the phase-factored Jost solutions, with the identities checked in sympy, and is cross-checked against the S1 code at 24 points (relative difference ≤ 4e−36). Controls: four known zeros counted (P_B l = 0, 1, 2; P_A l = 1); a 1/0 pair at 0.0103κ₋ from a PS mode; the symmetry checked at l = 0 and 2; a recomputed 5% sample. Written steps are routine: the argument principle, the conjugation symmetry, and analyticity off E₊ ∪ E_c (P2-D1′) | problems/P2/src/p2/zsol.py, part3.py, zeq_identity.py; problems/P2/certificates/part3/; problems/P2/notes/s3-preregistration.md, s3-report.md; problems/P2/tests/test_part3.py | 2026-10-17 |

**P2-K6, evidence amendment:** append "Re-derived in S3 with an independent Wronskian code (the
phase-factored solutions of problems/P2/src/p2/zsol.py): a unique zero of W₀ at
−0.909380978994805455524933762622i κ₋ (± 3.2e−31), contraction 5.3e−15
(problems/P2/certificates/part3/krawczyk_PB_l0.json)".

## 8. Commits

| commit | subject |
|---|---|
| `89a2336b` | P2: pre-register S3 (part iii counts at P_B, l = 0–19) with the phase-factored W and Cauchy discs |
| `3e7fca05` | P2: S3 control C0 calls arbw with enough Frobenius terms and enough precision; record run commits |
| `c7d1bfbc` | P2: S3 control gate passes (C0–C3, C6); lower bound of abs(W) per piece from the rectangle |
| `44904fbb`, `a3927dfb`, `ffa08710`, `b216311b`, `9c5714ac`, `a3b12803` | the certificates as they landed (l = 19; 17–18; 15–16; 12–14; 8–11; 0–7) |
| `1abe9d71` | P2: S3 report and the C4 symmetry check; verdict PASS at P_B for l = 0–19 |
| `09013c78` | P2: FORMULATION v0.5 records the S3 part (iii) implementation and result |

## 9. Reproduction (cwd `problems/P2/src`)

```
python -m p2.zeq_identity                                # sympy Z1-Z4 (<1 s)
python -m p2.part3 controls --workers 4                  # C0-C3, C6 (~1.5 min wall)
nohup nice -n 10 python -m p2.part3 run --point PB --l 0-19 --workers 4 &   # resumable
python -m p2.part3 assemble --point PB --l 0-19          # certificates/part3/winding_PB_l*.json
python -m p2.part3 symmetry --workers 4                  # C4
python -m p2.part3 verify --point PB --l 0-19            # C7
python -m pytest problems/P2/tests/test_part3.py -q      # from the repo root
```
