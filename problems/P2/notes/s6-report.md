# P2 S6 report: the two-sided β at P_B (P2-6)

Status: S6 exit, 2026-10-27. Branch `p2/s6-beta`, from main at `300c22a7`.
- Plan: `notes/s6-preregistration.md` (`c98ca16d`), with Amendments A1–A6 (`f82114bf`, dated
  2026-10-27, before any production run).
- Written proof: `notes/s6-proof.md`. Lemma E is written out in full there (§2, per A1).
- Machine check: `src/p2/beta_assembly.py`, which writes `certificates/beta/assembly_beta_PB.json`.
- Tests: `tests/test_beta.py`.

## 1. Verdict

- **P2-6 at P_B holds, at WRITTEN_PROOF**, pending the independent review.
  - The forward problem on Ω₂₃ has exactly two resonances with Im σ ≥ −0.92κ₋, namely 0 and ω_*.
    Both are of degree l = 0, and each is a simple zero of W₀. In D1's terms, the QNMs of any degree
    in S_e are 0 and ω_*.
  - **β = −Im ω_*/κ₋ ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹]**, two-sided, of width
    6.4·10⁻³¹ ≤ 10⁻⁸. The ball is P2-K6's, reused unchanged (A4).
- **Lemma E holds for every l.** −iκ₊ is not a resonance, and every W_l has a simple pole there.
  The levels are CERTIFIED at l = 0, WRITTEN_PROOF and CERTIFIED for 1 ≤ l ≤ 19, and WRITTEN_PROOF for
  l ≥ 20 (§5). D1's QNMs and HV17's resonances therefore agree in S_e.
- **Every pre-registered expectation held.**
  - The windings are 0 (l = 0) and −1 (1 ≤ l ≤ 19).
  - Part (i) certified at h₀ = 2/41 (l₀ = 20), and part (ii′) certified for every l ≤ 19.
  - All KATs and controls behave as required.
  - No kill criterion fired, and neither the K4 nor the K5 contingency was used.
- **The level is WRITTEN_PROOF**, the minimum of the ingredients, as pre-registered (§6 there).
  FORMULATION §2's target, CERTIFIED, is not reached.
- **Compute.** About 6.3 core-h in Phase 2, and 6.8 core-h for S6 with Phase 1. The projection was
  9–11 core-h and the cap 48.

## 2. Commits (branch `p2/s6-beta`)

| commit | content |
|---|---|
| `ec568eb8` | Phase-1 probes (`src/p2/s6_probe.py`, `results/s6/`) |
| `c98ca16d` | the pre-registration (Phase-1 checkpoint) |
| `f82114bf` | Phase 2 start: Amendments A1–A6; FORMULATION v1.1 (approved) |
| `7168b2e4` | production code: `beta_common`, `beta_epole`, `beta_part1`, `beta_part2`, `beta_part3` |
| `b7947fa8` | certificates: Lemma E, part (ii′), part (i) (T, A, S, L, cover, controls); the L-row provenance record |
| `5ee4796f` | the assembly check `beta_assembly.py`; `tests/test_beta.py`; the band-run records |
| `5a78bbcb` | certificates: the 20 band counts, KATs and controls, the 5% replay, B7 |
| `6cb88bca` | the assembly certificate |
| head | `notes/s6-proof.md`, this report, and a test fix (below) |

The test fix concerns the P2_SLOW test that recomputes the S+ box [0.06, 0.08]Ω_P. It now compares
the float design margin, which is used only to choose the design, to 10⁻⁹ relative; it had differed
in the 13th digit. The certified record is compared exactly, and it is identical.

## 3. Certificates, with the commit each was computed at

All are in `certificates/beta/`. "Computed at" is the commit tag the certificate records, and none
is `+dirty`. Every one passes.

| certificate | content | computed at | committed in |
|---|---|---|---|
| `epole_PB.json` | Lemma E: sympy first step; W̃_l(ω₀) for l ≤ 19 and samples; V_l > 0; B5 | `7168b2e4` | `b7947fa8` |
| `part2/A_fixed_PB.json` | part (ii′), A at fixed h_l, l ≤ 19 | `7168b2e4` | `b7947fa8` |
| `part1/T_PB.json`, `A_PB.json` | part (i): T6 (10 boxes), A (ζ_A = 6/5) | `7168b2e4` | `b7947fa8` |
| `part1/S_PB_plus.json`, `S_PB_c.json` | part (i): S, both sides, extended grid | `7168b2e4` | `b7947fa8` |
| `part1/L_PB.json` | part (i): L, band rows, pole-free series (1280 boxes) | `7168b2e4` (178 rows recomputed clean, §13) | `b7947fa8` |
| `part1/cover_PB.json` | part (i): the cover at depth e | `7168b2e4` | `b7947fa8` |
| `part1/controls_PB.json` | NB3, NB4, NB5 | `7168b2e4` | `b7947fa8` |
| `part3/band_PB_l{0..19}.json` | the band counts | `7168b2e4` (run); assembled at `5ee4796f` | `5a78bbcb` |
| `part3/kat_controls_PB.json` | B1–B4, B6, NB1, NB2, the left halves | `5ee4796f` | `5a78bbcb` |
| `part3/replay_PB.json` | the 5% replay of the band pieces | `5ee4796f` | `5a78bbcb` |
| `part1/B7_PB.json` | B7, the float check (A6) | `5ee4796f` | `5a78bbcb` |
| `assembly_beta_PB.json` | the assembly check | `5a78bbcb` | `6cb88bca` |

The S2–S5 certificates read are unchanged and equal to their blobs on main `300c22a7`:
- `assembly/assembly_PB.json`;
- `params.json`;
- `part3/krawczyk_PB_l0.json`;
- `kat6_winding_krawczyk.json`;
- `part3/zeq_identity.json`;
- `part1/bulk_PB.json`.

## 4. Band windings per l

B_l = [−a′_l, a′_l] × [−b′, −b₀], with b′ = 21012395/2²⁹ (eκ₋ rounded up; b′ − eκ₋ = 3.9·10⁻¹⁰ absolute) and
b₀ = 22839559/2³⁰ (κ₋/2 rounded down). a′_l is part (ii′)'s threshold, identical by the assembly
check. **Every winding is as pre-registered: 0 at l = 0 and −1 for 1 ≤ l ≤ 19.** With the pole
bound P ≤ 1, these give Z = 1 at l = 0 (ω_*, simple) and Z = 0 otherwise, with P = 1 throughout.

| l | winding | (sum over H)/π | Z, P | pieces | a′_l/κ₋ | min abs(W_l) on H | K_max | CPU (s) |
|---|---|---|---|---|---|---|---|---|
| 0 | 0 | [+/- 8.69e-5] | 1, 1 | 41 | 12.912 | 0.0321 | 64 | 69 |
| 1 | −1 | [-1.000 +/- 1.09e-4] | 0, 1 | 43 | 13.650 | 0.351 | 80 | 84 |
| 2 | −1 | [-1.00 +/- 2.28e-4] | 0, 1 | 53 | 18.446 | 0.477 | 88 | 118 |
| 3 | −1 | [-1.000 +/- 1.81e-4] | 0, 1 | 63 | 23.241 | 0.6 | 88 | 160 |
| 4 | −1 | [-1.00 +/- 2.68e-4] | 0, 1 | 71 | 27.668 | 0.71 | 88 | 201 |
| 5 | −1 | [-1.00 +/- 2.05e-4] | 0, 1 | 81 | 32.464 | 0.817 | 88 | 254 |
| 6 | −1 | [-1.000 +/- 2.62e-4] | 0, 1 | 92 | 38.367 | 0.913 | 88 | 324 |
| 7 | −1 | [-1.00 +/- 2.90e-4] | 0, 1 | 101 | 42.425 | 0.986 | 88 | 381 |
| 8 | −1 | [-1.000 +/- 3.32e-4] | 0, 1 | 112 | 48.082 | 1.1 | 88 | 466 |
| 9 | −1 | [-1.000 +/- 4.59e-4] | 0, 1 | 123 | 53.738 | 1.13 | 96 | 554 |
| 10 | −1 | [-1.00 +/- 5.02e-4] | 0, 1 | 128 | 56.813 | 1.27 | 96 | 616 |
| 11 | −1 | [-1.000 +/- 5.37e-4] | 0, 1 | 141 | 62.223 | 1.31 | 96 | 729 |
| 12 | −1 | [-1.00 +/- 4.83e-4] | 0, 1 | 151 | 67.634 | 1.43 | 96 | 842 |
| 13 | −1 | [-1.000 +/- 5.16e-4] | 0, 1 | 163 | 73.045 | 1.5 | 96 | 979 |
| 14 | −1 | [-1.00 +/- 5.34e-4] | 0, 1 | 173 | 78.455 | 1.61 | 96 | 1110 |
| 15 | −1 | [-1.000 +/- 6.59e-4] | 0, 1 | 182 | 83.866 | 1.68 | 96 | 1253 |
| 16 | −1 | [-1.000 +/- 6.03e-4] | 0, 1 | 195 | 89.277 | 1.78 | 96 | 1429 |
| 17 | −1 | [-1.000 +/- 5.97e-4] | 0, 1 | 205 | 94.688 | 1.82 | 104 | 1597 |
| 18 | −1 | [-1.00 +/- 6.84e-4] | 0, 1 | 216 | 100.098 | 1.93 | 104 | 1792 |
| 19 | −1 | [-1.00 +/- 7.01e-4] | 0, 1 | 228 | 105.509 | 1.96 | 104 | 2004 |

- All 2562 planned segments certified at their planned discs: no splits, no K5 retry, and K ≤ 104
  against K_MAX = 160.
- The run: launched at `7168b2e4` (clean), 3756 s wall with 4 workers at nice 10, 14 963 s of
  segment CPU. Assembled at `5ee4796f` (clean).
- The smallest abs(W₀) on H is 0.032. It is attained near ω_*, which lies 0.0106κ₋ above B₀'s
  bottom edge.

## 5. Lemma E: status per l range

Lemma E: ω₀ = −iκ₊ (−0.8229iκ₋, the one point of E₊ ∪ E_c in S_e) is not a resonance of the
forward problem on Ω₂₃ for any l, and every W_l has a simple pole there, with
Res_{ω₀}W_l = i(2κ₊r₊ + L)/(2r₊²)·W̃_l. The proof is in `s6-proof.md` §2, with the normalisation of A1
(Z₊ = Σa_m t^m, a₀ = 1, t = r − r₊; Z2 = t + Σ_{m≥2} b_m t^m; Z_c(r_c) = 1).

| l | what holds | level |
|---|---|---|
| 0 | Arb: W̃₀ = +0.038285581814409 ± 4.4·10⁻¹⁶. The band count (winding 0 with ω_* inside) forces the pole again. The written argument does not apply, because V₀ = FF′/r changes sign | **CERTIFIED** (the Frobenius steps are routine; the first is checked in sympy) |
| 1–19 | written (positivity, §2.4) and Arb (W̃_l < 0, from −22.80 at l = 1 to −2.96·10⁵³ at l = 19); the band counts (winding −1) force the pole | **WRITTEN_PROOF and CERTIFIED** |
| ≥ 20 | written (§2.4). Part (i)'s L certificate at depth e separates the two solutions at ω₀ for every h ≤ h₀ (Re y₊ ≥ 0.1617 > 0 > −0.1201 ≥ Re y_c in the grid box containing ω₀; W5, W6′). Arb samples at l = 20, 25, 30, 40, 50, 60, 100 | **WRITTEN_PROOF** |

New in the write-up (A1), beyond the pre-registration's sketch:
- **The closed form.** W̃_l = −(1/r₊)∫_{r₊}^{r_c}(L + rF′)Z2 dr/r (`s6-proof.md` §2.4(f)). It follows
  from the same Wronskian identity, and it makes the sign for l ≥ 1 immediate (L + rF′ > 0 and
  Z2 > 0).
- **A float cross-check of the identification** (not used by the proof). Integrating Z2 with scipy
  reproduces the Arb W̃_l at l = 0, 1, 2, 5 to < 10⁻⁸ relative, including the sign of W̃₀
  (`tests/test_beta.py::test_lemma_E_closed_form`). So the normalisations, the r*-constants and
  the sign conventions of §2.4 match zsol's Wronskian.
- **The residue's normalisation.** KAT B5 tests exactly this formula, in this normalisation. The
  relative error of ε·W_l(ω₀ + ε) against the residue is 9.16·10⁻⁵, 1.42·10⁻⁵, 4.49·10⁻⁵ and
  6.76·10⁻⁵ at ε = 10⁻⁵κ₋ (l = 0, 1, 5, 19). It is exactly 100 times smaller at 10⁻⁷κ₋, and zsol and
  arbw agree to every printed digit.

**What Lemma E is needed for** (`s6-proof.md` §1, item 6). For D3's β (D1's QNMs, which exclude
E₊ ∪ E_c by definition), only the bound "at most a simple pole at ω₀" (§2.3) is needed: it turns the
windings into zero counts. The non-resonance half is what makes the resonance statement hold, and
D1's QNMs agree with HV17's resonances in S_e.

## 6. KATs and controls

| item | test | required | result | certificate |
|---|---|---|---|---|
| B1 | replay of P2-K6's Krawczyk ball (not regenerated; A4) | unique zero; as stored | certified; the identical enclosure −0.909380978994805455524933762622i ± 3.16·10⁻³¹ (contraction 5.3·10⁻¹⁵) | `part3/kat_controls_PB.json` |
| B2 | l = 0, [−0.03, 0.03] × [−0.93, −0.89] (ω_* only) | +1 | +1 | same |
| B3 | l = 0 and l = 1, [−0.03, 0.03] × [−0.86, −0.79] (ω₀ only) | −1 | −1, −1 | same |
| B4 | l = 0, [−0.03, 0.03] × [−0.93, −0.79] (ω₀ and ω_*) | 0 | 0 | same |
| B5 | ε·W_l(ω₀ + ε) against the residue, l = 0, 1, 5, 19, zsol and arbw | agreement to O(ε) | relative error ≤ 9.2·10⁻⁵ at ε = 10⁻⁵κ₋, 100 times smaller at 10⁻⁷κ₋; the codes agree | `epole_PB.json` |
| B6 | l = 0 overtone ≈ −1.855iκ₋; l = 1 PS mode ≈ 7.2149 − 1.7533iκ₋ (no pole inside) | +1 each | +1, +1 | `part3/kat_controls_PB.json` |
| B7 | float, 100 random band points, l ∈ [20, 1000] (A6) | no zero | no zero: none passes the A6 rule; 20 secant searches stopped inside the band, each with abs(W(stop))/abs(W(seed)) = 1.000001, i.e. without moving to a root | `part1/B7_PB.json` |
| NB1 | l = 0 band box with its bottom at −0.905κ₋ (above ω_*) | winding −1 | −1 | `part3/kat_controls_PB.json` |
| NB2 | an edge at Im ω_* to 10 digits | not certify | not certified ("piece enclosure contains 0 at the maximal depth") | same |
| NB3 | L's ω₀ box with S4's series | not certify | not certified ("affine form may vanish") | `part1/controls_PB.json` |
| NB4 | S4's T design at depth e on [0.8, 0.85]Ω_P | not certify | μ = −0.0475a* | same |
| NB5 | A at h₀ = 2, ζ_A = 5 | not certify | μ = −0.392a* | same |
| NB6 | the assembly check under mutations | fail on the named items | NB6a–NB6i all fail as required (§7) | `assembly_beta_PB.json` |

Also checked:
- **The symmetry.** The left-half sums at l = 0 and l = 2 overlap the stored right-half sums.
- **The 5% replay** (`part3/replay_PB.json`). 141 stored pieces (at least 3 per l) recomputed with
  their stored K all overlap: 827 s CPU, 211 s wall with 4 workers.

## 7. The assembly check

`python -m p2.beta_assembly` writes `certificates/beta/assembly_beta_PB.json` (`6cb88bca`, computed at
`5a78bbcb`). It runs 332 checks at 256 bits with exact rationals, and all pass.
- **Pins (A2).** The six S2–S5 inputs equal their blobs on main `300c22a7`, and P2-4's 54 pinned inputs
  are unchanged; 38 inputs are pinned in all.
- **C1.** S_e ∩ (E₊ ∪ E_c) = {ω₀}: 2κ₊ − eκ₋ = 0.7258κ₋ and κ_c − eκ₋ = 1.2751κ₋, with ω₀ 0.0971κ₋
  inside. b′ − eκ₋ = 3.9·10⁻¹⁰, and κ₋/2 − b₀ = 5.2·10⁻¹⁰.
- **C2.** P2-4's assembly passes, with 366 checks and its controls.
- **C3.** Lemma E: W̃_l ≠ 0 for every l ≤ 19, the positivity inequalities, B5, and the normalisation
  statement.
- **C4.** Part (ii′) for every l ≤ 19, at depth e, with W7, Jost decay and a′_l ≥ ζ_AΩ_P/h_l.
- **C5.** For every l ≤ 19:
  - a′_l is identical to C4's;
  - b′ ≥ eκ₋;
  - the stored pieces (2562 in all) trace the right half of ∂B_l exactly, lie in their discs, and
    clear E₊ ∪ E_c;
  - the re-sums give 0 and −1, consistent with C3.
- **C6.** Part (i): T, A, S± and L at depth e and h₀ = 2/41. The cover replays identically, L's rows
  are [−19/20, −11/20]κ₋, and the regime overlaps hold (S reaches 82.75κ₋ above L).
- **C7.** l₀ = 20.
- **C8.** P2-K6's ball is certified and consistent with KAT-6. ω_* lies in B₀, 0.0106κ₋ above its
  bottom edge, off E₊ ∪ E_c and below S_δ. β ≤ 0.90938098, consistent with P2-4. The width is
  6.3·10⁻³¹.
- **C9.** The KATs, controls, replay and B7 are as recorded.
- **C10.** The minimum level is WRITTEN_PROOF, and no `+dirty` tag appears.
- **Negative controls.** All nine fail on their named items:

  | control | mutation | fails on |
  |---|---|---|
  | NB6a | `band_PB_l7.json` removed | C5.l7.present, C5.all_l_present |
  | NB6b | a′₇ lowered by 2⁻²⁴ | C5.l7.a_identical |
  | NB6c | the l = 3 winding set to 0 | C5.l3.winding_value, C5.l3.resum (and the C3 consistency) |
  | NB6d | b′ → the dyadic just below eκ₋ (A5) | C5.l{0..19}.b_prime_ge_e |
  | NB6e | W̃₅ replaced by a ball containing 0 | C3.l5.Wt_nonzero, C5.l5.pole_forced_consistent_C3 |
  | NB6f | part (ii′) at l = 9 only to depth 0.705κ₋ | C4.l9.pass |
  | NB6g | h₀ = 2/43 | C7.no_l_gap (and 7 more) |
  | NB6h | an S5 input off its main blob | C1.pins_main_blobs |
  | NB6i | P2-K6's ball widened to 2·10⁻⁷ | C8.beta_width, C8.consistent_with_P2_4 |
- **Verify.** `python -m p2.beta_assembly verify` is identical, and so is `python -m p2.assembly verify`
  (S5).
- **Tests.**
  - The fast P2 suite: 102 passed, 21 skipped (66 s). That is S5's 90 + 17, plus `test_beta.py`'s
    12 fast tests and 4 slow ones.
  - `P2_SLOW=1` on `tests/test_beta.py`: 16 passed, after the test fix in §2.

## 8. Compute

Phase 2. All runs used at most 4 workers at nice 10, and the long runs were detached and
resumable.

| item | CPU (s) | wall | workers |
|---|---|---|---|
| Lemma E (`epole_PB.json`) | 1.2 | 1.3 s | 1 |
| part (ii′), A at fixed h (20 l) | 65 | 65 s | 1 |
| part (i): T, A, cover (twice), controls NB3–NB5 | 29 | < 1 min | 1 |
| part (i): S+ and Sc | 230 | 3 min | 1 + 1 |
| part (i): L (first pass, plus the recomputation of 178 rows) | 288 | 4 min + 50 s | 1 |
| part (iii): the band counts (2562 segments) | 14 963 | 62.6 min | 4 |
| part (iii): KATs, NB1, NB2, left halves, B1 | 390 | 106 s | 4 |
| part (iii): the 5% replay (141 pieces) | 826 | 211 s | 4 |
| B7 (100 float points, l ∈ [20, 1000]) | 5457 | 33 min | 4 |
| assembly runs, verifies, tests (fast suite, P2_SLOW `test_beta`), development checks | ≈ 300 (estimate) | — | 1 |
| **total, Phase 2** | **≈ 22 500 s ≈ 6.3 core-h** | ≈ 2.2 h of session | |

- Phase 1 used about 0.5 core-h, so S6 used about 6.8 core-h.
- The pre-registered projection was 9–11 core-h. Part (iii) took 4.2 core-h against the projected
  5–7.
- B7 took 1.5 core-h, more than the part-(i) float items' 1.0 core-h line. The large-l T and A
  points cost up to 290 s each in floats.
- No projection ever approached 40 core-h, so "main" was not messaged.

## 9. Ledger proposals (A3): paste-ready, not applied

The column order is `| id | problem | statement | level | evidence | updated |`. Nothing here is
applied. The owner enters rows after the independent review and the user's approval.

**Should Lemma E be its own row?** Yes; I propose a row **P2-E** (the id is the owner's choice).
- Lemma E is a statement about the forward problem at one point, independent of the band and of β.
- Its level differs by l range: CERTIFIED at l = 0; WRITTEN_PROOF (and, for l ≤ 19, CERTIFIED) for
  l ≥ 1. Folded into P2-6, that profile would be lost in P2-6's minimum.
- It complements P2-D1′, which excludes E₊ ∪ E_c. It is what makes D1's QNMs and HV17's resonances
  agree in S_e, and a later use of [HV17] with α up to −Im ω_* would cite exactly that.
- It can be reviewed on its own.

The alternative, folding it into P2-6 as an ingredient, is also consistent: P2-6's level is
WRITTEN_PROOF either way.

**P2-6 (new).**

```
| P2-6 | P2 | At P_B (ΛM² = 3/50, Q = (999/1000)Q_max): the forward problem on Ω₂₃ (HV17 §2.5) has exactly two resonances σ with Im σ ≥ −(23/25)κ₋ (= −0.547γ₀), namely σ = 0 and σ = ω_*, P2-K6's l = 0 zero, ω_*/κ₋ ∈ −0.909380978994805455524933762622i ± 3.2·10⁻³¹; both are of degree l = 0 and simple zeros of W₀. In D1's terms, the QNMs of any degree in that strip are 0 and ω_*. The one point of E₊ ∪ E_c in the strip, −iκ₊, is not a resonance (P2-E). Hence D3's gap is attained at ω_*: β = −Im ω_*/κ₋ ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹] (two-sided; width 6.4·10⁻³¹ ≤ 10⁻⁸), and D1's QNMs and HV17's resonances agree in the strip. Assembled from P2-4 (S_δ), the conjugation symmetry (Re σ < 0), part (i) at depth 23/25 for l ≥ 20, and for l ≤ 19 regime A at fixed h = 1/(l+½) (Re σ ≥ a′_l) with pole-aware argument-principle counts on the band boxes [−a′_l, a′_l] × [−(23/25)κ₋, −κ₋/2]; with P2-E, P2-D1′'s set clause and P2-K6. The resonance multiplicity at ω_* (below −κ₊, outside P2-D1′'s window) is not claimed. Not P_A | WRITTEN_PROOF, the minimum of the ingredients: P2-4, P2-D1′ (set clause), part (i) at depth 23/25 (W1, W2, W5, W6′ written) and part (ii′) (W1, W2, W4) are WRITTEN_PROOF; the band counts, P2-E for l ≤ 19 and P2-K6's ball are CERTIFIED; FORMULATION §2's target CERTIFIED is not reached. Band counts (S3's Cauchy-disc machinery, every disc box clear of E₊ ∪ E_c; 2562 stored pieces; a 5% replay overlaps): winding 0 at l = 0 (ω_* and the simple pole at −iκ₊) and −1 for 1 ≤ l ≤ 19 (the pole, no zero). Part (i) at depth 23/25 with h₀ = 2/41 (l₀ = 20): L on the band rows with the pole-free series R̃ = (1+2ρ)R (W6′), S with an extended design grid, T redesigned, A at ζ_A = 6/5. Machine-checked bookkeeping (exact rationals and Arb, 256-bit; 332 checks; the six S2–S5 inputs pinned by sha256 and equal to their blobs on main 300c22a7; P2-4's 54 pinned inputs unchanged): a′_l identical between the band boxes and the A-at-fixed-h certificates; b′ − (23/25)κ₋ = 3.9e−10 > 0; the stored discs trace the right half of each band box exactly; the re-summed windings; part (i)'s cover replays identically; ω_* lies in B₀, 0.0106κ₋ above its bottom edge. Negative controls NB6a–NB6i fail as required (a band certificate removed; a′₇ lowered by 2⁻²⁴; a wrong winding; b′ just below (23/25)κ₋; and five more). KATs: B1 (P2-K6's ball replayed), B2 +1, B3 −1 (the pole; l = 0, 1), B4 0, B5 (the residue, with the independent S1 code), B6 +1 (two zeros below the edge); B7 (float, 100 band points with l ∈ [20, 1000]) finds no zero; controls NB1–NB5 fail as required. Independent review: pending | problems/P2/src/p2/beta_common.py, beta_epole.py, beta_part1.py, beta_part2.py, beta_part3.py, beta_assembly.py; problems/P2/certificates/beta/ (assembly_beta_PB.json); problems/P2/notes/s6-preregistration.md, s6-proof.md, s6-report.md; problems/P2/tests/test_beta.py | 2026-10-27 |
```

**P2-E (new; Lemma E).**

```
| P2-E | P2 | At P_B, the point ω₀ = −iκ₊ (≈ −0.8229iκ₋; the one point of E₊ ∪ E_c with Im ω ≥ −(23/25)κ₋) is not a resonance of the forward problem on Ω₂₃ (HV17 §2.5) in any sector l ≥ 0, and every W_l has a simple pole there: Res_{ω₀}W_l = i(2κ₊r₊ + L)/(2r₊²)·W̃_l ≠ 0. Here W̃_l is the Wronskian, in zsol's Jost-normalised form, of Z2 = t + Σ_{m≥2} b_m t^m and of the solution analytic at r_c. Z2 is the exponent-1 Frobenius solution at r₊; the exponent-0 one carries t·log t, since the first step reads 0·a₁ = 2κ₊r₊ + L. The normalisation is Z₊ = Σ a_m t^m, a₀ = 1, t = r − r₊. For l ≥ 1, W̃_l = −(1/r₊)∫_{r₊}^{r_c}(L + rF′)Z2 dr/r < 0 (V_l > 0 and Z2 > 0); W̃₀ = 0.038285581814409 ± 4.4·10⁻¹⁶. Hence D1's QNMs and HV17's resonances agree for Im ω ≥ −(23/25)κ₋. Not P_A | l = 0: CERTIFIED (Arb, 128-bit; routine written Frobenius steps, the first checked in sympy). 1 ≤ l ≤ 19: WRITTEN_PROOF (positivity, notes/s6-proof.md §2.4) and CERTIFIED (Arb W̃_l ≠ 0; P2-6's band counts force the pole). l ≥ 20: WRITTEN_PROOF; part (i)'s L certificate at depth 23/25 also separates the two solutions at ω₀ (W6′), and Arb samples at l = 20, 25, 30, 40, 50, 60, 100 are nonzero. KAT B5: ε·W_l(ω₀ + ε) tends to the residue as O(ε), with zsol and with the independent S1 code arbw (l = 0, 1, 5, 19). Float check, not used: the closed form reproduces the Arb W̃_l (l = 0, 1, 2, 5) to 10⁻⁸. No literature is cited. Independent review: pending | problems/P2/src/p2/beta_epole.py; problems/P2/certificates/beta/epole_PB.json; problems/P2/notes/s6-proof.md §2; problems/P2/tests/test_beta.py | 2026-10-27 |
```

**Notes for existing rows** (exact text, each appended to the end of the row's level column):
- **P2-K6** (needed: its statement calls the two-sided enclosure "the stretch goal P2-6"):
  ` Note 2026-10-27 (S6): P2-6 reuses this row's S3 ball (certificates/part3/krawczyk_PB_l0.json) unchanged; with it, β ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹] two-sided, which supersedes the one-sided bound for β (the zero itself is unchanged).`
- **P2-4** (needed, for the same pointer):
  ` Note 2026-10-27 (S6): the two-sided enclosure is P2-6, β ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹]; the bound 0.55 ≤ β ≤ 0.90938098 stated here stays true.`
- **P2-D1′** (optional; nothing in it changes):
  ` Note 2026-10-27 (S6): at P_B the one point of E₊ ∪ E_c with Im ω ≥ −(23/25)κ₋, −iκ₊, is not a resonance in any sector, and every W_l has a simple pole there (P2-E); so in that strip the resonances are exactly the zeros of the W_l (used by P2-6). No multiplicity is claimed below −κ₊.`

## 10. FORMULATION: a proposed exit entry (not applied)

The coordinator approved the v1.1 entry (S6 start), and it is applied (`f82114bf`). No exit entry
was authorised, so this one is a proposal in the style of v1.0 (S5 exit). It is for the owner to
apply, after the review, together with the header status line.

```
- **v1.2 (⟨date⟩, S6 exit; `notes/s6-proof.md`, `notes/s6-report.md`).** P2-6 at P_B is assembled.
  - **P2-6 at P_B: WRITTEN_PROOF**, the minimum of its ingredients. §2's target, CERTIFIED, is not
    reached.
    - The forward problem on Ω₂₃ has exactly two resonances with Im σ ≥ −0.92κ₋: 0 and ω_*. Each
      is a simple zero of W₀.
    - β = −Im ω_*/κ₋ ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹]. The frozen target (width
      ≤ 10⁻⁸) is met by P2-K6's ball, reused unchanged.
    - The bookkeeping is machine-checked, with negative controls:
      `certificates/beta/assembly_beta_PB.json` (`python -m p2.beta_assembly`).
  - **Lemma E** (row P2-E): −iκ₊ is not a resonance for any l, and every W_l has a simple pole
    there. So D1's QNMs and HV17's resonances agree in S_e, and D3's α is HV17's α₀ there.
  - **Unchanged:** the edge, the target, D1, D1′, D3, D7, the KATs of §5 (S6's B1–B7 are additive,
    v1.1) and the scope fences of §7.
  - **Open:** P2-7 (P_A).
```

Header status line, proposed: "S6 assembled P2-6 at P_B at WRITTEN_PROOF (`notes/s6-proof.md`;
Revision log v1.2); P2-7 (P_A) is still a target." 

## 11. What S6 does not establish

- **No resonance multiplicity at ω_*.** Only the radial order m₀(ω_*) = 1 is certified. ω_* lies
  below −κ₊, outside P2-D1′'s multiplicity window.
- **Nothing below −0.92κ₋.** The B6 KAT boxes (the l = 0 overtone near −1.855iκ₋ and the l = 1 PS
  mode near 7.2149 − 1.7533iκ₋) test the method; they are not claims.
- **Not P_A.** P2-7 is still a target, and a two-sided β at P_A is not targeted.
- **No new conversion.** With β known, [HV17] could be applied with any α < −Im ω_* = 0.9094κ₋ (also
  < γ₀). That would sharpen P2-5's H^{1.05−0} towards H^{1.409−0}. S6 does not claim it, and P2-5 is
  unchanged.
- **No CERTIFIED level.** P2-6 is WRITTEN_PROOF, below FORMULATION §2's target. CERTIFIED would need
  P2-0's written steps, W1/W2/W4/W5/W6′, D1′ and Lemma E's l ≥ 1 argument to be machine-checked.
- **Float evidence is not part of the proof:** the closed-form check of W̃_l and KAT B7.
- **The scope fences of FORMULATION §7 are unchanged:**
  - massless scalar on fixed RNdS only;
  - no nonlinear statement. β > 7/9 is now certified at P_B, but [CGNS18]'s other hypotheses are
    not addressed;
  - nothing about generic solutions being no smoother than H^{1/2+β};
  - no other parameters or neighbourhoods.
- **No ledger, FORMULATION-status or .tex change on this branch** beyond the approved v1.1 entry. §§9
  and 10 are proposals.

## 12. For the reviewer: what to attack first

1. **Lemma E at a point of E₊** (`s6-proof.md` §2.2). D1′'s characterisation of resonant states
   (smooth, H^s above threshold) is used at σ = ω₀ ∈ E₊, where D1′'s Lemma A does not apply. Attack:
   - the threshold 3/2 at r₊ against t·ln t ∈ H^{3/2−0};
   - the non-integer exponent κ₊/κ_c at r_c;
   - whether the log coefficient c = (2κ₊r₊ + L)/Δ₁ is right (it is what makes the exponent-0
     solution inadmissible).
2. **The positivity argument** (§2.4).
   - The limits W[Φ_c, Ψ] → 0, which uses κ₊ < κ_c, and W[Φ̃_c, Ψ] → 2κ₊r_ce^{κ₊c_c}.
   - Convergence of ∫V_lΦ₂Ψ at both ends.
   - The step from B ≠ 0 to W̃_l ≠ 0.
   - The closed form (E-int) is a compact target; its float check is in the tests.
3. **The counts with a pole** (§5).
   - P ≤ 1 (§2.3: only D₁ vanishes near ω₀).
   - The disc clearance near ω₀, which lies 0.0971κ₋ above the bottom edge and 0.32κ₋ below the top
     edge.
   - The l = 0 accounting: ω_* lies only 0.0106κ₋ above B₀'s bottom edge.
4. **Part (ii′)** (§4). Is S4's uniform-in-h bound valid at h = h₀ = h_l exactly, with W7 at depth e?
   Are the a′_l identities exact?
5. **Part (i) at depth e** (§3).
   - W6′: R̃ at ω₀ is the regular solution, and Re(n + 2ρ) > 0 for n ≥ 2 in S_e.
   - The extended S grid and the T redesign are design choices only, but they should be checked.
   - The cover inequalities at depth e.
6. **The seams.** The band against S_δ (−κ₋/2 against −0.55κ₋); b′ against eκ₋ (3.9·10⁻¹⁰
   absolute); l = 19 against l = 20 at h₀ = 2/41; Re ω = a′_l, held by both the band box and
   part (ii′).

## 13. Deviations, +dirty rows, denials

**Deviations from the pre-registration** (all additive, or provenance only):
- **178 L rows recomputed** at a clean commit. The first L run was launched at `7168b2e4` and
  imported its code at start. From 09:32 an untracked draft of `src/p2/beta_assembly.py`, which the L
  run does not import, made `git status` of `src/p2` non-empty, so 178 of the 1280 per-row commit
  tags read `7168b2e4+dirty`.
  - Those rows were moved to `results/beta/part1/L.w0.first-pass-dirty.jsonl` (kept, git-ignored).
  - The draft was moved out of `src/p2`, and the rows were recomputed at `7168b2e4` clean.
  - The recomputed leaves are identical (`results/beta/part1/L_rerun_clean_20261027.json`).
  - The cover was re-run clean, and `L_PB.json` and `cover_PB.json` carry only `7168b2e4`.
- **The assembly's negative controls** cover the four registered NB6 mutations (with A5's "just
  below"). They add five more: Lemma E at l = 5, part (ii′) at l = 9, h₀ = 2/43, an S5 pin, and a
  wide β ball.
- **The closed form (E-int)** and its float test are new. They are not used by the proof.
- **FORMULATION.** The v1.1 entry was approved and applied at `f82114bf`. The exit entry in §10 is a
  proposal only, since no exit entry was authorised.
- **No design fallback was needed.** T used no fallback design, every ζ_A(l) certified at its float
  value, S and L passed at the registered settings, and no K4/K5 contingency was used.

**+dirty rows (A2).** No certificate on this branch carries a `+dirty` commit tag: the assembly's C10
lists none. The only `+dirty` tags ever written are the 178 first-pass L rows above. They were
superseded by clean recomputations with identical results, and they are kept only in a git-ignored
results file.

**Denials.** No permission was denied during S6. Twice the harness blocked a plain `sleep`
(its rule against chained sleeps); I waited with its recommended monitors and condition loops
instead. Nothing else was affected.

**Kill criteria.** None fired:
- K1: W̃_l certified nonzero for every l ≤ 19, every band count consistent with one simple pole, and
  the written argument stands;
- K2: l₀ = 20, and the compute is far below 40 core-h;
- K3: no zero of any W_l in S_e other than 0 and ω_*, and ω₀ is not a resonance;
- K4: part (i) certified at h₀ = 2/41;
- K5: no band segment failed, so no retry was needed.
