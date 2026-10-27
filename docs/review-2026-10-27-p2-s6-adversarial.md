# Independent adversarial review: P2 stage S6 (the two-sided β at P_B, P2-6, and the exceptional point −iκ₊, P2-E)

Date: 2026-10-27. Reviewer: a separate agent with no part in authoring P2 S6.

Scope: branch `p2/s6-beta` at `172eb48c`, from main `300c22a7`. The S6 commits (2026-10-27, times −07:00):
- `ec568eb8` (09:12:19): the Phase-1 probes (`src/p2/s6_probe.py`, `results/s6/`);
- `c98ca16d` (09:16:17): the pre-registration `notes/s6-preregistration.md`;
- `f82114bf` (09:21:53): Amendments A1–A6 and FORMULATION v1.1;
- `7168b2e4` (09:28:06): the production code (`beta_common`, `beta_epole`, `beta_part1`, `beta_part2`, `beta_part3`);
- `b7947fa8` (09:37:17): the Lemma E, part (ii′) and part (i) certificates, and the L-row provenance record;
- `5ee4796f` (10:39:14): `beta_assembly.py`, `tests/test_beta.py`, the band-run records;
- `5a78bbcb` (11:18:48): the 20 band counts, the KATs and controls, the 5% replay, B7;
- `6cb88bca` (11:19:13): the assembly certificate;
- `172eb48c` (11:25:02): `notes/s6-proof.md`, `notes/s6-report.md`, and a test fix.

The review worked in its own worktree (`.claude/worktrees/p2-s6-review`, branch `p2/s6-review`, at `172eb48c`)
and changed no tracked file. It adds this document and `problems/P2/certificates/review_s6/` (scripts,
outputs and a README with the replay commands). For the `+dirty` check it read, without writing, two
git-ignored files in the author's worktree. No permission was denied. Once the harness blocked a
foreground `sleep` (its rule against chained sleeps); the review waited with a condition loop instead.

**Sources.** Fetched as TeX from arXiv on 2026-10-27 and read by the reviewer:
- **[HV17]** arXiv:1512.08004**v2** (posted 2020-05-27; TeX date line "December 25, 2015. Final revision: July
  12, 2017."), `cauchyhor-r.tex`, sha256 `a0d3302f…733cf2be` (the file the S5 review read). Read: §2.2
  (the radial sets, β_j and κ_j, ℓ. 836–870), Def 2.6 and the radial-point Proposition (ℓ. 1026–1067), §2.4
  through Lemma 2.14 (ℓ. 1200–1270), §2.5 (ℓ. 1313–1345), and the forms remark (ℓ. 1549).
- **[Vas13]** arXiv:1012.4391**v2** (updated 2011-11-09; TeX date line "Revised: November 1, 2011"),
  `kerr-ah-b.tex`, sha256 `b7ead254…b83bc18`. Read: §2.1–2.6 (the assumptions, the radial-point
  estimates and their proof, complex absorption, the global estimates, the Fredholm theorem
  `thm:classical-absorb` (ℓ. 1929–1941), Remark `rem:inv-indep-of-s`), and §6.1–6.3 (the Kerr–de Sitter geometry and β±, ℓ. 5613–6000).
- 16 passages used below are checked verbatim, with line numbers, by `review_s6/quotes_check.py`.

**Overall verdict: P2-6 at P_B is proven as stated, at WRITTEN_PROOF. Lemma E (P2-E) is proven as stated,
but its level must be stated by part: the simple pole of W_l at −iκ₊ (W̃_l ≠ 0) is CERTIFIED for l ≤ 19
and WRITTEN_PROOF for l ≥ 20, and the statement that −iκ₊ is not a resonance is at P2-D1′'s level,
WRITTEN_PROOF, for every l including l = 0 (S1). Nothing is blocking.**
- **The decisive question is settled in the author's favour.** At σ = ω₀ = −iκ₊ the localized forward
  problem is an index-0 analytic Fredholm family on spaces of any order above threshold, and [HV17] allows
  arbitrarily high regularity at the event and cosmological horizons. So a pole is exactly a nonzero
  kernel, kernel elements are smooth, and they separate into spherical harmonics. The exponent-0
  Frobenius solution at r₊ carries c·t·log t with c = (2κ₊r₊ + L)/Δ′(r₊) ≠ 0, so it is not even C¹. At
  r_c the second exponent is κ₊/κ_c = 0.3749, not an integer. Hence ω₀ is a resonance in sector l iff
  W̃_l = 0. There is no generalized resonant state without a resonant state, no pole without a kernel,
  and no state supported beyond the horizons.
- **Lemma E checks out independently.**
  - sympy re-derives the phase-factored equations, D_m, the first step, the log coefficient, the
    residue and every identity behind the closed form.
  - The reviewer's Arb code (no p2) recomputes W̃_l for l = 0, 1, 2, 3, 5, 10, 19, 20, 50, 100.
    W̃₀ = 0.0382855818144094394 ± 2·10⁻²² overlaps the author's ball; every sign is as stated.
  - ε·W_l(ω₀ + ε) tends to the stated residue with error exactly O(ε), and W₀ changes sign on the
    imaginary axis across P2-K6's ball.
  - The positivity argument for l ≥ 1 and the closed form (E-int) are correct.
- **The cover of S_e is complete for every l, and the counts mean what the author says.**
  - Winding 0 at l = 0 with ω_* inside forces Z = P = 1; winding −1 forces P = 1 and Z = 0.
  - The reviewer's own float W_l reproduces all 20 band windings, the symmetry halving, the pole at
    −iκ₊ for l = 0, 1, 7, 19, and winding −1 on l = 20 and l = 21 boxes that reach part (i)'s regime A.
  - Every seam is re-checked with independent code (215 checks; a 40 000-point sampled cover of S_e for l
    up to 10⁶).
- **Integrity.**
  - The plan and its amendments precede all production code and certificates, and were followed. The
    registered designs were used unchanged.
  - The `+dirty` episode is provenance only: 178/178 first-pass leaves equal the committed ones, and
    every committed row is clean.
  - Every part-(i), part-(ii′) and Lemma E certificate replays identically, including all 1280 L leaves
    and all 90 S leaves. 72 further band pieces, none in the author's sample, overlap. 36 mutations of
    the assembly's inputs are each caught on the right item.
- **Should-fix:** S1 (the level of P2-E's non-resonance statement) and S2 (the ledger and FORMULATION
  texts: the review, the band-box edges, W̃_l's formula, and the multiplicity caveat in v1.2). Nits W1–W10.

## 1. Verdicts

| item | verdict |
|---|---|
| A1. Plan and amendments before code; plan followed; §13 deviations | YES. `c98ca16d` and `f82114bf` precede `7168b2e4` and every certificate; the plan body is byte-identical after the amendments (only "(none)" was replaced). Cover, designs (ζ_A(l) equal to the registered table, one try each; T without fallback; S and L at the registered settings), KATs B1–B7, controls NB1–NB6 (NB6 with A5) and kill criteria are as registered; none fired. §13 is provenance or additive only. A5 corrects a real typo (the registered "just above" could not fail). A6 was followed |
| A2. Clean commits; replays; `+dirty` | YES. Every S6 certificate records a clean commit (`7168b2e4`, `5ee4796f` or `5a78bbcb`). `beta_assembly verify` and `p2.assembly verify`: identical. Recomputed identically: Lemma E, part (ii′) (20 l), T, A, cover, NB3–NB5, all 1280 L leaves, all 45 + 45 S leaves (float margins included). 72 fresh band pieces overlap. The 178 `+dirty` first-pass L rows have leaves identical to the committed ones; the final 1280 rows are all tagged `7168b2e4` clean |
| A3. FORMULATION v1.1 | ACCURATE. The diff touches only the header, a bracketed pointer in §2 P2-6 and the log entry; no definition, point, δ, KAT, target or edge changes. Item 1's correction (winding 1 on the full box, 0 on the band box) is right. Wording nit W1 |
| B4. Resonance at a point of E₊ | CORRECT. Index-0 analytic Fredholm family; pole ⟺ kernel; kernel smooth above threshold (any order allowed at r₊, r_c); thresholds 3/2 at r₊ and ½ + κ₊/κ_c = 0.875 at r_c match the explicit Frobenius regularity; t·log t excluded; no other mechanism (§4.B4). Citation nit W6 |
| B5. Log coefficient, pole order, residue | CONFIRMED in sympy and in Arb (residue error O(ε) at ε = 10⁻³, 10⁻⁵, 10⁻⁷ κ₋) |
| B6. Positivity argument (l ≥ 1) and the closed form | CORRECT, every step (§4.B6) |
| B7. W̃₀ and further W̃_l | CONFIRMED rigorously by reviewer code at l = 0, 1, 2, 3, 5, 10, 19, 20, 50, 100; simple pole confirmed |
| C8. Pole-aware counts | CORRECT. P ≤ 1; clearance ≥ 1.517·√2ρ from every point of E₊ ∪ E_c (all n); halving valid with the pole on the axis; float recount reproduces every winding |
| C9. Part (ii′) | CORRECT. At h = h_l exactly the bound needs no uniformity; W7 holds (arctan bound ≤ 0.0711 < 1/5); Jost decay holds; a′_l identities exact and a′_l the least dyadic ≥ ζ_AΩ_P/h_l |
| C10. Part (i) at depth 23/25 | SOUND. W6′ correct; Re(n + 2ρ) ≥ n − 1.118; every depth-dependent input reads `pt.beta`; the negative Jost margin is used nowhere; all certificates replay identically |
| C11. Seams; the assembly check | COMPLETE and FAITHFUL for the bookkeeping. By design it trusts each ingredient's own replay (an increment changed within ¼ passes) and lists, rather than fails on, `+dirty` tags (W5) |
| C12. β | CORRECT: α attained at ω_*, β ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹], width 6.4·10⁻³¹. Multiplicity not overclaimed in the rows; the v1.2 proposal needs the caveat (S2) |
| D13. Proposed rows and notes | ACCURATE apart from S1, S2 and W2–W4. P2-E as its own row: reasonable |
| D14. Proposed FORMULATION v1.2 | NEEDS S1/S2 (level split, multiplicity caveat, the review) |
| D15. Statements that become inaccurate | P2-K6's and P2-4's "stretch goal P2-6" (notes as proposed); FORMULATION's header; the .tex v1.3 entry's "(d) Not established: … a two-sided enclosure of β" needs a new revision entry, not an edit (§8) |
| E16. Tests | Fast: 102 passed, 21 skipped (67 s). `P2_SLOW=1`: 123 passed, 1 pre-existing warning (732 s). Both assembly verifies identical |

## 2. Findings by severity

### Blocking

None. Checked without finding a flaw that defeats P2-6 or P2-E:
- the notion of resonance at σ = −iκ₊ against [HV17] v2 and [Vas13] v2 (§4.B4);
- the Frobenius algebra, the log coefficient, the pole order and the residue (sympy; Arb; §4.B5);
- the positivity argument and the closed form, step by step (§4.B6);
- W̃_l(−iκ₊), rigorously, with independent code (§4.B7);
- the argument principle with the pole, the clearance, the halving, and independent float recounts (§5.C8);
- part (ii′) at fixed h and part (i) at depth 23/25, replayed in full (§5.C9–C10);
- every seam of the cover, a sampled cover, and 36 mutations of the assembly check (§5.C11);
- β and its multiplicity statements (§5.C12).

### Should-fix before the merge

**S1. P2-E's level: the non-resonance of −iκ₊ is at P2-D1′'s level, also at l = 0.**
(`s6-report.md` §1, §5 and §9 (P2-E row, level column); `s6-proof.md` §2.6 and §7; the v1.2 proposal; the
label "Lemma E, l = 0 and l ≤ 19 (Arb)": CERTIFIED in `beta_assembly.INGREDIENTS`. The frozen
pre-registration O1 says the same; leave it, and correct the exit texts.)
- **What the machine proves.**
  - For l ≤ 19, W̃_l(ω₀) ≠ 0 in Arb. With the routine Frobenius algebra of §2.3 this gives "W_l has a
    simple pole at ω₀ with nonzero residue".
  - The band windings force the same pole by a second machine route.
  - That statement is CERTIFIED for l ≤ 19.
- **What it does not prove.** "ω₀ is not a resonance of the forward problem on Ω₂₃" uses the written
  criterion of §2.2:
  - resonant states are smooth ([Vas13] §2.6, [HV17] Lemma 2.14's proof);
  - they separate into spherical harmonics;
  - the smooth solutions at r₊ are the multiples of Z2;
  - hence a resonance ⟺ W̃_l = 0.

  That is P2-D1′'s framework applied at a point D1′ excludes. Its level is P2-D1′'s, WRITTEN_PROOF. The
  same convention is already in the ledger: P2-K6 is "CERTIFIED for the zeros … Resonance identification at
  P2-D1′'s level".
- The README's CERTIFIED needs the written lemmas to be routine. The identification with [HV17]'s
  resonances is the step FORMULATION v0.2 flagged as needing new foundations; it is not routine.
- **Effect.** None on P2-6 (WRITTEN_PROOF either way). The P2-E row as a whole is WRITTEN_PROOF.
- **Fix.** Use the rows of §8, whose level columns separate the pole statement from the non-resonance
  statement. Relabelling the INGREDIENTS entry is optional: the minimum is unchanged, and relabelling would
  require re-running the assembly.
- **Failure scenario.** A later stage cites P2-E at l = 0 as a CERTIFIED fact about [HV17]'s resonances,
  for instance to apply [HV17] with α up to −Im ω_*, without knowing that it rests on the written D1′
  framework.

**S2. The ledger and FORMULATION texts need the review's results and three precision fixes.**
- **P2-6 row.**
  - Replace "Independent review: pending" by the review summary.
  - The band boxes' Im edges are dyadics, b′ = 21012395/2²⁹ ≥ (23/25)κ₋ and b₀ = 22839559/2³⁰ ≤ κ₋/2. They
    are not [−(23/25)κ₋, −κ₋/2].
  - Give P2-E's level as in S1, and qualify B7 (W4).
- **P2-E row.**
  - Give W̃_l's formula instead of "zsol's Jost-normalised form": a ledger row must be readable without
    the code.
  - The level as in S1.
- **FORMULATION v1.2 proposal.**
  - "Each is a simple zero of W₀" → add "(radial order 1; the resonance multiplicity of ω_*, which lies
    below −κ₊, is not claimed)".
  - Add S1's level split and the review.
  - Say that P2-5 is unchanged, and give the header line.
- The corrected texts are in §8.
- **Failure scenario.** "Simple" in the FORMULATION entry is read as a simple resonance at ω_* (pole of
  order 1, rank 1, as HV17 uses the word at 0). That is outside P2-D1′'s multiplicity window and not
  proved.

### Wording and nits

- **W1. FORMULATION v1.1 (applied).**
  - Item 1 writes the band box as [−a′₀, a′₀] × [−0.92κ₋, −κ₋/2]. The edges are the dyadics b′ and b₀ of
    S2.
  - Item 4's "Lemma E shows …" was a plan at v1.1's date.
  - Both are acceptable as dated text; v1.2 can be exact.
- **W2. P2-6 row.**
  - Write "the resonances … with Im σ ≥ −(23/25)κ₋ are exactly σ = 0 and σ = ω_*" (points, no
    multiplicity) for "has exactly two resonances".
  - Add "(real part ± 6.6·10⁻³⁵)" after ω_*'s ball, as in P2-K6.
- **W3. P2-E row.** "Hence D1's QNMs and HV17's resonances agree" also uses P2-D1′'s set clause. Say so.
- **W4. B7 is weak evidence.**
  - At 20 of the 100 points the secant search returned a point within 10⁻⁶ (relative) of its seed
    (|W(stop)|/|W(seed)| = 1.000001).
  - At 69 it did not converge.
  - At 11 it stopped at non-roots well below the band (Im ≈ −2.9 to −6.5 κ₋, |W| ratio 0.11–0.61).
  - It is correctly marked float evidence and is not used. A row should not list "finds no zero" without
    that qualifier, or should drop it. The reviewer's float counts at l = 20, 21 (§5.C8) are a stronger
    float check of the same region.
- **W5. `beta_assembly`.**
  - C10 lists `+dirty` tags but does not fail on them (mutation B2). That is what A2 asked for; the report's
    "no `+dirty` tag" is a listed fact, not one of the 332 checks.
  - `C6.L_rows_reach` compares two constants and does not read `L_PB.json`; `C6.L_band_rows` does.
  - The S5 pole bound reused by C5 handles n ≥ 3 by 3 min(κ₊, κ_c) + Im c > √2ρ, which ignores Re c. It
    rejects a far disc inflated by 5/4 (B3). The stored pieces pass it with margin ≥ 1.095, and the
    reviewer's all-n check passes with ≥ 1.517.
  - All harmless.
- **W6. The threshold citation (`s6-proof.md` §2.2).**
  - Cite [HV17]'s radial-point statement for "s > ½ − Im σ/κ": s > 1/2 + β_jα (ℓ. 1065–1067), with
    β_j = 2/|μ_*′(r_j)| (ℓ. 850) and κ_j = β_j⁻¹ (ℓ. 857).
  - Do not cite [Vas13]'s half-plane. [Vas13] v2 states it as Im σ > β⁻¹(k − 1 − 2s) (ℓ. 1940), i.e.
    s > (k − 1)/2 − (β/2)Im σ. The sign condition in its own radial estimate is
    s − (k−1)/2 + β Im σ − δ > 0 (ℓ. 1627), i.e. s > (k − 1)/2 − β Im σ, twice the shift. With its Kerr–de
    Sitter β± (ℓ. 5980), 2r²/|Δ′| = 1/κ at a = 0, only the second form matches the regularity of the second
    Frobenius solution. The published version was not checked. P2-D1′'s note quotes the v2 statement.
  - Immaterial here: the argument needs only kernel elements that are smooth on a space of some order above
    threshold, and [HV17] allows arbitrarily high regularity at r₂, r₃ (ℓ. 1549).
- **W7. `s6-report.md` §11.** "β > 7/9 is now certified at P_B" → "established (at WRITTEN_PROOF)".
- **W8. `beta_part2.float_zeta_A`.**
  - It also requires E₀ ≥ 0 in floats, which the registered rule ("float margin ≥ 0.02a* and W7") does
    not state.
  - Harmless: the Arb certification requires E₀ ≥ 0 anyway, and every ζ_A(l) equals the registered probe
    table at the first try.
- **W9. Authoring order cannot be checked.**
  - The 1474 lines of production code were committed 6 minutes after the amendments.
  - Nothing could favour the verdict post hoc. The registered designs were used unchanged, and the
    decisive outcomes (the windings, W̃_l ≠ 0) were registered expectations.
- **W10. Pre-registration O4 (frozen; no edit).** "ã₁ … with value −d₁ at h = 0" should read 1/d₁ =
  1/Δ′(r₊). In general ã₁(ω₀) = h²(2κ₊r₊ + L)/Δ′(r₊) (sympy E9). It is nonzero either way, and the proof note
  does not repeat the value.

## 3. Integrity (scope A)

### A1. The plan against execution

| registered (`c98ca16d`, with A1–A6 at `f82114bf`) | executed | verdict |
|---|---|---|
| §1 target and reduction: (a) the band, (b) Lemma E, (c) ω_* a resonance; set clause only; multiplicity at ω_* not claimed | `s6-proof.md` §1, §6, §8 | ✓ |
| O1 Lemma E: written argument for l ≥ 1; Arb W̃_l for l ≤ 19 and samples; residue cross-check; sympy first step | `epole_PB.json`; `s6-proof.md` §2 (A1's write-up: identifications, both expansions with error terms, convergence of ∫VΦ₂Ψ, why l = 0 is excluded, the normalisation) | ✓ (level: S1) |
| O3 part (ii′): A at fixed h_l; ζ_A rule; float table | `A_fixed_PB.json`: ζ_A(l) = 21/4, 37/20, 3/2, 27/20, 5/4, 6/5 (l = 5, 6), 23/20 (7–9), 11/10 (10–19), each at the first Arb try | ✓ (W8) |
| O4 part (i): T6 (s1 = 12, φ_∞ = −1, taper (3, 6)) on [0.8, 1.3]; A at 6/5; S with the extended grid; L with R̃ in every box, band rows | `T_PB.json` (10 boxes, no fallback, μ/a* 0.334 … 0.239), `A_PB.json` (1.972), `S_PB_{plus,c}.json` (45 leaves each), `L_PB.json` (1280 boxes, none refined) | ✓ |
| O5 part (iii): B_l with b′ = eκ₋ rounded up and b₀ = κ₋/2 rounded down to 2⁻³⁰; S3's discs; expected windings 0, −1 | 20 certificates, 2562 pieces, all certified at the planned discs | ✓ |
| §3 cover and its inequalities; the assembly with sha256 pins | `cover_PB.json` (11 inequalities, 6 tilings); `assembly_beta_PB.json` (332 checks, 38 pins) | ✓ |
| §4 KATs B1–B7, controls NB1–NB6 | `kat_controls_PB.json`, `controls_PB.json`, `epole_PB.json` (B5), `B7_PB.json`, `assembly_beta_PB.json` (NB6a–i) | ✓ (W4) |
| §5 kill criteria K1–K5 | none fired; no contingency used | ✓ |
| §7 budget (cap 48 core-h, ≤ 4 workers at nice 10) | about 6.8 core-h | ✓ |
| A2 pins and commit tags; A3 ledger proposals; A4 P2-K6 ball reused unchanged; A5 NB6 typo; A6 B7 | all done; the Krawczyk file equals its main blob | ✓ |

The pre-registration file at `f82114bf` differs from `c98ca16d` only in replacing "(none)" under
"Amendments" by A1–A6, and it is unchanged afterwards. `git diff --name-status 300c22a7 172eb48c --
problems/P2` shows only added files plus `FORMULATION.md`; no S2–S5 source, test or certificate changed,
and nothing outside `problems/P2` changed.

**§13 deviations.** All are provenance or additions:
- the 178 recomputed L rows (A2 below);
- five extra assembly controls;
- the closed form (E-int) with its float test;
- the v1.1 entry (approved), with v1.2 only proposed;
- no fallback used.

A5 is a real correction, not a convenience. As registered ("b′ just above eκ₋"), NB6d could not fail.

### A2. The certificates and the `+dirty` episode

- **Commit tags.** All 33 S6 certificates record clean commits. `computed_at` is `7168b2e4` for the band
  counts, which were launched at 16:36:22Z from that commit (launch record); `assembled_at` is `5ee4796f`.
- **Pins.** The six S2–S5 inputs equal their blobs on main `300c22a7`, P2-4's 54 pinned inputs are
  unchanged, and the 38 inputs of the S6 assembly equal their recorded sha256. The reviewer recomputed all
  of this with its own script.
- **Replays.** All identical; §6 has the table.
  - `beta_assembly verify` and `p2.assembly verify` print "identical".
  - Recomputed from the author's library functions and compared field by field (timing and commit
    excluded): Lemma E, part (ii′), T, A, the cover, NB3–NB5, all 1280 L grid boxes, and all S leaves on
    both sides. The S leaves match exactly, float design margins included.
- **The `+dirty` episode.**
  - The first L run started at 09:28:12 from `7168b2e4` clean. From 09:32 an untracked draft of
    `beta_assembly.py` in `src/p2` made `git_commit()` (which inspects `src/p2` only) append `+dirty` to 178
    rows.
  - The reviewer compared the author's git-ignored `L.w0.first-pass-dirty.jsonl` with `L_PB.json`, reading
    only: 178/178 leaves are identical.
  - The final `L.w0.jsonl` has 1280 unique rows, all tagged `7168b2e4` clean.
  - The reviewer's own full recomputation of the 1280 boxes at the review head reproduces every committed
    leaf. The modules L uses are unchanged since `7168b2e4`, so the recomputation was at a clean tree, and
    the episode is provenance only.
- **Band-run provenance.** The band run's pool processes imported their code once, at 09:36:22, after the
  draft had been moved out (the launch record is clean). Every module it uses is unchanged at the head.
  The author's 141-piece replay and the reviewer's 72-piece replay (§6) both overlap.

### A3. FORMULATION v1.1

- The diff against main touches only three places:
  - the header status;
  - one bracketed pointer in §2 P2-6;
  - the v1.1 log entry.
- D1–D7, §1.3's points, δ, the KATs of §5, the scope fences of §7, the P2-6 target (width ≤ 10⁻⁸) and the
  edge 0.92κ₋ are unchanged.
- Item 1 is a correct correction. The full box at depth 0.92 holds the zeros 0 and ω_* and the simple pole
  ω₀, so W₀'s winding there is 1, not 2; on the band box it is 0.
- Item 2 changes a method, not a definition.
- Item 5 keeps the target (CERTIFIED) and records the expected level (WRITTEN_PROOF). That is not a change
  of target.
- Nit W1.

## 4. Lemma E (scope B)

### B4. Is the author's notion of resonance right at σ = ω₀ ∈ E₊?

**Setting.** [HV17]'s spectral hypothesis concerns "the forward problem in Ω₂₃" (ℓ. 1313), P2-D1′'s notion
(v0.2). The forward problem on Ω₂₃ = {t₂₃ ≥ 0, r₂ − 2δ ≤ r ≤ r₃ + 2δ} carries no condition at the
final boundaries beyond the horizons.

**The threshold, and the smoothness of states.**
- [HV17]'s radial-point statement (ℓ. 1065–1067) propagates regularity into the radial set above
  1/2 + β_jα, where β_j = 2/|μ′(r_j)| (ℓ. 850) = 1/κ_j (ℓ. 857) and α is the weight; a mode e^{−iσt} has
  α = −Im σ.
- So the threshold is s > ½ − Im σ/κ_j. At σ = ω₀ this is 3/2 at r₊ and ½ + κ₊/κ_c = 0.8749 at r_c.
- These match the explicit regularity of the second Frobenius solutions:
  - at r₊, c·t·log t ∈ H^{3/2−0} ∖ H^{3/2} (its Fourier transform decays like |ξ|⁻², so it lies in H^s
    iff s < 3/2);
  - at r_c, (r_c − r)^{0.3749}·g̃ ∈ H^{0.875−0}.
- Above threshold the kernel is finite-dimensional and smooth: [Vas13] ℓ. 1888, "its elements are in
  C^∞(X)"; [HV17] ℓ. 1264, "u is smooth at r = r₂ and r = r₃ (where the a priori regularity exceeds the
  threshold value)".
- [HV17] may take arbitrarily high regularity at r₂ and r₃ (ℓ. 1549). So any s > 3/2 can be used near ω₀,
  and the exact value of the threshold is not even needed: t·log t has an unbounded derivative, so it is
  not smooth.
- On [Vas13]'s normalisation of its half-plane, see W6.

**No other mechanism.**
- **Fredholm theory.** On a fixed space of order s > 3/2 the family σ ↦ P̂₂₃(σ) is analytic (polynomial
  in σ) and Fredholm on a half-plane containing ω₀. It is invertible for Im σ ≫ 1, so its index is 0 and
  the inverse is meromorphic ([Vas13] ℓ. 1944; the poles do not depend on s, ℓ. 1968). At ω₀ the inverse is
  analytic iff P̂₂₃(ω₀) is invertible, iff (index 0) its kernel is trivial. So there is no pole without a
  kernel element.
- **Generalized resonant states.** These form Jordan chains, and every chain begins with a resonant state.
- **The E-points are not singular for P̂.** E₊ ∪ E_c are artefacts of normalising R₊ by a₀ = 1; P̂₂₃(σ)
  itself has no singularity there.
- **States beyond the horizons.**
  - A smooth kernel element vanishing on (r₊, r_c) vanishes beyond r₊ and beyond r_c. On a two-sided
    neighbourhood of r₊, the smooth solutions of each sector equation are B·Z2 with one B, matched by
    Z2′(r₊) = 1; so B = 0. At r_c, Z_c(r_c) = 1 and the other solution is not smooth.
  - Beyond the horizons r is a time function, and the Cauchy problem in r is unique.
- **Separation, and the choice of time function.**
  - P̂₂₃(σ) commutes with rotations, and each smooth state has smooth Y_lm components.
  - Near r₊, u = e^{−iσT}R = e^{−iσ(T+x) + iσc₊}·Z with T + x smooth (T̃(r₊) = −1). So u is smooth iff Z
    is, and Z is a multiple of Z2.
  - [HV17]'s t₂₃ differs from t* by a smooth function of r (P2-D1′ §4.1), at any σ.

**Conclusion.** ω₀ is a resonance of degree l iff W̃_l = 0, exactly as `s6-proof.md` §2.2 says. The author's
characterisation is the right notion for the localized forward problem.

### B5. The first step, the log term, the pole order and the residue

`sympy_lemmaE.py` checks every item, all true:
- **E1.** R = e^{−isω(x−c)}Z with x′ = r²/Δ turns (R) into
  ΔZ″ + (Δ′ − 2isωr²)Z′ − (2isωr + L)Z = 0, for s = ±1.
- **E2.** The coefficient of t^{m−1} in L[t^m] is m(Δ₁(m − 1) + B₀) = Δ₁m(m − iω/κ₊) (m = 1–5). The t⁰
  coefficient is D₁a₁ + C₀, with D₁(ω₀) = 0 and C₀(ω₀) = −(2κ₊r₊ + L). D_m(ω₀) = Δ₁m(m − 1), and the
  exponents at ω₀ are {0, 1}.
- **E3.** For Z = 1 + c·Z2·log t + d₁t, the non-log t⁰ coefficient is Δ₁c + C₀ = 2cκ₊r₊² − 2κ₊r₊ − L.
  So c = (2κ₊r₊ + L)/(2κ₊r₊²) ≠ 0, the author's "Δ₁c + C₀ = 0".
- **E4.** D₁ = −2ir₊²(ω − ω₀), and lim (ω − ω₀)a₁ = i(2κ₊r₊ + L)/(2r₊²).

**The pole order** (`s6-proof.md` §2.3).
- In S_e only D₁ vanishes near ω₀: D_m vanishes only at −imκ₊, and 2κ₊ = 1.6458κ₋ > eκ₋.
- Z₊ = U + a₁V, with U and V generated from (1, 0) and (0, 1) by the m ≥ 2 recurrence. Both are analytic
  in ω and converge for |t| < r₊ − r₋, and V(·; ω₀) = Z2.
- So (ω − ω₀)Z₊ is analytic, with value (residue)·Z2 at ω₀. Z_c is analytic, since ω₀ ∉ E_c.
- W_l is a bilinear form in (Z₊, Z_c) at fixed r, so it has at most a simple pole, with residue
  i(2κ₊r₊ + L)/(2r₊²)·W̃_l.

**Numerically**, with the reviewer's Arb code:

| l | rel. error at ε = 10⁻³κ₋ | 10⁻⁵κ₋ | 10⁻⁷κ₋ | ratio per factor 100 |
|---|---|---|---|---|
| 0 | 9.1446·10⁻³ | 9.1619·10⁻⁵ | 9.1621·10⁻⁷ | 99.8, 100.0 |
| 1 | 1.4217·10⁻³ | 1.4229·10⁻⁵ | 1.4229·10⁻⁷ | 99.9, 100.0 |
| 5 | 4.4840·10⁻³ | 4.4921·10⁻⁵ | 4.4922·10⁻⁷ | 99.8, 100.0 |
| 19 | 6.7381·10⁻³ | 6.7551·10⁻⁵ | 6.7553·10⁻⁷ | 99.7, 100.0 |

The 10⁻⁵ column equals the author's B5 values (9.16, 1.42, 4.49, 6.76 ·10⁻⁵). The normalisation of the
residue formula is therefore the one used.

### B6. The positivity argument for l ≥ 1 (`s6-proof.md` §2.4)

Each step was re-derived:
- **(a) V_l > 0 for L ≥ 2.**
  - E8 checks rF′ = (2/r)(M − Q²/r) − (2Λ/3)r².
  - The reviewer's P_B gives r₊ − Q²/M = 0.0745651 and 2 − (2Λ/3)r_c² = 0.6223715.
  - So L + rF′ > L − 1.378 > 0.
- **(b) The identifications.**
  - Φ₂ = r·e^{−κ₊(x−c₊)}Z2 ∝ rR^{(2)}, and its hyperboloidal form is Z2 times the smooth positive factor
    e^{−κ₊(T+x−c₊)}.
  - Φ_c = r·e^{κ₊(x−c_c)}Z_c is smooth at r_c, and Φ̃_c = r·e^{κ₊(x−c_c)}Z̃_c is the other branch.
  - At ω₀ the smooth solution at r_c grows like e^{κ₊x} and Φ̃_c decays. The signs and branches are as the
    note has them.
- **(c) The asymptotics.**
  - τ = e^{2κ₊(x−c₊)} = t(1 + O(t)) and σ_c = e^{−2κ_c(x−c_c)} = s(1 + O(s)), both analytic: A₊ = 1/(2κ₊),
    A_c = −1/(2κ_c), and the other terms of r* are analytic at the horizon.
  - So G(0) = r₊ and H(0) = H̃(0) = r_c (lim s/σ_c = 1), and the relative errors are O(τ) and O(σ_c).
- **(d) Φ₂ > 0.** This is the usual convexity argument: Φ₂,xx = (V_l + κ₊²)Φ₂ > 0 while Φ₂ > 0.
- **(e) The Wronskian identity.**
  - W[Φ₂, Ψ]′ = −V_lΦ₂Ψ, and W[Φ₂, Ψ](−∞) = 0, since it is O(e^{4κ₊x}).
  - The integrand is O(e^{4κ₊y}) at −∞ and O(e^{−2(κ_c−κ₊)y}) at +∞. It is integrable because
    κ_c − κ₊ = 0.058375 > 0.
  - E7 checks the limits W[Φ_c, Ψ] = 2κ_cσ_cH′e^{κ₊(2x−c_c)} → 0 and W[Φ̃_c, Ψ] → 2κ₊r_ce^{κ₊c_c}.
  - Hence B·2κ₊r_ce^{κ₊c_c} = −∫V_lΦ₂Ψ < 0 for l ≥ 1. So B ≠ 0 and Φ₂ ∦ Φ_c.
- **From B ≠ 0 to W̃_l ≠ 0.**
  - W̃_l is r₊⁻¹r_c⁻¹ times the Wronskian Δ(R^{(2)}R_c′ − R^{(2)′}R_c), up to the nonzero factor
    e^{κ₊(c₊−c_c)}.
  - It vanishes iff the two solutions are dependent, and E5 checks the Wronskian form.
- **(f) The closed form.**
  - E6 checks W[rR₁, rR₂] = Δ(R₁R₂′ − R₁′R₂) and V_lΦ₂Ψ dx = e^{κ₊c₊}(L + rF′)Z2 dr/r, and that the two
    relations give W̃_l = −I/r₊.
  - E7 checks W[Φ̃_c, Φ_c] → 2κ₊r_c².
  - The integrand is O(t) at r₊ and O(s^{−κ₊/κ_c}) at r_c, which is integrable.
  - W̃_l < 0 for l ≥ 1, and every W̃_l the reviewer computed for l ≥ 1 is negative.
- **l = 0.** V₀ = FF′/r changes sign, so the argument does not apply, as the note says. W̃₀ > 0 is
  certified, and its sign agrees with (E-int).

### B7. W̃_l(ω₀), independently and rigorously

`indep_lemmaE.py` is reviewer code on python-flint Arb only, at 200 bits:
- P_B is computed from its definitions by interval Newton (r_ext, then the four roots of Δ).
- Z2 comes from the exponent-1 Frobenius series at r₊ with B₀ = 0 exactly (only D₁ is modified), with a
  verified geometric-majorant tail.
- Z_c comes from the r_c series, again with its tail.
- Both are continued by validated Taylor steps to r_m = (r₊ + r_c)/2, each step with q·h ≤ 0.45 and a
  majorant tail.

| l | reviewer's W̃_l(−iκ₊) | author's ball | overlap |
|---|---|---|---|
| 0 | +0.038285581814409439436 ± 2.0·10⁻²² | +0.038285581814409 ± 4.4·10⁻¹⁶ | ✓ |
| 1 | −22.800094349715433981 ± 3.0·10⁻¹⁹ | −22.800094349715 | ✓ |
| 2 | −20790.693939591746910 ± 4.0·10⁻¹⁶ | −20790.693939592 | ✓ |
| 3 | −1.6146662780385774·10⁷ | −1.6146662780386·10⁷ | ✓ |
| 5 | −9.3233687215157748·10¹² | −9.3233687215158·10¹² | ✓ |
| 10 | −2.5336088082824510·10²⁷ | −2.5336088082825·10²⁷ | ✓ |
| 19 | −2.9600527340170672·10⁵³ | −2.9600527340171·10⁵³ | ✓ |
| 20, 50, 100 | −2.35236885869879·10⁵⁶, −3.27848569511787·10¹⁴³, −9.05049025036792·10²⁸⁸ | samples, nonzero | ✓ |

Further checks:
- **The pole.** The simple pole and its residue are confirmed in B5.
- **ω_*.** W₀(−iyκ₋) is real on the imaginary axis, since W(−ω̄) = conj W(ω). With the reviewer's code it is
  −1.4999·10⁻³⁰ just above P2-K6's ball and +4.4916·10⁻³⁰ just below it (y at the centre ∓ 2 radii). So
  there is a certified zero inside the ball's Im interval, independent of `zsol`, `arbw` and `winding`.
- **P_B.** The reviewer's κ₋, κ₊, κ_c, r_P, Ω_P and γ₀ overlap `params.json`.

## 5. P2-6 (scope C)

### C8. The pole-aware counts

**The argument principle.**
- W_l is meromorphic on a neighbourhood of B_l, and its only possible pole is ω₀. E₊ ∪ E_c meets the region
  Im ≥ −b′ only in ω₀, since b′ < 2κ₊ and b′ < κ_c.
- The pole has order ≤ 1 (B5), and there is no zero or pole on ∂B_l, by the certified pieces and the
  clearance.
- So winding = Z − P with P ∈ {0, 1}.
  - **l = 0.** ω_* ∈ int B₀ gives Z ≥ 1; winding 0 gives Z = P ≤ 1; so Z = P = 1.
  - **1 ≤ l ≤ 19.** Winding −1 gives P = Z + 1, so P = 1 and Z = 0.
- Both readings are forced, and the author's accounting is right.

**Clearance.**
- `indep_assembly.py` checks every one of the 2562 disc boxes against every point −inκ₊ and −inκ_c, all n.
  The minimum of dist/(√2ρ) is 1.517, above the planner's 1.5.
- ω₀ is interior: b₀ = 0.4999999878κ₋ < κ₊ = 0.8228895κ₋ < b′ = 0.9200000091κ₋.
- ω_* lies 0.010619κ₋ above B₀'s bottom edge, and |Re ω_*| ≤ 6.6·10⁻³⁵κ₋.

**The halving.**
- B_l is symmetric about the imaginary axis. Reflecting H gives the left half traversed backwards, along
  which arg W changes by −Δ_H. Traversed forwards it is +Δ_H, so the winding is Δ_H/π.
- The pole and ω_* on the axis do not matter: only boundary values enter.
- The author also computed the left halves at l = 0 and 2, and the reviewer's float run computes both halves
  at every l.

**Independent recount** (`indep_float.py`, reviewer's own vectorised W_l in floats).
- The float W_l agrees with the reviewer's 1500-bit Arb evaluator to ≤ 6·10⁻¹⁴ relative at 9 points, up to
  Re ω = 131κ₋ and l = 20.
- The winding is tracked with adaptive bisection until every step turns by less than 0.3 rad.

| check | result |
|---|---|
| full boundary of B_l, l = 0 | 0.000 (827 points; min abs W ≈ 0.046) |
| full boundary of B_l, l = 1, …, 19 | −1.000 each (819–1809 points) |
| right half H, /π, every l | equal to the full winding |
| [−0.03, 0.03] × [−0.86, −0.79]κ₋ (ω₀ only), l = 0, 1, 7, 19 | −1 each |
| [−0.03, 0.03] × [−0.93, −0.89]κ₋ (ω_* only), l = 0 | +1 |
| [−0.03, 0.03] × [−0.93, −0.79]κ₋ (both), l = 0 | 0 |
| ω_* by Newton | −0.909380978994807iκ₋ |
| l = 20, 21 on [−1.3(l+½)Ω_P, 1.3(l+½)Ω_P] × [−b′, −b₀] (into part (i)'s regime A) | −1, −1 (the pole, no zero; min abs W ≈ 2.2) |

The certified counts are the author's. The reviewer's 72-piece band replay (§6) also recomputed with the
author's `cauchy_segment` the 8 pieces nearest ω₀ and ω_* (l = 0, 1), and all overlap.

### C9. Part (ii′): regime A at fixed h_l

- **Validity at h = h_l.**
  - At h = h₀ exactly, hE = E₀/h₀ + E₁ + h₀E₂ ≥ F. No uniformity in h is used; the certificate also has
    E₀ ≥ 0.
  - The splitting V_l = U/h² + W̃ is exact at h = h_l, including l = 0 (h = 2, L = 0).
  - S2's contour and W1, W2 and W4 do not depend on l.
- **Monotonicity and decay.**
  - W7 at h_l: arctan(eκ₋h_l/(ζ_AΩ_P)) ≤ 0.0711 (at l = 0) < 1/5 = −φ₀, recomputed independently for every
    l.
  - Jost decay: ζ_AΩ_P/h_l ≥ 12.9κ₋ > eκ₋.
  - Each stored ν_low reaches −eκ₋.
- **The a′_l identities.**
  - a′_l is the same rational in `A_fixed_PB.json` and in every `band_PB_l{l}.json`.
  - It is the least 2⁻³⁰-dyadic ≥ ζ_A(l)Ω_P/h_l, recomputed with the reviewer's Ω_P.
- **Replay.** All 20 rows replay identically.

### C10. Part (i) at depth 23/25 for l ≥ 20

- **The depth enters only through `pt.beta`.** `deep_point()` sets pt.beta = (23/25)κ₋. Every
  depth-dependent S4 input reads it:
  - `throat.Params` (S's ν range) and its float design;
  - `part1_T.box_arbs` (T and A);
  - `monotonicity_A` (W7).

  L uses S6's own band rows [−19/20, −11/20]κ₋. No hard-coded 11/20 reaches an S6 regime.
- **W6′.**
  - In S_e, Re(n + 2ρ₊) = n + Im ω/κ₊ ≥ n − 1.118 > 0 for n ≥ 2, and Re(n + 2ρ_c) ≥ n − 0.419 > 0 for
    n ≥ 1.
  - In the τ-series the n = 1 denominator is A₂(1 + 2ρ). R̃ = (1 + 2ρ)R cancels it exactly, giving
    ã₁ = −(A₃ρ(ρ − 1) + B₂ρ + C₁)/A₂, and the n ≥ 2 recurrence is homogeneous.
  - At ω₀ (ρ = −½), R̃ is the exponent-½ (Z2-type, smooth) solution, provided ã₁ ≠ 0. sympy (E9) gives
    ã₁ = (1 + 2κ₊r₊h² − h²/4)/Δ′(r₊) = h²(2κ₊r₊ + L)/Δ′(r₊) > 0 for every h, the scaled log coefficient
    (the pre-registration's "−d₁ at h = 0" should read 1/d₁; W10).
  - The certificate's finite y₊ in the ω₀ box (Re y₊ ≥ 0.1617, Re y_c ≤ −0.1201) shows R̃ ≢ 0 there.
  - Off ω₀, R̃ is a nonzero multiple of R₊, so W5 applies unchanged.
- **The negative Jost margin κ₊ − eκ₋ < 0 is used nowhere.**
  - W2's Frobenius representation near r₊ needs only ω ∉ E₊, and S, T and A have Re ω ≥ 6κ₋.
  - W2's decay needs only Re ω + Im ω > 0.
  - The S closing box has Re ω ≥ Ω_S = 6κ₋.
  - L uses Frobenius series. S4's "n + 2ρ ≠ 0 since Im ω > −κ_h" (W6) is the only use in S4, and W6′
    replaces it.
- **The cover inequalities at depth e**, recomputed: 6κ₋ > eκ₋; 0.8Ω_P/h₀ > eκ₋ and 1.2Ω_P/h₀ > eκ₋;
  0.9Ω_P/h₀ − 8κ₋ = 82.75κ₋; T tiles [0.8, 1.3]; A at 1.2 < 1.3; S tiles [0, 0.9] on both sides; L's rows
  cover [−0.95, −0.55] ⊃ [−e, −0.55].
- **The T redesign and the S grid** are design choices. Their certificates are produced by S4's checkers
  (reviewed at S4) and replay identically. The controls show the redesign is needed: NB4, S4's T at depth e,
  gives −0.0475a*.

### C11. The seams, and the assembly check

**The seams** (`indep_assembly.py`, independent code):

| seam | check | value |
|---|---|---|
| band against S_δ | −b₀ = −0.4999999878κ₋ > −0.55κ₋ | overlap [−0.55, −0.5] (κ₋/2 − b₀ = 5.18·10⁻¹⁰) |
| b′ against eκ₋ | b′ the least 2⁻³⁰-dyadic ≥ eκ₋ | b′ − eκ₋ = 3.8790·10⁻¹⁰, with eκ₋ enclosed to about 10⁻⁶⁰ |
| l = 19 against l = 20 | l₀ = ⌈1/h₀ − ½⌉ = 20; 2/41 = h₂₀ ≤ h₀ < h₁₉ | no l missed |
| Re ω = a′_l | the band box is closed at a′_l; A covers Re ω ≥ ζ_AΩ_P/h_l ≤ a′_l | overlap |
| Im ω = 0, upper half-plane | inside S_δ (P2-4, via P2-0) | — |
| E_c and the higher E₊ points | κ_c − eκ₋ = 1.2751κ₋; 2κ₊ − eκ₋ = 0.7258κ₋ | outside S_e |
| sampled cover | 40 000 random (l, ω) in S_e, l ≤ 10⁶, both signs of Re ω, seams over-sampled, each assigned to S_δ, a band box, A at fixed h, L, S, T or A | none uncovered |

**Is `beta_assembly` faithful?** Yes, for the bookkeeping:
- the pins;
- the E-points;
- b′, b₀ and a′_l;
- the traced paths, the discs and their clearance;
- the exact re-sums and the stored windings with the pole correction;
- part (i)'s regimes, depths, tilings and cover replay;
- l₀;
- ω_* in B₀;
- β's width;
- the KAT and control records.

By design it does not re-verify an ingredient's analytic content. The ingredients' replays do that, and the
review re-ran them (§6).
- **Mutations** (`mutations.py`). M1–M36 are each caught on the named items:
  - a removed, swapped or truncated piece;
  - an increment off by 2π;
  - a certificate in the wrong slot;
  - the band top moved below S_δ;
  - a disc inflated onto −iκ₊;
  - an uncertified band;
  - a winding of −1 at l = 0;
  - a disc shrunk below its segment;
  - part (ii′) certified at h₀ = 2/41 instead of h_l, without W7, with a′₆ below ζ_AΩ_P/h₆, or with a row
    missing;
  - Lemma E's W̃₀ not separated from 0, a failed B5 row, or negative positivity inputs;
  - a T box removed or certified only to −κ₋/2;
  - A without W7;
  - L's rows ending at −0.9κ₋, or a negative L margin;
  - an S job missing, S_c at depth 11/20, or S_c missing;
  - a false cover inequality;
  - ω_* moved below the band, or inside it (caught by the KAT-6 consistency check);
  - Z4 false;
  - failed KAT, replay or NB-control records;
  - P2-4 not passing, or one of its inputs changed.
- **What it does not see, by design** (W5):
  - B1: an increment changed within the ¼ tolerance passes;
  - B2: a `+dirty` tag is listed, not failed;
  - B3: a far disc inflated by 5/4 is rejected only by the conservative n ≥ 3 bound.

### C12. β

- **D3's α.** Every QNM of every degree in S_e is 0 or ω_*, so α = −Im ω_*, and the infimum is attained.
  The QNMs are D1's: zeros of W_l off E₊ ∪ E_c, so ω₀ is excluded by definition.
  - For β alone, P2-E's non-resonance half is therefore not needed. Its pole half (order ≤ 1) is: it turns
    the windings into zero counts. `s6-proof.md` §1 item 6 says this correctly.
- **ω_* is a resonance.** By P2-D1′'s set clause, which needs no multiplicity: ω_* ∉ E₊ ∪ E_c, since
  0.9094 ∉ {0.8229, 1.6458, 2.1951, …}.
- **ω_* has degree l = 0 only.** For l ≤ 19, ω_* lies in B_l and Z = 0; for l ≥ 20 it lies in L's band rows.
- **The enclosure.** P2-K6's S3 ball, reused unchanged and equal to its main blob, gives
  β ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹]. The width is 6.3·10⁻³¹ ≤ 10⁻⁸. The reviewer's
  certified sign change confirms the zero inside its Im interval.
- **Multiplicity.**
  - "Simple" in the rows means radial order 1, and the rows say that the resonance multiplicity at ω_* is
    not claimed.
  - The v1.2 proposal's "each is a simple zero of W₀" needs the same caveat (S2).
  - Nothing else overclaims. P2-5 is unchanged, and the stronger conversion (α up to −Im ω_*) is noted and
    not claimed.

## 6. Independent checks, replays and tests (scope E)

All under `problems/P2/certificates/review_s6/` (README there), at the review head, `nice -n 10`.

| check | what | result |
|---|---|---|
| `indep_lemmaE.py` → `indep_lemmaE_out.json` | reviewer Arb code: W̃_l at 10 values of l, the residue at 3 values of ε, the ω_* sign change, P_B against `params.json` | pass (0.4 s) |
| `sympy_lemmaE.py` → `sympy_lemmaE_out.json` | E1–E9 (§4.B5, B6; §5.C10) | all true |
| `indep_float.py` → `indep_float_out.json` | float recount of all 20 band boxes (full and half), 6 small boxes, ω_*, l = 20, 21 seam boxes; accuracy against 1500-bit Arb | pass (51 s) |
| `indep_assembly.py` → `indep_assembly_out.json` | seams, pieces, clearance (all n), re-sums, W7, cover inequalities, 40 000-point sampled cover | 215/215 |
| `mutations.py` → `mutations_out.json` | M1–M36 against `p2.beta_assembly.check`; B1–B3 by design | all caught |
| `quotes_check.py` → `quotes_check_out.json` | 16 [HV17]/[Vas13] passages verbatim, with lines | 16/16; sha256 match |
| `dirty_episode_check.py` → `dirty_episode_check_out.json` | the 178 `+dirty` rows against `L_PB.json` (read-only) | 178/178 identical |
| `replay_author.py epole, part2, TA, controls, L, S+, Sc` | the author's certificates recomputed with the author's code, outputs only in `review_s6/` | all identical (L 1280/1280, S 45/45 and 45/45 exactly) |
| `replay_author.py band 20261027 64` | 72 band pieces (8 nearest ω₀/ω_*, 64 random outside the author's 5% sample), stored K | 72/72 overlap |
| `beta_assembly verify`; `assembly verify` | the author's replays | identical; identical |
| `pytest problems/P2/tests` | fast | 102 passed, 21 skipped (67 s) |
| `P2_SLOW=1 pytest problems/P2/tests` | slow | 123 passed, 1 pre-existing warning (S2's `m3` overflow), 732 s |

No run changed a tracked file: `git status --short` showed only the two new paths before each commit.

## 7. Compute and tool notes

| item | CPU |
|---|---|
| reviewer code (Lemma E Arb, sympy, float recounts and accuracy checks, seams, mutations, quotes, dirty check) | ≈ 0.03 core-h (≈ 100 s) |
| author-code replays: Lemma E, part (ii′), T/A/cover, controls, L (1280 boxes), S both sides | ≈ 0.16 core-h (≈ 565 s) |
| fresh band sample (72 pieces, 2 workers) | ≈ 0.10 core-h (377 s) |
| assembly verifies; fast pytest | ≈ 0.02 core-h (≈ 70 s) |
| `P2_SLOW=1` pytest | ≈ 0.20 core-h (727 s) |
| **total** | **≈ 0.51 core-h** (cap 4) |

All runs were at nice 10, with at most three of the reviewer's processes at once (for example the L replay,
the S replays and the quick replays in parallel, or the 2-worker band sample with the fast suite). No other
process was touched. No permission was denied. The harness blocked one foreground `sleep` (its rule against
chained sleeps), and the review waited with a condition loop instead.

## 8. The ledger

**Level of P2-6: WRITTEN_PROOF.**
- The minimum of the ingredients. P2-4, P2-D1′'s set clause, P2-E's non-resonance statement, part (i) at
  depth 23/25 (W1, W2, W5, W6′ written) and part (ii′) (W1, W2, W4 written) are WRITTEN_PROOF.
- The band counts, P2-E's pole statement for l ≤ 19 and P2-K6's ball are CERTIFIED.
- FORMULATION §2's target, CERTIFIED, is not reached.

**Level of P2-E: WRITTEN_PROOF**, by part:
- the simple pole (W̃_l ≠ 0) is CERTIFIED for l ≤ 19 and WRITTEN_PROOF for l ≥ 20;
- the non-resonance of −iκ₊ is at P2-D1′'s level for every l (S1).

**Should Lemma E be its own row?**
- Yes, as the author proposes. It is a statement about one point, independent of the band.
- It complements P2-D1′, which excludes E₊ ∪ E_c, and its level profile differs from P2-6's.
- A later use of [HV17] with α up to −Im ω_* would cite exactly it.
- The id is the owner's choice.

**P2-6 (new row; the author's text with S1, S2 and W2–W4 applied):**

| P2-6 | P2 | At P_B (ΛM² = 3/50, Q = (999/1000)Q_max): the resonances of the forward problem on Ω₂₃ (HV17 §2.5) with Im σ ≥ −(23/25)κ₋ (= −0.547γ₀) are exactly σ = 0 and σ = ω_*, P2-K6's l = 0 zero, ω_*/κ₋ ∈ −0.909380978994805455524933762622i ± 3.2·10⁻³¹ (real part ± 6.6·10⁻³⁵). Both are of degree l = 0 only and simple zeros of W₀ (radial order 1). In D1's terms, the QNMs of any degree in that strip are 0 and ω_*. The one point of E₊ ∪ E_c in the strip, −iκ₊, is not a resonance (P2-E). Hence D3's gap is attained at ω_*: β = −Im ω_*/κ₋ ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹] (two-sided; width 6.4·10⁻³¹ ≤ 10⁻⁸), and D1's QNMs and HV17's resonances agree in the strip. Assembled from P2-4 (S_δ), the conjugation symmetry (Re σ < 0), part (i) at depth 23/25 for l ≥ 20, and for l ≤ 19 regime A at fixed h = 1/(l+½) (Re σ ≥ a′_l) with pole-aware argument-principle counts on the band boxes B_l = [−a′_l, a′_l] × [−b′, −b₀] (b′ = 21012395/2²⁹ ≥ (23/25)κ₋, b₀ = 22839559/2³⁰ ≤ κ₋/2); with P2-E, P2-D1′'s set clause and P2-K6. The resonance multiplicity at ω_* (below −κ₊, outside P2-D1′'s window) is not claimed. Not P_A | WRITTEN_PROOF, the minimum of the ingredients: P2-4, P2-D1′ (set clause), P2-E's non-resonance statement, part (i) at depth 23/25 (W1, W2, W5, W6′ written) and part (ii′) (W1, W2, W4 written) are WRITTEN_PROOF; the band counts, P2-E's pole statement for l ≤ 19 (W̃_l ≠ 0) and P2-K6's ball are CERTIFIED; FORMULATION §2's target CERTIFIED is not reached. Band counts (S3's Cauchy-disc machinery, every disc box clear of E₊ ∪ E_c; 2562 stored pieces; a 5% replay overlaps): winding 0 at l = 0 (ω_* and the simple pole at −iκ₊) and −1 for 1 ≤ l ≤ 19 (the pole, no zero). Part (i) at depth 23/25 with h₀ = 2/41 (l₀ = 20): L on the band rows with the pole-free series R̃ = (1+2ρ)R (W6′), S with an extended design grid, T redesigned, A at ζ_A = 6/5. Machine-checked bookkeeping (exact rationals and Arb, 256-bit; 332 checks; the six S2–S5 inputs pinned by sha256 and equal to their blobs on main 300c22a7; P2-4's 54 pinned inputs unchanged): a′_l identical between the band boxes and the A-at-fixed-h certificates; b′ − (23/25)κ₋ = 3.9e−10 > 0; the stored discs trace the right half of each band box exactly; the re-summed windings; part (i)'s cover replays identically; ω_* lies in B₀, 0.0106κ₋ above its bottom edge. Negative controls NB6a–NB6i fail as required (a band certificate removed; a′₇ lowered by 2⁻²⁴; a wrong winding; b′ just below (23/25)κ₋; and five more). KATs: B1 (P2-K6's ball replayed), B2 +1, B3 −1 (the pole; l = 0, 1), B4 0, B5 (the residue, with the independent S1 code), B6 +1 (two zeros below the edge); controls NB1–NB5 fail as required; B7 (float, not used) found no zero, but most of its secant searches did not converge. Independent review 2026-10-27: proven as stated. Reviewer code sharing nothing with p2 recomputes Lemma E rigorously (W̃_l(−iκ₊) for l = 0, 1, 2, 3, 5, 10, 19, 20, 50, 100, all overlapping the stored balls; ε·W_l(−iκ₊ + ε) tends to the residue as O(ε); a certified sign change of W₀ on the imaginary axis across P2-K6's ball), reproduces every band winding, the symmetry halving, the pole and the l = 20, 21 seam in floats, and re-checks every seam (215 checks; a 40 000-point sampled cover of S_e for l up to 10⁶); 36 further mutations of the assembly's inputs are each caught; every part-(i), part-(ii′) and Lemma E certificate replays identically (all 1280 L leaves, all 90 S leaves), and 72 further band pieces were recomputed | problems/P2/src/p2/beta_common.py, beta_epole.py, beta_part1.py, beta_part2.py, beta_part3.py, beta_assembly.py; problems/P2/certificates/beta/ (assembly_beta_PB.json); problems/P2/notes/s6-preregistration.md, s6-proof.md, s6-report.md; problems/P2/tests/test_beta.py; docs/review-2026-10-27-p2-s6-adversarial.md; problems/P2/certificates/review_s6/ | 2026-10-27 |

**P2-E (new row; S1, S2 and W3 applied):**

| P2-E | P2 | At P_B, the point ω₀ = −iκ₊ (≈ −0.8229iκ₋; the one point of E₊ ∪ E_c with Im ω ≥ −(23/25)κ₋) is not a resonance of the forward problem on Ω₂₃ (HV17 §2.5) in any sector l ≥ 0, and every W_l has a simple pole there: Res_{ω₀}W_l = i(2κ₊r₊ + L)/(2r₊²)·W̃_l ≠ 0, where W̃_l := [Δ(Z2·Z_c′ − Z2′·Z_c) + 2iω₀r²·Z2·Z_c]/(r₊r_c) at any r ∈ (r₊, r_c). Here Z2 = t + Σ_{m≥2} b_m t^m (t = r − r₊) is the exponent-1 Frobenius solution at r₊ of ΔZ″ + (Δ′ − 2iωr²)Z′ − (2iωr + L)Z = 0, Z_c is the solution of ΔZ″ + (Δ′ + 2iωr²)Z′ + (2iωr − L)Z = 0 analytic at r_c with Z_c(r_c) = 1, and W_l has the same form with Z₊ = Σ a_m t^m, a₀ = 1. The exponent-0 solution at r₊ carries t·log t (the first step reads 0·a₁ = 2κ₊r₊ + L), so the smooth solutions at r₊ are the multiples of Z2, and ω₀ is a resonance in sector l iff W̃_l = 0. For l ≥ 1, W̃_l = −(1/r₊)∫_{r₊}^{r_c}(L + rF′)Z2 dr/r < 0 (V_l > 0 and Z2 > 0); W̃₀ = 0.038285581814409 ± 4.4·10⁻¹⁶. With P2-D1′'s set clause, D1's QNMs and HV17's resonances agree for Im ω ≥ −(23/25)κ₋. Not P_A | WRITTEN_PROOF. The pole (W̃_l ≠ 0): CERTIFIED for l ≤ 19 (Arb, 128-bit; routine Frobenius steps, the first checked in sympy; P2-6's band windings force it again), WRITTEN_PROOF for l ≥ 20 (positivity, notes/s6-proof.md §2.4; part (i)'s L certificate at depth 23/25 also separates the two solutions at ω₀, by W5 and W6′; Arb samples at l = 20, 25, 30, 40, 50, 60, 100 are nonzero). The non-resonance of ω₀, every l: at P2-D1′'s level (WRITTEN_PROOF), through the written criterion of notes/s6-proof.md §2.2 (resonant states are smooth above the radial-point threshold, 3/2 at r₊ and ½ + κ₊/κ_c at r_c for σ = ω₀). KAT B5: ε·W_l(ω₀ + ε) tends to the residue as O(ε), with zsol and with the independent S1 code arbw (l = 0, 1, 5, 19). Float check, not used: the closed form reproduces the Arb W̃_l (l = 0, 1, 2, 5) to 10⁻⁸. No literature is cited. Independent review 2026-10-27: proven as stated, at the levels given here. The notion of resonance at ω₀ was checked against HV17 v2 and Vas13 v2; the reviewer's sympy re-derives the phase-factored equations, the first step, the log coefficient (2κ₊r₊ + L)/Δ′(r₊), the residue and the identities behind the closed form; reviewer Arb code sharing nothing with p2 recomputes W̃_l for l = 0, 1, 2, 3, 5, 10, 19, 20, 50, 100 (all overlap the stored balls; W̃₀ = 0.0382855818144094394 ± 2·10⁻²²) and the residue as O(ε) at ε = 10⁻³, 10⁻⁵, 10⁻⁷ κ₋ | problems/P2/src/p2/beta_epole.py; problems/P2/certificates/beta/epole_PB.json; problems/P2/notes/s6-proof.md §2; problems/P2/tests/test_beta.py; docs/review-2026-10-27-p2-s6-adversarial.md; problems/P2/certificates/review_s6/ | 2026-10-27 |

**Notes for existing rows.** These are the author's (`s6-report.md` §9), checked and accurate. Each is
appended to the end of the row's level column.
- **P2-K6** (needed): ` Note 2026-10-27 (S6): P2-6 reuses this row's S3 ball (certificates/part3/krawczyk_PB_l0.json) unchanged; with it, β ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹] two-sided, which supersedes the one-sided bound for β (the zero itself is unchanged).`
- **P2-4** (needed): ` Note 2026-10-27 (S6): the two-sided enclosure is P2-6, β ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹]; the bound 0.55 ≤ β ≤ 0.90938098 stated here stays true.`
- **P2-D1′** (optional; S1's level made explicit): ` Note 2026-10-27 (S6): at P_B the one point of E₊ ∪ E_c with Im ω ≥ −(23/25)κ₋, −iκ₊, is not a resonance in any sector, and every W_l has a simple pole there (P2-E; the non-resonance at this row's level); so in that strip the resonances are exactly the zeros of the W_l (used by P2-6). No multiplicity is claimed below −κ₊.`

No other row changes. P2-5 stays as it is; S6 does not claim the stronger conversion.

**FORMULATION v1.1 (applied).** No correction is needed (W1 is optional).

**FORMULATION v1.2 (the author's proposal with S1, S2 applied; for the owner, after the user's approval):**

```
- **v1.2 (⟨date⟩, S6 exit; `notes/s6-proof.md`, `notes/s6-report.md`; independent review
  `docs/review-2026-10-27-p2-s6-adversarial.md`).** P2-6 at P_B is assembled.
  - **P2-6 at P_B: WRITTEN_PROOF**, the minimum of its ingredients. §2's target, CERTIFIED, is not
    reached.
    - The resonances of the forward problem on Ω₂₃ with Im σ ≥ −0.92κ₋ are exactly 0 and ω_*, both of
      degree l = 0 and simple zeros of W₀ (radial order 1; the resonance multiplicity of ω_*, which lies
      below −κ₊, is not claimed).
    - β = −Im ω_*/κ₋ ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹]. The frozen target (width
      ≤ 10⁻⁸) is met by P2-K6's ball, reused unchanged.
    - The bookkeeping is machine-checked, with negative controls:
      `certificates/beta/assembly_beta_PB.json` (`python -m p2.beta_assembly`).
  - **Lemma E** (row P2-E): −iκ₊ is not a resonance for any l, and every W_l has a simple pole there.
    The pole is CERTIFIED for l ≤ 19 and WRITTEN_PROOF for l ≥ 20; the non-resonance is at P2-D1′'s
    level. So D1's QNMs and HV17's resonances agree in S_e, and D3's α is HV17's α₀ there.
  - **The independent review** (2026-10-27) found P2-6 and P2-E proven as stated, with nothing
    blocking; its corrections (the level of P2-E's non-resonance statement, and the ledger and
    FORMULATION texts) are applied.
  - **Unchanged:** the edge, the target, D1, D1′, D3, D7, the KATs of §5 (S6's B1–B7 are additive,
    v1.1), P2-5 (no stronger conversion is claimed) and the scope fences of §7.
  - **Open:** P2-7 (P_A).
```

The header status line, when the rows are entered: "S6 assembled P2-6 at P_B at WRITTEN_PROOF
(`notes/s6-proof.md`; Revision log v1.2), and P2-E records the exceptional point −iκ₊; P2-7 (P_A) is still a
target." The header's list of ledger rows then gains P2-6 and P2-E.

**The problems document** (`docs/problems/open-problems-theoretical-cosmology-2026.tex`; the user approves
.tex changes, and the author drafts).
- Once P2-6 is entered, the v1.3 entry's "(d) Not established: … a two-sided enclosure of β (the stretch
  goal)" is out of date. v1.3 is a dated revision-log entry and should not be edited.
- A new entry (v1.4, P2 only) should state:
  - (a) at the primary point, for every ℓ, the only QNMs with Im ω ≥ −0.92κ₋ (≈ 0.55γ₀) are 0 and the
    ℓ = 0 near-extremal mode ω_*, so β = −Im ω_*/κ₋ ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹];
    this is the stretch goal of v1.2;
  - (b) the point −iκ₊, where the radial Wronskian has a pole, is not a resonance, so the two notions of
    QNM agree in that strip;
  - (c) the level: a written proof with machine-checked certificates, independently reviewed
    2026-10-27, not a fully machine-checked certificate;
  - (d) not established: the resonance multiplicity of ω_*; anything below −0.92κ₋; the second point;
    the Einstein–Maxwell problem; a conversion stronger than v1.3 (c);
  - that it supersedes v1.3 (d) for the two-sided enclosure.
- Optional. The entry point's parenthesis "(an interval-Newton enclosure of the ℓ=0 near-extremal mode)"
  describes the stretch goal more narrowly than it turned out to be: the band exclusion and the
  exceptional point were needed too. v1.4 could say so; the pre-registered entry point itself should not
  change.

### Before the S6 branch is merged into main

1. **S1.** Correct P2-E's levels in `s6-report.md` §1, §5 and §9, and in `s6-proof.md` §2.6 and §7: the
   pole is CERTIFIED for l ≤ 19; the non-resonance is at P2-D1′'s level. Relabelling the INGREDIENTS entry
   is optional; it requires re-running the assembly.
2. **S2.** Use the P2-6 and P2-E rows above, or apply their changes to `s6-report.md` §9; replace §10 by
   the v1.2 text above.
3. **W1–W10** as convenient. W4, W6 and W7 are one-line edits in the notes.
4. The owner enters P2-6 and P2-E and the notes for P2-K6 and P2-4 (P2-D1′ optional), and applies v1.2 and
   the header.
5. The .tex v1.4 entry, drafted by the author and approved by the user.
6. The user approves the merge into main, and the coordinator runs it.
