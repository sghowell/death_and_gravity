# P2 S6 pre-registration: the two-sided β at P_B (P2-6)

Status: written 2026-10-27 on branch `p2/s6-beta` (from main at `300c22a7`), after the Phase-1
probes and before any production run. The user approved S6 on 2026-10-27. **Checkpoint:** S6 stops
here for the coordinator's go.

The Phase-1 probes are committed at `ec568eb8`:
- code: `src/p2/s6_probe.py`;
- outputs: `results/s6/*.json`.

They are probes, not certificates, and no ledger text cites them. Any later change to anything fixed
here goes in a dated **Amendments** section at the end, never by silent edit.

**Notation.**
- e := 23/25 = 0.92, the strip edge of FORMULATION P2-6.
- S_e := {Im ω ≥ −eκ₋}.
- The band: B := {−eκ₋ ≤ Im ω ≤ −κ₋/2}. It overlaps S_δ = {Im ω ≥ −0.55κ₋}, which P2-4 covers.
- ω₀ := −iκ₊ = −0.822889544i κ₋, the first point of E₊.
- ω_* := the l = 0 zero of P2-K6, ω_*/κ₋ ∈ −0.909380978994805455524933762622i ± 3.2·10⁻³¹
  (`certificates/part3/krawczyk_PB_l0.json`).
- h = 1/(l+½), h₀ = 2/41.
- e/(γ₀/κ₋) = 0.5466, so the band reaches 0.547γ₀.

## 1. Target and its reduction

**P2-6 at P_B (target).**
- Every resonance σ ∉ {0, ω_*} of the forward problem on Ω₂₃, of any degree, has Im σ < −eκ₋.
- ω_* is a resonance.

Hence α = −Im ω_* (the infimum is attained) and

  **β = −Im ω_*/κ₋ ∈ [0.909380978994805455524933762622 ± 3.2·10⁻³¹]**,

which is two-sided, of width 6.4·10⁻³¹ ≤ 10⁻⁸. No new enclosure is computed: the width is that of
P2-K6's Krawczyk ball.

**Reduction.**
- **S_δ is P2-4** (WRITTEN_PROOF): no zero of any W_l in S_δ∖{0}.
- **Only one E-point in the strip.** S_e ∩ (E₊ ∪ E_c) = {ω₀}, because e < 2κ₊/κ₋ = 1.6458 and
  e < κ_c/κ₋ = 2.1951.
- **The set clause of P2-D1′** holds for every σ ∉ E₊ ∪ E_c, at any depth: σ is a resonance iff
  W_l(σ) = 0 for some l.
- So P2-6 follows from three facts:
  - **(a) the band:** for every l, W_l has no zero in B, except ω_* when l = 0;
  - **(b) Lemma E:** ω₀ is not a resonance, for any l;
  - **(c) ω_* is a resonance:** P2-K6 gives a zero of W₀, ω_* ∉ E₊ ∪ E_c, so D1′'s set clause
    applies.
- **Multiplicity at ω_*.** The radial order m₀(ω_*) = 1 is certified (Krawczyk, and the band count
  below). The resonance multiplicity at ω_* lies outside D1′'s window and is **not claimed**.

## 2. Foundations

### O1. The E₊ point ω₀ = −iκ₊ (Lemma E)

**Lemma E.** For every l ≥ 0, ω₀ is not a resonance of the forward problem on Ω₂₃, and W_l has a
simple pole at ω₀ with nonzero residue.

*Argument.* This uses zsol's Z₊-equation, ΔZ″ + (Δ′ − 2iωr²)Z′ − (2iωr + L)Z = 0, with
Z₊ = Σ a_m t^m, t = r − r₊ and a₀ = 1.

1. **Exponents and the log.** At r₊ the exponents are 0 and iω/κ₊, which are {0, 1} at ω₀.
   - The first step is D₁a₁ + C₀ = 0, with D₁ = Δ′(r₊)(1 − iω/κ₊) and C₀ = −(2iωr₊ + L).
   - At ω₀, D₁ = 0 but C₀ = −(2κ₊r₊ + L) ≠ 0 (2κ₊r₊ > 0, L ≥ 0).
   - So the exponent-0 solution carries t·log t. It lies in H^{3/2−0} only, below the threshold
     s > ½ − Im σ/κ₊ = 3/2 at r₊.
   - The solutions smooth at r₊ are therefore the multiples of Z2 = t + Σ_{m≥2} a_m t^m (exponent
     1; D_m = Δ′(r₊)m(m − 1) ≠ 0 for m ≥ 2).
   - At r_c the exponents are {0, κ₊/κ_c = 0.375}, so only Z_c is smooth there.
   - Hence ω₀ is a resonance in sector l **iff W̃_l := W(Z2, Z_c) = 0**. This is zsol's
     Jost-normalised Wronskian form at ω₀.
2. **The residue.** Near ω₀ only D₁ vanishes (D_m ≠ 0 for m ≥ 2 on S_e), so Z₊ has at most a simple
   pole:
   - Z₊ = [i(2κ₊r₊ + L)/(2r₊²)]·Z2/(ω − ω₀) + O(1);
   - Res_{ω₀} W_l = i(2κ₊r₊ + L)/(2r₊²)·W̃_l.

   So W_l has a simple pole at ω₀ iff ω₀ is not a resonance. Otherwise ω₀ is a removable point or a
   zero.
3. **l ≥ 1 (written).**
   - At ω₀ the radial equation is real: Φ″ = (V_l + κ₊²)Φ in r*.
   - V_l = F(L + rF′)/r² > 0 on (r₊, r_c) for L ≥ 2. This is because
     rF′ = (2/r)(M − Q²/r) − (2Λ/3)r² > −(2Λ/3)r_c² > −2.
     - Arb: r₊ − Q²/M = 0.0746 > 0, and 2 − (2Λ/3)r_c² = 0.622 > 0.
     - Output: `results/s6/E_PB.json`, `V_l_positive_for_l_ge_1`.
   - Φ₂ := rR^{(2)} (∝ Z2) behaves like c·e^{κ₊r*}, c > 0, at −∞. It stays positive and increasing,
     because Φ″ > 0 while Φ > 0.
   - With Ψ = e^{κ₊r*}:
     - W(Φ₂, Ψ)′ = −VΦ₂Ψ < 0, and W(Φ₂, Ψ)(−∞) = 0;
     - at +∞, W(Φ_c, Ψ) → 0 (κ₊ < κ_c) and W(Φ̃_c, Ψ) → 2κ₊ for Φ̃_c ~ e^{−κ₊r*};
     - so the coefficient of Φ̃_c in Φ₂ is −(2κ₊)⁻¹∫VΦ₂Ψ dr* < 0.

     Hence Φ₂ ∦ Φ_c and W̃_l ≠ 0.
4. **l = 0 (Arb), with l ≤ 19 as a cross-check.** W̃_l(ω₀) is computed with the exponent-1
   Frobenius series at r₊, zsol's Taylor continuation to r_m = r_P, and zsol's Z_c.
   - Probe values (128 bits):
     - W̃₀ = 0.03828558181 ± 9·10⁻¹⁵;
     - W̃₁ = −22.80009435;
     - W̃₂ = −20790.69394;
     - … to W̃₁₉ = −2.96005·10⁵³.
   - All 20 are certified nonzero, with real balls (the imaginary parts contain 0, as they must at
     real ω²).
   - Samples l = 20, 25, …, 60 and l = 100 are also nonzero.
   - The residue formula is cross-checked from the other side: ε·W_l(ω₀ + ε), with zsol at
     ε = 10⁻⁵κ₋ and 10⁻⁷κ₋ for l = 0, 1, 5, 19, converges to it with an error of O(ε).
5. **First step in sympy.** a₁ = −C₀/D₁, D₁(ω₀) = 0 and C₀(ω₀) = −(2κ₊r₊ + L) (`E_PB.json`,
   `sympy_first_step`).

**Two independent confirmations in Phase 2.**
- For l ≤ 19, the band counts (O5) force a simple pole at ω₀, hence W̃_l ≠ 0.
- For l ≥ 20, regime L's certificate at the box containing ω₀ (O4) certifies directly that the
  solution regular at r₊ (Z2-type at ω₀) is not proportional to the one regular at r_c.

**Verdict: ω₀ is not a resonance for any l.** The stop condition (β ≤ κ₊/κ₋) does not occur.
- Consequence: D1's QNMs (zeros of W_l off E₊ ∪ E_c) and HV17's resonances agree in S_e, so D3's
  α is HV17's α₀ there.
- No definitional change is needed.
- **Level.** l = 0: CERTIFIED (Arb, with routine Frobenius steps). l ≥ 1: WRITTEN_PROOF
  (positivity), and CERTIFIED for l ≤ 19.
- **Literature.** No source is cited for this point. The argument above is self-contained, and no
  paper was read and verified for it.

### O2. P2-D1′'s multiplicity window

S6 uses only D1′'s **set clause** (σ ∉ E₊ ∪ E_c, any depth) together with Lemma E at ω₀.
- β is an infimum over a set, so the multiplicity clause (Im σ > −κ₊) is not used.
- S6 claims no resonance multiplicity at ω_*, nor at any point below −κ₊.
- The windings count zeros of W_l with radial order. They are used only to show that there is no
  zero other than ω_* and that ω_* is a simple zero of W₀.

### O3. Part (ii) below −κ₊

**Why Liouville–Green fails.** On the real line, the Jost solution e^{−iωr*} is dominant at
r* → −∞ when Im ω < 0. Volterra iteration converges only while the potential decays faster, which is
Im ω > −κ₊.

**The replacement.** S4's regime A on S2's rotated contour, **at fixed h_l** = 1/(l+½) for each
l ≤ 19.
- A's bound E₀/h₀ + E₁ + h₀min(0, E₂) with h₀ := h_l is uniform for h ≤ h_l.
- It applies on Re ω ≥ a′_l := ζ_A(l)Ω_P/h_l, ν ∈ [−eκ₋, 0], monotone in Re ω above (S4 W7).
- On the rotated contour the Jost solutions are subdominant whenever Re ω + Im ω > 0 (S4 W2). The
  depth below κ₊ therefore costs nothing.

**Probe** (`results/s6/DESIGN_PB_e23_25.json`).
- Float thresholds (F > 0.02a*, W7) at e = 0.92:

| l | 0 | 1 | 2 | 3 | 4 | 5–6 | 7–9 | 10–19 |
|---|---|---|---|---|---|---|---|---|
| ζ_A/Ω_P | 21/4 | 37/20 | 3/2 | 27/20 | 5/4 | 6/5 | 23/20 | 11/10 |
| a′_l/κ₋ | 12.91 | 13.65 | 18.45 | 23.24 | 27.67 | 32.46, 38.37 | 42.42–53.74 | 56.81–105.51 |

- Arb confirmations at these ζ_A pass with W7 for l = 0, 1, 2, 5, 10 and 19. The margins are
  0.035, 0.039, 0.089, 0.053, 0.047 and 0.60 a*.
- a′_l is close to S3's P2-2 threshold a_l (Σa′_l = 1.04·Σa_l), except at l = 0–2.

**Level.** WRITTEN_PROOF, through S4's written W1, W2 and W4.

### O4. Part (i) at depth 0.547γ₀, l ≥ 20, h₀ = 2/41

Probes at e = 0.92 (`results/s6/*_e23_25.json`):

| regime | S4 design at e | S6 design at e |
|---|---|---|
| M3 windows (S2, per l) | m = 0.294, 0.344, 0.338 μ_HO at l = 20, 50, 100 (0.528, 0.568, 0.560 at d_δ) | (reference only) |
| T | fails on ζ ∈ [0.80, 0.95]Ω_P (−0.047 … −0.0009 a*) and above 1.70 | **redesigned**: s1 = 12, φ_∞ = −1, taper (3, 6), the rest as S4's T. Float 0.53 … 0.02 a* on [0.6, 1.7]; Arb 0.334, 0.308, 0.239 a* on [0.8, 0.85], [0.9, 0.95], [1.25, 1.3]. **T on [0.8, 1.3].** |
| A | ζ_A = 3/2: 6.63 a* | **ζ_A = 6/5**: Arb 1.972 a*, W7 holds |
| S (both sides, closing box) | c: 44/44 and closing pass; +: 43/44 and closing pass ([0.06, 0.08]Ω_P fails, −0.00027) | **design grid extended** (below): [0.06, 0.08] certifies at 0.0308 |
| L, band [0, 8] × [−19/20, −11/20] | 1272/1280 pass; the 8 failures are the boxes within 0.13κ₋ of ω₀ | **pole-free series R̃ = (1 + 2ρ)R in every box**: all 8 pass (Re y₊ ≥ 0.161), including the box containing ω₀; a 63-box random sample also passes (Re y₊ ≥ 0.128, Re y_c ≤ −0.079) |

**The S grid extension.** The S4 grid plus:
- w ∈ {1/2, 3/2, 5/2};
- y_a ∈ {1/40, 3/10, 3/5};
- φ_∞ ∈ {−π/5, −π/6};
- phase centre ∈ {−3/2, 1/2};
- phase width ∈ {3/4, 2}.

**E₊ inside L (W6 → W6′).** S4's W6 used n + 2ρ ≠ 0 "since Im ω > −κ_h". In S_e:
- n + 2ρ = 0 only for n = 1, ω = ω₀, since Re(n + 2ρ) ≥ n − eκ₋/κ₊ > 0 for n ≥ 2.
- R̃ has ã₀ = 1 + 2ρ and ã₁ = −(A₃ρ(ρ−1) + B₂ρ + C₁)/A₂. The factor is cancelled exactly; it is
  the scaled log coefficient and is nonzero, with value −d₁ at h = 0.
- R̃ is the solution regular at r₊ for every ω in the box:
  - off ω₀ it is a nonzero multiple of R₊;
  - at ω₀ it is the exponent-½ solution.
- y₊ = R̃′/R̃ does not depend on the normalisation, and W5 (Riccati invariance) is unchanged.

**The negative Jost margin κ₊ − eκ₋.** It is used nowhere:
- S, T and A use rotated contours, whose tails need only Re ω + Im ω > 0;
- L uses Frobenius series, not Jost Volterra.

S4's cover inequality "d_δ < κ₊ < κ_c" is replaced by "S_e ∩ (E₊ ∪ E_c) = {ω₀}, handled by W6′ and
Lemma E".

**Verdict: l₀ = 20 holds up in every regime at the probe level.**
**Level.** WRITTEN_PROOF, as for P2-1: W1, W2, W5 and W6′ are written.

### O5. Part (iii) for l ≤ 19: pole-aware band counts

**The boxes.** B_l = [−a′_l, a′_l] × [−b′, −b₀], where:
- b′ = eκ₋ rounded up to 2⁻³⁰;
- b₀ = κ₋/2 rounded down to 2⁻³⁰.

So B_l overlaps S5's R′_l (bottom −0.55000001κ₋) and P2-2's region (down to −0.55κ₋).

**The count.**
- S3's disc machinery (`part3.py`) evaluates W_l on ∂B_l. Every disc box clears E₊ ∪ E_c.
- Inside B_l, W_l is meromorphic with at most one pole, at ω₀, of order ≤ 1 (O1, step 2).
- The argument principle gives winding = Z − P with P ∈ {0, 1}.
- Expected:
  - **l = 0: winding 0**, with Z = 1 (ω_*, which lies in B₀ by Krawczyk) and P = 1;
  - **1 ≤ l ≤ 19: winding −1**, with P = 1 and Z = 0.
- These windings force P = 1, i.e. the residue at ω₀ is nonzero. That confirms Lemma E for l ≤ 19
  at the CERTIFIED level.
- The symmetry W_l(−ω̄) = conj W_l(ω) still halves the path (right half, /π), since ω₀ lies on the
  axis.

**Probes.**
- **Planned segments.** S3's planner gives 2562 segments (S3 itself: 2351). Of these, 41 are at
  l = 0 and 228 at l = 19.
- **Timing samples.**
  - l = 0: 1–2 s per segment near ω₀ and ω_*. The bottom edge passes 0.0106κ₋ below ω_* and
    certifies in one piece.
  - l = 19 near the PS mode: about 15 s per segment. S3's l = 19 average was 8.6 s.
- **KAT boxes** (`results/s6/KATBOX_PB.json`), full-boundary windings:

| l | box (κ₋ units) | contains | winding |
|---|---|---|---|
| 0 | [−0.03, 0.03] × [−0.86, −0.79] | ω₀ (pole) | **−1** |
| 1 | the same | ω₀ (pole) | **−1** |
| 0 | [−0.03, 0.03] × [−0.93, −0.79] | ω₀ and ω_* | **0** |
| 0 | [−0.03, 0.03] × [−0.93, −0.89] | ω_* | **+1** |

**Level.** CERTIFIED counts. The written steps are routine: the argument principle for a meromorphic
function, the symmetry, and the at-most-simple pole.

**The edge.** A shallower edge gains nothing material: e = 0.915 gives 0.5436γ₀ against 0.5466γ₀.
Every regime passes at 0.92, and the l = 0 box at 0.92 is well conditioned. **The edge stays 0.92**,
with no change to the frozen P2-6.

## 3. The cover of S_e for every l (Phase 2)

| region | l | method | certificate (new files only) | expected |
|---|---|---|---|---|
| S_δ = {Im ω ≥ −0.55κ₋} | all | P2-4 (S5) | existing, unchanged | the zero at 0 (l = 0) only |
| B, 0 ≤ Re ω ≤ 8κ₋ | ≥ 20 | L, band rows, R̃ in every box | `certificates/beta/part1/L_PB.json` | no zero; ω₀ handled by W6′ |
| −eκ₋ ≤ ν ≤ 0, 6κ₋ ≤ Re ω ≤ 0.9Ω_P/h | ≥ 20 | S, extended grid | `.../part1/S_PB_{plus,c}.json` | no zero |
| −eκ₋ ≤ ν ≤ 0, Re ω ∈ [0.8, 1.3]Ω_P/h | ≥ 20 | T, redesigned (10 ζ-boxes) | `.../part1/T_PB.json` | no zero |
| −eκ₋ ≤ ν ≤ 0, Re ω ≥ 1.2Ω_P/h | ≥ 20 | A, ζ_A = 6/5, W7 | `.../part1/A_PB.json` | no zero |
| B, abs(Re ω) ≤ a′_l | ≤ 19 | pole-aware band count | `.../part3/band_PB_l{l}.json` | winding 0 (l = 0), −1 (l ≥ 1) |
| −eκ₋ ≤ ν ≤ 0, Re ω ≥ a′_l | ≤ 19 | A at fixed h_l | `.../part2/A_fixed_PB.json` | no zero |
| Re ω < 0 | all | symmetry | `zeq_identity` Z4 | — |
| ω₀ | all | Lemma E | `.../epole_PB.json` plus the written argument | not a resonance |

- **l₀ = 20.**
- **Cover inequalities** (exact or Arb, in `certificates/beta/part1/cover_PB.json`):
  - 6 < 8;
  - ζ_TL = 0.8 < ζ_SU = 0.9;
  - ζ_A = 1.2 < ζ_TU = 1.3;
  - 0.9Ω_P/h₀ > 8κ₋;
  - Re ω + Im ω > 0 in S, T and A (6 > 0.92; 0.8Ω_P/h₀ and 1.2Ω_P/h₀ > 0.92κ₋);
  - eκ₋ < 2κ₊ and eκ₋ < κ_c (only ω₀ in S_e);
  - the l-ranges 0–19 and ≥ 20.
- **The rule for ζ_A(l), l ≤ 19.** Take the smallest multiple of 1/20 in [1.1, 12] (Ω_P units)
  with float margin ≥ 0.02a* and W7 (the O3 table). Certify it in Arb at h₀ = h_l.
  - If the Arb check fails, raise ζ_A by 1/20, up to 4 times.
  - Then a′_l := ζ_A(l)Ω_P/h_l, rounded up to 2⁻³⁰.
- **The rule for T.** One design (O4) for all ten ζ-boxes, with S4's target a*/10 and bisection
  depth 12.
  - If a box fails, T may use s1 ∈ {10, 14} or φ_∞ ∈ {−5/4, −3/4} for that box (the float grid of
    the probe), chosen by float maximin.
- **The rules for S and L.**
  - S: S4's procedure with the extended grid.
  - L: S4's procedure with R̃ everywhere and quadrisection down to 1/80.
- **The assembly check** (`src/p2/beta_assembly.py`, like `p2.assembly`). It reads every new
  certificate and the S5 ones, with sha256 pins, and checks:
  - the cover for every l;
  - the band counts with the pole correction;
  - the a′_l identities between the boxes and the A certificates;
  - Lemma E for l ≤ 19;
  - the Krawczyk enclosure and β.

  It writes `certificates/beta/assembly_beta_PB.json` and has negative controls.

## 4. Known-answer tests and negative controls (Phase 2)

**KATs** (each must give the stated result):

| KAT | test | expected |
|---|---|---|
| B1 | replay of P2-K6's Krawczyk ball (`krawczyk_PB_l0.json`) | unique zero; as stored |
| B2 | box around ω_* only, l = 0, with its top edge just above ω_* (e.g. −0.89κ₋) | +1 |
| B3 | box around ω₀ only, l = 0 and l = 1 | −1 (the pole is detected) |
| B4 | box around ω₀ and ω_*, l = 0 | 0 |
| B5 | residue from both sides: ε·W_l(ω₀ + ε) → i(2κ₊r₊ + L)/(2r₊²)·W̃_l, for l = 0, 1, 5, 19, with zsol and with the independent S1 code `arbw` | agreement to O(ε) |
| B6 | boxes around two zeros below the edge, with no pole inside: the l = 0 NE overtone ≈ −1.855iκ₋, and the l = 1 PS mode ≈ 7.2149 − 1.7533iκ₋ (S0 floats) | +1 each (the edge excludes them, not the method) |
| B7 | float cross-check (not used by the proof): \|W_l\| at 100 random (l, ω), l ∈ [20, 10³], ω in the band | no zero |

**Negative controls** (each must fail, or give the "wrong" count):

| control | test | must |
|---|---|---|
| NB1 | the l = 0 band box with its bottom edge at −0.905κ₋ (above ω_*) | give winding −1 |
| NB2 | a band edge through ω_* (Im ω = Im ω_* to 10 digits) | not certify (a piece enclosure contains 0) |
| NB3 | L's box containing ω₀ with S4's (non-pole-free) series | not certify |
| NB4 | T with S4's design at e on [0.8, 0.85] | not certify |
| NB5 | A at h₀ = 2 with ζ_A = 5 (below 21/4) | not certify |
| NB6 | the assembly check with one band certificate removed, a′₇ lowered by 2⁻²⁴, a winding set to 0 at l = 3, or b′ just above eκ₋ | fail on the named item (as `p2.assembly`'s N-controls) |

## 5. Kill criteria

- **K1 (E₊).** Stop and report if any of these happens:
  - the Phase-2 Arb certificate of W̃_l(ω₀) fails to exclude 0 for some l ≤ 19;
  - a band count for l ≤ 19 is inconsistent with exactly one simple pole at ω₀;
  - the written argument for l ≥ 1 is found wrong.
- **K2.** Stop and message "main" if the projected l₀ exceeds 100 or the projected compute exceeds
  40 core-h. The current projection is l₀ = 20 and about 10 core-h (§7).
- **K3.** Stop and report immediately on any certified zero of some W_l in S_e other than 0 and ω_*
  (both l = 0), or on ω₀ being a resonance for some l. Either would change β.
- **K4 (part (i)).** If part (i) cannot be certified at h₀ = 2/41 within the procedures of §3:
  - fall back to l₀ = 30, with band boxes and A at fixed h for 20 ≤ l < 30 (about 3 more core-h);
  - if that also fails, stop.
- **K5 (part (iii)).** If a band segment is not certifiable at ρ_min after one registered retry
  (ρ_min halved, K_MAX 160 → 240), stop and report. No unregistered change is made.

## 6. Level expected for P2-6

**WRITTEN_PROOF.** This is the minimum of the ingredients:

| ingredient | level |
|---|---|
| P2-4 | WRITTEN_PROOF |
| part (i) at e | WRITTEN_PROOF (W1, W2, W5, W6′) |
| part (ii′), A at fixed h | WRITTEN_PROOF (W1, W2, W4) |
| Lemma E for l ≥ 20 | WRITTEN_PROOF |
| D1′, set clause | WRITTEN_PROOF |
| part (iii) band counts | CERTIFIED |
| Lemma E for l ≤ 19 | CERTIFIED |
| ω_*'s Krawczyk ball | CERTIFIED |

§2's target for P2-6 is CERTIFIED. That would need P2-0's written steps, W1/W2/W6′ and D1′ to be
machine-checked, as for P2-4. **The expected level is therefore below the target.**

## 7. Budget and compute projection

The cap is 48 core-h, with at most 4 workers at nice 10. Production runs are detached and resumable
(JSONL plus logs, as in S3/S4). Commits are made at least daily and at each milestone.

| item | core-h (projected) |
|---|---|
| Lemma E certificates (l ≤ 19, and samples), sympy identity | 0.01 |
| part (i): L band (1280 boxes), S both sides with the extended grid, T, A, cover | 0.3 |
| part (i) controls NB3–NB5 and float check B7 | 1.0 |
| part (ii′): A at fixed h, l ≤ 19 | 0.1 |
| part (iii): band boxes l ≤ 19 (2562 planned segments, deeper edge) | 5–7 |
| part (iii): left halves l = 0, 2 (symmetry), 5% replay, KAT boxes B2–B4 and B6, NB1–NB2 | 1.2 |
| assembly check, tests (fast and `P2_SLOW`) | 0.4 |
| design contingency (T, S) | 1.0 |
| **total** | **≈ 9–11** |

**Phase 1 used about 0.5 core-h** (cap 4):
- M3 103 s;
- T/A 143 s;
- S 241 s;
- L 241 s;
- the first A-at-fixed-h scan 540 s (stopped; replaced by the float scan and Arb confirmations);
- timing and KAT boxes 200 s;
- design records 70 s;
- float studies.

Wall time: about 1 session for parts (i), (ii′) and Lemma E, about 2 h with 4 workers for
part (iii), and 1 session for the assembly, the proof note and the report.

## 8. FORMULATION: proposed changes (not applied; for the user)

Nothing in FORMULATION is edited in Phase 1. On approval, these go in a dated revision-log entry
(v1.1) at the start of Phase 2. A bracketed pointer is added to §2 P2-6; the frozen text is not
rewritten.

1. **P2-6's parenthesis "(where the l = 0 winding is 2)" is a correction.** At depth 0.92 the full
   l = 0 box holds two zeros (0 and ω_*) and **one simple pole (ω₀)**. The winding of W₀ is 1 on the
   full box and 0 on the band box. The v0.2 note already flagged this.
2. **P2-6's route "P2-1…P2-3 run with the strip edge at −0.92κ₋" changes method.**
   - S_δ is P2-4.
   - The band is covered by part (i) at e (T redesigned, S grid extended, L pole-free).
   - **Part (ii) is regime A at fixed h. This replaces Liouville–Green, which cannot reach
     Im ω < −κ₊.**
   - Part (iii) is pole-aware band counts.
   - Lemma E is new.

   The target is unchanged.
3. **New, additive.**
   - The band boxes B_l (D7 unchanged).
   - Lemma E.
   - The S6 KATs B1–B7.
4. **Unchanged.**
   - The edge e = 0.92.
   - The target β ∈ [β_lo, β_hi] with width ≤ 10⁻⁸. It is realised by P2-K6's ball, width
     6.4·10⁻³¹.
   - D1, D1′, D3 (Lemma E shows that the D1 and HV17 notions agree in S_e).
   - The scope fences.
5. **Level.** §2 states CERTIFIED as the target. The expected level is WRITTEN_PROOF (§6).

## 9. Phase-1 probe record (`results/s6/`)

| file | content |
|---|---|
| `E_PB.json` | sympy first step; positivity bound; W̃_l(ω₀) for l ≤ 19 and samples to 100; residue cross-check |
| `M3_PB_e23_25.json` | S2 windows at e for l = 20, 50, 100 |
| `TA_PB_e23_25.json` | S4 T (24 boxes) and A at e, S4 designs |
| `S_PB_{plus,c}_e23_25.json` | S4 S at e, S4 grid |
| `L_PB_e23_25.json`, `L_PB_e23_25_pole.json` | S4 L on the band; the box containing ω₀ with R̃ |
| `DESIGN_PB_e23_25.json` | T redesign (float and Arb), A at 6/5, A at fixed h (float table and Arb), S grid extension, L with R̃ on a sample |
| `P3_timing_PB.json` | band-box plans and timing samples (l = 0, 19) |
| `KATBOX_PB.json` | the four pole-aware KAT boxes |

## Amendments

All six are dated 2026-10-27 and made before any production run. A1–A4 come from the coordinator's
checkpoint review and Phase-2 go of 2026-10-27. A5 and A6 are the author's.

- **A1, the Lemma E write-up.** The exit proof note writes out the argument for l ≥ 1 in full:
  - the identification Φ₂ ∝ rR^{(2)} ∝ (hyperboloidal factor)·Z2, and of Φ_c, Φ̃_c with Z_c and the
    other branch at r_c;
  - both asymptotic expansions, with their error terms;
  - convergence of ∫VΦ₂Ψ dr* at both ends;
  - why l = 0 is excluded: V₀ changes sign.

  It also states the residue formula's normalisation, so that B5 tests exactly that formula:
  Z₊ = Σ_{m≥0} a_m t^m with a₀ = 1 and t = r − r₊ (zsol's `frob_z`), and Z2 = t + Σ_{m≥2} b_m t^m.
- **A2, the pins.**
  - `beta_assembly` pins every S5 certificate it reads by sha256, and checks each one against its
    blob on main `300c22a7` (as `p2.assembly` did).
  - Every new certificate records the commit it was computed at, clean or `+dirty`.
  - The stage exit lists any `+dirty` rows.
- **A3, ledger proposals at the exit.**
  - The exit proposes a P2-6 row in the house format.
  - It says whether Lemma E should be its own row or part of P2-6.
  - If an existing row (P2-D1′, P2-K6, P2-4) needs a note, it gives the exact text.
  - Nothing is applied.
- **A4, the frozen part of P2-K6.** β's enclosure reuses `certificates/part3/krawczyk_PB_l0.json`
  unchanged. B1 replays it into a new file and does not regenerate it.
- **A5, a typo in NB6.** "b′ just above eκ₋" should read "b′ just below eκ₋". That mutation must
  fail the check b′ ≥ eκ₋.
- **A6, how B7 is done.** l is drawn log-uniformly from [20, 1000] with seed 20261027. Re ω is taken
  in a regime's range, as in S4's K4, and Im ω uniformly in [−eκ₋, −0.55κ₋].
  - At each point the float W_l is evaluated, and a secant search is seeded there
    (`qnm_float.Solver.find`, at most 12 steps).
  - A stop counts as a root only under S4's Amendment-6 rule: |W(root)| < 10⁻⁸|W(seed)|, plus a
    zsol evaluation at 128 bits.
  - Expected: no root in the band.
  - This is float evidence and is not used by the proof.

## Dated note (2026-10-27, after the S6 review)

An erratum to O4 (review W10). It changes nothing in the plan, and O4's text above is left as it was
registered.
- O4 says that R̃'s coefficient ã₁ at ω₀ is "nonzero, with value −d₁ at h = 0". The value should read
  1/d₁ = 1/Δ′(r₊).
- In general ã₁(ω₀) = h²(2κ₊r₊ + L)/Δ′(r₊) = (1 + 2κ₊r₊h² − h²/4)/Δ′(r₊), which is positive for every
  h (the review's sympy check E9).
- The coefficient is nonzero either way, which is all that W6′ uses. `notes/s6-proof.md` §3 does not
  repeat the value.
