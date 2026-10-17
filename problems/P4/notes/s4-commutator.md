# S4-3: Lemma D closure and the commutator strip exclusion (b2)

> **Status banner, 2026-10-16 (recovery stage).  This banner overrides every
> "Status" or "Results" statement below that no stored artifact backs.**
>
> Certified, by artifacts produced 2026-10-16 from the committed code (commit
> `14da7d1e` sources):
> - the exact identities E1, E2, E3, E5 (§4): all hold.  E2 returns pmin = 14
>   with ok = True, i.e. val(Z∘u) ≥ 16 once A → 1.  Artifact:
>   `results/theorem_c/parts/checks.json`.
> - Lemma D (§1, Task 1).  `results/theorem_c/parts/lemma_d_certified.json`:
>   passes with centre balls that enclose the certified (a*, μ*).
>   `parts/lemma_d_midpoint.json` reproduces the `tc_energy.json` record
>   exactly, but that run's centre series sits on the a3 midpoint, which does not
>   enclose (a*, μ*).  See the §1 correction.
>
> Not certified: the commutator strip exclusion (Task 2).  Of the 576 region
> pieces (centre 30; tube 241 + band 18 + deep 129; gap 107; sonic 51), exactly
> one has a stored result: `parts/band/band-000.json`, the driver smoke test
> (pass, both phases, replayed bit for bit).  `results/theorem_c/tc_commutator.json` was never
> assembled.  The September region outputs were lost with `/tmp` in the
> 2026-09-07 reboot, and most of them came from earlier code versions.  The cloud
> production round was never retrieved.  The per-region numbers in §5 and §7–§9
> are September observations from the agent's transcript, not certificates.
> Inventory, cost model and re-run plan: `notes/s4-3-recovery.md`.  Driver:
> `scripts/s43_parts.py`.
>
> Update, 2026-10-16 evening (`s4-3-recovery.md` §11):
> - `_centre_expansion` and `tc_monotone.run_centre` now build the centre
>   series on the certified (a*, μ*) balls (`p4.validated.ec_centre`).
> - On that code, Lemma D reproduces `parts/lemma_d_certified.json` exactly
>   and P4-3 still passes.
> - The τ-quadratic-from-τ = 14 centre variant fails.
> - 310 region pieces are being computed; `s43_parts.py status` gives the
>   live count.

[Superseded 2026-10-16.  This status paragraph was written on 2026-09-01,
before the region runs finished, and no artifact supports it.]
Status: Task 1 complete (P4-8 now unconditional).  Task 2 certified with
b1 = 14 on tau in [14, 224] over the whole interval, and on tau > 224
everywhere except the near-sonic band xi <= XI_OPEN (value recorded in the
certificate JSON, `scope.xi_open_max`; ~0.048), which is a precisely
diagnosed open item (Section 8).  Code: `src/p4/validated/tc_commutator.py`
(+ `tc_tm.py`); certificate: `results/theorem_c/tc_commutator.json`.

## 1. Task 1 — Lemma D closed exactly

The deep-tail conditions (a)-(d) of `s4-energy.md` par. 4 on x <= -4.5 are now
certified unconditionally.  The obstruction was junk in the t^0/t^1 layers of
the centre reduction; it is removed by the exact layer identity

  (E2)  Z := QNr - (9/4) Dp^3 SDr^3 D2r^2 Bcr^2  has no centre layers below
        t^16 and no odd layers (fmpq_mpoly, centre grading p = -eN + 2eW + eV,
        A -> 1), while val(QD o u) = 14,

so q2, q1, q0 - 9/4 and the remainder carry honest t^2 with no junk slots.
One Taylor-model run over y = e^x in (0, e^{-4.5}] (1100 subranges, ratio
1.008, 62 s) certifies the hat bounds:

  q2_hat >= 1.9618 > 0,  mB = +3.475,  mC = +15.348,
  pole gaps |y/A|, |y F_N| >= 0.9995 away from 1.

Recorded in `results/theorem_c/tc_energy.json` under `"lemma_d"`
(`all_ok = true`); replay: `uv run python -m p4.validated.tc_commutator
lemma_d`.  P4-8 loses its conditionality (draft row in Section 9).

**Correction, 2026-10-16.**  (a) The replay command above prints its result
and writes nothing.  The `tc_energy.json` block was merged by hand from the
scratchpad file `lemmad_rep.json` (transcript L706).  (b) The re-run on
2026-10-16 reproduces that record exactly: ok, 1100/1100 subranges, the same
five minima, 62.6 s.  Artifact: `parts/lemma_d_midpoint.json`.  (c) That record
does not cover the EC background, however.  `_centre_expansion` builds the
centre series from the midpoint (a_m, μ_m) of `results/a3_midpoint.json`
(radius ≤ 1e-41) and omits the stored Krawczyk offset `K_minus_m`.  The
certified values are a* = a_m + 7.031e-12 ± 8.87e-16 and μ* = μ_m + 2.3757e-9 ±
4.61e-14, the balls that `modecount.certified_centre` and every tube build use.
With those balls Lemma D passes: q2 ≥ 1.9618364, mB ≥ 3.4752428,
mC ≥ 15.348292, pole gaps ≥ 0.99952 and 0.99999, 61.5 s.  Artifact:
`results/theorem_c/parts/lemma_d_certified.json`.  Replay:
`python problems/P4/scripts/s43_parts.py lemma-d --centre-balls certified`.
Lemma D is certified by that artifact.  The P4-8 ledger update is for Sean to
make (draft in `s4-3-recovery.md` §10).  `tc_monotone.run_centre` had the same
midpoint construction.  Both are fixed in the code since 2026-10-16 evening
(`p4.validated.ec_centre`; `s4-3-recovery.md` §11).

## 2. The commuted operator (m = 1)

Write the reduced eigenvalue ODE as chi'' + A0 chi' + B0 chi = 0 with
A0 = nT/dT, B0 = nB/dB (the S4-2 value polynomials; B0 = -U, nB = -U6).
One commutation step (chi1 := chi' satisfies the once-commuted equation
chi1'' + A1 chi1' + B1 chi1 = 0 with)

  A1 = A0 - (log B0)',   B1 = B0 + A0' - A0 (log B0)',

is legal wherever B0 != 0; on the strip this holds by the B0-quadratic
dominance (Section 5, thU < 1).  m = 2 is structurally dead: B1 has zeros on
the centre layer curve d(tau t)^2 ~ 6, so its log-derivative blows up.

## 3. The energy identity and the pointwise condition

Let kappa = eps + i tau (eps in [0, 6], tau >= b1; tau > 0 suffices by
conjugation symmetry), c = -1/32, and let t-frak = tt be the real kappa-lead
of T.  With the weight g = exp int (Re A1 - c tt) and the combination
Re - (tau/c) Im of the chi1-energy identity:

  0 = BT - int g |chi1'|^2 + int g K |chi1|^2 - int g nu {chi1},
  K = P - (tau/c) Q - (f' + rho f)/2,   P + iQ = B1,
  nu = Im A1 + tau tt,  rho = Re A1 - c tt,
  f = (c + tau^2/c) tt - (tau/c) nu,
  BT = Re[(1 + i tau/c) g chi1-bar chi1'] + (1/2) g f |chi1|^2.

Young's inequality on the nu-term gives: if

  S := K + nu^2/4 < 0  pointwise on the interval,

then chi1 == 0, hence chi == 0: no eigenvalue at that kappa.  (Toy-verified
numerically to 5e-10 before certification.)

Boundary terms: at the centre, g ~ e^{(3-2m)x} = e^{x} kills BT as
x -> -infty against the L^2-normalized decay; at the sonic point,
g ~ xi^{2m - (1-sigma0)(1-eps+|c|)} with the m = 1 inner-layer constraint
(1-sigma0)(1+|c|) < 1, i.e. |c| < 0.099 (satisfied by c = -1/32), so BT -> 0
at xi -> 0 on the eigenfunction branch.  Function space: on the strip
Re sigma(kappa) < 0 is automatic since Re kappa >= 0 > -sigma0/(1-sigma0)
= -0.0991, so the embedding needs no extra argument.  The kappa-poles of the
coefficients sit at real kappa in [kappa-bar, 1.89]: certified per piece via
|A|, |F_N| < b1 (real-pole bound), distance >= b1 - 1.89 >= 12.1 from the
excluded region.  Eta-elimination is legal since the Bcr cubic's roots obey
the Cauchy bound Rbar < b1 on every piece (certified; worst Rbar ~ 6.5).
[2026-10-17: "certified per piece" holds for the pieces carried by S4-3
parts.  A piece carried by a substituted S4-4 `layer_t14` part (`s43_parts.py
merge --substitute-layer-t14`, `s4-4-layer-engine.md` section 14) has no Rbar
or real-pole bound; there dT != 0 at every strip kappa of the piece is
certified instead, and since dT = Dp SDr D2r Bcr it gives Bcr != 0 and
D2r != 0 (no pole) at those kappa directly.  The real-pole bound is vacuous
in the strip anyway (the roots A and -F_N of D2r are real); review
`docs/review-2026-10-17-p4-layer-engine-adversarial.md` section 1.]

## 4. Certified structure constants and input raising

At the sonic point (exact identities in Q(V0, sqrt3)):

  (E1)  U6_9 * (-Dp) = (3 - V^2) * [kappa^7 coeff of denU],  Dp^2 | U6_9;
  (E3)  D(u0) = PW(u0) = PV(u0) = 0  (simple sonic vanishing: du, Delta~);
  (E5)  U6_j(u0) = 0 for ALL j = 0..9  (U6 = U denU, denU ~ xi^3,
        U ~ xi^{-2}),  and  Aer_3(u0) = Ber_3(u0) = 0.

E1/E3/E5 are named checks (`check_u6_top_identity`, `check_sonic_zeros`,
`check_sonic_vanishing`) run at merge time.  They justify
`sonic_raise_inputs`: dB jets raised by 3 Laurent orders, dT by 1, ALL nB
kdiv inputs by 1 (the top U6_9 by 2).  With raised inputs the B0-quadratic
kdiv is junk-free with certified hats

  u2 = -0.3951/xi, u1 = -0.909894/xi^2, u0 = +0.909894/xi^2 (u-hats),
  w0-hat = -sigma0/xi = -0.0901060/xi.

Two more containment-checked raises at the closing domain: nW-top by 1
(A1 has a simple pole while dA1-top = dT nB dB tops has val exactly 5), and
the G-remainder slots to val 4 (Section 6).

## 5. Certification architecture (per piece, two phases)

Each region piece builds atom kappa-polynomials of Laurent-Taylor models
(`tc_tm.py`: TM/LC/JLC, CAP = 80, recentred Taylor-shift evaluation), forms
the pairs (nT, dT, nB, dB), commutes once, and certifies S < 0 for all
eps in [0, 6] (13 graded cells) and tau >= b1 = 14 in two phases:

- direct phase, tau in [b1, tstar], tstar = 16 b1 = 224: adaptive ball
  evaluation of S on (ell-strip) x (eps, tau) boxes through value chains
  (A0, B0 quadratic + remainder jets, log-derivatives at depth 2).  On
  domains touching ell = 0 the chains carry (hull, ell-offset) pairs, so the
  sign of the hat certifies the open interval (0, hi].
- tau-quadratic phase, tau >= tstar: S <= K2 tau^2 + R1 tau + K0 + R0 with
  K2, K0 explicit LC-hats from the kdiv quotients (B1-quad b2, b1, b0;
  A0-lead tt; A1-quotient w0) and R1/R0 remainder slacks certified by
  `sup_ratio_boxes`: factored denominator lowers (|dT|, |dB| via atom hats,
  D2r through its exact linear factors (kappa - A)(kappa + F_N), |U6| >=
  (1 - thU)|quad||dB| after certifying thU < 1), kappa-boxes geometric in
  tau ending in one pseudo-box w = 1/kappa in [0, 1/T], and per-name Laurent
  allowances (net offsets).  The vertex bound closes sup_{tau >= tstar}.

Region assembly (all with b1 = 14, eps_max = 6).  [Correction 2026-10-16.
The region list below is the September design.  The pieces actually needed are
the 241 tube_full steps with x - h >= -3.05, minus the 4 h = 0.0075 steps in
[-0.169, -0.139]; plus 18 band re-steps (h = 0.002) and 129 deep re-steps
(h = 0.0125, x in [-4.56, -2.95]).  The 33 replaced steps fail their chain TM
inverses.  In the September per-step logs the gap tube's tau-quadratic phase
failed on steps 0-88 (xi <= 0.0357) and passed on steps 89-106, not "for
x <= -0.048".  The full list is in `results/theorem_c/parts/plan.json` and
`s4-3-recovery.md` §4.]

  centre   y = e^x in (0, e^{-4.5}], 30 pieces (geometric + closing piece);
  tube     x in [-4.5, -0.05], 274 certified ODE steps (tube_full);
  gap tube x in [-0.055, -0.0085], 107 new certified steps (Tube.build from
           x0 = -0.0085, same V0 interval as the main tube); both phases
           pass only for x <= -0.048 -- the steps nearer the sonic point
           are certified in the direct phase and their tau > tstar tail is
           open (the tube pieces have no Laurent structure, so the b0
           xi^{-3} split of Section 6 is not available to them);
  sonic    xi = -x in (0, 0.0091]: 50 Laurent pieces (geometric ratio 1.2
           from 1e-6) + the closing piece [0, 1e-6], each on its own domain
           [0, hi].  (The recentred quasi-step machinery is not used in the
           final assembly: its Taylor-shift hulls blow up below xi ~ 0.02,
           and the gap tube covers [0.0085, 0.055] instead.)

## 6. The sonic tau-quadratic: what certifies and what obstructs

On sonic pieces the explicit b0-quotient carries a GENUINE xi^{-3} hat: in
the kappa -> infty expansion b0 = u0 + t0' - t0 (log u2)' - tt lambda1 with
lambda1 = [u1' - (log u2)' u1]/u2, giving b0^{(-3)} = tt-hat u1-hat/u2-hat
= 0.91 * 2.303 = +2.10 / xi^3.  At any fixed kappa this mass is cancelled by
the remainder: the ell^3 layers of remB1 equal -b0^{(-3)} times the ell^6
layers of dB1.  The certified split keeps the cancellation:

  b0 = b0h + b0l (exact slot split),  B1 = b2 k^2 + b1 k + b0h + G/dB1,
  G := b0l dB1 + remB1  (slots containment-raised to val 4),
  Im(G/dB1) = Im(remB1/dB1)  (b0l is real and kappa-free).

Slack routing: R0 += DG ell^{-2} (Re part), R1 += DG/|c| ell^{-2} (Im part,
tau-linear); nu- and nu'-remainders use the real-rational tail lemma

  |F(eps + i tau)| <= (sup_{t >= tau} t^2 |F'|) / tau  for real-coefficient
  F with F(i infty) = 0, F' the kappa-derivative (gap-2 numerators
  NW1 = remW' dA1 - remW dA1', and the double combination for nu'),

whose sups are certified by the same box machinery (third output slot).
Den lowers at boxes where the quad hull self-cancels use the sign-aligned
branches |quadN| >= (l2/sqrt2) ell + cI |w| (min over valid bookkeepings).

This closes every Laurent piece's direct phase and all centre/tube pieces'
both phases.  It does NOT close the sonic tau > tstar tail: see Section 8.

## 7. Results (certificate summary) [SUPERSEDED 2026-10-16]

**Correction, 2026-10-16.**  The summary below was written on 2026-09-01 as
the expected outcome, before the runs finished.  `tc_commutator.json` was never
written, so no line except the first is backed by an artifact.  Status of each
line, from the transcript (details in `s4-3-recovery.md` §4):

- checks: re-verified 2026-10-16, all hold (`parts/checks.json`).  E2 returns
  pmin = 14, which gives val >= 16 after A -> 1.
- centre: not certified.  The 30/30 run (worst -5.99, 1022 s) started at
  09-01T04:18Z with code that had no direct phase: the tau-quadratic covered
  tau >= 14 directly, and `direct_S_boxes` was added at 04:32Z.  No run of the
  committed code finished (one was killed after 26 h wall / 25.3 h CPU).  The
  centre series was also the midpoint one (§1 correction).
- tube: not certified.  Steps 184-273 never finished locally.  33 steps failed
  their TM inverses and were re-stepped (band, deep).  12 steps in
  x in [-1.79, -1.54] failed the direct phase at the old strip-depth cap
  (sdepth < 3).  With the fix (sdepth < 9), two of the 16 cloud chunks are
  known to have passed: A1 (8 steps, worst -63145) and one that was never
  identified.  Neither was retrieved.  The other steps ran code older than the
  committed version.
- gap tube: not certified.  The tau-quadratic phase failed on steps 0-88
  (xi <= 0.0357).  A degrade run passed the direct phase on all of them, and
  steps 89-106 passed both phases.  The outputs are lost.
- sonic: not certified.  51/51 passed the direct phase in September (worst
  -5.9e-5; tau > 224 open), with code of 09-01T09:45Z.  The outputs are lost.

Original text, kept for reference:

  checks   E1, E2 (val >= 16), E3, E5 all pass exactly.
  centre   30/30 pieces, both phases; worst scaled margin -5.99.
  tube     274/274 steps, both phases; worst margin -1.7e5 (unscaled),
           hardest step x = -0.05: tau-quad -7.1e7, direct -0.17 (scaled).
  gap tube 107/107 steps, both phases (worst margins recorded in the JSON).
  sonic    51/51 pieces direct phase (worst hat margins ~ -0.003 .. -39);
           tau > 224 open on these pieces (Section 8).

Replay [superseded 2026-10-16 by `scripts/s43_parts.py`, which gives
per-piece resumable parts and a coverage-checked merge; see
`s4-3-recovery.md` §5]: `uv run python -m p4.validated.tc_commutator
checks|centre|sonic|tube <tube.json> [a b]` then `... merge <parts>` writes
`results/theorem_c/tc_commutator.json` (region parts can be sliced for
worker parallelism; the shipped JSON lists every piece).

## 8. The open item: near-sonic band x (tau > 224)

On the 51 sonic pieces (xi <= 0.0091) and on the gap-tube / main-tube steps
with xi below ~0.048 only tau in [14, 224] is certified.
The tau > tstar quadratic organization provably cannot close there: every
quotient/remainder split of B1 (and of nu, nu') separates pairs whose
cancellation is (ell tau)-coupled.  The truth is clean and scale-invariant:

  S ell^2 / tau^2 = F(s, eps) + o(1),  s = ell tau,

with F independent of ell to three digits already at ell = 1e-4 vs 1e-5
(float check, `corner3.py`): F(0+) = -0.90, F(1) = -3.21, F(infty) = -15.5
(= the certified K2 hat) at eps = 0, and F <= -41 at eps = 3.  So the
statement is true with a uniform margin >= 0.9 tau^2 ell^{-2}; certifying it
needs a JOINT (w = 1/kappa, ell) engine — a 2D Taylor model in (s, ell) or
(w, ell) evaluating the whole S-expression, exactly as the ell = 0-safe
direct chains do per box, but with the tau -> infty normalization built in.
This is bounded, mechanical work (~150 lines on top of `tc_tm.py`), flagged
for S4-4.  Two smaller alternatives: (i) an indicial/Frobenius endpoint
argument at the regular singular point replacing pointwise S < 0 on a fixed
sonic neighbourhood for large tau; (ii) a large-|kappa| perturbative
exclusion (GHJS Prop 3.14 style) for tau > 224, which would also make the
whole tau-quadratic phase redundant.

[2026-10-16: nothing in this paragraph is certified yet.  It states the
September expectation, and the September data put the open band at
xi <= 0.0357, not ~0.048 (see §5).]
Consequences for the b2 target: the certified exclusion is
  tau in [14, 224] on the full strip  (b1 = 14 achieved there), plus
  tau > 224 for xi > XI_OPEN (~0.048) pointwise.
No S4-4 tau-bridge below 224 is needed; the S4-4 hand-off is the single
large-tau sonic-sliver item above.

## 9. Draft CLAIMS rows (do not edit CLAIMS.md from this stage)

[Correction 2026-10-16.  The P4-8 update now rests on
`results/theorem_c/parts/lemma_d_certified.json` (centre balls enclosing
(a*, mu*)), not on the midpoint record in tc_energy.json (§1).  The P4-9 draft
is superseded: no certificate exists, and its piece counts are the September
design, not the 576-piece plan.  The current draft text is in
`s4-3-recovery.md` §10.]

  P4-8  (update)  Deep-tail energy conditions x <= -4.5 certified
        unconditionally (Lemma D closed via the exact layer identity E2 +
        one TM run; tc_energy.json all_ok = true).  Drops the previous
        "conditional on Lemma D" qualifier.

  P4-9  (new)  Commutator strip exclusion: on Re kappa in [0, 6], the
        once-commuted energy certificate (m = 1, c = -1/32) excludes
        eigenvalues with |Im kappa| in [14, 224] (full interval), and with
        |Im kappa| > 224 away from the near-sonic band xi <= ~0.048.
        Certificate: results/theorem_c/tc_commutator.json (E1/E2/E3/E5
        exact identities; 30 centre + 274 tube + 107 gap-tube pieces both
        phases; 51 sonic pieces direct phase).  Open: sliver x (tau > 224),
        see notes/s4-commutator.md Section 8.

## 10. S4-4 / S4-5 hand-off

- S4-4 (bridge stage): the only remaining task for full b1 = 14 strip
  exclusion is tau > 224 on the near-sonic band xi <= ~0.048 (Section 8;
  recommended: the large-kappa perturbative estimate, which would subsume
  all tau > 224 work, or the 2D (w, ell) engine on the band).  [2026-10-16:
  reassessed in `s4-3-recovery.md` §8.  The large-kappa matching route is not
  recommended for S4-4; the energy-framework layer engine is.  The September
  data put the band at xi <= 0.0357.]
- S4-5 (endgame): with the strip clear, ker = span{g}, S Delta-tilde != 0
  and the gauge argument proceed as planned (theorem-c-plan.md par. 6); the
  eta-elimination root bound (Rbar < 14) and the real-pole margins needed
  there are certified per piece in the JSON.  [2026-10-16: not yet.  These
  bounds are checked inside each piece run (`analyze_piece` raises if they
  fail), so they become certified piece by piece as parts land.
  2026-10-17: for pieces carried by a substituted `layer_t14` part the endgame
  uses dT != 0 at strip kappa instead (every polynomial it needs nonzero
  there divides dT up to powers of Dp: D2 = Dp^3 D2r, Bc = Dp^2 Bcr,
  S Delta-tilde = Dp SDr); Rbar and the real-pole margins are not recorded
  for those pieces and not needed (review 2026-10-17 section 1).]
- The dec-slots (DWdec etc.) remain in the JSON for the non-closing pieces;
  on sonic pieces the tail-lemma sups (KW/KWp third slots) replace them.

## 11. Engineering notes

- `sup_ratio_boxes` branch bookkeeping: branch-R replaces the quad lower
  lQ ell^{oQ} by lQR ell^{oQ+1} (net DECREASES by the number of lQ powers);
  an early version had the sign wrong — fixed; only the sonic closing piece
  ever used those branches, and its certificate was regenerated after.
- The direct phase's offset-aware chains (Section 5) subsume the previous
  scl-multiplication; margins are reported at hat level.
- The recentred sonic path (`sonic_components_recentred`) is kept in the
  code but unused by `run_sonic` (Taylor-shift hulls degrade below
  xi ~ 0.02); the gap tube supersedes it.
- Tube extension: `Tube.build(V0, 1e-16, x0=-0.0085, x_c=-0.055, K=28,
  hmax=0.004, prec=384)`, 107 steps, 134 s; saved as tube_gap.json next to
  the certificate inputs.  [2026-10-16: that file was lost.  All four tube
  caches are rebuilt by `scripts/s43_parts.py build-caches` into
  `problems/P4/data/cache/` (gap: 107 steps, 121 s), pinned by content hash in
  `results/theorem_c/parts/caches/`.]  `run_tube(..., degrade=True)` applies the same
  direct-phase fallback as `run_sonic` on the band.
- The centre certificate in the JSON is the original run (sound earlier
  code path; the region never touches the branch bookkeeping); a
  re-validation with the final code was started and aborted for CPU.
  [2026-10-16: no JSON exists.  That run (09-01T04:18Z) predates the direct
  phase and used the midpoint centre series, so it is not a run of the
  committed code (§7 correction).]
- Deep-tube diagnosis (x in [-4.56, -1.5]): the direct-phase failures there
  were pinned to the (tau ~ 14, eps ~ 0) corner on the first ell-strips: the
  chain-hull slop scales with the STRIP width (atom variation amplified
  ~1e3x through the nu = Im A1 + tau*tt cancellation), while (eps, tau) box
  splitting cannot shrink it -- so the old strip-split cap (sdepth < 3)
  pinned the hull at ~ +7% of scale against a true margin of -9%
  (float-verified: S < 0 on the whole band, no B0 zeros; min |B0| sits at
  the tau = 14 edge and decays like e^{2x}).  Raising the cap to sdepth < 9
  certifies the worst cell at -8.1e-4.  The h = 0.05 steps at x <= -3.05
  additionally fail their chain TM inverses outright (|q| hull >= 1) and are
  re-stepped at h = 0.0125 (tube_deep.json, 129 steps via extend_centre).
  Partial old-parameterization runs (tchunk/dchunk/tband2, gap_w0-4) are
  kept in the scratchpad as diagnostics only; they are superseded by the
  cloud re-runs with the strip-depth fix.  [2026-10-16: the scratchpad was
  deleted in the 2026-09-07 reboot.  The cloud results were not retrieved.
  See `s4-3-recovery.md` §2 and §6.]
