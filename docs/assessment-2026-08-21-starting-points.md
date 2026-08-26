# Starting-point assessment for the ten-problem portfolio

> **Status after the first working session (2026-08-25/26).** P9(a) has its
> first certified theorem: for the class C(G, L=1.5) (373-node grid), Δ = 4,
> frozen DESI DR2 BAO + Pantheon+ + Planck r_d (lower 2σ edge), every class
> member within Δχ² = 4 of the class best fit has H₀ ≤ 69.8417, and an explicit
> certified member reaches 69.5313 (so the true maximum is in [69.53, 69.84]).
> Certified feasible points exist for L ∈ {1, 2, 3, 5, 10} and Δ ∈ {1, 9}; the
> corresponding certified upper bounds are queued (`scripts/run_p9_curve.sh`).
> Getting there required three relaxation redesigns (`problems/P9/notes/
> relaxation-log.md`). P4 has the Evans–Coleman profile and spectrum
> reproduced (κ₁ = 2.8105525488, exactly one non-gauge unstable mode in the
> KHA box, GHJS-type monotonicity holds in V_z variables) and the sonic-point
> and centre expansions certified in Arb (Theorem A ingredients A1–A2); the
> validated shooting step (A3) is in progress. Ledger: `CLAIMS.md`.

Date: 2026-08-21. Scope: which entry points in
`docs/problems/open-problems-theoretical-cosmology-2026.tex` to attack first,
and what has to exist before any attack is credible. Literature claims below
were re-checked against arXiv/journals on 2026-08-21 (seven independent
search passes: P1, P2, P4, P5, P8, P9, plus a 20-entry citation audit).
P3, P6, P7, P10 were assessed from the document and prior knowledge only and
should get the same treatment before they are scheduled.

## 1. Recommendation in one paragraph

Start two tracks in parallel, on different verification stacks:
**P9(a)** (certified upper bound on H0 over an explicit late-time H(z) class;
convex duality against frozen DESI DR2 + SN + r_d; 4–8 weeks) and
**P4(b)** (computer-assisted existence of the Evans–Coleman radiation-fluid
self-similar solution plus a validated "exactly one unstable mode" count;
2+ quarters, existence first). P9 is the fastest route to a certified
artifact and builds the frozen-data/convex-dual stack; P4 is the best shot at
a theorem cosmologists would quote (a rigorous PBH critical exponent) and now
has a near-exact template (Guo–Hadžić–Jang–Schrecker, arXiv:2509.12435, Sep
2025, public VNODE-LP code). Add **P8** (mechanized Horndeski no-go + one
certified stable bounce) as the cheap third track that stands up the
computer-algebra/interval stack against known answers. Do **not** start P2
as written (its entry point does not deliver the stated theorem; see §3),
and do not start P5 until its definitional problems are fixed. P1, P3, P6,
P10 are long-tier and need formulation work first.

Before any of that: the repo is not a git repository, the problems document
has a placeholder author and several status claims that are already out of
date, and the only verification tool installed on this machine is Lean 4.
§4 lists the pre-work.

## 2. Ranking

| Rank | Entry point | Verdict after lit check | Effort to first certified result | Why this rank |
|---|---|---|---|---|
| 1 | P9(a) certified H0 bound | Open as methodology; physics answer (~69) is folklore | 4–8 weeks | Cheapest certificate; public data; convex (with one relaxation); builds frozen-data stack |
| 2 | P4(b) Evans–Coleman CSS existence + 1 unstable mode | Open; direct CAP template exists | 1 quarter (existence), 2–3 quarters (+spectrum) | Highest payoff/novelty among near-term; rigorous PBH exponent |
| 3 | P8 no-go mechanization + certified bounce | No-go solid under stated hypotheses; EFT-level boundary known by humans | days–2 wks (Lean core), 1–2 months (CAS derivation + certified example) | Low risk; stands up CAS + interval stack; real open piece (matter-coupled DHOST) comes later |
| 4 | P5 CET on tilted Bianchi A | Status claim wrong in letter; CET has no canonical Petrov-I definition | 2–4 wks numerical scan once definition fixed | Cheap counterexample paper possible, but only after definitional pre-work |
| 5 | P2 β>1/2 for near-extremal RNdS | Entry point as phrased does not yield the theorem | multi-quarter after redesign | Needs resonance-free strip uniform in ℓ (Hintz 2025 Conj. 1.5, open) |
| 6 | P1 period-3 cycle certificate | Qualitatively done (LHWG 2011); certified version open and hard | long | Non-hyperbolic passage; linearization fails at the 3-cycle (Buchner 2025) |
| — | P3, P6, P7, P10 | Not re-checked | — | See §3 notes; P10 entry point is not well-posed as stated |

## 3. Per-problem findings (deltas vs. the document)

### P9 — certified H0 bound (start)
- Prior art is all posterior/parametric: Bernal–Verde–Riess 1607.05617;
  Lemos et al. 1806.06781 (H0=68.42±0.88); Knox–Millea 1908.03663 (β≡c/(r_d H0)
  → H0≈69.0±0.95); Efstathiou 2103.08723 (the z<0.05 "hockey-stick" loophole:
  true H0 unconstrained at 70.5±3.6); Keeley–Shafieloo 2206.08440; Jiang et
  al. 2408.02365; Zhou et al. 2506.23556 (gradient-based "max H0" = 69.09±0.30,
  no certificate); Sabogal et al. 2607.21244 "shape wall" (GP, H0≲68.6–69.7).
  Nobody has framed it as a convex program with a dual certificate. Novel as
  method, not as number.
- Convexity: in u=1/H, D_M=c∫u, D_H=c·u, r_d scaling, SN offset M_B,
  smoothness/monotonicity constraints and the objective min u(0) are all
  linear/convex. Non-convex pieces: SN χ² in 5log10 D_L (concave log → use a
  rigorous outer relaxation: secant/tangent bounds on log over a data-derived
  box; essentially tight), D_V=(zD_M²D_H)^{1/3} for the BGS bin (drop or
  sandwich), r_d prior (outer 1-D loop). Relaxations only loosen, so the upper
  bound stays certified.
- Data (all public, machine-readable): DESI DR2 BAO means+cov
  github.com/CobayaSampler/bao_data/tree/master/desi_bao_dr2 (13 rows);
  Pantheon+ github.com/PantheonPlusSH0ES/DataRelease; DES-SN5YR
  github.com/des-science/DES-SN5YR; Union3 github.com/rubind/union3_release;
  r_d: Planck 147.09±0.26 Mpc, or CMB-free BBN+ω_m (Schöneberg et al.).
- Decisions that define the theorem: (1) the function class (derivative bound
  on ln H or u, z-range, behaviour above z=2.33 and below z≈0.01 — this alone
  decides whether the bound is ~69 or ~73); (2) feasibility threshold
  (χ²_min+Δ; Wilks is not rigorous over a function class); (3) continuous class
  vs. stated finite-dimensional spline class vs. grid with error term; (4) full
  SN covariance vs. diagonal boxes.
- Effort: 2–3 wks freeze data + cvxpy LP/SOCP + dual recovery; 2–4 wks exact
  rational/interval dual verification, SN relaxation, sensitivity.

### P4 — critical collapse (start, formulation first)
- Doc's "every rigorous statement for the fluid case is open" is too strong:
  Guo–Hadžić–Jang 2112.10826 (Ann. PDE 2023) prove existence of the
  relativistic Larson–Penston CSS solution for small ε=p/ρ (not ε=1/3);
  Guo–Hadžić–Jang–Schrecker 2509.12435 prove mode + nonlinear radial stability
  of the Newtonian LP profile with a CAP (energy identity for Re λ>1,
  high-order energy method for |Im λ|≫1, VNODE-LP validated ODE integration
  with hand-built Taylor series at origin and sonic point for the compact strip;
  code github.com/mrischrecker/Larson-Penston-Stability).
- Nothing rigorous exists for Evans–Coleman, for ε=1/3, or for any
  one-unstable-mode CSS solution (not even Newtonian Hunter (a)). Entry point
  is open. Also relevant: Singh–Zheng 2605.16095 / Zheng 2605.16235 (May 2026,
  stability of Christodoulou CSS scalar-field naked singularities), Li–Zhu
  2505.20766.
- Mode-count facts to encode in the formulation (Harada review gr-qc/0302004;
  Koike–Hara–Adachi gr-qc/9503007 γ=0.3558019; Maison gr-qc/9504008):
  Evans–Coleman = GR "Hunter (a)": exists for 0<k≤1, exactly one analytic
  unstable mode for all k, sonic point is a saddle for k≲0.41 (so at 1/3).
  Ori–Piran/GRLP: zero unstable modes, analytic only for k≲0.036. Hunter
  (b),(c),(d): 2,3,4 modes.
- Obstacles: (1) existence at ε=1/3 is non-perturbative (validated Frobenius
  at centre and sonic point, Krawczyk/Newton–Kantorovich shooting in the sonic
  point location + metric constant, continuation past the sonic point);
  (2) spectrum: non-self-adjoint coupled ODE system with singular endpoints;
  fix the function space (exclude the kink sector), identify gauge modes (time
  translation λ=1, lapse rescaling); a priori confinement of eigenvalues to a
  compact region — GHJS's large-|Im λ| argument uses LP monotonicity that EC
  (oscillatory) lacks, so this needs new analysis; then argument-principle count
  on an Evans/Wronskian function; (3) universality proper needs nonlinear
  codim-1 stability in GR, out of scope for the entry point.
- Tooling: VNODE-LP (closest template), CAPD, Arb/Arblib.jl, INTLAB,
  IntervalArithmetic.jl/RadiiPolynomial.jl; Reiterer–Trubowitz C source;
  Chen–Hou 2305.05660 and Buckmaster–Cao-Labora–Gómez-Serrano 2208.09445 as
  CAP templates for sonic-point Taylor expansions.

### P8 — bounces (start as infrastructure track)
- No-go holds under hypotheses that must be stated: flat FLRW, a>ε,
  G_T,F_T>ε asymptotically, Θ finite, unitary gauge regular. Inside Horndeski
  the "strong gravity in the past" loophole (Ageeva–Evseev–Melichev–Rubakov
  1810.00465, 2003.01202; Ageeva–Petrov–Rubakov 2104.13412, 2207.04071) gives
  stable bounces when ∫aF_T dt converges; a mechanized statement must carry
  the divergent-integral assumption. Rigorous formulation: Mironov, Universe
  5:52 (2019).
- EFT-level boundary already known: Creminelli–Pirtskhalava–Santoni–Trincherini
  1610.04207 (⁽³⁾R δN with sign-changing coefficient is necessary and
  sufficient); Ye–Piao 1901.02202 (GR + a single DHOST operator). Genuinely
  open: theorem-level classification over covariant DHOST operator subsets,
  minimality, matter coupling (Mironov–Rubakov–Volkova 2005.12626, 2011.14912
  exceptional subclass is necessary-only; An et al. 2501.09985 two-field DHOST
  bounce), UV/positivity (none exists).
- Part (a) (semiclassical singularity theorems): no Fewster–Kontou sequel since
  2108.12668; Graf–Kontou–Ohanyan–Schinnerl 2209.04347; Engelhardt–Nagar
  2605.05326. Open: reference-state-free QSEI, null QEI, non-minimal coupling.
- Tooling: field uses xAct/xPand (Mathematica) and cadabra2; no independent
  verification pair exists. No Lean/rigorous-numerics work anywhere; PhysLean
  has only basic FLRW.
- Effort: Lean core of the no-go given the F_S identity as axiom: days–2 wks;
  CAS derivation of the identity with cross-check: 1–2 months; certified bounce
  example (reproduce KMSV 1705.06626 / MRV 1807.08361 / Ye–Piao with interval
  sign certificates): 2–4 wks; new matter-coupled DHOST classification:
  3–6 months, research-grade.

### P2 — strong cosmic censorship (redesign before starting)
- The entry point as phrased ("validated enclosure of the dominant QNMs") does
  not prove β>1/2. That needs a resonance-free strip {Im σ > −κ₋/2−δ}∖{0}
  for all ℓ and all Re σ, i.e. a spectral-gap lower bound. Missing inputs,
  none available with explicit constants: (1) uniform large-ℓ exclusion
  (barrier-top/normally-hyperbolic-trapping results are asymptotic;
  re-derive via Olver-type explicit WKB error bounds); (2) explicit
  large-|Re σ| exclusion per ℓ; (3) only then validated argument-principle
  enclosure on the compact remainder.
- Hintz 2504.01734 (Apr 2025) rigorously pins the near-extremal family on
  |σ|≲κ_C — and corrects the CCDHJ near-extremal formula the document quotes
  (ω≈−i(ℓ+n+1)κ₋ is right only for ℓ=0; the general form is
  −i(λ_ℓ⁺+n)κ_C). His Conjecture 1.5 (all shallow QNMs are near-horizon ones)
  is exactly the missing gap statement and is stated as open; no 2025–26 paper
  closes it.
- "No rigorous verdict for the Λ>0 interior problem in any formulation" is
  overstated: interior conditional on exterior decay is settled
  (Hintz–Vasy 1512.08004, Costa–Franzen 1607.01018, CGNS 1707.08975, Rossetti
  2309.14420); rough-data Christodoulou-SCC is a theorem
  (Dafermos–Shlapentokh-Rothman 1805.08764). What remains open is the smooth
  exterior gap bound.
- Redesigned entry point: fixed explicit subextremal parameters (numerically
  β≈0.6–1), three-part theorem (explicit large-ℓ, explicit large-|Re σ|,
  validated compact remainder). Tooling angle: the RNdS radial ODE is
  Fuchsian, so Mezzarobba's rigorous D-finite connection (ore_algebra-analytic,
  1607.01967) plus Arb can compute connection coefficients/Wronskians with
  certified error. No certified-QNM package exists.

### P5 — gravitational entropy (definitional pre-work first)
- CET (1303.5612) is defined via the Bonilla–Senovilla square root of
  Bel–Robinson, unique only for Petrov D and N; CET themselves say other types
  are "not unique, further study necessary". Sarma–Nájera–Sussman 2605.20611
  (May 2026) give the first Petrov-I extension and it is non-unique (three
  factorizations, free couplings). Generic Bianchi class A is Petrov I, so
  "certified verdict on CET for Bianchi class A" is not well-posed until a
  factorization, trace/gauge fixing, frame vector z^a, congruence (fluid vs
  normal — they differ under tilt), time parameter and constants are fixed.
- Non-monotonicity is already documented (LTB decaying modes, collapsing
  Petrov-D stages, LRS Bianchi I with some sources: Sussman–Larena 1310.7632,
  1503.04589; Chakraborty–Guha–Goswami 1912.01414); the doc's "no known
  counterexample" is wrong in letter. Accurate: no tilted-Bianchi computation
  exists; no better alternative; no no-go theorem for observer-dependent
  Bel–Robinson ansätze (Pelavas–Lake gr-qc/9811085 only kills
  observer-independent dimensionless scalars, but note Bianchi equilibria are
  self-similar).
- Tilted class A dynamical systems exist in the literature (II: gr-qc/0008037,
  1004.3661; VI₀ gr-qc/0403040; VII₀ gr-qc/0509032; VIII gr-qc/0512070;
  IX incomplete); no public code.
- Also: the statement "non-decreasing along Einstein evolution" needs a
  slicing convention before it is a theorem statement.
- Effort: numerical scan 2–4 wks after definitions; certified monotonicity on
  compact type II plausible (SOS/CAD); VII₀/VIII/IX likely false.

### P1 — BKL (long tier; fix status text now)
- Status text is partly outdated: Béguin–Dutilleul, CMP 399 (2023) 737–927
  prove a positive-4-dim-Lebesgue-measure set of Bianchi VIII/IX data shadowing
  Kasner-map chains (Pesin theory); Reiterer–Trubowitz 1005.4908 cover a
  full-measure set of itineraries on the circle (not of initial data); Brehm
  1606.08058: a.e. solution converges to the Mixmaster attractor. "No
  full-measure statement about realized Kasner statistics" remains accurate.
  Suggested rewording: "codimension-one (LHWG, Béguin) and positive-measure
  (Béguin–Dutilleul 2023) families; full measure open."
- Entry point: LHWG 1004.1989 Thm 4.2 already gives a codim-1 Lipschitz stable
  manifold of the period-3 cycle with explicit-form but non-numeric constants;
  "open neighbourhood contracts" is the wrong phrasing (backward time expands
  x_u, x_s, contracts x_c) — the right target is an explicit cone-condition /
  graph-transform certificate with numeric ε₀, δ, rates. Buchner 2503.02664
  shows Takens linearization fails at the 3-cycle (resonance λ₂=λ₁+λ₃), so the
  CAP cannot lean on normal forms; passage near the line of equilibria defeats
  Lohner-type integrators. Only CAP precedent on Bianchi-type ODEs:
  Church–Hénot–Lappicy–Lessard–Sprink 2203.03763 (Hořava–Lifshitz, not GR;
  Julia RadiiPolynomial.jl).

### P3, P6, P7, P10 — not re-checked; notes
- P10: "for polynomial phase functions the saddle set is an algebraic variety"
  holds for FLRW (after the exact q-integral) but not for triaxial Bianchi IX,
  where the β± path integral is not Gaussian and saddles are not algebraic;
  the biaxial case was analysed by Janssen–Halliwell–Hertog (1904.11602). The
  "certified atlas" needs a definition of the integration domain before it is
  a target.
- P6: the positivity framework to adopt is probably the dS Källén–Lehmann /
  non-perturbative bootstrap line (Hogervorst–Penedones–Salehi Vaziri,
  Di Pietro–Gorbenko–Komatsu, Loparco et al.); axioms must be pinned first.
- P7(a): HKTT-type no-gos are scaling inequalities; mechanizing them is cheap
  and could be done in Lean, but the physics content sits in the 10D reduction
  that Lean will not check; value is mainly methodological. P7(c) certified
  flux enumeration is the substantive target.
- P3: known partial Λ>0 uniqueness proofs (Borghini–Mazzieri, Ambrozio) use
  non-polynomial weights in V; a polynomial-ansatz impossibility result may
  miss the real proof. Needs an invariant-theory canonicalizer for curvature
  polynomials (xAct/Invar-class) before search is meaningful.

## 4. Pre-work checklist

### 4.1 Pre-registration hygiene (days)
> Status 2026-08-21 (later the same day): the citation corrections, status-text
> fixes, and entry-point rewrites below have been applied to the `.tex` as
> revision 1.1 (see its Revision log); the PDF was rebuilt with tectonic and the
> original v1.0 PDF preserved at `docs/problems/revisions/`. Still open: git
> init/tag, author name, external timestamp, re-audit of P3/P6/P7.
- `git init`, commit, tag `problems-v1.0`; the document claims to freeze
  statements, which is meaningless without version control. Consider an
  external timestamp (arXiv/Zenodo) if the pre-registration is meant to be
  legible to outsiders.
- Fill `[Author Name]`.
- Citation audit (20 newest entries): 0 fabricated. Corrections: hubert2025
  print is Nature 651, 607–613 (2026) (online Nov 2025); chenhou2022 title ends
  "... with smooth data I: Analysis" (Part II = 2305.05660). Enrichments:
  nexus2026 authors Tsoukalas et al.; desi2025 → PRD 112, 083515; freedman2025
  → ApJ 985, 203 and label 68.8±1.8±1.3 as "JWST-only TRGB" (CCHP headline is
  70.4±1.2±1.3±0.7); cai2026 → RAA 26, 084011; poulin2023 full subtitle;
  tadpole2020 → JHEP 11 (2021) 223.
- Status-text fixes (the doc invites a revision log; start it): P1 measure
  status; P2 near-extremal formula (Hintz 2025) and "no rigorous verdict"
  overstatement; P4 fluid-case overstatement; P5 "no known counterexample";
  P8 hypotheses of the no-go.
- Entry-point rewrites: P2 (gap theorem, not enclosure), P1 (cone-condition
  certificate, not "contracting neighbourhood"), P5 (definition first), P10
  (integration domain first).

### 4.2 Per-problem formulation documents (1–2 weeks each, human-led)
One `docs/problems/P<n>-formulation.md` per pursued entry point, fixing:
exact statement, spaces/boundary conditions/gauge, the class or data being
frozen (with hashes), the certificate format and what checks it, acceptance
criteria, and the known-answer regression tests. This is the defence against
statement drift and the numerics-to-theorem gap the document itself flags.
P9 and P4 first.

### 4.3 Verification infrastructure (currently: Lean 4 only)
Installed: Lean 4 (elan 4.28/4.31). Missing: Julia, Mathematica, Sage/Maxima/
cadabra2, TeX, Arb/FLINT; the default `python3` is another project's venv.
- Project-local `uv` environment: numpy/scipy, cvxpy + a conic solver
  (Clarabel/SCS), python-flint (Arb), mpmath, sympy, astropy; exact rational
  dual verification via `fractions` or flint.
- Validated numerics: Julia (IntervalArithmetic.jl, RadiiPolynomial.jl,
  Arblib.jl) and/or C++ (CAPD, VNODE-LP — the latter is what the P4 template
  uses).
- Computer algebra pair for P8/P5/P3: cadabra2 + SymPy (free) vs. xAct
  (Mathematica licence decision).
- SOS/SDP with rational rounding: later (P6/P7/P3).
- TeX (for the problems doc and papers).

### 4.4 Harness: reuse Empiricist
`~/dev/empiricist` already implements the harness this document assumes:
epistemic ledger (HEURISTIC→CONJECTURED→VERIFIED_N→CERTIFIED→FORMALIZED,
promotion only with machine evidence), blake3 CAS, sandboxed executor, Claude
Code transport with the model never holding a shell, a kernel-anchored Lean
FORMALIZED gate, cost accounting. Decide early whether death_and_gravity is a
new Empiricist *domain* (as P5 FT-FBQC was) or a fork. Lessons carried over:
SEARCH generations cost ~$8–9 each; these cosmology entry points are
verification-heavy rather than search-heavy, so the ledger/provenance/
certificate pieces matter first and the SEARCH loop matters later (P5, P7, P8
search phases).

### 4.5 Concrete first-week actions
1. git init + tag; fix author; apply citation and status corrections; open
   revision log.
2. Write P9 formulation doc: class, threshold, discretization; download and
   hash DESI DR2 BAO, Pantheon+ (and DES-SN5YR/Union3 as alternates), fix r_d
   prior(s).
3. Write P4 formulation doc: similarity variables, sonic-point saddle
   structure at k=1/3, gauge modes, spectral domain; port the GHJS VNODE-LP
   setup as the starting skeleton.
4. Stand up `uv` env + Julia; decide CAS stack; decide Empiricist reuse.
