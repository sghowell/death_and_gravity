# RATE4-D8.WITNESS: retained-source obstruction and the next input boundary

Date: 2026-09-22. This executes the parent-evidence milestone admitted by
the [RATE4-D8 conditional construction](assessment-2026-09-21-p8-rate4-d8-candidate.md).

**Outcome: the retained-source route is insufficient. No physical witness
has been established.** The specified classical action, scoped quantum
prescriptions and four rate targets do not imply the proposed tail or
pole/source enclosures. The conditional construction remains valid, but
its hypotheses have not become verified properties of a quantum parent.
Original RATE4, all frozen P8 calculations and M/V/G/B/R statuses are
unchanged. P8 remains **OPEN**.

The new result is stronger than finding one missing Taylor coefficient:
even any finite weighted-degree S/U coefficient prefix, together
with the same four first-order rates and known first-order cuts, does not
bound the remaining analytic tail in the formal EFT comparison class.
Analyticity without a magnitude bound does not supply the Cauchy witness.
The pole/source inputs are independent even when the analytic tail is
held fixed. These are data-sufficiency results, **not** unitary UV
counterexamples, a failure of a verified D8 member, or an exclusion of
every possible common parent.

## 1. Trace the proposed inputs to their actual owners

The audit reads the current source chain rather than identifying a
classical Lagrangian with a complete physical quantum prescription.

| Retained source | What it actually supplies | What cannot be imported from it |
|---|---|---|
| [S238 full covariant parent](../problems/P8/s6/continuation/s6_238/FORMULATION.md) | A specified classical action, tree matching and scoped classical clock-tube bounds | Independent finite-gravity matching values or a norm of its unknown quantum tail |
| [S239 matter prescription](../problems/P8/s6/continuation/s6_239/FORMULATION.md) | Complete formal limiting-action first-loop matter amplitude, with specified flat finite conditions | Commutation of finite-gravity quantization and decoupling, complete curved counterfunctional or omitted-loop bounds |
| [S240 QG2 state/profile](../problems/P8/s6/continuation/s6_240/FORMULATION.md) | The actual H8A420-SLE-PRE-T0 Gaussian heavy state, fixed subtraction scheme and reference stress profile | A complete interacting light/metric state, global quantitative response or same-parent nonlinear bounce |
| [S294 matching](../problems/P8/s6/continuation/s6_294/notes/matching.md) | Known mixed-loop representative and an explicit finite curvature/residue distinction | A value for c_RH, additional heavy mass/residue matching or higher-EFT coefficients |
| [S297 source](../problems/P8/s6/continuation/s6_297/notes/source.md) and [matching](../problems/P8/s6/continuation/s6_297/notes/matching.md) | Finite-order equality of the complete relevant matter jets/graphs and the stated inclusive assembly | New finite gravitational matching or a common interacting quantum parent |
| [S302 matching](../problems/P8/s6/continuation/s6_302/notes/matching.md) | A definite calculable minimal metric-loop representative | The independent finite alpha, beta and delta_kappa inputs; its API explicitly requires them |
| [S336 curved comparison](../problems/P8/s6/continuation/s6_336/notes/matching.md) | A real radiative curvature direction invisible to retained flat first-order matching | Its physical coefficient from the known flat or soft data |

S297's [formulation](../problems/P8/s6/continuation/s6_297/FORMULATION.md)
also makes the order/domain limits explicit. The source manifest pins all
ten directly reviewed files, including its [contour boundary](../problems/P8/s6/continuation/s6_297/notes/forward.md).

This is not a claim that no state or renormalization information exists:
S239 and S240 supply substantial, specific data. Their scope is the reason
they cannot silently fill the new inputs. Similarly, the smoothness of a
background-field profile is neither a proof nor a refutation of holomorphy
in the different variables S and U. The decisive issue here is the
unfixed finite quantum matching data, not that distinction in regularity.

## 2. A tail obstruction beyond any finite coefficient prefix

Use mass-one units, `S=s^2+t^2+u^2`, `U=stu`, and the same polynomial as
the earlier remainder audit:

```text
R5=3259SU-18645S^2+2463528U-152808S+21316224.
```

It vanishes at all four retained RATE4 centers. For any finite weighted
degree D whose coefficients have been determined, choose an integer

```text
k>=3, 2k>D,             W_k(S,U)=S^k R5(S,U).
```

Every monomial of W_k has weighted degree at least 2k, where the weights
of S and U are two and three. Thus adding `zeta W_k`:

- preserves every coefficient through weighted degree D, including all
  four D8 fitted coefficients;
- vanishes at the four matching points without refitting those coefficients;
- is a real crossing-symmetric polynomial, hence is entire in S and U
  and introduces no pole or branch cut;
- can be realized as a formal finite local four-scalar counterterm by the
  same labelled-slot construction, now with at most `4k+10` derivatives.

Take this comparison at order hbar. Its insertion into a loop or into the
first real-radiation correction costs a further order, so the retained
first-order absorptive cuts and rate-conversion term are unchanged.
Multiplication by the existing V(X) preserves the reference-clock action
jets through order 1023. These statements do not bound an arbitrary
nearby curved tube or higher-order effects.

If the formal gravitational decoupling bookkeeping is also to be retained,
one can multiply this quantum-order comparison by `(kappa0/kappa)^2`.
It has the displayed value at the original kappa=kappa0 and vanishes in
that formal limit with canonical classical inputs held fixed. This does
not prove interchange of a quantum limit or an all-order decoupling theorem.

### The exact norm and forward projection

Since k>=3, all five words lie in the D8 tail. Their absolute norm is

```text
N64[zeta W_k] = |zeta| 64^(2k) B_R5,
B_R5=21316224+152808*64^2+2463528*64^3
                  +18645*64^4+3259*64^5
    =4458582098560.
```

On the forward line `S=8+2v^2`, `U=0`, the previous exact identity gives

```text
R5=18900480-902256v^2-74580v^4,
L W_k=8^k(4725120k-902256)>0,                 k>=3.
```

The positivity follows already at k=3 and then increases in the linear
factor. Therefore choosing a formal comparison parameter
`zeta=+/-8lambda/[8^k(4725120k-902256)]` gives shifts `+/-8lambda`, with
the same finite coefficient prefix and the same four first-order rates.
The sample subtraction does not change this shift because W_k is zero
at each sample. This shows that the retained information does not bound
even the **weaker direct projection**, not only the particular N64 norm.

These controls deliberately violate the proposed D8 norm. They do not
show that a member already proved to satisfy that norm violates its
consequences. Nor are arbitrarily large finite-order polynomials assumed
to satisfy a physical high-energy growth bound, exact unitarity or a UV
completion. Those additional restrictions could remove this freedom;
they have not been supplied for the candidate under test.

The all-D statement follows from the displayed construction, not from a
finite numerical test. The diagnostic independently checks prefixes
`D=4,8,16,32,128`, both signs, all four sample equalities and their exact
forward coefficients.

## 3. Why the optional analytic witness does not supply its own bound

The previous double-Cauchy implication remains correct:

```text
A holomorphic near |S|<=128^2, |U|<=128^3,
sup|A|<=23000lambda
    => N64,D8 <=(29/336)*23000lambda <2000lambda.
```

But holomorphy is not a bound on `sup|A|`. The polynomial W_k is entire,
and its value at the boundary point `(S,U)=(128^2,0)` is nonzero:

```text
W_k(128^2,0)=128^(2k)
             [-18645*128^4-152808*128^2+21316224].
```

Scaling it makes that value, and hence the supremum, arbitrarily large
without changing the four rates or the finite prefix above. Knowing a
large heavy mass or the locations of known singularities cannot rule out
such a polynomial contribution. A dispersion reconstruction would need
the subtraction/growth or contour information that controls it.

Two qualifications prevent overinterpreting this result:

1. The R=128 witness is only sufficient. A large supremum at that radius
   does not prove that the defining radius-64 coefficient norm fails.
2. A tail norm above its proposed cap need not make its decision projection
   large. The diagnostic separately normalizes W_k to `N64=4000lambda`
   and obtains a projection below `lambda/1000` for its five test prefixes.
   The distinct `+/-8lambda` controls in section 2 establish the direct
   projection obstruction; this is not inferred from norm failure alone.

No such control is inserted into the physical candidate. They distinguish
missing logical implications from genuine bounds or model exclusions.

## 4. The pole/source bounds are independently unidentified

Let `P4 f=f-b4 M^-1[f(q_i)]`, with `b4=(1,S,U,S^2)`. For any one of the
declared complementary shapes Q, J and T, the formal deformation

```text
Delta F=delta_r P4 Q+delta_eta P4 J+delta_D P4 T
```

vanishes at all four centers. Its compensating local polynomial is
uniquely available from the same invertible M. In the declared
pole-plus-polynomial decomposition it adds **no E_tail term**. Thus a
tail bound does not supply the separate r, eta or D_N enclosure.

They remain distinct modulo the polynomial span. For example, at fixed
t=-2, s=z, u=6-z, their formal heavy Laurent coordinates are

```text
Res_(z=n) Q=-1,
coeff_(z-n)^(-2) J=1,          Res_(z=n) J=0.
```

At s=8, t approaching zero, u=-4-t, the T pole has residue 34. Local
polynomials cannot remove these Laurent coordinates. These are identities
of the specified first-order shapes, **not** a construction of an
unstable-particle physical pole or an admissible gravitational observable.

Even their conditional forward projections are not exactly null. Define
`K_Q=L P4 Q`, `K_J=L P4 J` and `K_T=L P4 T`, using the earlier conditional
functional with L(T)=0. Exact algebra gives

```text
K_Q<0, K_J<0 for all n>=32,
K_J=-d K_Q/dn,
lim_(n->infinity) n^6 K_Q=-2255640/3259,
lim_(n->infinity) n^7 K_J=-13533840/3259,
K_T=536987/2346480>0.
```

The half-line signs follow from shifted-polynomial coefficient positivity,
and independent Fraction solves check the sample cancellation and
projection at n=32,64,1024 and the original n. Small sensitivities at the
original n do not bound an independent unrestricted coefficient.

For any baseline coordinate, the two shifts `+/-2 cap` cannot both fit
an interval of radius cap, while the compensating local terms keep all
four rates fixed. This is an information requirement, not a recommendation
to choose either shift as a physical value. Physical resonance and
Newton/source matching, or a justified bound on them, is still needed.

## 5. Fixing these flat bounds still would not fix the curved parent

The already identified identity `H=(n+2)Q-3` yields the exact compensation

```text
delta_B H-(n+2)delta_B Q+3delta_B=0.
```

It leaves the complete on-shell four-scalar comparison identically zero,
not just the four sample values. In the retained convention
`delta_B=2g delta_c_RH/kappa`; a compensating finite heavy-residue term
and quartic constant can therefore keep the **combined** r fixed while
changing the individual curved c_RH. The D8 flat tail and combined-pole
bounds do not distinguish this direction. The covariant lifts differ.

S336 supplies another, radiative curvature example at order hbar whose
flat four-point effect starts at order hbar squared. Its chi is not the
known S347 loop coefficient. These observations identify the remaining
curved/source data; they do not prove that every such term affects a
homogeneous FLRW background. B needs the actual relevant causal/state
projection. The inherited S240 heavy Gaussian state is retained, not
replaced by a new state choice or mistaken for the complete interacting
light/metric state.

## 6. Could general amplitude bounds replace a microscopic input?

They can in principle constrain a comparison class more strongly than
locality and finite data alone. The polynomial controls above are not a
no-go for that research route. A focused primary-source check found no
automatic application supplying the numerical D8 witness:

- [Guerrieri and Sever, section II](https://arxiv.org/html/2106.10257v2)
  impose fixed-transfer twice-subtracted dispersion relations with massive
  thresholds and nonperturbative partial-wave unitarity. The current
  four-dimensional gravity amplitude does not inherit that analytic
  problem by deleting its graviton pole or using an inclusive rate.
- [Haering and Zhiboedov](https://arxiv.org/pdf/2202.08280), sections I-II
  and VII, make the semiclassical large-impact-parameter and other
  amplitude assumptions explicit. Their analysis leaves the properly
  defined four-dimensional infrared-finite observable outside its scope.
  It is not a ready-made R=128 supremum for this candidate.
- [Tokuda, Aoki and Hirano](https://arxiv.org/abs/2007.15009) retain
  Regge-dependent finite corrections in gravitational positivity.
  Their model-dependent UV data cannot be replaced by the smallness of
  a known low-energy loop.

This is an applicability assessment, not an exhaustive exclusion of all
possible bounds. A successful new theorem with its hypotheses verified
for the actual observable could supply the missing input. It must include
the needed absolute constants and subtractions; it cannot assume the
candidate's desired norm as part of the proof of that same norm.

## 7. Decision and intervention boundary

The **unchanged-source witness acquisition attempt stops with a precise
input obstruction**. More finite known-loop terms, a larger finite
Taylor prefix, or replaying the same diagnostics do not close it.
The broader witness objective remains unfulfilled; the conditional D8
candidate has not been excluded or physically matched.

Further progress requires new justified information, for example:

1. A concrete common quantum-parent matching prescription sufficient to
   determine the decision-relevant finite data and bound its remainder.
   A controlled finite-cutoff construction may be considered; a cutoff
   plus unspecified bare/finite coefficients is not sufficient.
2. Independent physical pole/source matching and a certified analytic or
   direct-projection bound, or a new amplitude-bound argument with all
   hypotheses and constants established for this four-dimensional theory.

Neither route requires an all-orders UV theory merely to state a finite
necessary test. Neither has been furnished by the current definition.
Any finite-cutoff parent must control the actual domain/contour used,
not assert asymptotic properties that its definition does not supply.

**Recommendation:** undertake a bounded feasibility study of a concrete
quantum matching prescription, retaining the original frame, targets and
known state data. Its first output should be a specified source of the
missing bounds or an explicit incompatibility/remaining input, before
another matching solver is built. This is a new parent-construction
direction; adopting new finite matching conditions or a new quantum
model requires an explicit decision. It is not a request for the user to
invent zeta, c_RH, chi or a convenient Cauchy supremum.

No new parent, cutoff, finite zero, normalization target or relaxed P8
finish line is adopted by this audit. The original project remains open.

## 8. Reproduction and evidence level

```sh
.venv/bin/python -B scripts/p8_rate4_d8_witness.py
.venv/bin/ruff check scripts/p8_rate4_d8_witness.py
.venv/bin/ruff format --check scripts/p8_rate4_d8_witness.py
```

The [read-only diagnostic](../scripts/p8_rate4_d8_witness.py) protects
311 inherited input files. It checks fifteen forward-coefficient
identities/controls, twenty tail sample nulls, five analytic-domain-without-
supremum controls, five large-norm/small-projection controls, 48 nuisance
sample nulls, twelve independent
Fraction/SymPy projections, two half-line sign proofs, six cap-shift
controls, four Laurent identities, the exact curved/residue compensation
and ten rejected invalid inputs. The arbitrary-prefix statement is a
written all-k proof, not an extrapolation from the five finite tests.

Evidence level: **VERIFIED_N** for the scoped source/identifiability audit,
not a new physical matching certificate. See the
[validation receipt](validation/p8-rate4-d8-witness-2026-09-22.json).
No full 115,487-test publication replay or new frozen packet is involved.
