# RATE4.PARENT: compare complete input routes before adopting a new candidate

Date: 2026-09-21 (PDT). This is the user-requested continuation of
[RATE4.REMAINDER](assessment-2026-09-21-p8-rate4-remainder.md).
**Outcome: input designs compared; no complete parent supplied or adopted.**
Physical M/V/G/B/R and original P8 remain OPEN. The approved RATE4 v1
conditions and all frozen calculations are unchanged.

**Recommendation:** start with the minimum-change **RATE4-D8, parent-bounded**
proposal: keep the four approved targets, fit the complete local polynomial
through eight derivatives, and require independent bounds on the pole/source
directions and the entire omitted analytic tail. It has an explicit
sufficient error allocation below `lambda/2`. This is a new bounded-EFT
restriction, not a conclusion about the unrestricted RATE4 family.

The alternative **RATE5-D10** adds one off-angle target and fits through
ten derivatives. It determines the first previously invisible contact
direction instead of putting its highest-degree coefficient in the tail
bound. Both designs fit the same stated error allocation; the fifth rate
is useful if that coefficient should remain unrestricted by the tail
hypothesis, but is not mandatory once the stronger D8 bound is supplied.
Neither proposal or its missing bounds is adopted here.

Eight low-energy rates were also tested. They do determine the enlarged
eight-shape first-order ansatz, with a stable target-coefficient projector.
However, at sample precision `lambda/100` their individual heavy-mass
coordinate has an independent-error-box radius between `10^1168` and
`10^1169`. This cannot supply a controlled heavy-sector or curved parent
from that precision. It is not a contradiction in exact eight-rate data.

The distinction is important: this comparison produces a usable input
contract, not the physical data, a full hard-loop calculation, a
massless-gravity dispersion theorem or a quantum bounce.

## 1. Routes compared and the common-parent requirement

| Route | What it fixes | Assessment |
|---|---|---|
| More RATE4 loops, same finite data | Known nonlocal contributions | Does not constrain the R5 contact or independent curved data. Not a matching-resolution route. |
| Eight low-energy inclusive rates | Eight flat shapes including the heavy-mass insertion | Algebraically identifiable, but individual heavy coordinates are extraordinarily ill-conditioned. Requires additional control for the proposed parent. |
| Four retained rates plus parent bounds (D8) | Four local coefficients; all terms starting with SU belong in the bounded tail | Recommended minimum-change design, conditional on a stronger tail input. No new rate target; not yet physical matching. |
| Five rates plus parent bounds (D10) | Five local coefficients; SU is fitted and the bounded tail starts later | Explicit off-angle alternative if the stronger D8 tail condition is not appropriate. Not yet physical matching. |
| A specified microscopic or finite-cutoff quantum parent | Potentially all matching and state-dependent data | Would supply the strongest common input, but no matching parent with this spectrum, frame and bounce is currently supplied. A regulator and unspecified bare coefficients are not that construction. |
| Minimal-subtraction zeros or coefficient naturalness alone | A reference representative or a prior | Neither is a certified physical error bound. Not adopted. |

A complete parent packet must identify a **single action and finite
counterfunctional in the original physical matter frame**, its matching
reference, spectrum and domain, and the state/causal response used for B.
The flat and curved calculations cannot independently choose those data.
This does not add a requirement to construct an all-orders UV theory.

The local/nonlocal distinction follows the usual EFT separation explained
by [Donoghue and Holstein, section II](https://arxiv.org/html/1506.00946):
calculable nonlocal effects do not determine unknown local matching
coefficients. A naturalness-based statistical truncation interval would
be a different evidence level: [Furnstahl et al.](https://arxiv.org/abs/1506.01343)
explicitly require coefficient priors. No such interval is promoted here
to a deterministic bound.

## 2. Why eight rates do not automatically reconstruct a healthy parent

In mass-one units take

```text
S=s^2+t^2+u^2, U=stu,
H=sum_a(a+2)/(n-a), J=sum_a 1/(n-a)^2,
T=sum_cyclic(2-2a-bc)/a.
```

The full eight-shape ansatz considered is
`(1,S,U,S^2,SU,T,H,J)`. Its analytic/pole-subtracted target row is

```text
l8=(0,2,0,32,0,0,2(n+2)/(n-2)^3,6/(n-2)^4).
```

Besides the original four points, use
`(10,-1,-5), (12,-1,-7), (12,-2,-6), (16,-1,-11)`.
All lie in the retained physical nonforward window. The evaluation
matrix is nonsingular for `n>=32`. The diagnostic proves this by an exact
rank-six constant system and a two-by-two heavy Schur complement, not
floating-point inversion of an almost degenerate matrix.

For `n>=1024`, its coefficient weights `w8=l8 M8^-1` have
`sum|w8_i|<1`. Thus the target projection itself is stable even at the
original enormous n. That does not control each recovered coefficient.
At the original parameters, an independent absolute sample-error box
of radius `lambda/100` gives the exact coordinate-radius bounds

```text
10^774  < radius(B_H) < 10^775,
10^1168 < radius(eta_J) < 10^1169.
```

The radii are `epsilon*sum|row(M8^-1)|`; opposite corners attain them.
If `eta_J=-g^2 delta_n`, ensuring the **uncertainty radius** alone fits
inside the cap `g^2*n/4` would require relative sample precision
`epsilon/lambda` between `10^-982` and `10^-981`. The inferred central
value would still need a separate check. These numbers are conditioning
bounds, not actual enormous physical coefficients or uncertainty in exact
stipulated targets. Correlated analytic input can be better than an
independent box, but such correlations must be established.

This is why the eight-rate construction is not recommended as a
low-energy-only determination of the full parent. It does not rule out
eight exact conditions with suitable parent control. No heavy pole is
being treated as an external stable particle: the retained heavy field
can decay, and genuine resonance matching needs its decay/line-shape
prescription. See the retained S233/S234 qualifications and
[Beneke's unstable-particle EFT review](https://arxiv.org/abs/1501.07370).
The review supplies context, not the missing resonance data for this model.

## 3. Two bounded local designs and their exact conditioning

### Four retained rates, complete eight-derivative polynomial

With the pole/source directions treated as independent bounded inputs,
fit `(1,S,U,S^2)` using the original four RATE4 centers and targets.
These are **not** the old `(S,H,1,T)` interpolation coordinates. The
new bounds on the remaining directions are essential; this is not a
relabeling of the already underdetermined unrestricted RATE4 family.

For `l4=(0,2,0,32)`, the exact projector is

```text
w4=(-3579/13036,5084/9777,-3509/13036,232/9777),
sum|w4_i|=3544/3259.
```

The matrix determinant is `10011648`, and the maximum coefficient error
amplification is `85231/3259<27`, independent of n. The SU coefficient is
not fitted: it must be included in the parent-provided tail norm below.
In particular, the R5 contact remains invisible if its coefficient is
unbounded. The new tail condition bounds it; it does not claim the four
rates alone have suddenly determined it.

### Five rates, complete ten-derivative polynomial

Keep the original four centers and propose adding

```text
q5=(10,-2,-4),       x_det=1/256.
```

This has a different angle from the four `t=u` centers. Its full unexpanded
Born amplitude obeys `232<a5/lambda<233`; all five amplitudes are positive
and below `650lambda`. The proposed fifth target extends the existing
matter-calibrated rule `R5_target=1+2 Re L_m(q5)/a5`, with the same complete
S239 reference and the same detector and regulator order. If adopted,
its formal higher target coefficients would be zero as in RATE4 v1;
only first-order matching is analyzed here. It is not measured data or
an old-parent prediction. **No fifth condition is adopted by this document.**

Separate the complete non-matter hard correction as

```text
F(q)=K(q)+p5(q)+r Q(q)+eta J(q)+D_N T(q)+E_tail(q),
Q=sum_a 1/(n-a),       p5=c0+cS S+cU U+cSS S^2+cSU S U.
```

`K` is the whole specified calculable hard reference, with S293/S294/S302
and only the S305 increment counted once. The complete S239 matter loop
remains separate. Here `r` is the **combined** heavy-residue direction,
including any equivalent H direction using `H=(n+2)Q-3`; the absorbed
constant belongs in p5. This does not identify the individual curved
`R H_heavy` coupling. The five local coefficients are solved for fixed
remaining inputs, which are bounded rather than set to zero.

Let `M5` evaluate `(1,S,U,S^2,SU)` at the five points and
`l5=(0,2,0,32,0)`. The exact weights are

```text
w5 = (-47677/221072, 15256/41451, -18538/96719,
       9941/663216, 18797/773752),
sum|w5_i| = 630347/773752 < 1.
```

Consequently `L p5=w5 p5(q_i)`. The complete local R5 ambiguity from the
remainder audit is now in the fitted span and is no longer invisible:
its value at q5 is `-37140096`. The maximum absolute row sum of `M5^-1`
is `1906116/96719<20`, independent of n. This bounds the error amplification
of **every** fitted local coefficient in these stated mass-one units,
not just their target projection. It is not a bound on their central
values without complete physical residuals.

The corresponding subtraction is `P5 f=f-b5 M5^-1 f(q_i)`. Under a
finite reference shift in the five polynomial directions the coordinates
shift oppositely and the total hard amplitude is unchanged. Changes in
the tail or pole/source reference must also transform their input bounds;
one cannot keep a bound numerically fixed while changing what it bounds.

## 4. Explicit sufficient bounds on the unfitted directions

The following are proposed **parent-input acceptance conditions** in this
same reference, not consequences of the existing normalization data:

```text
|r| <= g^2/4,    |eta| <= g^2 n/4,    |D_N| <= 1/(4kappa).
```

In the linear insertion convention the second is equivalent to
`|delta_n|<=n/4`, while the last corresponds to `|delta_kappa|<=kappa/4`.
A physical parent must justify the mapping from its pole/source data to
these residual coordinates, including the already retained loop terms.
They are not independently measured parameters already available here.
They do not imply that all higher orders are small. Other enclosures can
be used; these are explicit sufficient, not necessary, choices.

Why are these relatively broad bounds enough for the **flat projection**?
The expansion of Q through channel degree five is in the p5 span:

```text
Q=sum_{k=0}^5 (sum_a a^k)/n^(k+1)
       +sum_a a^6/[n^6(n-a)].
```

The power sums through degree five are polynomials in `(1,S,U,S^2,SU)`
when `s+t+u=4`, so P5 kills that entire truncation. For `|a|<=16`, n>=32,
the complete Q and J remainder bounds are

```text
B_Q=3*16^6/[n^6(n-16)],
B_J=3/n^2 * (16/n)^6 * [7-6(16/n)]/[1-16/n]^2.
```

They also bound the analytic forward circle `|v|=1` since its channels
have modulus at most 3. Cauchy's estimate and the sample weights give

```text
|L P5(r Q+eta J)| <= (1+sum|w5|)(|r| B_Q+|eta| B_J).
```

For the Newton direction, conditional on the specified full pole
prescription with `L T=0`, the exact projection is
`L P5 T=26901281/79585920`. At the original parameters the three displayed
parent caps therefore imply a combined nuisance error below
`10^-200*lambda`. This is a first-order finite-direction bound. It is not
a resummed resonance theorem, a Newton measurement, or an admissible full
massless-gravity functional.

For the D8 design, the same proof truncates Q at channel degree four:
replace `16^6/n^6` by `16^5/n^5` in B_Q, and replace
`(16/n)^6[7-6(16/n)]` by `(16/n)^5[6-5(16/n)]` in B_J. Use w4 instead
of w5. Its Newton projection is `536987/2346480`. The original-parameter
combined nuisance bound remains below `10^-200*lambda` under the same
caps. No claim that a bounded first-order insertion controls all higher
insertions is needed or made.

## 5. A controlled tail without declaring higher operators absent

For D10, require a bound on the residual analytic matching tail, in the fixed
common reference after all known light cuts, pole/source directions and
the low-degree polynomial have been separated:

```text
E_tail(S,U)=sum_{2i+3j>=6} c_ij S^i U^j,
N64=sum_{2i+3j>=6} |c_ij| 64^(2i+3j) <= B.
```

This is an absolute summability/analyticity input about the **remaining
local matching data**, not about the full scattering amplitude. Its scale
64 is a chosen norm radius in mass-one Mandelstam units, not a derived
physical cutoff or the heavy mass. The known light/graviton nonanalytic
functions are retained explicitly in K; they cannot satisfy this norm by
being silently Taylor-expanded through their cuts. EFT power counting,
an apparent gap, or a few small terms do not themselves prove N64 finite
or bound B. A bound proved on a different domain could be converted with
its actual constants instead.

Let `x_i=S(q_i)/64^2` and `y_i=U(q_i)/64^3`. They lie in [0,1).
Each omitted monomial has normalized sample magnitude at most
`max(x_i^3,x_i^2 y_i,y_i^2)`: split into `i>=3`, `i=2,j>=1`, or `j>=2`.
On the forward line U=0 and S=8+2v^2. The normalized v^2 coefficient
of S^i, i>=3, is largest at i=3, with value `384/64^6`; successive ratios
are at most `1/384`. Absolute summability permits termwise extraction.
It follows, for **all** omitted terms, that

```text
|L P5 E_tail| <= C_tail B,
C_tail=384/64^6 + sum_i |w5_i| max(x_i^3,x_i^2 y_i,y_i^2)
      =217375713/6490702217216 < 1/20000.
```

For example, the **additional condition** `B<=2000lambda` would allocate
less than `lambda/10` to this entire analytic tail. It does not set any
of its coefficients to zero or prescribe their signs. No such bound has
been derived from the approved RATE4 data. The general B-dependent formula,
not this optional sufficient allocation, is the result of the calculation.

For D8, replace `2i+3j>=6` by `2i+3j>=5`; the additional omitted word is
SU. The sample majorant becomes `max(x_i^3,x_i y_i,y_i^2)`. The forward
majorant is unchanged because SU vanishes on the forward line. Using the
four weights gives the complete bound

```text
|L P4 E_tail| <= C_tail,D8 N64,D8,
C_tail,D8=332221287/6998649208832 < 1/20000.
```

Thus `N64,D8<=2000lambda` also gives a tail error below `lambda/10`.
Despite the same numerical cap, it is a **stronger input**: it additionally
restricts the SU coefficient. The unbounded R5 comparison from the previous
audit is not invalidated; arbitrary scaling of it violates this new norm
condition. Neither N64 bound has been established for RATE4 v1. If a parent
cannot support the D8 tail restriction, the off-angle fifth condition can
instead fit SU, while all later terms still require the D10 tail bound.

## 6. Complete flat error contract and a sufficient allocation

For either proposal, with its corresponding weights and tail norm,
errors in the subtracted projection obey

```text
epsilon_total <= sum_i |w_d,i| epsilon_sample,i
    + epsilon_known,L + epsilon_higher,L
    + C_tail B + epsilon_pole/source.
```

Each sample error must include its physical rate and target uncertainty,
complete hard-reference conversion, relevant real radiation, higher-order
rate terms, bin/regulator errors if used, and known-hard subtraction error.
The separate higher-L term concerns omitted orders in the analytic test
functional, not just an amplitude square or its values at the samples.
If the exact defining rates are stipulated, their measurement error is
zero, but the calculation/conversion/omitted-order errors are not thereby
zero. Small real interference does not bound an imaginary part or a square.

One sufficient **conditional allocation**, using `a_i<650lambda`, is:

| Required input | Sufficient bound |
|---|---|
| Combined normalized-rate/conversion/omitted-rate uncertainty at each center | `<=1/10000` |
| Complete known-hard sample subtraction error at each center | `<=lambda/100` |
| Complete known-hard L evaluation error | `<=lambda/10` |
| All omitted orders in L | `<=lambda/4` |
| Analytic tail norm N64 | `<=2000lambda` |
| Pole/source inputs | The three stated caps, implying much less than `lambda/1000` |

Then `epsilon_sample <=(13/400+1/100)lambda`. For D10, using `sum|w5|<1`,

```text
epsilon_total < (13/400+1/100+1/10+1/4+1/10+1/1000)lambda
              = (987/2000)lambda < lambda/2.
```

This fits the existing absolute lambda target with margin. It does not
supply any missing physical input in that table, give the central matched
coefficient, or prove its sign. The known retained first-order radiation
bound remains available at q5 because it is in the same window; it is
not a bound on the higher-order entries. The full V/G pole/cut/infrared
prescription and contour allowance must still be justified independently.

For the four-rate D8 alternative, the sample norm is slightly above one.
With the **same** allocation table, now including SU in N64, the bound is

```text
epsilon_total < [(3544/3259)(17/400)+451/1000]lambda
              = (1620429/3259000)lambda < lambda/2.
```

It therefore meets the target without adding a fifth rate, conditional on
that stronger parent input. Neither bound is a physically established
error bar: the missing inputs are precisely the ones listed in the table.

## 7. What a genuinely complete same-parent packet must add

The recommendation does not turn four or five flat numbers into a curved action.
The complete input packet would contain these distinct items:

1. **Identity and domain:** one action, spectrum, physical matter frame,
   finite counterfunctional, subtraction convention, and controlled domain.
2. **Flat inputs:** the retained four targets (plus a fifth only if adopted),
   the pole/source residual enclosures, tail bound, and complete physical
   errors in the preceding contract. Parent evidence must support the
   bounds; the arithmetic interface cannot authenticate their provenance.
3. **Covariant matching:** the finite lift of the fitted polynomial and
   independent curved/source directions with a domain bound. The flat
   heavy-residue combination does not determine c_RH separately. S336's
   real radiative contact is an explicit remaining direction, not exhausted
   by four, five or eight scalar rates or by the known S347 loop coefficient.
   A matched physical radiative/source observable could constrain it, but
   one such datum is not automatically a bound on every curved operator.
4. **Causal/state control for B:** the prepared quantum state, canonical
   map including the retained boundary transport, local and nonlocal
   response/stress errors, and the mode/volume/domain limits actually used.
   In-out scattering data are not substituted for a state-dependent
   causal response. The corrected finite S275 hybrid is not this missing
   global, regulator-controlled construction.
5. **Admissible V/G test:** the physical observable, specified subtractions,
   and actual contour/infrared error needed by the chosen necessary test.
   A finite inclusive rate is not itself a crossing-analytic amplitude.
   If using V, establish its actual quantum decoupling limit with the
   canonical inputs held fixed. The proposed finite-coordinate caps do
   not imply that every additional interaction vanishes in that limit.

The low-degree contact part has explicit covariant representatives using
the same labelled-slot construction and V(X) as in the remainder audit.
That demonstrates the local algebra is realizable as a formal finite EFT
counterfunctional. Its zero clock jets do not bound an arbitrary nearby
tube, fix the curvature-only kernel, or construct the required state.
No finite coefficient, tail, curved kernel or higher-order error has been
chosen zero in this comparison.

## 8. Decision, approval boundary and reproducibility

The comparison favors **retaining four rates with explicit D8 parent
bounds first**, since that input design already meets the error target.
It does not make the unrestricted RATE4 family MATCHED. If the stronger
tail restriction on SU is inappropriate, RATE5-D10 supplies an off-angle
condition for that coefficient while retaining a bounded later tail.
Neither design uses low-energy data to infer unbounded heavy coordinates.

The next substantive choice is whether to develop RATE4-D8 as a separately
named bounded-EFT proposal, use its five-rate alternative, or seek a
concrete microscopic parent. Adopting either set of bounding assumptions
(and the fifth target if applicable) is a model choice. User approval
would authorize that choice, not prove its physical realization or P8
closure. The complete source and error obligations above remain visible.

The [read-only diagnostic](../scripts/p8_rate4_parent.py) runs with:

```sh
.venv/bin/python -B scripts/p8_rate4_parent.py
.venv/bin/ruff check scripts/p8_rate4_parent.py
.venv/bin/ruff format --check scripts/p8_rate4_parent.py
```

It pins 309 inherited files before/after execution. Original SymPy and
independent Fraction solves agree on 33 projection weights; eighteen synthetic
fitted coordinates, exact half-line determinant/sign/norm proofs, eleven
power-sum identities, the infinite-tail generating identities, 329 omitted
monomial checks, 48 sample-error corners and fourteen invalid/incomplete
packet rejections pass. The infinite-tail and Cauchy arguments are written
proofs above, not conclusions inferred solely from finitely many samples.
The preliminary point search was exploratory; no globally optimal design
claim is made. All final numbers use exact arithmetic at the stated points.

Evidence level: **VERIFIED_N** for this comparison and conditional input
map, not a new frozen physical matching certificate. See the
[validation receipt](validation/p8-rate4-parent-2026-09-21.json).
No full 115,487-test regression was rerun for these root diagnostics.
