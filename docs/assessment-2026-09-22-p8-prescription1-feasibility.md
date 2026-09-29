# PRESCRIPTION-1: cutoff, projected matching and common-parent feasibility

Date: 2026-09-22. The user authorized the feasibility study of all three
avenues recommended after the
[RATE4-D8.WITNESS audit](assessment-2026-09-22-p8-rate4-d8-witness.md).
This authorizes investigation, not adoption of new physical matching
conditions, a cutoff, a state or a replacement parent.

**Outcome: the feasibility screen is complete; no complete common quantum
parent is established.** Two explicit identifications of a finite regulator
with the original physical local field fail necessary tests. There are also
positive results: a controlled direct projection of a selected UV shell,
a one-moment bound for a restricted positive-spectral matching class, and
a quantitative high-momentum readout bound for the actual retained heavy
state. None supplies all the missing finite-gravity matching data.

The work changes neither the original RATE4 targets nor the D8 hypotheses.
No original P8 gate is promoted. M/V/G/B/R and P8 remain **OPEN**; completed
scoped P8(a) and the linear classification are preserved.

## 1. Decisions, rather than another conditional parent definition

| Avenue examined | Result | Consequence |
|---|---|---|
| Finite proper-time covariance treated as an exact physical local scalar propagator | **EXCLUDED ANSATZ** already in its free, gravity-decoupled seed | Do not adopt this as the quantum parent; it remains usable as a computational regulator |
| Sharp spatial cutoff field identified with the original exact local field | **Fails that locality identification** | A coarse-grained Hamiltonian description is possible, but needs matching to the local theory |
| Wilsonian regulator with a separately specified physical boundary action | **Not excluded; incomplete input** | Regulator choice and RG running do not determine the physical boundary coefficients |
| Direct projection of the selected light-bubble UV shell | **Bound obtained** | Demonstrates a controlled matching mechanism, not a bound on every sector or the unknown boundary action |
| Positive scalar-exchange spectral model | **Restricted bound mechanism obtained** | A physical spectral moment can replace a much stronger coefficient norm in this class; generic gravity/curved matching is not in the class |
| Finite-gravity positivity observable | **Still not admitted for a verdict** | Neither a matching projector nor a parametric finite-detector error is an admissible complete gravity test |
| Original heavy-state high-momentum readout | **Bound obtained in the existing reference domain** | No state reset is needed for this estimate; interacting light/metric and nonlinear/global B obligations remain |

The negative results below concern explicitly tested regulator ansatzes,
not Wilsonian EFT in general, all possible quantizations or the existence
of an allowed P8 parent. A complete all-orders UV theory is not added to
the project's finish line.

## 2. Test a concrete covariant cutoff before promoting it to physics

### 2.1 Finite proper-time free covariance

The first trial identifies the gravity-decoupled physical scalar's exact
Euclidean covariance with

```text
C_Lambda(x)=integral_(Lambda^-2)^infinity exp[-tau(x+m^2)] d tau
           =exp[-(x+m^2)/Lambda^2]/(x+m^2),  x=p_E^2.
```

Here Lambda is finite and m>0. This test is upstream of any choice of
finite interaction counterterms. It is not a claim about a gauge-dependent
graviton propagator.

For an ordinary positive-spectral scalar covariance,

```text
C(x)=integral d rho(z)/(x+z),   d rho>=0,
d[x C(x)]/dx=integral z d rho(z)/(x+z)^2 >=0.
```

Put L=Lambda^2 in the trial. Direct differentiation instead gives

```text
d[x C_Lambda(x)]/dx at x=L
    =-exp[-(L+m^2)/L] L/(L+m^2)^2 <0.
```

Thus this covariance does not have the required positive spectral
representation. The ordinary free pole remains, but the additional entire
term is not a harmless local contact term.

An independent position-space argument avoids relying on unsubtracted
spectral notation. At fixed spatial momentum its time kernel is

```text
c(t)=integral_(Lambda^-2)^infinity
       exp[-tau omega^2-t^2/(4tau)] d tau/sqrt(4pi tau).
```

It is smooth and even, with c(0)>0, c'(0)=0 and c''(0)<0. Reflection
positivity would require the two-time matrix with entries c(t_i+t_j)
to be positive semidefinite. For t_1=t, t_2=2t,

```text
c(2t)c(4t)-c(3t)^2=c(0)c''(0)t^2+O(t^4)<0
```

for sufficiently small positive t. Smooth positive-time test functions
approximating the two time probes, and spatial wave packets near the
chosen momentum, retain the strict negative direction by continuity.
This excludes the proposed exact local-field Gaussian covariance.

These are the study's direct checks. As an independent literature check,
[Handrack and Salmhofer, appendix B](https://arxiv.org/html/2608.10685v1#A2),
also treat the **scalar** exponential-regulator covariance and construct
a reflection-positivity violation. Their factor omits the positive constant
exp(-m^2/Lambda^2), which cannot change the sign. We do not apply a theorem
about rational propagators to an exponential one.

This does **not** forbid using proper time to calculate a regulated
effective action. Regulated intermediate fields need not themselves be
physical local observables. What fails is promoting this particular
finite-regulator covariance, without its missing physical completion,
to the exact parent assumed by the S6 contract.

### 2.2 Spatial Hamiltonian cutoff

A spatial momentum cutoff can retain a positive Hilbert space for its
retained modes, but it does not leave the original local field intact.
Its equal-time canonical commutator uses

```text
delta_Lambda(r)=[sin(Lambda r)-Lambda r cos(Lambda r)]/(2pi^2 r^3)
```

instead of a spatial delta distribution. In particular
`delta_Lambda(pi/Lambda)=Lambda^3/(2pi^4)>0` at nonzero separation.
The time derivative at t=0 of its free scalar commutator is proportional
to this kernel. Consequently it is nonzero for some arbitrarily small
positive times at that fixed spacelike separation. Identifying this
band-limited field with the original exact local field fails microcausality.

This is not a no-go for coarse-grained observables or Hamiltonian
regularization. It means that locality, boosts, Ward identities and the
physical observable require a matching/reconstruction argument. The
existing finite-band bounce estimates do not provide that argument.

### 2.3 What a genuine Wilsonian prescription still needs

A calculational cutoff and a physical Wilsonian boundary action are
different data. At a fixed loop order write schematically

```text
F_phys=F_loops,Lambda+A_Lambda,
D F_phys=D F_loops,Lambda+D A_Lambda.
```

Here D is the matched projection defined below, not a positive dispersion
functional. Changing cutoff or subtraction reference requires compensating
changes in A_Lambda. Small regulator dependence of the first term does
not bound the second term's integration constant. For example
`D F_loops,Lambda=b/Lambda^10` and
`D A_Lambda=C-b/Lambda^10` give the cutoff-independent result C for
any supplied C. The diagnostic checks this identity for three distinct
formal C values; it does not assert their physical realizability.

Setting every remaining renormalized coupling to zero at a fixed scale
would be **new physical boundary data**, not merely a subtraction convention.
It can be proposed as a separately defined quantum EFT, but must specify
the complete relevant covariant functional and its frame, normalization
and state prescription. Zero values in a few flat amplitudes do not specify
that functional. No such physical boundary condition is adopted here.
The distinction agrees with the local/nonlocal EFT discussion in
[Donoghue and Holstein, section II](https://arxiv.org/html/1506.00946#S2).

## 3. Direct matching: a worked UV-shell bound

### 3.1 Projection and actual source of the example

Keep mass-one units, the four RATE4 centers, `b4=(1,S,U,S^2)`, and

```text
D f = L f - sum_i w_i f(q_i),
L f = coefficient of v^2 in f(2+v,0,2-v),
w=(-3579/13036,5084/9777,-3509/13036,232/9777),
A=1+sum|w_i|=6803/3259.
```

The retained amplitude normalization is lambda=10^-600; no rate target
or physical parameter is changed.

D annihilates every symmetric channel polynomial through degree four.
This is the projected effect of refitting the four local coordinates,
not a statement that the unfitted amplitude itself is small.

The [S235 graph formula](../problems/P8/s6/continuation/s6_235/notes/diagrams.md)
contains the light-bubble term `C_tree^2 sum B(a)/(32pi^2)`.
[S239](../problems/P8/s6/continuation/s6_239/FORMULATION.md) and
[S297](../problems/P8/s6/continuation/s6_297/notes/source.md) establish the
specified inheritance of that sector in the retained first-order matter
calculation. This is **only the C_tree^2 term**, not the entire bubble,
triangle/box sum, gravity correction or physical matching remainder.

For a proper-time determinant calculation, its omitted short-time piece,
after subtracting a channel-independent constant, is the entire function

```text
H_Lambda(a)=integral_0^1 dx integral_0^(Lambda^-2) d tau/tau
            exp(-m^2 tau)[exp(a x(1-x)tau)-1].
```

Derive the separation in the Euclidean/subthreshold domain and continue
the long-time amplitude on its original Feynman branch. H_Lambda is
entire and adds no cut; an above-threshold long-time integral is not
treated as an ordinary convergent positive integral.

This is a **total proper-time determinant prescription** for the selected
loop. It is not asserted to equal a bubble formed by separately cutting
the proper time of each propagator. Nor is it a single completed quantum
construction with the rejected physical covariance of section 2.

### 3.2 Infinite remainder controlled before any numerical truncation

For |a|<=16 and Lambda>=2, subtract the first four powers of a. Using
`x(1-x)<=1/4` and the exponential remainder gives

```text
|H_Lambda(a)-Taylor_4 H_Lambda(a)|
 <= exp(4/Lambda^2) 4^5/(600 Lambda^10).
```

The forward circle |v|=1 has channels of modulus at most three and is
inside the same bound. Cauchy's coefficient estimate, the three channels
and the four sample weights therefore give

```text
|D [beta sum_a H_Lambda(a)]|
 <= |beta| 3 A exp(4/Lambda^2) 4^5/(600 Lambda^10)
 < (2612352/81475)|beta|/Lambda^10
 < 33|beta|/Lambda^10.
```

The physical prefactor for this selected example is
`0<beta=C_tree^2/(32pi^2)<1`, verified from the original parameters.
For a zero prefactor both sides vanish; the strict comparisons apply to
the nonzero case.
At the **test scale**, not adopted physical cutoff, Lambda=10^100, the
bound is below `10^-398 lambda`. The same proof applies to any
nonnegative loop mass squared. This proves that this particular cutoff
contribution is controllable; it does not make the independent physical
boundary matching small.

There is also an exact sign check. At m^2=1,

```text
D sum_a H_Lambda(a)
 =-18797/(45169740 Lambda^10)+O(Lambda^-12).
```

The coefficient is derived from `D sum a^5=-2255640/3259` and the exact
parameter integrals. To bound the remainder after proper-time order N,
expand both exponentials in
`exp[(-m^2+a x(1-x))tau]-exp(-m^2 tau)`. If Lambda^2>=m^2+4, their
combined exponential remainder, the tau integral and the same projector
estimate give

```text
E_N <=18 A [(m^2+4)/Lambda^2]^(N+1)/[(N+1)(N+1)!].
```

The stated sign is not inferred from big-O: at m^2=1, N=5, the
order-six remainder bound at Lambda=1000 is only
`1571493/4812032` of the leading magnitude and this ratio decreases as
Lambda^-2. Thus the omitted short-time projection is strictly negative
for positive beta and every Lambda>=1000. The cutoff-minus-full
amplitude difference has the **opposite** sign in this convention.

Sixteen independently integrated series coefficients and four finite
series enclosures check the algebra. The infinite-remainder argument
above, rather than those finite checks, is what controls the tail.

## 4. A possible source of direct bounds: positive scalar spectra

To investigate an actual matching mechanism rather than assume a general
Taylor norm, restrict provisionally to scalar-exchange amplitudes

```text
F_spec(s,t,u)=integral_(z>=32) d rho(z) sum_a 1/(z-a),
d rho>=0,       M6=integral d rho(z)/z^6 < infinity.
```

For an unbounded continuous spectrum, understand this expression with
five channel Taylor subtractions and a separately specified polynomial
in `span(1,S,U,S^2)`. M6 alone need not make the unsubtracted integral
converge. D annihilates that polynomial, and the subtracted integral
converges uniformly on the sample/forward domain by the bound below.
Finite scalar sectors need no such convergence qualification.

This is a restrictive model class, **not** a spectral representation
theorem for the unknown P8 four-point remainder. It describes, for
example, leading exchange of explicit positive-kinetic heavy scalars
with residues h_j^2, for which M6=sum h_j^2/z_j^6 is fixed by their action.
A heavy state above the two-light threshold is understood in its
perturbative exchange expansion; no exact stable-heavy-particle
asymptotic state is inferred.

The exact identity

```text
1/(z-a)=sum_(k=0)^4 a^k/z^(k+1)+a^5/[z^5(z-a)]
```

and `|a|<=16`, `z>=32` imply, after D kills the first sum,

```text
|D F_spec| <= 6 A 16^5 M6 = (42800775168/3259) M6.
```

Thus one physically supplied moment can suffice:
`M6<=lambda/[24 A 16^5]` allocates at most lambda/4 to this sector.
This is weaker information than a norm of every analytic coefficient,
but it still must be supplied by a model, not chosen from the desired
answer. Other-spin contributions, local contact data and gravitational
counterterms are not bounded by this formula.

For a finite scalar sector, its elementary stability check is constructive.
With positive masses squared z_j and mu, take

```text
V=mu phi^2/2
  +sum_j z_j[H_j-h_j phi^2/(2z_j)]^2/2
  +lambda_stab phi^4/24,   lambda_stab>=0.
```

Its standard local expansion gives the exchanges and a definite quartic
contact. The diagnostic checks the square completion and positive-spectral
fixtures, including a fixture below the lambda/4 sector allocation.
Those fixture couplings are **not** P8 inputs or a proposed retuning.
Reimposing the four RATE4 targets would require further matching; its
effects cannot be ignored. No all-energy UV completion or common bounce
follows from this elementary matter-sector example.

[Alviani, Falkowski and Marinellis](https://arxiv.org/html/2507.11426v1)
give related examples of explicit heavy-field matching, but their
massless-scalar spectrum and curvature operators differ from ours.
Their discussion also retains unresolved UV-sensitive contributions
in a vector-matter example. We use this as a methodological comparison,
not an import of coefficients or a completion of the H8A420 model.

**Decision:** this is an admissible source-of-bounds mechanism to consider
in a genuinely specified parent. It is not enough to add another small
heavy sector while leaving the original independent local/gravity
counterfunctional unchanged. Such an addition would not resolve M.

## 5. Gravity compatibility: the projector is not the positivity theorem

The fitted projection D is a reconstruction device. In fact, for a
positive-residue scalar exchange with mass squared z>=32,

```text
L Q_z=2/(z-2)^3>0,           D Q_z=K_Q(z)<0.
```

The half-line sign of K_Q follows from the earlier exact rational
identity and is replayed here. Therefore imposing `D F>=0` as a
standalone unitarity test would already reject an ordinary positive
scalar-exchange contribution. A small direct matching error must instead
be combined with the full reference amplitude and an independently
admissible functional. Neither D nor the four physical rates supplies
that functional automatically.

The primary-source check of
[Bellazzini et al., sections 4-5](https://arxiv.org/html/2512.13780v2#S5)
retains both the all-spin/all-energy kernel condition and a finite-coupling
remainder outside its detector scaling limit. Its fixed-resolution result
does not assign a numerical upper bound to that remainder for this massive
parent. This agrees with the existing
[S278 applicability audit](../problems/P8/s6/continuation/s6_278/notes/dispersion.md).
We neither change the detector resolution 1/256 nor import a massless
example's constants. The original fixed-negative-transfer Regge route
also still needs its actual residue/trajectory or finite-contour input.

For a candidate to pass the next gravity admission check it must provide:

- the same physical massive-channel observable and its regulator-removal
  order, including a quantitative conversion from the stipulated rates;
- a justified positive spectral/partial-wave functional or an appropriate
  finite-contour replacement, with the required boundary information;
- an absolute gravitational and omitted-order error in the units of the
  final decision coefficient, not an unevaluated O(G) term.

No all-spin functional or gravitational allowance is proved by this
feasibility study. G remains UNTESTED and V remains OPEN.

## 6. Bounce compatibility: an actual same-state tail bound

### 6.1 Keep the prescribed state and all local terms

[S240](../problems/P8/s6/continuation/s6_240/FORMULATION.md) supplies the
actual heavy Gaussian state H8A420-SLE-PRE-T0, not merely a class of
possible states. Its [bound proof](../problems/P8/s6/continuation/s6_240/notes/bounds.md)
controls the renormalized nonlocal mode remainder, for energy and pressure
and each time derivative through order five, by

```text
10^301 (n+p^2/16)^(-5/2),    t in [-1,1], p>=0.
```

Keep the full fixed local heat-counteraction and vacuum terms. Split
only this already-subtracted mode integral at a **fixed comoving** P>0.
Since `(n+p^2/16)^(-5/2)<=1024/p^5`, its omitted normalized readout obeys

```text
|tail^(j)|/kappa
 <=10^301/(2pi^2 kappa) integral_P^infinity p^2(n+p^2/16)^(-5/2) dp
 <=256*10^301/(pi^2 kappa P^2)
 <29*10^301/(kappa P^2),       0<=j<=5.
```

The exact radial integral, before the measure, is

```text
64/(3n) [1-P^3/(P^2+16n)^(3/2)].
```

At the original kappa=10^800 and the test split P=10^100, the uniform
normalized bound is below `10^-697`. This is a bound on an **existing
state's readout**, not a naturalness estimate or a new Gaussian choice.
It lies well below the existing 10^-399 reference-profile budget; that
stress budget is distinct from the amplitude's lambda budget.

The split does not delete modes from the state or replace its CCR.
It is a certified way to evaluate the original integral. Since it is
fixed in comoving momentum, differentiation creates no moving-boundary
terms. A time-dependent physical cutoff P(t)=a(t)Lambda would require
those terms and its own matching proof. The raw unsubtracted stress
cannot be truncated using this bound.

### 6.2 What remains genuinely independent

The bound concerns the free heavy reference on [-1,1], not its complete
interacting response or a global-time theorem. It does not bound the
light/metric state, the full causal kernel or the nonlinear physical
volume. The existing linear light/matter/tensor preparations and
[S264 covariance estimates](../problems/P8/s6/continuation/s6_264/FORMULATION.md)
are retained; they are not a completed interacting constrained state.

Flat matching also still leaves the curvature/residue compensation
identified in the witness audit. A flat spectral construction cannot
silently fix that covariant information. And the
[S265 Gaussian coefficient obstruction](../problems/P8/s6/continuation/s6_265/notes/integrability.md)
prevents promoting unrestricted substitution of the linear Gaussian lapse
into every nonlinear coefficient to a quantum construction. It does not
exclude full-channel cancellations or every possible constrained
quantization. The
[S277 scope boundary](../problems/P8/s6/continuation/s6_277/notes/scope.md)
likewise does not identify different classical fold solutions with the
original prepared quantum bounce.

**Decision:** the heavy reference tail passes this narrowly defined
compatibility check. A cutoff alone supplies none of the remaining
interacting-domain, curved matching or global B estimates.

## 7. Next admission boundary

No examined construction is ready to be adopted as the complete common
parent. The regulator-as-physical-local-field proposals should not receive
more matching computation. Computational regularization remains available,
with an explicitly matched physical boundary action.

The promising surviving mechanism is **local, positive-spectrum matching
with action-determined correlations**, targeting the necessary projection
rather than an unnecessarily strong norm. Its next admission requirement
is a concrete model that determines or bounds the *remaining* independent
local/gravity and relevant curved terms. The scalar-sector construction
above demonstrates only one part of that mechanism. An added negligible
sector with the old arbitrary counterfunctional would fail admission.

Possible new boundary conditions can be proposed, including a complete
fixed-scale quantum-EFT definition. They must be presented as new physical
assumptions for review, with their effects on the original targets,
operator support and prepared-state matching stated. They are not
predictions extracted from the current classical action. This study does
not adopt a finite-zero boundary or require the user to invent a Wilson
coefficient.

Until such a model is specified, keep the present H8A420/RATE4-D8 route
**conditional and parked for full matching**, rather than claim that the
selected shell or reference-state tail has closed its input gap. The
feasibility authorization has been used; it need not be requested again.
Adopting a replacement quantum model remains a separate decision.

## 8. Reproduction and evidence

Run the [read-only diagnostic](../scripts/p8_prescription1.py):

```sh
.venv/bin/python -B scripts/p8_prescription1.py
.venv/bin/ruff check scripts/p8_prescription1.py
.venv/bin/ruff format --check scripts/p8_prescription1.py
```

It protects 318 inherited/direct source files and checks five cutoff
identities, thirteen independent channel projections, sixteen parameter-
integrated shell coefficients, four shell enclosures, six boundary-
compensation controls, three positive-spectral fixtures, a half-line sign
proof, a stable-potential square completion, the state-tail antiderivative,
four tail-scaling identities and seventeen invalid inputs. Infinite-tail,
reflection-positivity and microcausality conclusions use the written
arguments, not extrapolation from a finite grid.

Evidence level is **VERIFIED_N** with written proofs in the stated scopes,
not FORMALIZED or a complete physical matching certificate. See the
[validation receipt](validation/p8-prescription1-2026-09-22.json).
No frozen certificate or 115,487-test publication replay is changed.
