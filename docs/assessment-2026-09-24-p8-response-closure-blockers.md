# PRESCRIPTION-2 and gravity: continuum extension and closure blockers

Date: 2026-09-24. Continuation of the authorized tracks after the
[constrained-source packet](assessment-2026-09-24-p8-source-response-pole-match.md).

**Outcome: blocked for physical completion with the retained inputs.**
The source-squared noise kernels exist in the continuum and their
time-retarded responses admit local extension families. But the retained
normalizations do not select a unique physical finite response: two
explicit covariant functionals preserve the first-order flat matching
data and reference mean while changing independent response directions.
The current gravity route likewise lacks a quantitative same-observable
remainder bound. This is a specified-input obstruction, not exclusion of
the parent, every gravity route, or original P8.

No new finite boundary condition, parent, state, physical cutoff, detector,
RATE5 target or weakened completion criterion has been adopted. M/V/G/B/R
and original P8 remain OPEN; the physical gravity verdict is UNTESTED.

## 1. Admission and stopping condition

This checkpoint addresses the common-reference input to M/B and the
physical gravitational allowance in G. The intended positive outcome was
a unique renormalized source response with usable bounds and a justified
finite-gravity allowance. The alternative decisive outcome was an
explicit obstruction showing which additional physical input is needed.

The previous finite-regulator vertices, conditional pole matching and
small known-sector contributions could not distinguish these outcomes.
Here the nonzero rank of the finite-response ambiguity supplies the
obstruction. For gravity the comparison margin is still the original
lambda=10^-600; a representative allocated absolute error is lambda/100,
not a new necessary finish line or an assumed physical bound.

Apply the living plan's CONDITIONAL ONLY stop rule: do not expand more
known-sector calculations as if they could determine missing independent
finite data. A separately specified physical boundary candidate can
reopen the calculation. The remaining computations are listed below,
but their completion alone cannot resolve the demonstrated ambiguity.

## 2. Continuum source-squared distributions

### Inherited state, not a new preparation

The scalar state is the full coupled two-mode preparation from
[S251](../problems/P8/s6/continuation/s6_251/notes/ultraviolet.md), with
the original momentum transports retained. The Proca state is the
all-order reference state of
[S176](../problems/P8/s6/continuation/s6_176/notes/state.md), not a
finite-order WKB vacuum. Work on compact subslabs of (-1/2,1/2), with
compact smooth spacetime tests. No unrestricted scalar shift potential
is required by the N, P and G rows used here.

S251 supplies two strictly positive characteristic frequencies c_s/a
and 1/a, all-order oriented parametrices and smooth low-momentum
contributions. In the original prepared phase variables, position
amplitudes have order k^(-1/2), while uncleaned momenta have order at most
k^(3/2). These statements include both chart transitions and the full
mixed covariance. They are not smallness estimates.

Use the previous packet's exact rows

```text
n=N y=[Theta P+ell E ps-(2E q+3 T_profile)v]/(2J),
P=pv+3ell sigma,  G=div(Pi_W)/(kappa*a^3).
```

Both n and P have amplitudes of order at most k^(3/2). For the canonical
Proca coordinate A=sqrt(kappa*zeta) W, the longitudinal kinetic factor is

```text
K_L=a/(1+zeta*k^2/a^2),
pi_A=sqrt(K_L)*(v_L_dot-d_L*v_L).
```

Thus pi_A has order k^(-1/2), and G has order k^(1/2). Its leading
adiabatic covariance symbol is

```text
zeta*k^2*K_L*sqrt(k^2/a^2+1/zeta)/(2*kappa*a^6)
    ~ k/(2*kappa*a^4).
```

The displayed symbol divided by k^2 tends to
sqrt(zeta)/(2*kappa*a^5) at zero momentum. This last limit is a
normalization control, not identification of the actual infrared state
with the leading adiabatic approximation. Infrared regularity follows
from the original state's smooth finite-momentum construction.

### Ordered noise and time-retarded extension

Every contracted Wightman factor has the same strictly negative
first-leg time orientation, including the Proca factor. Their wavefront
covectors cannot sum to zero. The scalar-vector fish and the two-scalar
longitudinal bubble therefore define ordered distribution products on
the whole doubled slab. Real and imaginary parts are taken after the
ordered products, not by forming an unsupported product W*conjugate(W).
The centered Hermitian quadratic insertions have positive noise by the
same-state smoothing/GNS limit used in
[S252](../problems/P8/s6/continuation/s6_252/notes/wick.md).

Source-specific power counting sharpens the general S252 bound. On a
double internal momentum shell of size rho, the six-dimensional measure
contributes rho^6. The scalar-vector amplitude contributes at most
rho^(3+1), and the two-scalar amplitude at most rho^(3+3). Consequently

| Kernel | Before external differential rows | Complete prepared-phase kernel |
|---|---:|---:|
| Scalar-vector fish | scaling degree <=10 | <=14 |
| Two-scalar longitudinal bubble | scaling degree <=12 | <=14 |

Each external scalar row has at most two spatial derivatives and each
external G row one. Derivatives falling on smooth coefficients do not
increase these bounds. Unequal momenta and all cross entries remain;
there is no diagonal-mode replacement in the convolution.

Here is the distributional, rather than pointwise, estimate. Use center
coordinates Z and four normal coordinates Y=(t-s,x-y) to the diagonal.
For a test supported on a normal annulus, a dyadic shell obeys the S252
bound C_N rho^omega (1+epsilon*rho)^(-N), now with omega as in the table.
Away from equal times, integrate in the normal time using the positive
sum of internal frequencies. Near equal times the annulus is spatially
separated; integrate in a high internal momentum. Split the shell to
avoid differentiating angular symbols at zero momentum. Tangential
phase derivatives cost only epsilon*rho, absorbed by increasing N.
The same argument applies to theta(t-s) times the commutator off the
diagonal. Its wavefront closure at the diagonal is conormal. Low
momenta and sufficiently deep parametrix remainders are smoother.

Summing above and below rho=epsilon^-1 gives scaling degree at most
omega, uniformly for the required tangential derivatives. These are
the hypotheses needed for normal-coordinate extension, not an
assumption that this constrained two-speed system is Klein-Gordon.
The finite-scaling-degree extension results used are
[Brunetti and Fredenhagen, Theorems 5.3 and 6.9](https://arxiv.org/pdf/math-ph/9903028).
Their general distribution theorems apply after the preceding checks;
their scalar-field construction is not imported as a DHOST theorem.

For example choose smooth compact w(Y)=1 near zero and subtract

```text
W_10 phi(Z,Y)=phi(Z,Y)
  -w(Y)*sum_{|alpha|<=10} Y^alpha/alpha! * partial_Y^alpha phi(Z,0).
```

The remainder vanishes through normal degree ten. Codimension four and
the bound fourteen leave a positive power gain
4+11-14=1, sufficient also with harmless logarithms. A shrinking-hole
limit applied to this remainder gives an extension; assigning the
removed jets completes it on all tests. A smooth partition of unity
patches the local construction. Differences of such extensions have
the form

```text
sum_{|alpha|<=10} c_alpha(Z) partial_Y^alpha delta^4(Y).
```

Diagonal additions preserve time-retarded support. Before the external
rows the corresponding caps are six and eight. The respective four-
normal-coordinate Taylor spaces have dimensions 210, 495 and 1001.
These are coarse jet dimensions, NOT numbers of independent covariant
physical operators or evidence that every allowed derivative occurs.

This completes the scoped continuum existence/ambiguity argument. It
does not choose the c_alpha, prove a single covariant variational
counterfunctional or its Ward identities, evaluate the response, or
bound it for B. The local quartic contacts contain separately
renormalized coincident composites such as <G^2> and <nP>; connected
noise positivity does not set their physical finite parts to zero.
The source-linear mixed block and other reference sectors remain.

## 3. Explicit physical finite-data obstruction

### Two regular covariant witnesses

Write X=-g^(mu nu) partial_mu u partial_nu u, and F(W)=dW. Consider
the following possible DIFFERENCES between quantum boundary completions:

```text
DeltaGamma = hbar*kappa*integral sqrt(-g) [cN*O_N+cG*O_G],
O_N = X^4*(X-1)^2,
O_G = zeta^2*X^4*[(nabla_nu u) nabla_mu F^{mu nu}(W)]^2.
```

They are smooth local covariant expressions with no inverse X. No
claim is made that the preceding loops generate precisely these
counterterms or that every extension coefficient has this form. These
are two allowed extra finite directions unless a physical boundary
specification excludes or fixes them. They are NOT adopted changes to
the frozen classical action.

At the flat vacuum u=0, partial u=0, W=0 their minimum field degrees are
eight and twelve. Their jets through degree four vanish, so they do not
change the retained masses, residues, cubic sources or four-point RATE4
data at first quantum order. An order-hbar vertex with d fields first
contributes to a four-external-leg graph at total loop grade at least
1+(d-4)/2, namely three and five here. This lower bound allows the
original additional vertices of valence at least three and does not
infer nonzero values for any particular contraction. It is not an
all-orders equivalence of the full amplitudes.

At the reference u=t, X=1, W=0, both densities and all their first
functional variations vanish. Thus the same reference mean and its
fixed stress profiles are unchanged at this grade. Their second
variations do not vanish. With lapse N_lapse=1+e*n,
X=N_lapse^-2 and the actual physical volume, their normalized quadratic
action densities are

```text
4*cN*n^2,    cG*G^2.
```

The vector equality uses (nabla u).div(F)=G/zeta on the reference,
including the canonical normalization Pi_W=kappa*zeta*a*E_i.
Squaring removes the Fourier phase convention. The vanishing first
variations hold for compact general variations, not just homogeneous
lapse variations: O_N has a double zero at X=1 and the scalar squared
in O_G is zero at W=0.

### Rank two after the original constraints

At the formal tree bounce q=1, J0=243/160 and
F0=1199/800=J0-1/50. The original prepared lapse row is

```text
n=(80*v-4*ps)/243.
```

The Hessian response map on the independent normalized probes v and
longitudinal p_W (with G amplitude one) is

```text
diag(51200/59049, 2),     determinant=102400/59049 != 0.
```

Consequently the two coefficients cannot be inferred from the matching
conditions that both witnesses leave invariant. The scalar witness
shifts J to J+4*hbar*cN while F is unchanged. The actual two-mode
principal pencil from
[S241](../problems/P8/s6/continuation/s6_241/notes/clock.md) has roots
1 and F/(J+4*hbar*cN). At the tree bounce,

```text
d(c_s^2)/d(hbar) = -153472*cN/59049.
```

This is a response change, not a disposable constant or a forgotten
canonical boundary phase. Small opposite coefficients preserve the
local tree scalar gap in this control; no complete quantum stability
or global neighborhood theorem is inferred.

An independent frozen-coefficient longitudinal check at a=1 writes

```text
zeta_eff=zeta+2*hbar*cG*zeta^2*k^2,
L=zeta_eff*(Wdot-k*W0)^2/2+W0^2/2-W_L^2/2,
W0=zeta_eff*k*Wdot/(1+zeta_eff*k^2),
H=(1/zeta_eff+k^2)*p_W^2/2+W_L^2/2.
```

At first order deltaH=-hbar*cG*k^2*p_W^2, giving
delta(omega_L^2)=-2*hbar*cG*k^2 in this frozen-coefficient control.
This is at fixed finite momentum, not a stationary dispersion law on
the evolving curved clock, a high-k front-velocity statement, or a
complete healthy UV completion.

### Physical boundary data versus scheme transport

A completed common reference contributes some coordinates rN,rG; the
physical response depends on r+c. A reference shift r->r+eta must be
accompanied by c->c-eta. The rank-two response is then unchanged.
Resetting c to zero after changing reference instead changes it by
the displayed response matrix times eta.

Thus covariance, the original four rates, reference means and positive
noise do not determine these physical finite directions. A specific
complete quantum theory could fix them; the retained inputs have not
done so. The statement is not that two theories with every physical
observable fixed can nevertheless differ.

## 4. Gravity: an asymptotic remainder is not a physical error budget

Rechecking [Bellazzini et al., *Positivity with Long-Range Interactions*,
version 2, sections 4-5](https://arxiv.org/html/2512.13780v2) does not
remove the gap. The finite-detector argument has a finite-coupling
remainder of order G M^2 (notably Eqs. 5.9-5.17); it does not provide a
numerical constant and validity domain for this project's exact
observable. The paper also distinguishes finite-resolution approximate
bounds from its scaling-limit unitarity argument. No assertion that
our fixed detector makes the error order one is needed or made.

The previous conditional pole match fixes C(0)/j'(0)=pi/(2*kappa)
only under its explicit Regge/dispersion assumptions. It still does
not give finite-angle shape bounds or the residual contour/observable
conversion. A finite-energy alternative is admissible, but its actual
moments and finite-contour errors must be bounded for the same
observable; an all-energy UV theorem is not mandatory by definition.

Absorb fixed units into C_err and write an absolute remainder as
C_err/kappa. With the original kappa=10^800 and lambda=10^-600, the
representative allocation requires

```text
C_err/kappa < lambda/100  =>  C_err < 10^198,
```

AND a proved validity domain containing the actual parameters. This
large threshold may be encouraging, but an unspecified constant has no
certified upper bound. For an exact logical control, f(g)=g and
f(g)=(2*kappa*lambda)*g are both entire O(g) functions with fixed
constants. At g=1/kappa one is below the allocation and the other is
2*lambda, above it. These are not claimed physical UV amplitudes or
counterexamples; they show why the asymptotic notation alone cannot
certify the inequality, even without an asymptotic-radius issue.

The current route needs a certificate containing:

1. The same physical observable and its low-cut subtraction definition.
2. A bound on conversion to the positive reference comparison.
3. Finite-angle or finite-contour spectral control in that definition.
4. A quantitative finite-coupling remainder constant.
5. A validity domain containing the actual parameters and detector.

None is supplied by assigning an order-one Regge factor or by counting
the small power of 1/kappa. No missing input is assigned a value here.
This blocks certification through the retained-input route, not all
possible future gravity proofs.

## 5. Decision needed, and what it would not solve

**Recommended decision:** authorize the separately named
RATE4-COVZERO-1 finite-boundary candidate described in the
[initial PRESCRIPTION-2 proposal](assessment-2026-09-22-p8-prescription2-gravity.md).
Its zero-complement rule would be additional physical boundary data
relative to a fully specified common reference, not a prediction of
the four old rates. Do not request arbitrary numerical cN or cG from
the user.

Adoption must remain conditional on spelling out that reference:
dimensional/evanescent continuation, gauge fixing, constrained measure,
source contacts, causal/state prescription and common finite
conversions. The existing reference is not yet executable as a full
covariant prescription. Later changes of reference must transport the
whole boundary functional, not reset its complementary coefficients.

If authorized, complete and audit this reference and evaluate the
remaining same-state kernels and contacts for that named candidate.
The source-linear mixed block, source-independent sectors, quantitative
curved bounds and higher-order errors still need work. The witnesses
above would then be excluded or fixed by the new boundary rule, not
proved numerically zero in every convention.

Keep gravity as an independent certificate obligation: authorizing a
finite-boundary candidate does not authorize a Regge envelope, spectral
budget or physical allowance by assumption. A candidate-specific
physical construction or another proved same-observable theorem is
needed. Neither this authorization nor a successful first-order
renormalization alone would close original P8.

If the original retained physical inputs are to remain the entire
specification, these two routes stop as CONDITIONAL ONLY. Further
finite-regulator or negligible known-sector refinements do not remove
the obstruction. This is the point for user direction, not another
unqualified instruction to continue the same calculation.

## 6. Reproducibility and scope

The new [diagnostic](../scripts/p8_response_closure.py) pins 494 inherited
files. It checks the refined orders, normal Taylor projectors, Proca
normalization limits, covariant field degrees and clock jets, the
rank-two reduced response, principal scalar sensitivity, independent
temporal-vector Legendre solve, finite-reference transport and exact
gravity allocation controls. Nine malformed inputs are rejected.
These are exact algebraic controls accompanying the written continuum
argument, not a machine-checked theorem of quantum gravity.

The [validation receipt](validation/p8-response-closure-2026-09-24.json)
records all twelve focused root diagnostic replays, lint/format and
preservation checks. The frozen 115,487-test suite is not rerun for this
root-level analytical checkpoint. A first exact-budget control exposed
Python floating-point underflow in 1/10^800; using Fraction(1,kappa)
corrected the implementation before the final successful runs.

Prior milestone files and unrelated P4/P9 work are preserved. Frozen
P8 sources and the three tracked precursor diagnostics are unchanged.
The new work is local, uncommitted and unpushed. No original completion
gate is promoted and no background research is claimed after handoff.
