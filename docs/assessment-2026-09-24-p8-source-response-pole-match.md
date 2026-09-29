# PRESCRIPTION-2: constrained source vertices and gravity pole matching

Date: 2026-09-24. Continuation of the authorized two tracks after the
[reference/Regge packet](assessment-2026-09-24-p8-reference-regge.md).
This is a new calculation in the existing parent, not adoption of the
proposed finite-zero prescription. M/V/G/B/R and original P8 remain open;
the physical gravity verdict remains UNTESTED.

**New result:** the actual source-tagged cubic and quartic Hamiltonian
vertices are derived after the temporal-vector and lapse constraints,
with both original canonical momentum shifts. They determine the
source-squared one-loop scalar and longitudinal-vector response kernels
at a common finite regulator. The quartic contact is not optional and
has no fixed sign. Renormalized continuum values and bounds are not yet
supplied. On gravity, a conditional forward Regge normalization is fixed
by the actual graviton pole, but finite-angle shape and remainder bounds
are still missing.

## 1. Inherited results, and the perturbative scope

The previous packet's heavy-source filtration is already stronger in
[S250](../problems/P8/s6/continuation/s6_250/notes/graphs.md). Likewise,
[S253](../problems/P8/s6/continuation/s6_253/notes/coupled.md) already solves
the joint Gaussian constraint/measure problem, and
[S257](../problems/P8/s6/continuation/s6_257/notes/measure.md) proves the
finite-regulator nonlinear second-class determinant cancellation on a
regular single-root branch. These are inherited, not new discoveries.
They do not supply the common covariant continuum measure or its finite
conversion. The new work carries the nonlinear source through this
reduction to the interaction vertices required for the response.

Use the formal grade-zero clock root for vertices inside a one-loop
calculation. The original fixed stress profiles have loop grade one;
their insertions in these loops have higher total grade. This does not
delete the profiles from the physical action, change their held-fixed
variation rule, or reset the initial state. Using the unexpanded S253
current pivot in a kernel is a resummed notation, not proof of a live
nonlinear root with nonzero profiles or an all-orders error bound.

In particular, the off-clock unforced classical datum in S257 is NOT
substituted for the original quantum-reference clock. The required
common-reference/curved counterfunctional work remains open.

## 2. Joint trace and temporal-vector elimination

Use [S257's full canonical Hamiltonian](../problems/P8/s6/continuation/s6_257/notes/canonical.md).
Let `a_tr=-M/3`, `M=R^(1/4)`, `U=R^(-3/4)`, `r=R-1`,
`c=-3 Hclock r/N`, `p=2 tr(pi_metric)/(3 V)`, and
`G=div(pi_W)/V` in its kappa-normalized notation. Introduce a marker xi
in the ENTIRE source square, `U[T-xi(r K+c)]^2/2`; xi is bookkeeping
and the physical value remains one. R, M, U and all source-independent
coefficients are held unchanged.

Starting from `L=a_tr K^2+B K+U[T-xi(r K+c)]^2/2`, the two stationary
equations for `p K-L-T G` give

```text
K=(p-B-xi r G)/(2 a_tr),
T=xi(r K+c)-G/U,
H_red/N=(p-B-xi r G)^2/(4 a_tr)+G^2/(2U)-xi c G+H_other.
```

The new diagnostic solves these two equations independently. The other
Hamiltonian terms and the spatial momentum constraint are source-marker
independent in these canonical variables. They still determine the
actual linear lapse and the source-independent response.

At the reference, `a=(1+u^2)^2`, `h=(1+u^2)^3`, `H=4u/(1+u^2)`,
`ell=1/(10a^3)`, `R=1`, `R_N=-2/h`, and `B=0`.
Write the scalar variables as `y=(v,sigma,pv,ps)`, where the momenta are
the physical density momenta divided by `kappa*a^3`.
The [original boundary transport](../problems/P8/s6/continuation/s6_275/notes/boundary.md)
is

```text
pv_raw=pv-18H v+3ell sigma,    ps_raw=ps+3ell v.
```

Both shifts and the time generator are retained. The time generator is
source-marker independent, so it adds no xi-dependent cubic or quartic
vertex. This is the same prepared picture, not a freshly minimized state.

Define

```text
E=1-3/(2h),  Theta=H-u/(1+u^2)^4,  q=k^2/a^2,
P=pv+3ell sigma,
n=[Theta P+ell E ps-(2E q+3 T_profile)v]/(2J).
```

Here `J` and `T_profile` have exactly their S253 meaning; the tree-clock
coefficient uses `J=J_tree` and `T_profile=0`. In particular
`J_tree(0)=243/160`, including the primitive's nonzero differentiated
contribution. The script binds this value to the pinned S180 expression,
not to a guessed oscillator.

The xi-linear Hamiltonian term before lapse elimination is
`xi r [3N(p-B)/(2M)+3H] G`. Its two lower field degrees vanish. Expand
`p=(pbar+pv_raw/3) exp(-3v)` with `pbar=B-2H`. The boundary shift cancels
`9 pbar v` and leaves P. At arbitrary reference time the cubic is

```text
h3_tag1 = -xi G n (P-6Theta n)/h.
```

All displayed h's are densities before multiplication by `kappa*a^3`.
The finite clock jets come from the exact full-parent factorization in
[S253](../problems/P8/s6/continuation/s6_253/notes/source.md), not a global
replacement of R or F by their clock polynomials.

## 3. The lapse-generated quartic contact

The direct xi-squared term after temporal elimination is
`-3 xi^2 G^2 n^2/h^2`. It is not the entire reduced quartic. The lapse
Hamiltonian has quadratic pivot `-2J`. Eliminating its second-order
source-induced correction adds the Schur contact
`(partial_n h3)^2/(4J)`. Thus

```text
h4_tag2 = xi^2 G^2/h^2 * [ (P/2-6Theta n)^2/J - 3n^2 ].
```

For an independent check, introduce a field-order marker e and solve
the following lapse polynomial exactly:

```text
H(n)=-J n^2+L n-e xi G n(P-6Theta n)/h-3e^2 xi^2 G^2 n^2/h^2,
n_star=(L-e xi G P/h)/(2[J-6e xi Theta G/h+3e^2 xi^2 G^2/h^2]).
```

The e and e-squared coefficients after substitution reproduce h3 and
h4 with `n=L/(2J)`. Terms omitted from this polynomial cannot affect
the specified xi-squared quartic: source-independent nonlinear
constraint corrections first enter that coefficient at field degree
five, and the other quartic corrections are xi-linear. The nonlinear
spatial reduction is xi-independent; its linear dictionary suffices at
this source degree. No claim is made to have computed the xi-linear
quartic or all other interactions.

At the tree bounce,

```text
n=(q v-ps/20)/(2J), P=pv+3sigma/10, J=243/160,
h3_tag1=-xi G n P,
h4_tag2=xi^2 G^2[P^2/(4J)-3n^2].
```

The contact is positive for n=0, P nonzero and negative for P=0,
n nonzero. No sign, smallness, stability or parent exclusion follows
from it alone. Dropping the cross-boundary shift changes h3 by the
nonzero term `-3xi ell G n sigma/h` even at the bounce.

## 4. Same-state finite-regulator response kernels

Use the original reference product of the COUPLED two-scalar state and
the original Proca state, with the canonical transports of
[S253](../problems/P8/s6/continuation/s6_253/notes/ctp.md) and S275.
Do not factor the two scalar modes or erase their mixed covariances.
If z denotes unit-CCR scalar phase variables, then
`y=C_y^-1 z`, `C_y=sqrt(kappa) diag(1,1,a^3,a^3)`.
The Wightman kernel W_y is obtained by these time-dependent endpoint
maps from `S(u)[V0+i Omega/2]S(v)^T`, retaining the same V0.
W_G is the divergence-momentum kernel from the original Proca transport.

All products below are defined first at the SAME finite regulator, with
Weyl-ordered phase-space insertions. This is not an unproved choice of
continuum covariant subtraction. Let N and R_P denote the linear rows
for n and P, and set

```text
F=N^T R_P+R_P^T N-12Theta N^T N,
A=R_P/2-6Theta N,   K=A^T A/J-3N^T N,
D(u)=kappa*a(u)^3/h(u).
```

For an external scalar variation e define
`n_e=N e`, `P_e=R_P e` and the internal row
`L_e=n_e R_P+(P_e-12Theta n_e)N`.
The cubic variation is `-D G L_e y`. Its connected two-vertex kernel is

```text
C_scalar(e,f;x,y)=D(x)D(y) W_G(x,y)
                  [ L_e(x) W_y(x,y) L_f(y)^T ].
```

This is the mixed scalar-vector loop, not a source-free vector
determinant. In position space N is a differential row; its derivative
acts on the appropriate leg. In momentum space this is a CONVOLUTION:
the external, internal-scalar and internal-vector momenta are different.
The row function accepts a specified momentum, and does not silently
use the external momentum for every leg. Homogeneous external scalar
probes can still couple to nonzero internal spatial momenta.
The implemented polarized matrices use distinct left/right rows:
`F_lr=N_l^T R_r+R_l^T N_r-12Theta N_l^T N_r` and
`K_lr=A_l^T A_r/J-3N_l^T N_r`. They obey `F_lr=F_rl^T` and
`K_lr=K_rl^T`, but need not individually be symmetric at unequal
momenta. Six unequal-momentum controls check the external/internal
variation and reject the diagonal-mode shortcut. F and K in the
position-space kernels mean these leg-wise differential bilinears.

For a longitudinal-vector probe write `g_e=delta G`. The two-internal-
scalar bubble is

```text
C_long(e,f;x,y)=D(x)D(y) g_e(x)g_f(y)
               Tr[F(x) W_y(x,y) F(y) W_y(x,y)^T]/2.
```

The source-squared mixed scalar/vector block has zero Gaussian mean
because it leaves one internal vector and an odd scalar product. This
does NOT set the separate source-linear mixed response to zero.

For either connected kernel, `chi=-i<[V_e,V_f]>=2 Im C`; in the
inherited effective-action convention the nonlocal retarded piece is
`-theta(u-v) chi`. The ordinary transpose in the bubble is essential.
The two local terms in the Hamiltonian second variation are

```text
local_scalar(e,f;x)=2 kappa*a^3/h^2 * <G(x)^2>
                    [ (A e)(A f)/J -3 n_e n_f ],
local_long(e,f;x)=2 kappa*a^3/h^2 * g_e g_f Tr[K V_y(x,x)].
```

Their effective-action sign is minus. Include both before attempting
any continuum bound. TT and transverse-vector external legs have no
linear n, P or G and hence no xi-squared one-loop quadratic term from
these vertices. Other source-independent terms and the xi-linear
mixed block are not computed here.

The completeness claim is narrow but useful: for a source-squared
one-loop two-point graph, field-degree counting permits two cubic
source vertices or one xi-squared quartic. Higher vertices cannot enter
at this loop/field degree. For a pure-light one-point graph neither is
possible at one loop. These facts do not eliminate the original free
Proca/heavy metric responses or higher-order source effects.

Three correlated Gaussian fixtures independently compare the scalar
kernel with the full block-matrix Wick trace. An indexed four-scalar
contraction checks the longitudinal bubble's factor of one half, and
both local Hessians are checked. These fixtures test algebra, not replace
the physical state or evaluate its continuum loop integrals.

## 5. A subtraction-sensitive contact, not a cutoff prescription

As an independent ultraviolet control ONLY, freeze a flat metric and
take the massive longitudinal Proca vacuum, with canonical field
`A=sqrt(kappa*zeta) W`, `m^2=1/zeta` and `omega=sqrt(k^2+m^2)`.
The longitudinal kinetic factor is `m^2/(k^2+m^2)` and
`<Pi_A^2>=m^2/(2omega)`. Consequently

```text
<G_k G_-k>=k^2/(2kappa sqrt(k^2+m^2)),
I_G(Lambda)=1/(4pi^2 kappa) integral_0^Lambda k^4/sqrt(k^2+m^2) dk,
primitive=[k sqrt(k^2+m^2)(2k^2-3m^2)+3m^4 asinh(k/m)]/8,
I_G(Lambda) ~ Lambda^4/(16pi^2 kappa).
```

The derivative and asymptotic coefficient are checked exactly. This
control is not the original curved state's finite part, a physical
cutoff, or a proof that a divergence survives the combined fish/contact
and covariant-counterterm calculation. It demonstrates why setting the
local contact to zero by reference normal ordering is not a completed
common-reference prescription. Its finite conversion must be tracked
alongside the rest of the effective action.

## 6. What graviton pole matching actually fixes

Keep the original mass-one amplitude convention, `s0=2`, the same two
unequal denominators, detector `E=1/256` and trial split T at least 32.
For this paragraph assume a single leading even-signature Regge sector
on both cuts, with

```text
rho(s,t)=C_T(t) s^2 (s/T)^(j(t)-2)+residual,
j(0)=2,  j'(0)=alpha'>0.
```

Also require the appropriate subtracted dispersion relation, a
controlled complex contour, regular finite-window terms after the
specified pole/low-cut subtractions, and no extra unaccounted 1/t term
from the residual. These are substantive hypotheses, not conclusions
about the present physical observable. In particular, uncontrolled
massless-loop singularities do not meet them automatically.

The two-cut tail has singular part
`2 C_T(t)/[pi(2-j(t))]`. The original graviton pole in the quadratic
coefficient is `-1/(kappa t)`. Matching these coefficients gives

```text
C_T(0)/alpha' = pi/(2kappa).
```

This follows in our normalization, not by importing the numerical
coefficient from a massless amplitude. For an EXACT leading term the
finite part of its tail, after removing that pole, is also fixed as

```text
1/kappa * [ -C_T'(0)/C_T(0) + j''(0)/(2alpha') + alpha' D_T ],
D_T=log[T/(T-2)]+4/(T-2)+2/(T-2)^2.
```

Here D_T is the integral of `s^2/(s-2)^3-1/s` from T to infinity; the
crossed denominator agrees at t=0. A residual contributes its own
finite part, which is NOT assigned zero. Thus even exact pole matching
does not fix the finite subtraction. The script checks the Laurent
coefficient, the finite part and D_T's derivative and endpoint.

For the detector interval suppose additionally that physical bounds
establish `|C_T(-q^2)|<=B C_T(0)` and
`[2-j(-q^2)]/q^2>=a_rel alpha'>0` for EVERY `E<=q<=1`.
Then the preceding packet's tail estimate implies

```text
R_tail < 95 B/(kappa a_rel) + supplied_Regge_residual_error.
```

With the actual kappa and lambda this leading allowance is below
lambda/100 IF `B/a_rel<=10^196`. This sensitivity test supplies neither
B nor a_rel. Pole matching fixes a forward value, not these interval
bounds, the finite-window moment, the contour or the physical-observable
conversion. Moving the split also transports the residue:
`C_T2(t)=C_T1(t)(T2/T1)^(j(t)-2)`; it is not an automatic improvement.

The research cross-check was
[Noumi and Tokuda, arXiv:2212.08001v2](https://arxiv.org/abs/2212.08001),
especially sections2,4,6. Their massless-scalar analysis uses leading
Regge dominance and a finite complex-arc approximation. Its finite-energy
sum rules relate spectral moments to trajectory/residue data; they are
not identities between four RATE4 samples and a complete Regge envelope.
Their treatment of massive internal light loops does not establish
control of the full massless-graviton infrared problem here. We import
no benchmark slope, string scale or numerical positivity bound. The
algebraic two-moment inversion is checked only as a control, not applied
to unsupplied physical moments.

## 7. Closure-driven next work and validation

The next model calculation is now explicit: evaluate the two source-
squared response blocks and their local subtractions using the SAME
prepared propagators, compute the source-linear mixed block where
needed, and combine them with the remaining source-independent sectors
and common finite conversions. A continuum bound must keep the finite
boundary functional's transport; normal ordering cannot supply a new
zero condition without adoption. Neither the full nonlinear mean nor
the original B gate is established by the finite-regulator formulas.

On gravity, the remaining job is physical finite-window spectral data,
a uniform finite-angle Regge/contour estimate and observable-conversion
errors, or a different admissible proof of the same original gravity
functional. The conditional forward normalization narrows the task but
does not certify any of those hypotheses. No new model, detector,
cutoff, RATE5 condition, boundary zero or weaker finish line is adopted.

The read-only [diagnostic](../scripts/p8_source_response.py) pins 419
inherited report/source inputs. It checks the joint Legendre solve,
literal arbitrary-time source jet, original boundary omission control,
exact lapse Schur contact, primitive-sensitive bounce value, nine
time/momentum row controls, six unequal-momentum polarizations, three
correlated Wick controls, indexed
longitudinal bubble, both local blocks, ultraviolet primitive, conditional
pole/finite-part/transport identities and 27 invalid/missing input
rejections. These are exact algebra checks plus the written field-degree
argument, not a formal proof assistant or a continuum numerical bound.

The [receipt](validation/p8-source-response-2026-09-24.json) records the
focused replays and preservation checks. Frozen P8 packages, prior
milestone packets and unrelated P4/P9 edits remain unchanged. The full
115,487-test release suite is not rerun for this standalone continuation.
