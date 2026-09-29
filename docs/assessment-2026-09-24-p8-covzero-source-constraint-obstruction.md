# COVZERO coupled-source admission: the fixed-count reference does not extend

Date: 2026-09-24. Continuation after the user-authorized
[COVZERO boundary and dimensional seed](assessment-2026-09-24-p8-covzero-reference-seed.md).
The finite-boundary authorization remains valid. No classical source,
prepared state, physical cutoff or field domain is changed here.

**Result:** the retained nonlinear Proca source lifts the scalar/metric
kinetic null direction on generic tilted slices. On a purely spacelike
scalar gradient its added scalar-acceleration Hessian is
`(R-1)^2/X`, nonzero for X<0 and R!=1. The temporal Proca constraint
does not compensate for this loss. A written bound exhibits the effect
inside the ACTUAL retained parent's coefficient domain, arbitrarily
close to the vacuum in field/gradient values, not just in a toy R.

This rejects the proposed shortcut to **one exact fixed-count canonical
measure throughout that domain**. It does not invalidate the earlier
clock/unitary constraint calculation, the vacuum quadratic spectrum,
the pure-light dimensional seed or the legacy scale transport. Nor does
it by itself exclude a controlled perturbative, order-reduced EFT.
The latter requires its own validity and error argument; treating the
extra solutions as absent is not a consequence of COVZERO's finite
one-loop boundary data.

The common physical reference is not completed. The current exact
canonical-extension route has an action-level obstruction, preceding
the Ward and loop-counterfunctional calculations. A different classical
source completion, a restricted domain or an additional physical mode
prescription would be additional model data, not an authorized repair
already made in this packet. Original M/V/G/B/R/P8 remain open.

## 1. Literal action and principal variables

After the already exact 56-mode algebraic connection elimination, the
retained vector action is

```text
L_vector/kappa = -zeta F(W)^2/4 + (W-S)_mu (W-S)^mu/2,
S_mu = u_mu (R-1) [ -3 H_clock + 3 R_u/(4R) + Box(u)/X
                        + (-1/X^2+3R_X/(2RX)) u^mu u_munu u^nu ].
```

These are the pinned [S174 full source](../problems/P8/s6/continuation/s6_174/notes/source.md)
and [full affine reduction](../problems/P8/s6/continuation/s6_174/notes/affine.md),
with S238's entire new R substituted as required by that parent. Both
the linear source coupling and the source square remain. The mass term
is the physical-metric Lorentz contraction, not a positive Euclidean
square. The vector curl has no scalar acceleration or metric velocity
when W is the independent retained vector. Neither minimal scalar
matter field supplies such a velocity either.

Conventions: X=g_phys^(mu nu)u_mu u_nu, signature +---. R is the positive
tensor coefficient. There is a further prose sign correction to the
preceding seed note: the physical curvature density is
**-R*Ricci_P8/2**, not its displayed +R*Ricci_P8/2. This follows the
explicit [S177 Einstein convention](../problems/P8/s6/continuation/s6_177/notes/gravity.md).
The seed diagnostic and every kinetic identity already use the correct
inherited ADM density R*(KijKij-K^2)/2. No computed coefficient changes.

Use a local physical ADM frame and rotate the spatial scalar gradient
onto axis 1. Put

```text
u_mu=(a,b,0,0), X=a^2-b^2,
V=normal derivative of a, K_ij=metric extrinsic velocities,
Z=a V+b^2 K_11.
```

In the usual auxiliary-gradient first-order presentation, spatial
integrability fixes the time derivative of u_i by spatial derivatives
of u_0. Thus the principal scalar Hessian components are

```text
(u_00)_vel=V,
(u_0i)_vel=-K_ij u_j,
(u_ij)_vel=-a K_ij.
```

Their literal Lorentz contractions give

```text
(Box u)_vel=V+a K,
(u^mu u_munu u^nu)_vel=a Z,
[(u_munu u^nu)(u^murho u_rho)]_vel=Z^2.
```

The normalized pure-light principal density is exactly

```text
L0=R/2*(Kij Kij-K^2)-2R_X K Z
   +(R_X/X)*a*(V+a K)*Z
   +[-R_X/X-7R_X^2/(4R)+a^2 R_X^2/(RX)]*Z^2.
```

Let q be the velocity-dependent scalar multiplying u_mu in S_mu:

```text
q=(R-1)[(V+a K)/X+(-1/X^2+3R_X/(2RX))*a Z].
```

The R_u and H_clock source terms have no principal velocity. Setting
them aside for this Hessian does not delete them from the action or
replace the state. If c is the seven-component velocity gradient of q,
then the complete scalar/metric principal Hessian obeys

```text
H_coupled = H0 + X c c^T.
```

The terms linear in W do not change this velocity Hessian, but they are
essential in the separate temporal-constraint check below.

## 2. The source lifts the null direction

For a!=0, X!=0 normalize the pure-light null velocity by Z=1. With
l=b^2/X, its components are

```text
K11=R_X/R*(l-1/2), K22=K33=-R_X/(2R), off-diagonals=0,
V=[1-b^2 K11]/a.
```

Substitution into H0 gives zero. On this same direction the source is

```text
q_null=(R-1)*b^2/(a X^2)*[-1+3X R_X/(2R)].
```

This is generically nonzero at a tilted gradient. For an invariant
determinant proof, first change velocities from (Kij,V) to (Kij,Z),
whose determinant is a, and complete the pure-light metric square.
The six-dimensional DeWitt velocity Hessian has determinant -16 R^6.
Its remaining row is zero. A rank-one update then has determinant
`(-16 R^6)*X*q_null^2` in those coordinates. Transforming back yields

```text
det H_coupled = -4 R^4 (R-1)^2 b^4 (2R-3X R_X)^2 / X^3.
```

This is a closed-form identity, not an inference from numerical ranks.
The five independent rational controls in the diagnostic reproduce it
using both SymPy and Fraction elimination. They include aligned,
timelike tilted, spacelike tilted, purely spatial and R=1 controls.
The ranks are respectively 6,7,7,7,6; the pure-light rank is six in
all five controls. These generic-R fixtures are not substituted for the
physical parent. Section 4 gives a separate actual-parent argument.

For the branch a=0 do not divide by a. Directly from the action,

```text
X=-b^2, q=(R-1)V/X,
(H_coupled)_VV=(R-1)^2/X,
det H0_metric=-4 R^4 (2R+3b^2 R_X)^2.
```

There is no V-metric mixing on this branch. The scalar acceleration is
therefore a genuinely nondegenerate principal variable when R!=1 and
the displayed metric block is regular. Its negative entry for X<0 is
not, in isolation, a fully reduced spectral or instability-time result.
What is established is loss of the specific scalar primary constraint
that the proposed fixed-count reference would need to transport.

## 3. The temporal vector does not supply a compensating pair

On the rank-lifted branch use q itself as a velocity coordinate y;
this change is invertible precisely when q_null!=0. After completing
the pure-light metric square its block is independent of y. Write T
for the normal vector component, W for its spatial component along b,
and J=a T-b W. The part relevant to the temporal Hessian is

```text
L_mass=X y^2/2-J y+(T^2-W^2)/2,
p_y=X y-J,
H_mass=(p_y+J)^2/(2X)-(T^2-W^2)/2.
```

The other source terms are affine translations of y and introduce no
additional T^2. The complete Maxwell Legendre transform retains the
Gauss multiplier T G, which is also linear in T. Consequently

```text
partial_T^2 H = a^2/X-1 = b^2/X.
```

This is nonzero on the generic rank-lifted branch; it equals -1 at
a=0. Preservation of the temporal primary therefore solves T by its
usual second-class partner. It does not impose the missing scalar
primary/secondary pair. Dropping the Gauss term is unnecessary and
would be incorrect; the calculation retains it explicitly.

Ordinary diffeomorphism constraints persist. They must not be counted
a second time as a replacement for the separately lost scalar
degeneracy. Conversely, this local calculation is not a theorem about
every possible nonlocal boundary condition or EFT order-reduction.

## 4. An actual-parent spacelike domain, without numerical underflow

It is important to bind the obstruction to the literal fixed action.
At the anchor kappa0 and u=0 the old analytic R plus the S238 addition is

```text
R=1+(X-1)[T+(1-T)N X^2 exp(-N X^2)]
      +N X^2(1-X)^8 exp(-A X^2),
T=X^N/[X^N+(1-X)^N], N=1024, A=10^420.
```

The diagnostic extracts the six defining algebraic assignments from
the pinned original source without importing its ancestor computation.
The extra term is independently bound to the frozen S238 report.
At u=0 its germ is

```text
R-1=-7168 X^3+1024*(1052-10^420)X^4+O(X^5).
```

This germ is not used as a replacement action. To obtain an actual
nonzero interval, set X=-e, 0<e<=10^-430. Let B0=N e^2 exp(-N e^2).
Then exactly

```text
R-1=N e^2[(1+e)^8 exp(-A e^2)-(1+e)exp(-N e^2)]
     -(1+e)T(1-B0).
```

Use `1-z<=exp(-z)<=1`, binomial coefficients, 0<B0<1 and
`0<T<e^1024<e^4`. The bracket's lower bound is
`7e-A e^2(1+8e)>6e`. Subtracting the T term still gives
`R-1>6N e^3-2e^4>5N e^3`. Its upper bound is
`8e+N e^2(1+e)<10e`, so

```text
5N e^3 < R-1 < 10N e^3.
```

Differentiate the literal formula. The old-R derivative has magnitude
less than `e^4+N e^2+2N e^3+4N e<5N e`. The new-R derivative has
magnitude less than `4N e+16N e^2+4NA e^3<5N e`. Hence

```text
|R_X|<10N e,    2R-3X R_X>2-30N e^2>1,
1<R<6/5,       -1/4096<X<0,
(H_coupled)_VV=-(R-1)^2/e < -25N^2 e^5 < 0.
```

All side inequalities are exact rationals. Monotonicity extends the
endpoint checks to every e in the stated interval. No exponentially
small nonzero value is rounded to zero. Since e can be arbitrarily
small, these coefficient configurations approach the vacuum, also
in fixed canonical variables at the fixed finite kappa0.

This is an off-shell action/constraint-domain result. It does not
assert that each constant spacelike-gradient configuration solves the
coupled background equations, supply a physical pole frequency or
establish an instability on the actual prepared bounce. A common
off-shell canonical measure nevertheless cannot borrow a constraint
that is absent on this domain.

## 5. Consequences, surviving alternatives and authority

The old [clock result](../problems/P8/s6/continuation/s6_174/notes/constraints.md)
uses b=0 and the scalar-unitary chart. The source does not lift the
null direction there. At the exact vacuum R=1 the quadratic source
vanishes too. The original packet explicitly did not claim a full
constraint count on every background; this audit resolves part of
that previously open extension problem, not a failure of a stated
frozen theorem. The new seed likewise explicitly limited its
arbitrary-tilt degeneracy proof to the pure scalar/metric sector.

The literature distinguishes kinetic-rank loss from failure of
secondary constraints under matter coupling; our two separate tests
address both issues directly for the actual source. The generalized
Proca no-go theorem is not imported wholesale: its matter-coupling
assumptions differ from this explicitly scalar-dependent source.
[Garcia-Saenz](https://arxiv.org/pdf/2106.14960).

A tilted-slice rank change in a unitary-degenerate theory can represent
a nonpropagating, spatial-boundary-controlled mode on timelike scalar
backgrounds. That interpretation needs an actual equation and boundary
analysis, and is not automatically valid on spacelike gradients or at
du=0. The cited example does not prove such an interpretation for this
coupled model. [De Felice, Mukohyama and Takahashi](https://arxiv.org/pdf/2110.03194).

There are two distinct possible continuations:

- **Unchanged-action EFT route:** formulate a strictly perturbative,
  order-reduced common reference and prove its validity, omitted-error
  bounds and compatibility with the original clock/vacuum states. It
  must not resum the high-derivative interaction into an exact
  fixed-count propagator by assertion. No such bound, cutoff, new state
  or completed prescription is supplied here. Formal one-loop diagrams
  alone would not close physical M or B.
- **Classical source-completion route:** investigate a separately named
  scalar/vector/metric completion whose FULL kinetic and secondary
  constraints work on the required shared domain, preserving the stated
  flat targets and clock requirements where possible. Changing or
  deleting the retained source is additional physical model data, not a
  finite subtraction or an implementation of the authorization already
  granted to COVZERO.

Neither a clock-only restriction that loses the original vacuum nor a
newly chosen physical cutoff may silently repair the obstruction. A
first-order finite-boundary counterterm also cannot cancel an order-zero
constraint defect in the formal quantum expansion. No reference or
classical completion is adopted by this diagnostic.

**Recommendation:** do not proceed with the uniform exact canonical
extension as though this test had passed. Before changing the parent,
obtain authorization for a bounded, separately named source-completion
feasibility study; compare it with the unchanged-action EFT route and
require preservation tests before adopting any repair. The existing
COVZERO finite-boundary approval need not be repeated. The independent
gravity certificate remains missing as before. This is an obstruction
to this completion route, not a no-go for P8 or every EFT realization.

## 6. Validation

The [new read-only diagnostic](../scripts/p8_covzero_source_constraints.py)
pins 500 inherited/source files. It independently reconstructs the
covariant Hessian contractions, checks the source's null projection,
five coupled ranks and rational determinants, the a=0 branch, the
temporal-Proca Legendre Hessian, literal parent bindings and every
rational side inequality in the actual-domain proof. Four unsupported
measure-admission inputs are rejected. One initial structural SymPy
matrix equality was replaced with componentwise exact cancellation;
the derived formulas did not change.

The [validation receipt](validation/p8-covzero-source-constraints-2026-09-24.json)
records all fourteen focused replays, lint/format, local links, source
hashes and preservation checks. The 38 pre-existing non-status files
at the start of this admission check, including the completed seed
packet and unrelated P4/P9 changes, remain unchanged. Frozen P8 and
the original three tracked diagnostics are untouched. No full frozen
115,487-test replay is claimed. Work remains local, uncommitted and
unpushed; no background work continues after handoff.
