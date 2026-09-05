# S6.2 — Obstructions to the first matching routes

These are conditional classical/tree parent tests, not a UV-completion
classification. The explicit positive-lambda vacuum has a conventional
massive-mediator realization at leading order, but that route does not
extend perturbatively to the fixed D bounce under the hypotheses below.

## 1. Stable algebraic elimination has the wrong second X jet

Fix phi and write U=V-ZX/2. On a smooth stationary branch U_a=0, let
K_ab=U_ab be positive definite and j_a=Z_a/2. Differentiating stationarity
gives `K_ab*h_*X^b=j_a`. Since L_alg,a=0 on that branch, the envelope theorem
here is just the chain rule:

```
P_X=Z(phi,h_*)/2
P_XX=j^T K^(-1) j >= 0.
```

For any number of heavy scalars, positivity follows by setting
w=K^(-1)j: `j^T K^(-1)j=w^T K w>=0`. Zero is permitted when j=0. K needs
to be invertible at each point; uniform gap/derivative estimates would be
additional matching work, not a consequence of this identity.

The replay independently solves a two-heavy Gaussian model and a nonlinear
one-heavy example Z=z0+2jh+b h^2, V=V0+mu^2 h^2/2. In the latter,
`h_*=jX/(mu^2-bX)`, `K=mu^2-bX` and
`P_XX=j^2 mu^4/(mu^2-bX)^3`. The condition K>0 is essential.
The positive vacuum mediator is recovered from Z=Z0(1+g h):
`lambda=g^2/(8mu^2)` in the canonically normalized vacuum field.

For the pinned D action, direct differentiation instead gives

```
F_D,XX = -2*(45u^8+122u^6+96u^4-12u^2+5)/(1+u^2)^8 < 0,
45u^8+122u^6+96u^4-12u^2+5
  =45u^8+122u^6+96*(u^2-1/16)^2+37/8.
```

Thus its second X jet is negative at every finite u. The whole S6.1
extension equals that action on the tube, regardless of positive vacuum
lambda, Z or m. In particular F_D,XX(0,1)=-10. A matching decomposition
`F_target=P_alg+R` therefore needs `R_XX(0,1)<=-10`, so
`|R_XX(0,1)|>=10` in M=tau=1 units. In physical units the lower bound is
`10 M^2/tau^2`: at least the full magnitude of the target second jet.
Increasing M*tau does not reduce this relative mismatch. Exact matching,
or matching with error strictly smaller than that magnitude in this jet,
is impossible for the stipulated stable affine-X elimination branch.

This is **not** a positivity bound F_XX>=0 on DHOST. F by itself is not the
physical reduced scalar action. Heavy-derivative corrections, phi-dependent
heavy motion, basis/IBP changes or curvature mixing can feed the coefficient
and are included in R if this matching convention is used. They cannot be
called a perturbatively small remainder in this second-jet norm here.
Other power countings/parent classes have not been excluded.

### Physical heavy gap is not the stationary Hessian

The code retains a negative control to prevent an overstrong conclusion.
For canonical quadratic fields pi,h with mixing c*h*pi_dot and stationary
Hessian kappa, direct equations give

```
(omega^2-k^2)*(omega^2-k^2-kappa)-c^2*omega^2=0.
heavy_gap_squared(k=0)=kappa+c^2
c_low^2=kappa/(kappa+c^2).
```

kappa=-1, c^2=2 has a positive heavy frequency squared 1 but negative
stationary Hessian and a long-wavelength instability. For k^2=1/10, the
low root is `3/5-3*sqrt(5)/10<0`; both roots are positive for k^2>1.
It is not a healthy bounce or complete UV model. It shows why K>0 cannot
be replaced by a check of one heavy pole, nor inferred from positivity of
the two bare time-kinetic terms. A full spectrum and its domain are needed.

## 2. The minimal parent cannot bounce in the physical metric

The parent with constant M and positive field metric has

```
rho+p=G_AB*Phi_dot^A*Phi_dot^B >= 0
-2M^2 H_dot=rho+p.
```

The potential cancels from the null contraction; field-space curvature and
mixing cannot change the positive quadratic form. This is an exact statement
about homogeneous classical solutions, without integrating out a heavy mode.
At the selected bounce H_dot=2/tau^2, the required null stress is instead
`-4M^2/tau^2`. Any stress-matching remainder would need magnitude at least
`4M^2/tau^2` there. An exact minimal parent solution realizing the target
therefore does not exist. A bounce in a truncated description of that same
parent cannot be justified by a uniformly small stress-error estimate.

This does not address quantum expectation values of stress, nonminimal
gravity, higher derivatives or a different physical metric. It does not
assert that every unitary UV theory obeys classical NEC.

## 3. Quantifying the nonminimal escape

For `-Q(Phi)R/2`, Q>0, define T=G^J_AB*Phi_dot^A*Phi_dot^B. Do not assume
T>=0: a negative Jordan kinetic contribution can be offset by its conformal
mixing with gravity. After a time boundary, direct homogeneous variation gives

```
L = -3Qa*a_dot^2/N -3a^2*a_dot*Q_dot/N
    +a^3*T/(2N) -N*a^3*V
3QH^2+3H Q_dot-T/2-V=0
2QH_dot+3QH^2+2H Q_dot+Q_ddot+T/2-V=0
-2QH_dot=T+Q_ddot-H Q_dot.
```

The Ricci scalar used here is independently obtained from the flat-FLRW
metric with lapse, in P8's (+---) curvature convention. The Einstein-frame
field metric is
`G^E_AB=M^2 G^J_AB/Q+3M^2 Q_A Q_B/(2Q^2)` and is assumed positive definite.
Use q=Q/M^2, dT_E=sqrt(q)dt and a_E=sqrt(q)a. Direct differentiation gives

```
H_E=(H+q_dot/(2q))/sqrt(q)
dH_E/dT_E=(H_dot+q_ddot/(2q)-3q_dot^2/(4q^2)
                         -H*q_dot/(2q))/q <= 0.
```

The inequality follows from Einstein-frame NEC. It agrees identically with
the physical-frame field equations and the transformed kinetic energy,
including the Q_dot^2 term. No matter-frame equivalence is asserted: the
target scale factor is still the original physical a, not a_E.

At t=0 write q0=q, q1=tau*q_dot and q2=tau^2*q_ddot. The exact bounce needs

```
4*q0+q2-3*q1^2/(2*q0) <= 0.
```

If |q0-1|<=epsilon0<1, |q1|<=epsilon1 and |q2|<=epsilon2, the left side is
bounded below by

```
B=4*(1-epsilon0)-epsilon2-3*epsilon1^2/(2*(1-epsilon0)).
```

Indeed the expression with q2=-epsilon2 and q1^2=epsilon1^2 is increasing
in q0>0, with derivative `4+3*epsilon1^2/(2q0^2)>0`. Therefore **B>0 excludes
this near-constant conformal parent matching**. For a common error epsilon,

```
B=(2-epsilon)*(4-7epsilon)/(2*(1-epsilon)).
```

Any epsilon<4/7 is impossible. The necessary lower bound epsilon>=4/7 is
not sufficient, nor a constructed healthy solution at equality. Correlations
between Q and the fields can strengthen it. Merely permitting epsilon>=4/7
does not pass matching. This bound already allows an indefinite Jordan field
metric provided its Einstein counterpart is healthy; assuming a positive
Jordan metric would give a stronger but unnecessary restriction.

For example, ten-percent bounds on all three jets give B=209/60>0. Small
amplitude alone is not enough: arbitrarily rapid variation can evade that
error budget. If H_dot itself is matched only within
`|tau^2*H_dot-2|<=eta<2` with H(0)=0, replace 4 by `4-2eta`; the same argument
gives a sufficient exclusion test, not a claim for unrestricted geometry errors.

### Explicit nonminimal escape control

To check that this reasoning does not wrongly exclude every nonminimal
bounce, reconstruct a different background-only model with
`q=exp(-3u^2)`, `a=1+u^2`, phi=t, and define

```
P(u)=1+11u^2+15u^4+9u^6 > 0
H_E=-u*(1+3u^2)*exp(3u^2/2)/(tau*(1+u^2))
dH_E/dT_E=-exp(3u^2)*P(u)/(tau^2*(1+u^2)^2)
Z_E=2M^2*P(u)/(tau^2*(1+u^2)^2) > 0
V_E=M^2*(3H_E^2+dH_E/dT_E)
Q=M^2*q; G_J=q*Z_E-3M^2*q_dot^2/(2q); V_J=q^2*V_E.
```

Z_E is the Einstein field metric in the phi coordinate. Both Einstein
background equations, the scalar equation and the independently derived
physical-frame lapse/scale equations vanish exactly. The physical bounce
therefore coexists with a decreasing H_E and positive Einstein kinetic
energy. Its normalized Planck jets are (1,0,-6), outside the small-error
domain as required. Crucially q->0 in both tails: it violates P8's nonzero
Einstein tensor asymptotics and does not match the constant D Planck
coefficient. It is a countercheck, not a new accepted P8 witness; its
perturbations, vacuum and UV origin are not certified here.

## 4. What these failures do and do not decide

The conventional massive-scalar realization works as a vacuum matching
control, but the stable algebraic and minimally coupled routes fail on the
fixed bounce. A simple nonminimal scalar-curvature extension also cannot
match its constant Planck coefficient with perturbatively small C2 errors
in background units, under the healthy Einstein-frame parent hypotheses.

The S6.1 splice remains an algebraic family; these tests do not disprove
every possible parent for it. A viable next search must genuinely leave at
least one tested ansatz: leading derivative/curvature interactions, sizable
frame/threshold dynamics with an explicit matching map, quantum contributions,
or another witness. Each needs its own gap and observable analysis. No
finite-gravity positivity inequality, UV completion, M1 interaction result or
whole-row exclusion is earned here. The previous linear/tree results stand.

The [follow-up frame tests](frame-tests.md) further close the large-frame-
change escape within conventional Einstein/NEC parents once the required
tensor tails are imposed. They also derive the finite-time loss of
invertibility of any stipulated exact D-to-quartic-Horndeski map. Thus a
frame change alone is not an untested rescue for those seed classes.
