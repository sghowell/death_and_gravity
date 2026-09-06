# Proof of monotonicity and the sharp CD endpoint obstruction

This proof uses exactly the action, geometry and separate classical NEC
matter hypotheses in the [formulation](../FORMULATION.md). Define
`n_i=rho_i+p_i>=0`. Every time derivative and sign statement refers to
the fixed physical `g` metric; `Z` below is a diagnostic function, not a
new metric frame or a canonical perturbation variable.

## 1. Full arbitrary-lapse identity

The pinned lapse-first equations imply, at all times,

\[
-2G D_gH_g=n_g+(y-c)P,\qquad
-2F D_fH_f=n_f+\frac{c-y}{cy^3}P,
\]

where `D_g=N_g^{-1}d/dt`, `D_f=N_f^{-1}d/dt=D_g/c` and
`P=m4(beta1+2 beta2 y+beta3 y^2)`. The undivided Bianchi equation is

\[
P(y)B=0,\qquad B=N_g\dot b-N_f\dot a.
\]

The exact kinematic error is

\[
S:=H_f-H_g/y=\frac{B}{N_gN_fb}.
\]

Write `A=G+Fy^2>0`, `Z=H_g/sqrt(A)`. Before imposing a branch,

\[
\begin{split}
W_{\rm full}&=-2G D_gH_g-2Fcy^3D_fH_f,\\
W_{\rm dyn}&=-2A D_gH_g+2FyH_gD_g y,\\
W_{\rm full}-W_{\rm dyn}&=-2Fcy^3D_f S,\\
W_{\rm dyn}&=-2A^{3/2}D_gZ.
\end{split}\tag{1}
\]

The source combination is `N=n_g+cy^3n_f>=0`, where
`cy^3=N_fb^3/(N_ga^3)` is positive. If `E` denotes the sum of the
two normalized acceleration-minus-lapse equations weighted by `1,cy^3`,
the exact off-shell identity is

\[
2A^{3/2}D_gZ+N=E-2Fcy^3D_fS. \tag{2}
\]

On a parent solution `E=0`. Equation (2) shows explicitly which
differentiated branch error must also vanish before inferring the sign
of `D_g Z`. It is not permissible to drop `D_fS` solely because `P`
vanishes at a point.

On the open set `U={t:P(y(t))!=0}`, the Bianchi equation gives `B=0`
on a neighborhood of every point. Thus `S` vanishes identically there
and so does its derivative. Equation (2) becomes

\[
\boxed{D_gZ=-\frac{n_g+cy^3n_f}{2A^{3/2}}\le0
       \quad\hbox{on }U.}\tag{3}
\]

In particular no assumption that `H_g` is nonzero is needed. The
arbitrary-lapse and ratio terms in (1) are checked symbolically and by a
separate exact Fraction Laurent calculation using only the pinned
polynomial primitive. The nonzero omission polynomials detect dropping
`D_g y` or the factor `c`.

## 2. Covering every root point

First suppose `P` is not identically zero as a polynomial in `y`.
It has finitely many real roots (at most two here), though the time-dependent
root set can have arbitrarily many components. Let `R=I\U` on the
connected regular time interval `I`. Continuity of `y` makes `R`
relatively closed.

Every point of the interior of `R` has a connected open neighborhood
on which `P(y(t))=0`. A continuous map from that connected neighborhood
into the finite root set has constant image. Therefore `y` is locally
constant and `D_g y=0` there. The `g` null equation, without using the
`f` dynamics, gives

\[
\boxed{D_gZ=\frac{D_gH_g}{\sqrt A}
       =-\frac{n_g}{2G\sqrt A}\le0
       \quad\hbox{on }\operatorname{int}_I R.}\tag{4}
\]

Every remaining point of `R` is a relative boundary point and is
approached by points of `U`. The assumed smooth positive metric and
lapse functions make `H_g,y,Z` continuously differentiable and `D_gZ`
continuous. Passing to the limit in the inequality (3) gives
`D_gZ<=0` at every such boundary point.

This covers `I=U∪int_I R∪boundary_I R`. It assumes neither finitely
many branch switches nor a simple root, and does not ignore a boundary
set of positive measure. An isolated root, an accumulated root, a root
interval and a branch-switching point all obey the same derivative sign.
No integration over a discarded exceptional set is involved.

Since `N_g>0`, `dZ/dt=N_gD_gZ<=0`. The mean-value theorem on any
compact subinterval of `I` now gives `Z(t2)<=Z(t1)` for `t2>t1`.
No uniform lower/upper bound on all-time lapses or scale factors is
needed beyond regularity at every point and continuity on each compact
subinterval. This conclusion does not extend across a singular endpoint.

## 3. The identically zero polynomial is different

If `beta1=beta2=beta3=0`, the interaction null stress in each metric is
zero everywhere, regardless of `y(t)`. Its zero set is no longer finite,
so the locally constant-ratio argument must not be used. Instead,

\[
D_gH_g=-\frac{n_g}{2G}\le0. \tag{5}
\]

Thus `H_g` itself is nonincreasing. This case includes two independent
cosmological terms `beta0,beta4` and zero interaction.

There is an exact control showing why replacing (5) by `D_gZ<=0` would
be false. In reference mass units equal to one, take

\[
G=F=m4=1,\quad \beta=(3,0,0,0,3),\quad
a=e^t,\ b=e^{-t},\ N_g=N_f=1,\quad n_g=n_f=\rho_g=\rho_f=0.
\]

All four parent lapse/scale equations hold. Yet

\[
H_g=1,\quad y=e^{-2t},\quad
Z=\frac1{\sqrt{1+e^{-4t}}},\qquad
\dot Z=\frac{2e^{-4t}}{(1+e^{-4t})^{3/2}}>0.
\]

The control is an exact pair of decoupled de Sitter metrics; it is not
a failed numerical approximation. It does not contradict (5).

## 4. Sign preservation and all degeneracies

For nonzero polynomial `P`, `Z` has the same sign as `H_g` because
`sqrt(A)>0`. Its monotonicity therefore implies

\[
H_g(t_1)\le0\ \Longrightarrow\ H_g(t_2)\le0
       \quad(t_2>t_1).
\]

If the first sign is strict, the second is also strict, since
`Z(t2)<=Z(t1)<0`. The same implication follows directly from (5) for
the zero polynomial.

There can consequently be no physical contraction-to-expansion transition
within one connected regular interval. The argument does not refer to
`Hdot_g` at a crossing, the order of a zero, or an analytic expansion.
It covers an arbitrarily flat zero, a zero plateau and any other proposed
degenerate route from a negative to a positive Hubble value. Turnaround
from expansion to contraction is not excluded by this sign statement.

## 5. Sharp endpoint mismatch for CD

In physical `g` cosmic time, CD has

\[
H_{CD}(-\tau/2)=-h_*,\quad H_{CD}(\tau/2)=h_*,
\qquad h_*=\frac8{5\tau}>0.
\]

If both endpoint parent errors were strictly smaller than `h_*`,
the earlier parent Hubble value would be negative and the later one
positive. Section 4 rules this out, with no hypothesis on `Hdot`
or on the errors at intermediate times. Hence

\[
\boxed{\max\{|H_{\rm parent}(-\tau/2)+h_*|,
             |H_{\rm parent}(\tau/2)-h_*|\}\ge h_*.}\tag{6}
\]

The bound is sharp across the stipulated class. The old regular beta1
Minkowski member `a=b=N_g=N_f=1`, zero matter and
`beta=(-3,1,0,0,-1)` has `H_g=0` and satisfies all four equations.
Both endpoint errors then equal `h_*`. No inferred derivative
accuracy or large-mass limit is required for (6).

## 6. Interacting control: H itself need not be monotone

One must not inflate the theorem to `D_gH_g<=0` on every interacting
branch. A local exact background supplies a counterexample. Again in
reference mass units, set

\[
G=F=m4=1,\quad \beta=(-2,5/4,0,0,7/4),\quad y=1+t,
\]
\[
H=\sqrt{\frac{7y^2+5/y}{12}},\quad c=y+\frac1H,\quad
\rho_g=3H^2+2-\frac{15y}{4},\quad
n_g=-\frac{\rho_g'}{3H}
    =\frac{15-14y+5/y^2}{12H},\quad
\rho_f=p_f=0.
\]

Here prime means `d/dy=d/dt`. Let
`a(t)=exp(integral_0^t H(1+s) ds)`, `b=ya`, `N_g=1` and `N_f=c`.
Then `bdot=c adot` exactly. Both Friedmann and both null equations
hold identically, and the code also substitutes the exact jets into
all four pinned Euler expressions.

On `9/10<=y<=11/10`, `H,c,y` are positive and smooth. The numerator
of `n_g` is at least

\[
15-14(11/10)+5/(11/10)^2=\frac{2258}{605}>0.
\]

Thus this is a local strictly NEC-respecting `g` matter configuration,
not a forbidden-NEC jet. It admits a positive-kinetic canonical
reconstruction: define `phi_dot=sqrt(n_g)>0`, invert the resulting
monotone `phi(t)` locally, and set `V(phi(t))=rho_g-n_g/2`.
The exact source conservation equation and `phi_dot!=0` then imply
the canonical scalar equation. An additional free M1 field may stay
constant. This constructs a local background, not a vacuum, global
completion, or perturbative health result.

At `y=1`,

\[
H=1,\quad \dot H=3/8>0,\quad n_g=1/2>0,\quad
\dot Z=-\sqrt2/16<0.
\]

This control shows that the ratio derivative in (1) matters even on an
actual NEC solution.

## 7. Units and evidence boundary

`G,F,A` have mass dimension two, `H` dimension one, `Z` dimension
zero, and `n_i` dimension four; both sides of (3) have consistent
dimensions. For the exact controls a reference mass `mu0` was set to
one. Restore `G=F=mu0^2`, `m4=mu0^4`, `t_dimensionless=mu0*t`,
`H=mu0*h` and `n_i=mu0^4*nbar_i`. Then `Hdot` scales as `mu0^2`
and `Zdot` as `mu0`. Equation (6) is already in physical units.

The certificate replays the pinned input and exact new identities.
Coefficientwise Fraction checks and independent covariant/control tests
support the algebra; the closed-set, mean-value and canonical-reconstruction
steps are explicit written arguments rather than formalized PDE results.
The source audit makes no global novelty assertion.

No arbitrary-beta vacuum or perturbation-health theorem, parent matching
construction, radiative closure, interaction cutoff or UV verdict follows.
The hypotheses exclude non-flat or non-common/non-bidiagonal geometries,
other or singular square-root branches, singular-endpoint continuation,
time-dependent interaction coefficients, quantum/NEC-violating matter,
nonminimal/derivative/shared matter and changed physical frames.
General S6 matching and P8 therefore remain open.
