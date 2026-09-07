# S6.20: variable-coupling local background and tensor audit

This is a new, bounded `P8-S6.20.VARIABLE` result. It neither edits nor extends
the hypotheses of the frozen constant-coefficient tree theorem S6.17. The
adopted S6 matching contract is unchanged; S6.19 is not a dependency.

## Specified action and physical frame

Use the P8(b) convention `+---`,
`R_B=-6(DH+2H^2)`, and `S_EH=-M^2 integral sqrt|g| R_B/2`.
Both Einstein coefficients are the same positive constant `M^2`. The actual
physical metric is `g`, not a conformal or composite metric. The action is

\[
S=-\frac{M^2}{2}\int\sqrt{|g|}R_B[g]
  -\frac{M^2}{2}\int\sqrt{|f|}R_B[f]
  +\frac12\int\sqrt{|g|}\bigl[(\partial\phi)^2_g+(\partial\chi)^2_g\bigr]
  -2\int\sqrt{|g|}\sum_{n=0}^4\beta_n(\phi)e_n(\sqrt{g^{-1}f}).
\]

The square root is the regular positive diagonal FLRW root. All coefficients
are genuine functions of a dynamical canonical clock, not prescribed functions
of coordinate time. In particular, the clock stress is interaction-sourced.
The additional scalar `chi` is genuinely free and minimally coupled to the
same physical `g`. An optional free canonical scalar on `f` can be constant;
no such field or independent source is needed for the fixture.

## Exact local family

Let `M,tau>0`, `u=T/tau`, `d=1+u^2`, where `T` is physical `g` proper time.
On `|u|<=1/10`, choose either `c=1` or any fixed `2<c<=4`, and set

\[
N_g=1,\quad N_f=c,\quad a=d^2,\quad b=2d^{-2},\quad
y=b/a=2d^{-4},\quad h=4u/d,
\]
\[
\bar n=2(y^3/c-1)h',\quad
\chi_{,T}=\frac{M}{10\tau a^3},\quad
\bar k=\bar n-\frac{1}{100a^6},\quad
\phi=M\int_0^u\sqrt{\bar k(v,c)}\,dv.
\]

A prime here means `d/du`. Define the dimensionless coefficient profiles

\[
b_1=\frac{y^3h'}{c(c-y)},\quad b_2=b_3=0,\quad
b_4=\frac{3h^2}{2c^2}-\frac{b_1}{y^3},\quad
b_0=\frac{3h^2-\bar n/2}{2}-3b_1y.
\]

The action coefficients are
`beta_n(phi)=(M^2/tau^2) b_n(u(phi/M),c)` on the resulting clock range.
The inverse exists analytically for each finite admitted `c`. Uniform strict
kinetic bounds are `kbar>3599/100` for `c=1` and `kbar>449/100` for `2<c<=4`.
These are canonical clock bounds, not a reduced coupled-scalar stability test.
There is no uniform coefficient bound as `c` approaches `2` from above.

The scale factor and free-`chi` trajectory agree locally with the normalized
CD/M1 trajectory. The interacting clock, parent action, extra metric and
physical perturbations are not identified with the original DHOST clock or
C/D operators. A background coincidence is not the original S6 matching map.

## Certified tensor result

Literal lapse/scale and scalar variations certify all background equations.
Literal elementary symmetric polynomials and a two-dimensional square root
certify the two-TT and transverse-shift quadratic coefficients. The complete
time-dependent canonical two-TT map retains normalization, moving weights,
the mixed derivative operator, and its time boundary.

For `c=1` the transverse-shift coefficient is negative and the formal
high-momentum relative-vector kinetic coefficient becomes negative. This is
not a computed below-cutoff growth rate. For `2<c<=4` that kinetic coefficient
is positive, but the `f` tensor principal cone is wider than physical `g`.
Neither fact alone decides the light-only finite-band S6 contract.

The prospective hierarchy `c=2+epsilon^2`, `u=epsilon x` has, on every fixed
compact `x` interval and fixed finite dimensionless momentum, the exact
canonical operator limit

\[
l_{,xx}=0,\qquad Q_{,xx}+\frac{80}{1+8x^2}Q=0.
\]

The second adiabatic ratios of the algebraic mass and the actual canonical
heavy diagonal tend to `-1/5` for mass squared and `-1/10` for its positive
square root. Thus the growing algebraic mass does not produce a vanishing
adiabatic parameter by this scaling alone. The finite-interval propagator
convergence asserted here is only in the shrinking physical window
`|T|<=epsilon tau X`; it is not a fixed-CD-window estimate.

## Verification and exclusions

The source-hashed certificate replays S6.17 and pins the adopted S6 contract.
Exact symbolic identities are supplemented by an independently implemented
Fraction/Taylor and coefficient engine, actual-to-primary bridges, and a
separately authored covariant audit. The real-analysis continuation argument
is written in `notes/proof.md`; it is not proof-assistant formalization.

No full scalar/vector health, finite-time rolling spectral gap, retarded
physical-source matching error, complete-CD solution, loop closure, cutoff,
UV completion, positivity theorem, or original C/D operator match is claimed.
The completed scoped photon objective and frozen linear classification are
unchanged. Original P8 remains open.
