# The first omitted action on a physical-metric anisotropy

This S6.32 calculation concerns the **formal** order-six stationary action of
the unchanged VARIABLE parent. It is not a full covariant operator reduction,
the CD tensor perturbation equation, a controlled derivative expansion, or
an extraction of `Xi=A1-2F2_X`. The test metric below is deliberately off
shell. Every differentiation and integration by parts precedes evaluating
the coefficient functions at the clock center.

## 1. Fixed representative and clock

The frozen own-f construction uses signature `(+---)`, `EH=-M^2 R_B/2`,
`beta2=beta3=0`, and `r^3=-beta1/beta4>0`. Its metric correction is

\[
 f_0=r^2g,\qquad f_2=r^2 A\mathcal P,
 \qquad A=\frac{M^2r}{\beta _1},\quad
 \mathcal P_{\mu\nu}=R_{\mu\nu}[f_0]-\frac16f_{0\mu\nu}R_B[f_0].
\]

The indices of the displayed mixed `P` below are raised with the physical
`g`, not `f0`. To avoid confusion, `A` in this note is the inverse-potential
coefficient, **not** the DHOST function `A1`. The first omitted action is

\[
 S_{(6)}=\frac12 f_2 S_{\rm EH,f}'' f_2
       +\frac16S_{\rm pot}'''[f_2,f_2,f_2].
\]

The functional derivatives hold `g,theta,chi` fixed. This is equivalently
the coefficient of `epsilon^2` in the hidden Einstein action at
`f=f0+epsilon*f2`, plus the coefficient of `epsilon^3` in its potential.
The `f4` terms cancel from this action coefficient by leading stationarity;
this does not mean that `f4` or a remainder estimate is unnecessary for the
original equations.

Restrict the physical metric to

\[
 g=\operatorname{diag}(1,-e^{q(t)},-e^{-q(t)},-1),\qquad\sqrt{|g|}=1.
\]

Here `t` is physical-g proper time. The functions `r=r(theta(t))` and
`A=A(theta(t))` are arbitrary smooth compositions until the later center
restriction; write `ell=(log r)'`. The canonical free `chi` and its physical
metric are not changed. This metric has zero homogeneous volume expansion:
it is not obtained by setting the CD volume expansion to zero in its
equations. In particular, using the actual coefficient jets at `theta=T`
does not turn this off-shell unit-volume probe into the actual CD spacetime.

## 2. Literal coframes and both action boundaries

The diagonal mixed Schouten entries are

\[
 \begin{split}
 P_0&=-5q'^2/12-2\ell'+\ell^2,\\
 P_1&=-q''/2+q'^2/12-\ell q'-\ell^2,\\
 P_2&=+q''/2+q'^2/12+\ell q'-\ell^2,\\
 P_3&=q'^2/12-\ell^2.
 \end{split}
\]

These follow from the literal diagonal Ricci tensor and the conformal
Schouten relation, not an assumed tensor mass. Put `h_i=A P_i`. The exact
formal coframes are

\[
 N_f=r\sqrt{1+\epsilon h_0},\quad
 a_{f1}=r e^{q/2}\sqrt{1+\epsilon h_1},\quad
 a_{f2}=r e^{-q/2}\sqrt{1+\epsilon h_2},\quad
 a_{f3}=r\sqrt{1+\epsilon h_3}.
\]

For each fixed regular field jet a sufficiently small formal neighborhood
of `epsilon=0` has positive roots. This supplies no uniform interval in
`delta=c-2` for the full nonlinear stationary solution.

Let `b_i=ell+(q'/2,-q'/2,0)_i` for the three spatial directions. Expanding
the exact volume/lapse ratio gives `V_f/N_f=r^2(1+epsilon V1+epsilon^2 V2+...)`,
where

\[
 V_1=\tfrac12(\sum_{i=1}^3h_i-h_0),\qquad
 V_2=\tfrac18(\sum_{i=1}^3h_i-h_0)^2
       -\tfrac14(\sum_{i=1}^3h_i^2-h_0^2).
\]

The logarithmic rates are
`b_i+epsilon*h_i'/[2(1+epsilon*h_i)]`. Define

\[
 \begin{split}
 S_0&=3\ell^2-q'^2/4,\\
 S_1&=\tfrac12\sum_{i=1}^3h_i'(3\ell-b_i),\\
 S_2&=\tfrac14\sum_{1\le i<j\le3}h_i'h_j'
       -\tfrac12\sum_{i=1}^3h_i h_i'(3\ell-b_i).
 \end{split}
\]

The hidden Einstein coefficient, divided by `M^2`, is exactly
`-r^2(S2+V1 S1+V2 S0)` in the ADM representative. For arbitrary diagonal
scales, a direct curvature calculation gives

\[
 -\frac{M^2}{2}N_f V_f R_B[f]
 =-\frac{M^2V_f}{N_f}\sum_{i<j}(\log a_{fi})'(\log a_{fj})'
 +M^2\frac{d}{dt}\left[\frac{V_f}{N_f}\sum_i(\log a_{fi})'\right].
\]

Thus the order-two covariant-EH-to-ADM boundary, divided by `M^2`, is

\[
 F_{\rm EH,2}=r^2\left[3\ell V_2+\tfrac12V_1\sum_{i=1}^3h_i'
                         -\tfrac12\sum_{i=1}^3h_i h_i'\right].
\]

Using the **same** interaction normalization `-2 sum beta_n e_n`, its
third-order coefficient divided by `M^2` is

\[
 \frac{r^2}{24 A}
 \left[(\sum_{i=0}^3h_i)^3-6(\sum_{i=0}^3h_i)(\sum_{i=0}^3h_i^2)
                          +5\sum_{i=0}^3h_i^3\right].
\]

The relation `r beta1=M^2 r^2/A` was used, with no g/source equation.
The API retains the full Einstein and cubic expressions before extracting
their coefficient quadratic in `q`. The independently authored Fraction
test expands the literal square roots, coframes and potential; it does not
reuse this coefficient formula as its numerical oracle.

## 3. Quadratic action and full Euler operator

Write the quadratic result divided by `M^2` as
`L_raw=sum_{1<=i<=j<=3} a_ij q^(i)q^(j)`. A precise time normal form is

\[
 L_{\rm raw}=B_3q'''^2+B_2q''^2+B_1q'^2+F_{\rm time}',
\]
\[
 F_{\rm time}=\tfrac12a_{23}q''^2+a_{13}q'q''
                 +\tfrac12(a_{12}-a_{13}')q'^2.
\]

The full boundary relative to the covariant Einstein representative is
`M^2([F_EH,2]_(q^2)+F_time)`. Both terms are exposed in `derive()`. In
particular a term proportional to `q'q'''` produces a **second derivative
of its coefficient** in `B1`; it cannot be set to its center value first.

The exact normal coefficients are

\[
 B_3=\frac{r^2A^2}{16},\qquad
 B_2=-\frac{r^2A}{16}(7A\ell^2+2A\ell'+6A'\ell+A''),
\]
\[
 \begin{split}
 B_1=\frac{r^2}{16}[&23A^2\ell^4+4A^2\ell^2\ell'
 +8A^2\ell\ell''+8A^2\ell'^2\\
 &+24AA'\ell^3+24AA'\ell\ell'+2AA'\ell''
 +6AA''\ell^2+4AA''\ell'+2AA'''\ell\\
 &+8A'^2\ell^2+2A'^2\ell'+2A'A''\ell].
 \end{split}
\]

The contribution to the variational equation is

\[
 \mathcal E_{(6)}=M^2\{-2D_t(B_1q')+2D_t^2(B_2q'')-2D_t^3(B_3q''')\}.
\]

Its coefficients of `(q^(6),q^(5),q^(4),q''',q'',q')`, divided by `M^2`, are

\[
 (-2B_3,-6B_3',2B_2-6B_3'',4B_2'-2B_3''',2B_2''-2B_1,-2B_1').
\]

There is no `q` term; constant anisotropy is a spatial rescaling. The
highest symbol `-M^2 r^2 A^2/8` agrees with the constant-r own-f Schur
six-derivative symbol after the norm-two polarization dictionary for `q`.
This is a formal higher-operator coefficient, not a new resummed physical
mode or a ghost diagnosis.

## 4. Actual coefficient center, units and a nonzero operator limit

Only now pull back the unchanged coefficient functions with
`theta=T=tau*u`, `d=1+u^2`, `y=2/d^4`,

\[
 b_1=\frac{y^3h_u}{c(c-y)},\quad b_4=\frac{3h^2}{2c^2}-\frac{b_1}{y^3},
 \quad r^3=-b_1/b_4,\quad
 \bar A=A/\tau^2=r/b_1,\quad h=4u/d.
\]

The positive regular family has `2<c<=4`, `|u|<=1/10`; `c=2` is not an
action in this family. Limits refer to its cancelled rational coefficient
jets. Differentiating those functions, rather than assuming a constant
clock or lapse ratio, gives

\[
 r(0)=2,\quad \bar A(0)=\frac{c(c-2)}{16},\quad
 \bar A_{uu}(0)=\frac{13c^2-22c+8}{8},\quad
 (\log r)_{uu}(0)=-\frac{4(c+2)}c.
\]

Odd jets of `r` and `A` vanish at the center, but their even jets do not.
The positive canonical parent clock is the old `phi(theta)`; this scalar
reparametrization neither changes the physical metric nor varies the frozen
free `chi` coupling.

Put `delta=c-2`. If `B_n=tau^(2n-2) Bbar_n` before the overall `M^2`, the
center coefficients are

\[
 \bar B_3=\frac{c^2\delta^2}{1024},\quad
 \bar B_2=-\frac{c\delta(9c^2-22c+24)}{512},\quad
 \bar B_1=-\frac{\delta(c+2)(9c^2-22c+24)}{32}.
\]

All three tend to zero. This does **not** imply that their full Euler
operator tends to zero at the center. Write the physical coefficient of
`q^(n)(T)` as `M^2 tau^(n-2) e_n(c)`. All odd `e_n` and `e_0` vanish there,
whereas

\[
 \begin{split}
 e_6&=-c^2\delta^2/512,\\
 e_4&=-c\delta(75c^2-154c+120)/256,\\
 e_2&=-(349c^4-2220c^3+5780c^2-7872c+5056)/128\longrightarrow-2.
 \end{split}
\]

One way to see the last result without the large polynomial is to retain
the even jets in section 3. At the center,

\[
 2B_2''-2B_1=-\frac{r^2}{8}
 [A''^2+22AA''\ell'+26A^2\ell'^2+2A^2\ell'''+AA''''].
\]

As `delta->0+`, `A->0` but `A''->2`, so the square of the second inverse
coefficient derivative remains. The term is `2B2''`, not a nonvanishing
`B1`. Freezing the coefficients before differentiating would instead give
the incorrect limit zero; at `c=4` it even gives `+60` in place of `-104`.

There is an exact continuous lower bound on this **individual** Euler
contribution. In delta coordinates,

\[
 e_2=-2+\tfrac74\delta-\tfrac{209}{32}\delta^2
              -\tfrac{143}{32}\delta^3-\tfrac{349}{128}\delta^4.
\]

Consequently `e2<-793/400` for `0<delta<=1/100`. For a smooth compactly
supported test `q` equal to `u^2/2` near the center, this gives
`|E_(6)(0)|>(M^2/tau^2)*793/400`. Relative to `M_*^2=5 M^2`, its limiting
coefficient is `-2/5`. This is a local off-shell operator test on a specified
jet, not a mode/WKB estimate, a fixed-window evolution theorem, or a bound
on the complete omitted action.

## 5. Mandatory fourth-plus-sixth cancellation control

The same Schouten entries give the frozen four-derivative restriction

\[
 L_{(4)}^{(2)}=\kappa[q''^2/2+2\ell q'q''+\ell'q'^2]
 \simeq\kappa q''^2/2-\kappa'\ell q'^2,
 \quad \kappa=M^2r^2A/4.
\]

The exact additional boundary is `kappa ell q'^2`. Its full Euler
coefficient of `q''` is `kappa''+2kappa' ell`, which at the same center equals

\[
 M^2\frac{9c^2-22c+24}{8}\longrightarrow+2M^2.
\]

It cancels the order-one sixth-order contribution. The combined center
coefficient divided by `M^2` is exactly

\[
 \tfrac72\delta-\tfrac{173}{32}\delta^2
                    -\tfrac{143}{32}\delta^3-\tfrac{349}{128}\delta^4=O(\delta).
\]

Thus the individually non-small omitted term does **not** establish an
uncancelled order-one retained-plus-omitted mismatch. The opposite omission
controls are both retained: ignore coefficient derivatives and the S6 limit
is wrongly zero; ignore S4 and the combined limit is wrongly nonzero. The
observed cancellation is not an all-orders convergence result. Neither
`S_(8)` nor the exact nonlocal parent response has been bounded here.

## 6. What this does not extract

The Euler operator is unchanged by all recorded boundaries, so its tested
contribution is action-level evidence. It nevertheless is not the
normal-form covariant `F2_X` or `A1`, nor the full physical CD propagation
operator. Higher-curvature/curvature-derivative operators also contribute
to `q''` after their coefficient derivatives are taken. Assigning the
displayed low symbol to a DHOST function would discard precisely that
distinction.

For example the exact boundary

\[
 \nabla_\mu[W(\theta)R X\theta^\mu]
 =W_\theta R X^2+W X\theta^\mu\nabla_\mu R
  +W R X\Box\theta+2W R\theta^\mu\theta_{;\mu\nu}\theta^\nu
\]

shows why a clock-center coefficient projection is inadequate: `W=0` at
the center need not imply `W_theta=0`. A claimed covariant Xi extraction
requires a fixed independent higher-operator basis and its complete IBP
dictionary. This gate instead fixes the action representative, exposes
every time boundary, and tests the Euler operator without that identification.

No vacuum, all-scalar stability, finite-band source matching, rolling spectral
gap, loop bound, physical cutoff, ultraviolet verdict or P8 closure follows.
The next matching question is whether a specified full reduction and error
norm controls the complete series/nonlocal response on its actual physical
domain. The present calculation supplies both a failure of a naive small
coefficient argument and a concrete cancellation that such a computation
must reproduce.
