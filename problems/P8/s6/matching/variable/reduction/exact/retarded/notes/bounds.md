# Finite-parameter coefficient bounds for the own-f retarded problem

This proof uses the **exact** own-f homogeneous background already established
in frozen S6.33. The physical metric is the external off-shell probe `g=eta`,
with the same prescribed clock `theta=T`; neither its field equation nor the
clock equation is imposed. The hidden homogeneous tensor response is a
linearized own-f response, not a full nonlinear parent solution.

The background interval calculation is inherited explicitly from
`p8_exact_stationary.intervals`. This child does not claim to reproduce that
background enclosure independently. It independently proves the domain
inclusion, desingularized gap identity and polynomial inequalities, and the
final monotone rational coefficient bounds. All arithmetic in `bounds.py`
is exact integer/Fraction arithmetic.

## 1. The inherited domain

Use the physical-g proper time and dimensionless clocks

\[
u=T/\tau,\qquad x=u/\sqrt\delta,\qquad v=u^2=\delta x^2,
\qquad c=2+\delta.
\]

The domain here is

\[
-\frac14\le x\le0,\qquad 0<\delta\le\frac1{625}.                    \tag{1}
\]

It is contained in the frozen branch domain because

\[
0\le v\le\frac1{625}\left(\frac14\right)^2
=\frac1{10000},\qquad \frac1{625}<\frac1{100}.                       \tag{2}
\]

Consequently, the entire required background lies in the explicit S6.33
rectangle. Its unique fixed point has `11<zeta<13`. We call the frozen
Fraction engine on the rectangular **enclosure**

\[
v\in[0,1/10000],\quad\delta\in[0,1/625],\quad\zeta\in[11,13].
                                                                         \tag{3}
\]

The exact relation `v=delta*x^2` is not erroneously assumed for arbitrary
points of (3); (3) is a superset used for rigorous enclosure. The returned
lapse is the full implicit lapse, including the fixed-point derivative,
not a partial-derivative substitution.

The zero-`delta` face is only the regularized limiting extension of S6.33.
The literal action still requires `delta>0`. In contrast, `x=0` is an
ordinary allowed observation point of every literal positive-`delta`
background. `denominator_ratio()` rejects `delta=0` unless its explicitly
named `extension=True` option is requested.

## 2. Exact gap identity with a polynomial positivity proof

Let `d=1+v`. The same coefficient denominator and numerator as S6.33 are

\[
D=c-2/d^4=\delta+2(1-d^{-4}),\qquad
Q_{\rm profile}=\frac{32(1-v)}{cd^{14}},\qquad b_1=Q_{\rm profile}/D.
\]

The profile numerator `Q_profile` is unrelated to the hidden tensor
amplitude `Q` used below. Define

\[
P(v)=4+6v+4v^2+v^3.
\]

Direct polynomial multiplication gives

\[
d^4-1=vP(v),\qquad
4d^4-P(v)=10v+20v^2+15v^3+4v^4.                                    \tag{4}
\]

For every `v>=0`, therefore, `0<P(v)/d^4<=4`. Substituting `v=delta*x^2`
and dividing only for literal `delta>0` gives the exact identity

\[
\frac D\delta
=1+\frac{2x^2(4+6v+4v^2+v^3)}{(1+v)^4}.                            \tag{5}
\]

The right side is regular at `delta=0` and defines its limiting extension;
it does not make the original ratio `0/0` a literal coefficient. From
(1), (4), and (5),

\[
1\le D/\delta\le1+8x^2\le\frac32.                                 \tag{6}
\]

The module checks the complete coefficient identities in (4), not merely
values at sample points. In particular, every nonconstant coefficient of
`4d^4-P` is positive. This supplies the continuum inequality in (6).

## 3. Physical action, time derivatives, and coefficient normalization

The literal norm-two homogeneous tensor coordinates are

\[
g={\rm diag}(1,-e^q,-e^{-q},-1),\qquad
f={\rm diag}(N^2,-b^2e^Q,-b^2e^{-Q},-b^2).
\]

The S6.33 own-f action variation gives the quadratic tensor density

\[
\mathscr L_2=\frac{M^2}{4}q_T^2+\frac{K_f}{4}Q_T^2
-\frac\nu4(Q-q)^2,\qquad K_f=M^2b^3/N,\qquad \nu=2\beta_1b,
                                                                         \tag{7}
\]

where `beta1=(M^2/tau^2)*b1`. Prescribing `q`, rather than varying it as a
second dynamical tensor, gives

\[
\partial_T(K_f Q_T)+\nu(Q-q)=0.                                    \tag{8}
\]

Since `T=tau*sqrt(delta)*x`,

\[
\partial_T=(\tau\sqrt\delta)^{-1}\partial_x,\qquad
K_f=M^2 k_\delta,\qquad
s_\delta=\frac{\delta\tau^2\nu}{M^2}
=\frac{2Q_{\rm profile}b}{D/\delta}.
\]

Multiplying (8) by `tau^2*delta/M^2` yields the exact finite-parameter
divergence-form equation

\[
\boxed{\ (k_\delta Q_x)_x+s_\delta(Q-q)=0,\qquad
k_\delta=b^3/N.\ }                                                   \tag{9}
\]

No factor of `M`, `tau`, or `sqrt(delta)` is left in these dimensionless
coefficients. Equivalently, the action after the change of integration
variable has the overall factor `M^2/(4*tau*sqrt(delta))` multiplying
`q_x^2+k_delta*Q_x^2-s_delta*(Q-q)^2`.

The kinetic derivative has not been dropped. For example,
`(k_delta)_x=2*delta*x*(k_delta)_v`, with the **total** branch derivative in
`v`. Equation (9) retains its entire contribution to `(k_delta Q_x)_x`.
The retarded Green proof can work with this divergence form using only
positive upper/lower coefficient bounds; it does not require an unproved
uniform estimate on `N_v`, `k_x`, or a canonical pump term. S6.33's analytic
branch does supply the ordinary smoothness needed to define the equation.
Zero data `Q=Q_x=0` at the left endpoint are equivalent to `Q=Q_T=0` for
each admitted positive `delta`.

## 4. Monotone exact rational inequalities

The frozen engine on (3) returns positive enclosures of `b`, `b_cubed`,
`N`, and `Q_profile`. Use its rational `b_cubed` directly in `k=b^3/N`,
rather than recubing the wider root bracket `199/100<b<201/100`. Positivity
proves the monotone endpoint bounds

\[
\begin{aligned}
\frac{(b^3)_{\min}}{N_{\max}}
&\le k_\delta\le\frac{(b^3)_{\max}}{N_{\min}},\\
\frac{2(Q_{\rm profile})_{\min}b_{\min}}{3/2}
&\le s_\delta\le2(Q_{\rm profile})_{\max}b_{\max}.
\end{aligned}                                                        \tag{10}
\]

A separate Fraction endpoint calculation, with integer floor/ceiling
rounding directed outward, gives

\[
\frac{987657}{250000}\le k_\delta\le\frac{1008683}{250000},\qquad
\frac{42355819}{1000000}\le s_\delta\le\frac{1608}{25}.               \tag{11}
\]

The following independent integer comparisons turn (11) into the simpler
strict theorem bounds:

\[
\begin{aligned}
5(987657)&>19(250000),&5(1008683)&<21(250000),\\
42355819&>42(1000000),&1608&<65(25).
\end{aligned}
\]

Thus, uniformly on the entire physical domain (1),

\[
\boxed{\quad 19/5<k_\delta<21/5,\qquad42<s_\delta<65.\quad}          \tag{12}
\]

These are continuous bounds for the exact finite-`delta` branch; they are
not inferred from its limiting differential equation. `checks()` verifies
strict positive rational margins, the polynomial identities, domain
inclusion, and the direction of each final rounding. `calibration()` exposes
`k_min,k_max,s_min,s_max,delta_max,width` for the separate Green-kernel proof.
No certificate file or frozen source is written by these functions.

## 5. Limiting interpretation and exclusions

Joint regularity and the center values of S6.33 give `k_delta -> 4` and
`s_delta -> 64/(1+8x^2)` on compact inner intervals. Dividing the limiting
own-f equation by four therefore gives coefficient `16/(1+8x^2)`. The
coefficient `80/(1+8x^2)` belongs to the coupled relative tensor equation,
which also varies the physical metric. It is not a substitute for (9).

The input `q` is a prescribed physical-metric tensor perturbation, not an
identified external matter stress or an arbitrary physical-g source loading.
A profile fixed in `x` has physical duration proportional to
`tau*sqrt(delta)`. The coefficient bounds make no fixed physical
low-frequency-band, full-parent background, general EFT, or UV claim.
