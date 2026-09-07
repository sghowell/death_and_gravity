# Exact rational uniform branch enclosure

This note proves the continuous quantitative step for S6.33. Its inputs are
the same named VARIABLE coefficient profiles, composed with the original
clock, on the **off-shell** physical probe `g=eta`, `theta=T`. It does not
impose the physical-metric or clock equations. The proof constructs an exact
own-`f` homogeneous background, not a tensor Green inverse or a controlled
local effective action.

The independent engine is `src/p8_exact_stationary/intervals.py`. It imports
only standard-library exact rational arithmetic; it does not call SymPy,
Arb, the primary symbolic implementation, or a frozen ancestor. All intervals
below enclose the entire stated box, not a grid of sample points.

## 1. Regular variables and the closed proof box

Let `u=T/tau`, `v=u^2`, `c=2+delta`, and set

\[
0\le v\le10^{-4},\qquad 0\le\delta\le10^{-2},\qquad
11\le\zeta\le13.                                                    \tag{1}
\]

Define the following rational functions, with derivatives holding `c` fixed:

\[
\begin{aligned}
d&=1+v,& J&=cd^4(1-7v)+12v,\\
Q&=\frac{32(1-v)}{cd^{14}},&
R&=\frac{d^8J}{8c(1-v)},&D&=c-\frac2{d^4},\\
L&=R_v/R,&H&=D(J_v/J-6/d)-D_v,\\
F&=\frac23(v\zeta H-L),&
A&=1-vD\zeta,&b&=(A/R)^{1/3},\\
\mathcal T&=\frac{3bF^2}{2Q}.&&&
\end{aligned}                                                        \tag{2}
\]

The cube root is the positive real root. The letter `D` in this note is the
regular coefficient denominator in (2), **not** the derivative combination
`b1_v+b4_v*b^3` used in the unregularized equations. The symbol `R` in (2)
is an algebraic profile, not a spacetime curvature.

For the literal action, `beta_n=(M^2/tau^2)*b_n` and

\[
b_1=Q/D,\qquad b_4=-QR/D.                                           \tag{3}
\]

Equation (3) follows directly from the frozen `profile()` formulas, not from
a change of the action. On the admitted physical family `delta>0`, `D>=delta`
on (1), so these coefficients remain finite. At `delta=v=0`, however,
`D=0` and `Q>0`. The literal action is undefined there. Including that face in
the **regularized calculation** proves a uniform limiting enclosure; it
does not supply a new `c=2` physical theory.

## 2. Derivatives used by the independent engine

All required derivatives are evaluated by rational formulas:

\[
\begin{aligned}
J_v&=cd^3(-3-35v)+12,&J_{vv}&=cd^2(-44-140v),\\
Q_v/Q&=-1/(1-v)-14/d,&D_v&=8/d^5,&D_{vv}&=-40/d^6,\\
L&=8/d+J_v/J+1/(1-v),\\
L_v&=-8/d^2+J_{vv}/J-(J_v/J)^2+1/(1-v)^2,\\
H_v&=D_v(J_v/J-6/d)
 +D[J_{vv}/J-(J_v/J)^2+6/d^2]-D_{vv},\\
F_v&=\tfrac23[\zeta(H+vH_v)-L_v],&F_\zeta&=\tfrac23vH,\\
b_v&=\tfrac b3[-\zeta(D+vD_v)/A-L],&b_\zeta&=-bvD/(3A),\\
\mathcal T_v&=\frac3{2Q}
 [b_vF^2+2bFF_v-bF^2Q_v/Q],\\
\mathcal T_\zeta&=\frac3{2Q}[b_\zeta F^2+2bFF_\zeta].
\end{aligned}                                                        \tag{4}
\]

Here `b_v` and `T_v` are partial derivatives before solving for `zeta`.
On a fixed-point branch, implicit differentiation gives

\[
\zeta_v=\frac{\mathcal T_v}{1-\mathcal T_\zeta},\qquad
\frac{db}{dv}=b_v+b_\zeta\frac{\mathcal T_v}{1-\mathcal T_\zeta},\qquad
N=\frac2F\left(b_v+b_\zeta\frac{\mathcal T_v}{1-\mathcal T_\zeta}\right).
                                                                         \tag{5}
\]

The last term in (5) is retained. It is nonzero, for example, in the interior
of the positive-`v`, positive-`delta` box. A partial derivative in its place
would not recover the exact lapse.

## 3. Continuous rational arithmetic proof

For intervals `[a,b]` and `[c,d]`, sums use endpoint sums, products use the
minimum and maximum of the four endpoint products, and reciprocals use
`[1/b,1/a]` only when zero is excluded. Integer powers include zero as the
lower endpoint for an even power of a zero-crossing interval. Every endpoint
is an integer or `fractions.Fraction`; bool, binary float, infinity, NaN,
and implicit string conversions are rejected. Input sub-boxes must be
contained in (1).

First, rational evaluation proves positivity of `J,Q,R,A`. For the only
nonrational operation, evaluation of `A/R` proves

\[
\frac{7945004}{10^6}\le b^3\le\frac{8045632}{10^6},\qquad
(199/100)^3<\frac{7945004}{10^6},\qquad
\frac{8045632}{10^6}<(201/100)^3.                                     \tag{6}
\]

Monotonicity of the positive cube root therefore supplies
`199/100 < b < 201/100`. No approximate root evaluation enters the proof.
Substitution of this enclosure into (4)-(5), including all dependencies as
interval inclusions, gives the following outward rational bounds:

| Quantity | Lower endpoint | Upper endpoint |
| --- | --- | --- |
| `T` | `11864487/1000000` | `12164087/1000000` |
| `T_zeta` | `-5/1000000` | `1627/1000000` |
| `N` | `1981609/1000000` | `2042137/1000000` |

`calibration()` additionally reports bounds on the positive denominators,
the negative `F`, the negative total derivative `db/dv`, and each remaining
quantity needed for replay. To obtain concise output, a lower endpoint is
rounded **down** and an upper endpoint **up** to an integer multiple of
`10^-6`. The rounding itself is integer floor/ceiling arithmetic, also for
negative values. The engine checks that each concise interval contains its
underlying exact interval and that all advertised strict margins are
positive rational numbers.

Consequently, throughout (1),

\[
11<\mathcal T<13,\qquad |\mathcal T_\zeta|<1/100,\qquad
19/10<b,N<21/10,\qquad F<0,\qquad db/dv<0.                           \tag{7}
\]

The derivative bounds here are continuous-box estimates. They are neither
finite differences nor an inference from values at the endpoints.

## 4. Banach existence and analytic patching

For each `(v,delta)` in the closed parameter rectangle, (7) makes
`T(v,delta,.)` a strict self-map of the complete interval `[11,13]` and a
contraction with the uniform constant `1/100`. Banach's theorem therefore
gives a unique fixed point in that interval. Starting anywhere in the
interval gives the ordinary geometric contraction error bound; a numerical
iteration is not needed to prove existence.

All rational denominators in (2)-(4) stay nonzero on the box, and `A/R` stays
strictly positive. These functions, including its positive cube root, are
real analytic in a neighborhood of each point of the box. Moreover,
`1-T_zeta>0`. The analytic implicit-function theorem applies to
`zeta-T(v,delta,zeta)=0` at every fixed point. These local branches agree on
overlaps by the already-proved fixed-point uniqueness. They therefore
define one joint real-analytic branch on the rectangle, locally through its
boundary as well. Equation (5) is its actual derivative, so (7) proves a
positive lapse for this branch, rather than for an independent arbitrary
choice of lapse.

Composing with `v=u^2` gives a two-sided even analytic metric for
`|u|<=1/100` and every literal `0<delta<=1/100`. The joint desingularized
extension at `delta=0` is mathematical limiting information only.

## 5. Relation to the exact own-f Euler equations

For clarity the equivalence is recorded without suppressing the central
turning point. With the common dimensionless factor removed, the literal
own-f minisuperspace action is

\[
\mathscr L_f=-3bb_u^2/N-2b_1N-6b_1b-2b_4Nb^3.
\]

Its lapse and spatial Euler expressions are

\[
\begin{aligned}
C&=3bb_u^2/N^2-2b_1-2b_4b^3,\\
E&=3b_u^2/N+6bb_{uu}/N-6bb_uN_u/N^2-6b_1-6b_4Nb^2.
\end{aligned}
\]

Direct differentiation gives the undivided identity

\[
C'-(b_u/N)E=-2[b_{1,u}+b_{4,u}b^3-3b_1b_u/N].                       \tag{8}
\]

Holding `b` fixed in the coefficient derivatives, (2)-(3) imply

\[
b_{1,v}+b_{4,v}b^3=\tfrac32b_1F,
\qquad 1-Rb^3=vD\zeta.                                             \tag{9}
\]

These formulas make the original algebraic consistency equation exactly

\[
\frac{2vb}{3b_1^2}(b_{1,v}+b_{4,v}b^3)^2-b_1-b_4b^3
 =vQ(\mathcal T-\zeta).                                             \tag{10}
\]

On the constructed fixed point, (5) and (9) give the consistency bracket in
(8) equal to zero. The lapse equation is `C=2vQ(T-zeta)=0`. For nonzero `u`,
the last inequality in (7) makes `b_u=2u db/dv` nonzero. Equation (8) then
forces `E=0`. All expressions are smooth for each literal `delta>0`, so the
spatial equation extends to `u=0` by continuity. No division by `u` or
`b_u` is used to assert the equation at the center.

As exact center checks, at `v=0` one obtains

\[
\zeta_0=\frac{3(c+2)^2}{2c},\qquad b_0=2,\qquad
b_v(0)=-c(c+2),\qquad N_0=c^2/2.                                   \tag{11}
\]

The fixed-point uniqueness is a statement within the stated `zeta` box.
The separate regular fixed-`c` Euler/IFT argument identifies the local
solution through `b(0)=2` with this branch without assuming parity first.
No classification of other algebraic roots or remote singular branches is
claimed. In particular, possible simultaneous zero-derivative branches
outside this box have not been removed by declaration.

## 6. Verification and scope

`checks()` returns Boolean proof checks; every entry must be true.
`calibration()` returns exact enclosures and margins; `report()` only encodes
them as JSON-compatible rational strings and writes no files. Tests include
separate Fraction dual-number checks of the derivative formulas, rejection
of enlarged domains and inexact inputs, the genuine implicit-lapse term,
and an explicit distinction between the limiting and literal families.

This uniform background construction does not choose tensor homogeneous
data, a retarded/advanced/Feynman prescription, or a Green-function norm.
Nor does it prove the physical-`g`/clock field equations, a full parent or CD
background solution, a uniform derivative expansion, or a UV conclusion.
