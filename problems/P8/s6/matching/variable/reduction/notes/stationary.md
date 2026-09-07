# Source-preserving stationary reduction of VARIABLE

This note proves the formal own-
`f` expansion used by the separately reviewed S6.30 candidate. It neither
establishes convergence of that expansion nor promotes the original S6
matching gate. The physical metric is always `g`, and the original free
canonical `chi` action is retained literally. No frozen ancestor is changed.

## 1. Action, domain and derivative order

With the P8(b) convention `(+---)`,

\[
 R_B=-6(\dot H+2H^2),\qquad
 S=-\frac{M^2}{2}\int(\sqrt{|g|}R_B[g]+\sqrt{|f|}R_B[f])
 +\int\sqrt{|g|}\left\{\frac{X+Y}{2}
 -2[\beta _0(\phi)+\beta _1(\phi)e_1(S)+\beta _4(\phi)e_4(S)]\right\},
 \quad S=\sqrt{g^{-1}f}.
\]

Here `X=(grad phi)^2_g`, `Y=(grad chi)^2_g`, the clock has mass dimension one,
and `beta_n` have mass dimension four. `M,tau>0` are constants. The two
canonical scalars couple only to `g`; no source is applied to `f`.

Hold `g,phi,chi` fixed when solving the **own-f** equation. Its derivative-zero
stationary point is

\[
 f^{(0)}=r(\phi)^2g,\qquad r^3=-\frac{\beta _1}{\beta _4}>0.
\]

This is not a requirement that `g` solve its metric equation, and is not a
vacuum construction. `r` is distinct from the actual FLRW ratio `y` and the
actual lapse `c`. Choose a positive regular square-root branch and smooth
nonzero opposite-sign `beta1,beta4`. The full traceful algebraic inverse below
requires `r beta1 != 0`; its existence does not establish a bounded inverse
of a differential operator.

For the actual positive-lapse family in this gate,

\[
 u=T/\tau,\quad d=1+u^2,\quad a=d^2,\quad y=2d^{-4},\quad
 h=4u/d,\quad
 b_1=\frac{y^3h_u}{c(c-y)},\quad
 b_4=\frac{3h^2}{2c^2}-\frac{b_1}{y^3},\quad
 \beta_n=\frac{M^2}{\tau^2}b_n,
 \qquad 2<c\le4,\quad |u|\le1/10.
\]

The continuous root proof does not rely on a sampled vacuum or an unpinned
descendant. On this domain `19/10<y<=2`, `h_u>3`, and `c-y>0`. Put

\[
 \theta=\frac{3h^2(c-y)}{2ch_u}
 =\frac32\frac{4u^2}{1-u^2}\frac{c-y}{c}.
\]

Then `0<=theta<7/220<1`, since `4u^2/(1-u^2)<=4/99` and
`(c-y)/c<21/40`. Therefore `b1>0`,
`b4=-(b1/y^3)(1-theta)<0`, and
`r^3=y^3/(1-theta)` is positive analytic. The exact rational margins also give
`19/10<r<21/10`. Analytic inversion of the original positive canonical clock
is available for each finite admitted `c`. No action at `c=2` is defined by
the limiting formulas used below.

Derivative counting is formal: a light-field derivative has degree one;
curvature and a clock Hessian have degree two. Coefficient functions of
`phi` have degree zero. This is not a claim that `X` is small on the physical
clock tube, or that formal derivative order bounds an operator evaluated on
the rolling background. A scalar-only reparametrization cannot turn a missing
quantitative hierarchy into one.

## 2. Literal full Hessian and its inverse

Write `f=r^2(g+h)`, `H=g^{-1}h`, and `b=r beta1`. Expanding the actual matrix
square root and volume gives

\[
 \operatorname{tr}\sqrt{I+H}=4+\tfrac12[H]-\tfrac18[H^2]
 +\tfrac1{16}[H^3]+O(H^4),
\]
\[
 \sqrt{\det(I+H)}=1+\tfrac12[H]
 +\tfrac18([H]^2-2[H^2])
 +\tfrac1{48}([H]^3-6[H][H^2]+8[H^3])+O(H^4).
\]

Using `beta4*r^4=-b` only after these expansions, the first-order potential
vanishes and the quadratic/cubic densities divided by `sqrt|g|` are

\[
 L_{\rm pot,2}=-\frac b4([H^2]-[H]^2),\qquad
 L_{\rm pot,3}=\frac b{24}([H]^3-6[H][H^2]+5[H^3]).
\]

`potential_hessian()` uses all ten independent **covariant** components of
`h`, not only TT components. In an orthonormal Lorentz frame its Hessian has
determinant `3 b^10/16`. Equivalently, the map
`h -> h-g tr_g h` has eigenvalue one on the nine tracefree components and
minus three on the trace. The indefinite algebraic trace form is not a
physical propagating-ghost calculation.

Because the Einstein action is `-M^2 R_B/2`, its covariant-metric variation
has positive Einstein-tensor sign. The term linear in `h` from the hidden
Einstein action is

\[
 L_{f,2}^{(1)}=\frac{M^2r^2}{2}B^{\mu\nu}h_{\mu\nu},
 \qquad B_{\mu\nu}=G_{\mu\nu}[f^{(0)}],
\]

where both contractions and raising on the displayed `B` use `g`. Therefore

\[
 h^{(2)}_{\mu\nu}=\frac{M^2r}{\beta _1}
 \left(B_{\mu\nu}-\frac13g_{\mu\nu}\operatorname{tr}_g B\right),
\]
\[
 \boxed{\ f^{(2)}_{\mu\nu}=\frac{M^2r^3}{\beta _1}
 \left(R_{\mu\nu}[f^{(0)}]-\frac16 f^{(0)}_{\mu\nu}R_B[f^{(0)}]\right).\ }
\]

The tensor in parentheses is twice one common Schouten convention. The
trace inverse is `1/3`, not the massless source-projector factor `1/2`.
These statements are independently checked by a literal Fraction determinant
and all-component stationary equations.

## 3. Clock derivatives and the formal action

Let `s=(log r)_phi`, `t=s_phi`. With physical-g derivatives,

\[
 B_{\mu\nu}=G_{\mu\nu}[g]-2s\phi_{;\mu\nu}
 +2(s^2-t)\phi_\mu\phi_\nu
 +g_{\mu\nu}[2s\Box\phi+(2t+s^2)X].
\]

In particular the Schouten-type tensor in section 2 equals

\[
 R_{\mu\nu}[g]-\frac16g_{\mu\nu}R_B[g]
 -2\nabla_\mu\nabla_\nu\log r
 +2\partial_\mu\log r\partial_\nu\log r
 -g_{\mu\nu}(\partial\log r)^2.
\]

No `r_phi` or `r_phiphi` term has been dropped. The conformal scalar curvature
is `R_B[f0]=r^-2[R_B[g]-6 r^-1 Box r]`. After a physical-g integration by
parts, the zero-plus-two-derivative density is

\[
 L_{0+2}=-\frac{M^2(1+r^2)}2R_B
 +\left(\frac12-3M^2r_\phi^2\right)X+\frac Y2
 -2(\beta _0+3r\beta _1).
\]

The negative conformal clock contribution follows from the stated B-sector
curvature sign. It is not obtained from the other P8 track's curvature
dictionary. The literal two-lapse FLRW check retains the boundary
`3 M^2 d_T(a^3 r r_T)` and verifies the sign independently.

Stationarity makes the order-four value one half of the linear term:

\[
 \boxed{\ L_4=\kappa(\phi)
 \left(B_{\mu\nu}B^{\mu\nu}-\tfrac13(\operatorname{tr}_gB)^2\right),
 \qquad \kappa=\frac{M^4r^3}{4\beta _1}.\ }
\]

Every density here is divided by `sqrt|g|`. Equivalently its invariant measure
is `sqrt|f0| kappa (Ric[f0]^2-R_B[f0]^2/3)`.

For a useful IBP check put `w=d log r`, `H=Hess log r`, `h=Box log r` and
`v=w^2`. Defining `Q_g=Ric[g]^2-R_B[g]^2/3`, the exact conformal identity is

\[
 Q_{\rm hat}-Q_g=-4\nabla_\mu V^\mu,\qquad
 V^\mu=G^{\mu\nu}w_\nu+w^\mu h-H^{\mu\nu}w_\nu+v w^\mu.
\]

It uses the contracted Bianchi identity and
`div(w h-H.w)=h^2-H:H-Ric(w,w)`. Hence, modulo a recorded local boundary,

\[
 L_4\simeq\kappa Q_g+4\kappa_\phi\{
 sG^{\mu\nu}\phi_\mu\phi_\nu
 +s^2[X\Box\phi-\phi^\mu\phi_{;\mu\nu}\phi^\nu]+s^3X^2\}.
\]

This compact expression agrees with the separately implemented coefficient
dictionary. The scalar-tensor portion is quartic Horndeski at this formal
order (`A1=2 F2_X`, `A3=0`); the independent `kappa Q_g` must still be kept.
Although `Q_g=(C^2-E4)/2` in four dimensions, `kappa(phi) E4` is **not** a
boundary for variable `kappa`. No free-matter contact is generated while
`g` is kept literally; a derivative redefinition of `g` is a different
source dictionary and must account for the changed matter action.

## 4. Independent constant-r source calibration

With `r,beta1` constant, use one TT polarization of norm one and

\[
 S_T=-\frac18[(G D+\nu)h_g^2-2\nu h_gh_f+(F D+\nu)h_f^2]
 +\frac14 jh_g,\quad
 G=M^2,\ F=M^2r^2,\ \nu=2r\beta _1.
\]

Here `D` is the flat wave operator. The interaction's factor two is part of
the literal VARIABLE normalization. Define the positive canonical stress by
`delta S_chi=(1/2) integral sqrt|g| T_mn delta g^mn`, equivalently
`delta S_chi=-(1/2) integral sqrt|g| T^mn delta g_mn`. Therefore
`delta g_ij=-h_g e_ij`, with `e:e=1`, gives
`delta S_source=+Pi_TT h_g/2` and **`j=+2 Pi_TT`**, not `j=Pi_TT`.
The own-f equation gives

\[
 h_f=\frac{\nu}{\nu+F D}h_g,\qquad
 K_g=(G+F)D-\frac{F^2D^2}{\nu}
 +\frac{F^3D^3}{\nu(\nu+F D)}.
\]

Thus `c_C=F^2/(4 nu)=kappa/2` exactly, matching the pinned S6.3 physical-source
Schur result. This is an off-shell constant-coefficient calibration, not an
assertion that the frozen variable action has the required flat vacuum.
The first omitted flat quadratic symbol is `+F^3 D^3/nu^2` in the kernel,
or `-F^3/(2 nu^2)` multiplying `Ric Box Ric-R Box R/3` in a covariant
six-derivative action representative. Curvature-cubic terms and all
variable-clock completions are not determined by this flat quadratic check.
Solving the g equation for f instead is not this elimination. True massive
eigenmode integration also requires its induced source/source contact and
the physical-g observable map; it cannot silently replace this source chart.

## 5. Exact center values are not a uniform derivative hierarchy

For the actual background, let `ell=(log r)_u`. With `kappa_bar=kappa/(M^2 tau^2)`,

\[
 f^{(2)}_{TT}=4\bar\kappa[-2(h_u+\ell_u)-h^2+\ell^2],\qquad
 \frac{f^{(2)}_{ii}}{a^2}=4\bar\kappa(h+\ell)^2,
 \qquad \bar\kappa=-\frac1{4b_4}.
\]

The symbol `TT` in this section means the **time-time covariant component**,
not a transverse-traceless perturbation. The code calls it `f2_00`.
At `u=0`, `r=2`, `r_u=0`, `r_uu=-8-16/c`. Put `delta=c-2`. Then

\[
 f^{(0)}_{TT}=4,\quad f^{(2)}_{TT}=4\delta,\quad
 f^{\rm actual}_{TT}=c^2,\quad
 \mathcal R_{TT}:=f^{\rm actual}_{TT}-f^{(0)}_{TT}-f^{(2)}_{TT}
 =\delta^2\quad (u=0).
\]

The conformal stationary lapse is **2**, whereas the actual lapse is **c**.
The small center value does not control derivatives. Exact identities are

\[
 (f^{(2)}_{TT})_{uu}(0)
 =-\frac4c(27c^3-116c^2+196c-176)\longrightarrow64,
\]
\[
 (\mathcal R_{TT})_{uu}(0)
 =\frac{4(c-2)}c(27c^2-62c+80)\longrightarrow0,
\]
\[
 (\mathcal R_{TT})_{uuuu}(0)
 =\frac{48}{c^2}(559c^4-2648c^3+5452c^2-6560c+4448)
 \longrightarrow10752.
\]

Indeed `c^2[(R_TT)_uuuu-10752]` equals
`22272 delta+132288 delta^2+87552 delta^3+26832 delta^4`, positive throughout
the admitted positive branch. For `0<delta<=1/100`, the first displayed
second derivative is greater than 60 by an exact positive polynomial margin.

To track the actual canonical clock, put `varphi=phi/M`. Its center jets are

\[
 k_0=(\varphi_u)^2=\frac{6400-801c}{100c},\qquad
 k_{uu}(0)=\frac{1206}{25}-\frac{1920}{c}<0.
\]

For these even profiles, differentiating the **same physical-T component as
a scalar-valued function of the clock label** gives

\[
 (f^{(2)}_{TT})_{\varphi\varphi}(0)=\frac{(f^{(2)}_{TT})_{uu}(0)}{k_0}
 \longrightarrow\frac{6400}{2399},
\]
\[
 (\mathcal R_{TT})_{\varphi^4}(0)
 =\frac{(\mathcal R_{TT})_{u^4}}{k_0^2}
 -\frac{2(\mathcal R_{TT})_{uu}k_{uu}}{k_0^3}
 >\frac{107520000}{5755201}.
\]

The last strict inequality holds for every `2<c<=4`: `k0<2399/100`,
`R_uu>0`, and `k_uu<0`. These are **not** the tensor components
`f_varphi,varphi` after changing the spacetime time coordinate; that change
would introduce the extra factor `tau^2/kbar`. Derivatives with respect to
the dimensionful field `phi` carry `M^-2` or `M^-4`; time derivatives with
respect to `T` carry `tau^-2` or `tau^-4`.

There is therefore no C2-small first metric correction and no C4-small
retained metric defect in this prescribed physical chart as `delta->0+`.
Also `kappa_bar(0)=c delta/16->0` but
`kappa_bar,uu(0)=(9c^2-22c+24)/8->2`. A small coefficient value alone is not
a coefficient-and-derivatives bound.

On each **fixed compact inner interval** `u=sqrt(delta) x`, the cancelled
profile expressions are analytic near `(u^2,delta)=(0,0)`, and Taylor's theorem
gives

\[
 \mathcal R_{TT}=\delta^2(1+64x^2+448x^4)+O_L(\delta^3).
\]

This is consistent with, rather than an exception to, the derivative floor.
The local own-f principal inverse parameter is
`alpha_bar=r/(2b1)=2 kappa_bar/r^2`; its center jets are
`alpha_bar=c delta/32` and `alpha_bar,uu=(13c^2-22c+8)/16->1`. Thus its inner
leading term is `delta(1+8x^2)/16`, and `alpha_bar d_u^2` becomes
`(1+8x^2)d_x^2/16+O_L(delta)`: it does not tend to the zero operator by
making delta small. This is not a proof that every bounded-band inverse
fails; different function spaces, fixed-frequency bands or a controlled
nonadiabatic inverse require their own estimates.

## 6. What is omitted, and what matching still requires

Let `S0` be the algebraic potential and `S2` the hidden Einstein action, with
derivatives in the following formula taken at fixed `g,phi,chi`. The first
omitted action term is structurally

\[
 S_6=\tfrac12\,\delta f^{(2)}S_2''\delta f^{(2)}
 +\tfrac16S_0'''[\delta f^{(2)},\delta f^{(2)},\delta f^{(2)}].
\]

Terms involving `f4` cancel in this action coefficient by
`S0'' f2+S2'=0`; derivatives acting on `f2`, `r` and the variable coefficients
in `S2''` must not be deleted. It contains curvature-derivative/cubic and
clock-jet operators, not merely another number multiplying `C^2`. The
Fraction and generic stationary-series controls do not compute or bound
this full operator.

The separate physical-frame dictionary measures the CD defect by
`I=X(A1-2 F2_X)/G_T`. The retained scalar-tensor action gives `I=0`; at the
CD center the frozen witness gives `I=1` on its positive-X tube. Therefore
an action/coefficient matching assertion needs an explicitly nonzero
omitted-symbol contribution in that seminorm. In particular, at normalized
`X=1`, exact matching requires
`|Delta A1|+2|Delta F2_X| >= M_CD^2` after the stated Planck/clock dictionary.
This is a necessary remainder floor in that specified comparison, not a
claim that the actual remainder is bounded by a small number or that no
more general finite-band light EFT can match observables.

Although `f0+f2` suffices for the action through degree four, its own-f
equation has a degree-four residual. The chain rule for the clock equation
contains that residual multiplied by `partial f0/partial phi`; merely
substituting `f0+f2` into the original partial clock equation is not an
order-six full-system residual theorem. Such a claim would need the next
stationary solution and the appropriate source/boundary estimates.

An exact own-f elimination would solve a differential equation with specified
boundary or initial/homogeneous data. An algebraic Hessian inverse does not
certify that functional inverse, nor does it select a state. A retarded
solution map cannot automatically be treated as a symmetric single-copy
action. All these issues, loop corrections, finite-band matching errors,
the physical source map and the adopted V/G/B gates remain distinct.

The source/boundary distinction is already explicit in Hassan, Schmidt-May
and von Strauss, [arXiv:1303.6940, section 2.1, equations (7)–(9)](https://arxiv.org/html/1303.6940#S2.SS1).
Their principal higher-curvature construction instead solves the g equation
algebraically; its coefficients are not imported as our own-f/source-preserving
answer. The present factor calibration uses the repository's independently
checked S6.3 Schur calculation. This note is a formal local expansion and
a quantitative failure of a proposed delta-uniform jet hierarchy, not a
claim of novelty, a parent ghost, a quantum cone verdict or P8 closure.
