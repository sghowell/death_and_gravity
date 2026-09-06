# Exact composite-source branches and a local CD construction

All claims use the [formulation](../FORMULATION.md), including the changed
physical metric and shared matter coupling. Dots denote the arbitrary
coordinate derivative in sections 1–2. Sections 3–5 use the dimensionless
effective proper time `u=T/tau`; primes there mean `d/du`.

## 1. Vary first, retaining both lapses

Write `y=b/a`, `c=N_f/N_g`, `r=alpha+beta*y`, `s=alpha+beta*c`.
The homogeneous scalar action is

\[
L_\chi=\frac{a_e^3\dot\chi^2}{2N_e}-N_ea_e^3V,
\quad a_e=\alpha a+\beta b=ar,
\quad N_e=\alpha N_g+\beta N_f=N_gs.
\]

Let `rho=chidot^2/(2N_e^2)+V`, `p=chidot^2/(2N_e^2)-V`, and
`n=rho+p>=0`. The induced densities and pressures are

\[
\rho_g^m=\alpha r^3\rho,\quad p_g^m=\alpha sr^2p,
\qquad
\rho_f^m=\beta r^3\rho/y^3,\quad
p_f^m=\beta sr^2p/(cy^2). \tag{1}
\]

The different density/pressure weights follow from different lapse and
scale variations. Replacing either by a guessed separately conserved
perfect-fluid source would change the equations.

Define the interaction functions

\[
\begin{aligned}
U_g&=m4(\beta_0+3\beta_1y+3\beta_2y^2+\beta_3y^3),\\
U_f&=m4(\beta_4+3\beta_3/y+3\beta_2/y^2+\beta_1/y^3),\\
P&=m4(\beta_1+2\beta_2y+\beta_3y^2),\\
P_g&=-U_g+(y-c)P,\qquad
P_f=-U_f+(c-y)P/(cy^3).
\end{aligned}
\]

Literal variation of the lapse-retaining HR plus scalar Lagrangian gives

\[
\begin{aligned}
3GH_g^2&=U_g+\alpha r^3\rho,\\
3FH_f^2&=U_f+\beta r^3\rho/y^3,\\
G(2D_gH_g+3H_g^2)+P_g+\alpha sr^2p&=0,\\
F(2D_fH_f+3H_f^2)+P_f+\beta sr^2p/(cy^2)&=0.
\end{aligned}\tag{2}
\]

Here `D_g=N_g^{-1}d/dt`, `D_f=N_f^{-1}d/dt` and
`D_e=N_e^{-1}d/dt`. The scalar equation and its undivided energy identity
are

\[
E_\chi=D_e^2\chi+3H_eD_e\chi+V_\chi=0,
\quad C:=D_e\rho+3H_en=(D_e\chi)E_\chi.
\]

No division by `chidot` is required in deriving the background identities.
Define `B=N_g bdot-N_f adot` and `Q=P-alpha*beta*r^2*p`.
Off shell, the two total source balances are

\[
\begin{aligned}
D_g(U_g+\rho_g^m)+3H_g(U_g+P_g+\rho_g^m+p_g^m)
 &=\alpha r^3sC+\frac{3QB}{N_g^2a},\\
D_f(U_f+\rho_f^m)+3H_f(U_f+P_f+\rho_f^m+p_f^m)
 &=\frac{\beta r^3s}{cy^3}C-
   \frac{3QB}{N_gN_facy^3}.
\end{aligned}\tag{3}
\]

The lapse/scale Noether identities set these balances to zero on a
gravitational solution. With the scalar equation, (3) implies exactly
`QB=0`. Neither `P` alone nor a separately conserved induced matter
source gives the correct branch equation.

Subtracting each lapse equation from its scale equation gives

\[
\begin{aligned}
-2G D_gH_g&=\alpha r^3n+(y-c)Q,\\
-2F D_fH_f&=\beta r^3n/y^3+(c-y)Q/(cy^3),\\
-2G D_gH_g-2Fcy^3D_fH_f&=sr^3n\ge0.
\end{aligned}\tag{4}
\]

Equations (1)–(4) are checked by literal symbolic variation and a separate
Fraction Laurent engine. The latter includes the composite chain rule,
clears positive monomial denominators, substitutes the two affine sums,
and checks every polynomial coefficient. This is not interpolation.

## 2. Both branches and all nondegenerate bounce intersections

The kinematic error is

\[
S=H_f-H_g/y=\frac{B}{N_gN_fb},\qquad
H_e=H_g/r+\frac{\beta cy}{sr}S. \tag{5}
\]

On a connected open set where `B=0`, both `S` and its derivative vanish.
With `A=G+Fy^2` and `Z=H_g/sqrt(A)`, (4) then gives

\[
-2A^{3/2}D_gZ=sr^3n\ge0,\qquad H_e=H_g/r. \tag{6}
\]

Both lapse clocks preserve time orientation. Monotonicity of `Z` and
its positive relation to the sign of `H_e` exclude contraction-to-expansion
within that branch, including degenerate transitions. The symbolic audit
also retains the off-branch term `-2Fcy^3D_fS`; it is removed only on an
open branch, not at an isolated zero of `B`.

On `Q=0`, the individual equations instead give

\[
D_gH_g=-\frac{\alpha r^3n}{2G},\qquad
D_fH_f=-\frac{\beta r^3n}{2Fy^3}. \tag{7}
\]

They do not require `B=0` and do not force a composite Hubble sign.
In particular, the constant-polynomial root topology of S6.5 is not
applicable: `Q` now depends on the evolving scalar pressure.

For completeness consider any proposed *nondegenerate* effective bounce
with `H_e=0` and `B=0`. Equation (5) forces `H_g=H_f=0`. If `Q!=0`
at the point, continuity supplies an open dynamical neighborhood and

\[
D_eH_e=-\frac{r^2n}{2(G+Fy^2)}\le0.
\]

If `Q=0` at that point, differentiate the physical scale retaining both
lapses, then use (7). All terms involving first constituent velocities
vanish, yielding

\[
D_eH_e=-\frac{r^2n}{2s^2}
\left(\frac{\alpha^2}{G}+\frac{\beta^2c^2}{Fy^2}\right)\le0.\tag{8}
\]

These two cases exhaust `B=0`, including branch intersections and zero
pressure/polynomial specializations. A CD bounce has positive acceleration,
so necessarily `B!=0`. By continuity it has a whole pressure-branch
neighborhood `Q=0`, with opposite-sign constituent Hubbles. This does not
claim that arbitrary mixed branches cannot contain other degenerate
behavior; a universal composite no-bounce theorem is already false below.

The pressure equation is quadratic in `y`, but no leading coefficient,
discriminant, `p_y`, `H_e`, `Q` or `B` is divided out. The independent
free example crosses a pressure double root regularly.

## 3. An explicit free canonical bounce

Restore physical scales by `G=F=M^2`, `m4=M^2/tau^2`, `u=T/tau`,
`chi=M*varphi`, `rho=(M^2/tau^2)*rho_bar` and likewise for pressure.
Choose `alpha=beta=1`, `beta_n=(0,0,1,0,0)` and `V=0`. In this
section all Hubble variables `X,Y,h` are multiplied by `tau` and all
pressures are barred. Set

\[
\rho=p=\frac{2y}{(1+y)^2},\quad
X=\sqrt{\frac{y(2+5y)}3},\quad
Y=-\sqrt{\frac{5+2y}{3y^2}},\quad
D=(y+5)X-(5y+1)Y,
\]
\[
N_g=-\frac{(5y+1)Y}{D},\quad
N_f=\frac{(y+5)X}{D},\quad
y'=\frac{6y(1+y)XY}{D},\qquad y(0)=1.\tag{9}
\]

For every `y>0`, the roots and `D` are positive in the stated signs,
both lapses are positive, and `N_g+N_f=1`. The right side of (9) is
analytic near `y=1`, giving a regular local solution. Define

\[
a_e=\left[\frac{(1+y)^2}{4y}\right]^{1/6},\quad
a=\frac{a_e}{1+y},\quad b=ya,\quad
\varphi'=\frac{2\sqrt y}{1+y},\quad
h_e=\frac{(y-1)XY}{D}. \tag{10}
\]

Exact identities verify both Friedmann and both null equations in (2),(4),
the constituent scale kinematics, `a_e'=a_e h_e`, the pressure branch,
and the free scalar current `a_e^3 varphi'=1`. Equations (1)–(2) connect
these to the literal action variations; none of the four gravity equations
is dropped.

At `y=1`,

\[
N_g=N_f=\tfrac12,\quad y'=-\sqrt{7/3},\quad
\varphi'=1,\quad h_e=0,\quad h_e'=7/36>0.\tag{11}
\]

Thus physical `dH_e/dT=7/(36tau^2)` is strictly positive despite
positive scalar kinetic energy. Since `y'<0` and `XY<0`, (10) also
has the required negative-to-positive Hubble signs near the crossing.
At this instant the pressure factor is

\[
Q=2y-p(1+y)^2=-\tfrac12(y-1)^2\quad (p=1/2).
\]

The pressure has a maximum, `p_y(1)=0`, `p_yy(1)=-1/4`, while
`y'` is nonzero. Differentiating `h_e=-p_y y'/(6p)` gives
`h_e'(0)=-p_yy(1)y'(0)^2/(6p(1))=7/36`; the second audit checks
this last arithmetic using Fractions. A two-valued inverse `y(p)` would
obscure the smooth time-domain solution.

This free history is not CD even if one adjusts its positive time scale.
Keeping the actual nonlinear clock `d/du=y'(y)d/dy` gives
`h_e'''(0)=-73/216` and `a_e''''(0)=-97/432`. Hence

\[
a_e''''(0)-\tfrac32[a_e''(0)]^2=-9/32\ne0.
\]

Every normalized CD shape `(1+(u/k)^2)^2`, for any positive `k`, has
that combination equal to zero. Equivalently,
`h_e'''(0)+(3/2)[h_e'(0)]^2=-9/32` here and zero for CD. The
nonzero residual scales as `tau^-4` in physical time and cannot be
removed by any finite positive time scaling. This rejects this explicit
free example, not every possible composite free-scalar model.
No global hidden-curvature or EFT conclusion follows from (9).

## 4. A regular ODE that realizes CD exactly on a finite interval

Use the same gravitational parameters, but now construct a scalar potential.
Write `R=rho_bar`, `p(y)=2y/(1+y)^2`, `r=1+y`, and prescribe
`h(u)=4u/(1+u^2)`. Define

\[
X=\sqrt{y^2+r^3R/3},\quad
Y=-\sqrt{y^{-2}+r^3R/(3y^3)},\quad D=X-yY,
\]
\[
\boxed{
y'=\frac{yr[XY-h(X+Y)]}{D},\qquad
R'=-3h[R+p(y)],\qquad (y(0),R(0))=(1,1/2).
}\tag{12}
\]
\[
N_g=\frac{rh-yY}{D},\quad N_f=\frac{X-rh}{D},\quad
a_e=(1+u^2)^2,\quad a=a_e/r,\quad b=ya.\tag{13}
\]

The identities

\[
h-y'/r=N_gX,\qquad h+y'/(yr)=N_fY,\qquad N_g+N_f=1
\]

give the constituent and physical clocks/kinematics exactly. Both
Friedmann equations follow from the definitions of `X,Y`. With (12),
direct differentiation verifies both null equations as well. These
identities hold for an arbitrary symbolic prescribed value `h`; no
division by `h` occurs at the bounce. The prescribed CD function then
fixes the actual time dependence, not just a finite Taylor jet.

### Explicit positivity and continuation window

On the box

\[
|u|\le1/64,\quad 3/4\le y\le5/4,\quad 1/4\le R\le3/4,
\]

the following elementary bounds hold:

\[
\begin{gathered}
|h|\le1/16,\quad 7/4\le r\le9/4,\quad
0<p\le1/2,\quad R+p\ge59/108>1/2,\qquad R+p\le5/4,\\
775/768\le X^2\le1129/256<(9/4)^2,\quad
(-Y)^2\le535/108<(5/2)^2,\quad y(-Y)>1,\\
2<D\le43/8,\qquad
rh-yY>55/64,\quad X-rh>55/64.
\end{gathered}\tag{14}
\]

For `p<=1/2`, use the exact factorization
`1/2-p=(y-1)^2/[2(1+y)^2]`. The lower bound on `p` used in
`59/108` is the conservative `2(3/4)/(9/4)^2=8/27`.
The `X` bounds use the endpoints separately in `y^2+r^3R/3`.
For `-Y`, use `y^-2+(1+1/y)^3R/3`; and
`y^2Y^2=1+r^3R/(3y)>1` gives the useful weighted lower bound.

Consequently `N_g,N_f>55/344>1/7`; since they sum to one each is
also below `6/7`. The right sides of (12) obey

\[
|y'|\le\frac{17055}{2048}<9,\qquad
|R'|\le\frac{15}{64}<1/4.\tag{15}
\]

The first bound uses
`(5/4)(9/4)[(9/4)(5/2)+(1/16)(9/4+5/2)]/2`.
The second uses `3(1/16)(5/4)`. The vector field is analytic on an
open neighborhood of this compact box: both radicands and `D` are
strictly positive. Local ODE uniqueness/existence applies in both time
directions. Before any first exit and for `|u|<=1/64`, (15) gives

\[
|y-1|<9/64<1/4,\qquad |R-1/2|<1/256<1/4.
\]

Thus no exit is possible through a spatial face. Compactness and local
continuation extend the solution throughout the closed interval
`|u|<=1/64`, with strict positivity margins. The certificate computes
every rational margin in (14)–(15) and the no-exit estimates exactly.
The analytic ODE theorem and inequalities are the written existence
argument; a finite numerical integration is not used as its substitute.

### The single canonical potential and its scope

Define

\[
\varphi'=\sqrt{R+p}>0,\qquad
\overline V(\varphi(u))=\frac{R-p}{2},\qquad
\chi=M\varphi,\quad V(\chi)=\frac{M^2}{\tau^2}\overline V(\chi/M).
\tag{16}
\]

The positive derivative and analyticity give a unique inverse clock on
the traversed field interval and hence one analytic potential there.
The density and pressure obtained from (16) are exactly `R,p`.
Furthermore, with `n=R+p`,

\[
\varphi''+3h\varphi'+\overline V_{\varphi}
=\frac{R'+p_y y'}{2\sqrt n}+3h\sqrt n
 +\frac{R'-p_y y'}{2\sqrt n}=0.
\]

This also checks the field-space derivative: it is `V_u/varphi'`,
not `V_u` alone. No potential was fixed in advance. It need not be
nonnegative, and no unique extension outside this field interval is
asserted. The resulting local action/solution reproduces the CD physical
scale exactly on `|T|<=tau/64`, but is not the old free M1 action and
does not provide an all-time CD or quantum matching construction.

## 5. Evidence boundary

Every physical energy density in sections 3–4 carries `M^2/tau^2`,
every Hubble carries `tau^-1`, and each physical scalar velocity carries
`M/tau`; `G=F=M^2` and `m4=M^2/tau^2` therefore give dimensionally
consistent equations. Positive `M,tau` do not by themselves certify an
EFT cutoff hierarchy on the rolling scalar background.

The bare composite coupling generically has a BD-mode/cutoff limitation.
Absence of that mode in minisuperspace does not certify all propagating
scalar/vector/tensor modes, interactions or omitted quantum operators.
The [source audit](sources.md) discusses these distinctions and the
higher-source completion, which is not the action used in this proof.
This gate establishes background algebra and a local exact construction,
not full parent matching, positivity closure or completion of P8.
