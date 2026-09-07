# General HR potentials: a retained-order theorem and a genuine next-order freedom

This is the general stationary-potential component of S6.31.
It does not alter frozen S6.30 or supply a controlled functional inverse.
The degree-six curvature table below is a **potential-only** identity.
The correlated deformation in section 6, by contrast, determines an exact
**difference of the full formal degree-six action coefficients**. Neither
statement is a finite-band or finite-amplitude approximation theorem.

## 1. Action, physical frame, and the root being used

Use the P8(b) convention `+---`, `R_B=-6(Hdot+2H^2)`, and

\[
 S=-\frac{M_G^2}{2}\int\sqrt{|g|}R_B[g]
   -\frac{M_F^2}{2}\int\sqrt{|f|}R_B[f]
   +\frac12\int\sqrt{|g|}(X+Y)
   -2\int\sqrt{|g|}\sum_{n=0}^4\beta_n(\phi)e_n(\sqrt{g^{-1}f}).
\]

Here `X=(dphi)^2_g`, `Y=(dchi)^2_g`, both Einstein coefficients are finite
positive constants, and all beta functions are smooth on a clock interval.
The genuinely free canonical chi and the canonical clock are coupled to
the same physical g. There are no derivative interactions or independent
sources for f. A clock-only potential can be included in beta0. No metric
or matter redefinition is made. In code `MG2` and `MF2` already mean the
**squared** Planck masses; for example `MF2**2` means `M_F^4`.

The specified branch is a positive smooth proportional stationary root
`f0=r(phi)^2 g` of the *algebraic own-f potential equation*. Define

\[
 F_{\rm root}=\beta_1+3r\beta_2+3r^2\beta_3+r^3\beta_4=0,
 \quad P=\beta_1+2r\beta_2+r^2\beta_3,
 \quad Z=\beta_1+3r\beta_2+2r^2\beta_3.
\]

`F_root` is not the DHOST coefficient `F2`. It is also not the full
proportional-vacuum equation, which additionally uses the g equation,
curvature, and scalar conditions. The actual rolling lapse and scale-factor
ratio need not equal r. The theorem assumes a simple root. At a root,

\[
 \partial_r F_{\rm root}=-3P/r.
\]

Thus simplicity is exactly `P!=0`; the implicit-function theorem then gives
a local smooth root. Derivatives of r must be retained. For example
`r_phi=r (partial_phi F_root)/(3P)`, with the partial derivative taken at
fixed r. No sign of P follows from positive Einstein coefficients alone.
The inverse APIs accept either sign, without declaring either fixture a
healthy parent spectrum.

The exact rational evaluator is only an arithmetic interface: the symbolic
theorem is not restricted to rational coefficients. `polynomial_at` also
accepts off-root and P=0 data, without inversion. `inverse_at` additionally
requires `F_root=0`, `P!=0`, and positive r, MG2, MF2. Neither routine
accepts binary floats, booleans, nonfinite values, or unresolved symbols.

## 2. Literal full Hessian and cubic

Set `f=r^2(g+h)`, `H=g^{-1}h`, and let `[H^j]` denote mixed traces. These
are all ten covariant symmetric h components, not a tensor-only ansatz.
At fixed g and phi, expand `sqrt(I+epsilon H)` through cubic order. With
`C_j=binomial(4-j,n-j)` and out-of-range binomials zero,

\[
\begin{split}
 e_n(r\sqrt{I+\epsilon H})=r^n\{&\binom4n
 +\tfrac12\epsilon C_1[H]
 +\tfrac18\epsilon^2(C_2[H]^2-(C_1+C_2)[H^2])\\
 &+\tfrac1{48}\epsilon^3(C_3[H]^3-3(C_2+C_3)[H][H^2]
 +(3(C_1+C_2)+2C_3)[H^3])\}+O(\epsilon^4).
\end{split}
\]

Substitution into the literal `-2 sum beta_n e_n` gives the linear term
`-r F_root[H]`. Only after this expansion impose the root equation. The
constant, quadratic and cubic potential densities, divided by `sqrt|g|`, are

\[
 V_0=-2(\beta_0+3r\beta_1+3r^2\beta_2+r^3\beta_3),
 \qquad V_2=-\frac{rP}{4}([H^2]-[H]^2),
\]
\[
 V_3=\frac r{24}\{Z[H]^3-3(P+Z)[H][H^2]+(3P+2Z)[H^3]\}.
\]

The Hessian determinant in the ten independent lower-index entries at a
Lorentz frame is `3(rP)^10/16`. Its invariant trace-reversal map has the
usual `1/3` trace inverse, not `1/2`. The relative numerical factors follow
from the actual potential and Einstein normalizations. They are checked by
an independent Fraction matrix square-root series and all principal minors,
including shifts and unequal diagonal perturbations.

The P=0 exception is real, not a small-denominator convention. At r=1,
`(beta1,beta2,beta3,beta4)=(1,-2/3,1/3,0)` gives
`F_root=(r-1)^2`, `P=0`, `Z=-1/3`. The whole quadratic Hessian vanishes,
while `H=diag(1,1,1,0)` gives `V3=-1/12`. This is not a regular inverse
branch and is not classified for health. When beta1 through beta4 all vanish,
f is instead algebraically disconnected and cannot be eliminated this way.

## 3. Own-f inverse and the complete retained scalar normal form

Let `B_mu nu=G_mu nu[f0]`; its indices in the following contractions are
raised with g, not f0. The own-f Einstein action varies by
`+(M_F^2 r^2/2) B:h`. Completing its quadratic potential gives

\[
 h_2=\frac{M_F^2r}{P}\left(B-\frac13g\,\mathrm{Tr}_gB\right),
 \qquad \delta f^{(2)}=r^2h_2,
\]
\[
 L_4=\kappa\left(B:B-\frac13(\mathrm{Tr}_gB)^2\right),
 \qquad \kappa=\frac{M_F^4r^3}{4P}.
\]

Equivalently
`h2=(M_F^2 r/P)(Ric[f0]-f0 R_B[f0]/6)`. This is twice one common Schouten
normalization. h2 is the **relative** correction; the additive metric
correction contains another factor r squared. No g, phi or chi equation
has been used in this elimination. Coefficient derivatives enter through
the full own-f Einstein source, not a frozen constant-r substitute.

Write `s=(log r)_phi`, `t=s_phi`, `H_phi=Hess(phi)`, `v=dphi`,
`T=Box phi`, and `vHv=v^mu H_phi,mu nu v^nu`. The leading density is

\[
 L_{0+2}=-\frac{M_G^2+M_F^2r^2}{2}R_B
 +\left(\frac12-3M_F^2r_\phi^2\right)X+\frac Y2+V_0.
\]

The minus sign of the Jordan-frame clock correction follows by a weighted
conformal-curvature IBP. It is not itself a reduced-scalar stability test.
The full conformal tensor inserted in L4 is

\[
 B=G[g]-2sH_\phi+2(s^2-t)vv
       +g[2sT+(2t+s^2)X].
\]

The exact conformal identity for `Q_g=Ric[g]^2-R_B[g]^2/3` is
`Q_hat-Q_g=-4 div V`, where
`V^mu=G^{mu nu}w_nu+w^mu Box log r-(Hess log r)^{mu nu}w_nu+w^2 w^mu`
and `w=d log r`. Thus the weighted normal form retains `kappa Q_g` and
replaces its conformal remainder by

\[
 4\kappa_\phi\{sG(v,v)+s^2(XT-vHv)+s^3X^2\}.
\]

For example `j Ric(v,v) ~= j(L2-L1)+j_phi(XT-vHv)` and
`m vHv ~= -(m/2)XT-(m_phi/2)X^2` when j,m depend only on phi. These follow
from explicit divergences and do not use any field equation. The resulting
quadratic-Hessian normal form, with `kappa Q_g` **separately retained**, is

\[
 F_2=-\frac{M_G^2+M_F^2r^2}{2}-2s\kappa_\phi X,
 \quad A_1=-4s\kappa_\phi=2F_{2X},\quad A_2=-A_1,
 \quad A_3=A_4=A_5=0,
\]
\[
 K=6[(s^2+t)\kappa_\phi+s\kappa_{\phi\phi}]X,
 \quad F=V_0+(1/2-3M_F^2r_\phi^2)X+f_4X^2,
\]
\[
 f_4=2\{\kappa_\phi(2s^3+2st+s_{\phi\phi})
 +\kappa_{\phi\phi}(s^2+2t)+s\kappa_{\phi\phi\phi}\}.
\]

Constant kappa removes all these extra scalar fourth-order terms even if
r varies. Variable kappa does not: `kappa times Euler density` is not a
disposable boundary. The full retained action is **not** purely quartic
Horndeski because its independent curvature-square operator is still there.

Consequently `Xi=A1-2F2_X=0` for the retained scalar-tensor part, for every
admitted choice of the five smooth beta functions. Scalar-only clock changes
multiply Xi by the nonzero squared clock Jacobian, preserving its zero.
Where `G_T!=0`, `I=X Xi/G_T` is correspondingly zero in this specified
normal form. The frozen CD center has `Xi=M_CD^2` and `I=1` on its positive
X tube. Hence any claimed exact match in that same coefficient dictionary
requires `|Delta A1|+2|Delta F2_X|>=M_CD^2`; tolerance eta<1 in normalized
Xi requires at least `(1-eta)M_CD^2`. This is a necessary omitted-coefficient
floor, not a proved upper bound on the actual remainder.

## 4. All curvature-linear terms in the potential cubic

This section is only `V3[h2]`, not the full formal degree-six action. Put

\[
 \alpha=M_F^2r/P,\quad N=r\alpha^3/8,\quad
 A=\mathrm{Ric}[g]-gR_B/6,\quad
 D=-2sH_\phi+2(s^2-t)vv-gs^2X,
 \quad h_2=\alpha(A+D).
\]

With `a=Tr_g A`, `d=Tr_g D`, `d2=Tr_g(D^2)`, the coefficient linear in
physical curvature obtained from the literal Newton cubic is

\[
 \frac{L_{3,R}}N=Zad^2-(P+Z)\{ad_2+2d\,\mathrm{Tr}_g(AD)\}
 +(3P+2Z)\mathrm{Tr}_g(AD^2).
\]

Let `G=G[g]`, `L1=H_phi:H_phi`. The following is one clean **algebraic**
operator basis for this sector, before any curvature-Hessian IBP. A square
on H inside a contraction means mixed matrix multiplication.

| Operator | Coefficient in L3,R / N |
|---|---:|
| G:H_phi^2 | 4s^2(3P+2Z) |
| T G:H_phi | -8s^2(P+Z) |
| X G:H_phi | 4s[P s^2-2(P+Z)t] |
| G(v,H_phi v) | -8s(3P+2Z)(s^2-t) |
| T G(v,v) | 8s(P+Z)(s^2-t) |
| X G(v,v) | 4(s^2-t)[2(P+Z)s^2-Pt] |
| R_B(L1-T^2) | (4/3)(2P+Z)s^2 |
| R_B vHv | -(8/3)(2P+Z)s(s^2-t) |
| R_B X T | (4/3)s[(P+2Z)s^2-2(2P+Z)t] |
| R_B X^2 | (P+2Z)s^4-4P s^2t |

The module retains the first five and the three R-Hessian operators as a
named eight-operator complement. They are not deleted on a homogeneous
background, or because some coefficient has a zero at a point. Their full
four-dimensional Lorentzian matrix contraction is checked exactly against
the raw cubic expression.

At a point with s=0 the **raw density projection** is

\[
 L_{3,R}=C XG(v,v),\qquad C=4NPt^2=r\alpha^3Pt^2/2.
\]

Z cancels in this statement. It is not a normal-form theorem about the
full center coefficient: a coefficient vanishing at the point can have
nonzero derivatives in an IBP.

## 5. What a coefficient projection does and does not mean

If the eight-operator complement above is fixed and retained, one may reduce
just the two Hessian-free entries. For any `j(phi,X)`,

\[
 j\,\mathrm{Ric}(v,v)\simeq j(L_2-L_1)
 +j_\phi(XT-vHv)+2j_X(TvHv-vH_\phi^2v).
\]

In particular `C(phi) XG(v,v)` contributes
`F2=-C X^2/2`, `A1=-C X`, `A2=C X`, `A3=2C`, `A4=-2C`, and hence
`Xi=C X`, together with its lower-Hessian terms. Applying precisely this
projection to the two entries in section 4 gives

\[
 \Xi_{\rm projection}=4N\{P(s^4+s^2t+t^2)-2Zs^2t\}X.
\]

This is not a scalar-tensor invariant of the *full* higher-operator action.
For example Bianchi gives the exact relation

\[
 a(\phi)XG:H_\phi\simeq-a_\phi XG(v,v)-2aG(v,H_\phi v).
\]

For the actual table coefficient
`a=4Ns[P s^2-2(P+Z)t]`, s=0 gives `a=0` but
`a_phi=-8N(P+Z)t^2`. Eliminating that one operator while retaining the
changed complement changes the *assigned* center Xi from `4NPt^2 X` to
`4N(3P+2Z)t^2 X`. The action has not changed. Thus taking s=0 before IBP
does not determine a full center coefficient. The curvature-free cubic
Hessian operators, curvature-quadratic/cubic terms, and the Einstein second
variation also cannot silently be discarded when discussing full S6.

The fourth-order theorem in section 3 remains true. Sixth-order operators
may contribute to a declared matching remainder, but an isolated projection
does not show that they supply the required physical CD invariant, preserve
degeneracy, or yield a smaller subsequent remainder. In particular the
frozen CD inverse-X coefficient changes naive derivative homogeneity;
polynomial formal grading is not automatically an approximation on X near1.

## 6. A genuine new cubic freedom, not just basis ambiguity

Fix the smooth functions r,P and any smooth function eta(phi). The exact
coefficient deformation

\[
 (\Delta\beta_0,\Delta\beta_1,\Delta\beta_2,\Delta\beta_3,\Delta\beta_4)
 =(3r\eta,-2\eta,\eta/r,0,-\eta/r^3)
\]

preserves the root equation, P and V0 as **functions**, while `Z` increases
by eta. Hence r and all its derivatives, kappa and all its derivatives,
the complete degree<=4 action, and h2 are unchanged. This is a comparison
of named parent coefficients, not permission to change a frozen parent.

For clarity, write `S0` for the algebraic potential and `S2` for the own-f
Einstein action, with variations at fixed g,phi,chi. The degree-six
coefficient of stationary elimination is

\[
 S_6=\tfrac12 S_2''[f_2,f_2]+\tfrac16 S_0'''[f_2,f_2,f_2].
\]

Terms with f4 cancel because `S0'' f2+S2'=0`. The deformation preserves
S2, S0'', and f2, so its full formal degree-six change is exactly

\[
 \Delta L_6=\Delta V_3[h_2]
 =\frac{r\eta}{24}\{[h_2]^3-3[h_2][h_2^2]+2[h_2^3]\}
 =\frac{r\eta}{4}e_3(h_2).
\]

This is stronger than changing one projected coefficient. It does not
compute the common, eta-independent part of full S6 or any order-eight
remainder, nor establish an actual rolling solution for a deformed parent.

There is an explicit non-boundary control. Take r,P,eta,MF2 constant and
an **off-shell physical** Einstein metric `Ric=Lambda g`, with nonzero
Lambda. Its Schouten-type mixed tensor is `p I`, `p=Lambda/3`, so

\[
 \Delta L_6=r\eta\alpha^3\Lambda^3/27.
\]

A compact conformal variation `delta g=2sigma g` gives
`delta Schouten_mixed=-2sigma p I-2 Hess(sigma)_mixed` and

\[
 \delta(\sqrt{|g|}e_3)=\sqrt{|g|}\{-2\sigma e_3-6p^2\Box\sigma\}.
\]

Because p is constant, the Box term integrates to zero for compact sigma.
The covariant-metric Euler trace is therefore `-Delta L6`, nonzero when
`eta Lambda!=0`. This is not an inference from a nonzero density alone.
In the B-sign flat-FLRW de Sitter chart `Lambda=-3H^2`, the density is
`-r eta alpha^3 H^6`. The fixture need not solve the parent field equations;
it proves non-boundary character of this deformation, not generic
all-order independence, health, or a physical singularity conclusion.

## 7. Source boundary and the next discriminating test

The literal constant-beta HR action and elementary polynomials are given
in [Gording--Schmidt-May, arXiv:1807.05011v4, equations2.2--2.4](https://arxiv.org/pdf/1807.05011v4).
Their equations3.1--3.10 implement own-f elimination; section2.3 explicitly
distinguishes that procedure from eliminating only the massive eigenmode.
The smooth-beta and variable-r statements here are derived directly above,
not attributed to their constant-beta vacuum calculation.

[Hassan--Schmidt-May--von Strauss, arXiv:1303.6940v2, section2.1 equations7--9](https://arxiv.org/html/1303.6940#S2.SS1)
distinguishes own-f differential elimination from solving the g equation,
including the chain-rule term and boundary dependence. Its appendixA also
retains source terms when comparing higher-derivative descriptions. No
massive-gravity-only theorem, vacuum pole estimate, or scalar stability
claim is imported here. Frozen S6.30 supplies the local physical-frame,
clock and necessary-coefficient dictionary; none of its files are modified.

Arbitrary beta2/beta3 do not repair the retained-order CD mismatch. The
new discriminating test is a **specified full higher-operator/source
comparison**: include the actual degree-six Einstein contribution and
this genuine cubic freedom, fix the normal-form complement, and determine
whether the required CD remainder can coexist with a quantitatively
smaller omitted effect on a stated physical interval and source class.
The algebraic Hessian is not that differential inverse. Derivative or
nonminimal parent couplings, X-dependent interactions, hidden-metric
matter, a nonsimple root, or a redefined physical metric are outside the
theorem; none is an automatically healthy or source-preserving repair.
