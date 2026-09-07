# Independent signature, coefficient and global lower-order audit

This component reads [Aoki--Shimada, arXiv:1806.02589v2](https://arxiv.org/pdf/1806.02589v2)
directly, especially Eqs. (2.3), (2.13), (4.1)--(4.12). It checks the
source's metric normal form against the actual CD functions independently
of the primary lift dictionary. The full connection variation is a
separately authored calculation, not duplicated here.

The paper-lookup route did not provide indexed markdown, so the original
versioned PDF and its layout-preserving text were inspected. No later
version or secondary summary is substituted. The printed Eq. (4.8) is
distinguished explicitly from the corrected coefficient needed below.

## 1. A signature dictionary, not a change of physical matter frame

Use the frozen [P8 conventions](../../../../FORMULATION.md): signature
+---, scalar kinetic variables with no factor one half, and Einstein
density \(-M^2R_B/2\), with \(R_B=-6(\dot H+2H^2)\). The paper uses the
opposite metric sign, with positive curvature coefficient and negative
minimal scalar kinetic density. Its curvature-index convention agrees
with the repository's after the metric sign is changed.

Write \(g_{\rm p}=-g_{\rm r}\) as a convention translation, not a disformal
redefinition of the physical geometry. The Levi-Civita connection and
lower Ricci tensor are unchanged, while inverse metric and curvature
scalar change sign. With the same scalar and covariant Hessian,

\[
x=X_{\rm p}=-X_{\rm r},\qquad R_{\rm p}=-R_B,\qquad
\Box_{\rm p}\phi=-\Box_{\rm r}\phi.
\]

Counting inverse metrics, or directly contracting a general symmetric
Hessian and gradient, gives

| scalar | paper value divided by repository value |
|---|---:|
| \(L_1,L_2\) | +1 |
| \(L_3,L_4\) | -1 |
| \(L_5\) | +1 |
| \(v^\mu\phi_{\mu\nu}v^\nu\) | +1 |

Thus the paper's metric normal-form coefficients map as

\[
F_{2,\rm r}(u,X)=-f(u,-X),\quad
(A_1,A_2,A_3,A_4,A_5)_{\rm r}
=(\alpha_1,\alpha_2,-\alpha_3,-\alpha_4,\alpha_5)(u,-X).
\tag{1}
\]

The paper's \(-(\partial\chi)^2_{\rm p}/2\) is exactly the repository's
\(+Y/2\). This free scalar uses partial derivatives and does not source
the independent connection. No fermion or other connection-dependent
matter coupling is added. The Fraction tests use dense, traceful Hessians,
not a time-only background which could conceal a sign error.

## 2. Independent substitution in the principal source formulas

Work in the frozen CD coefficient units \(M=\tau=1\), with clock label u.
These normalizations do not identify differently dimensioned coupling
functions outside those units. Let

\[
h=(1+u^2)^3,\quad w=h-1-x,\quad
f=\frac{w}{2h},\quad p=\frac12\sqrt{\frac wh},\quad
c=\frac{2(p-f)}x,\quad F_4=\frac{1/2-f}{x^2}.
\tag{2}
\]

Here the source couplings \(f_1,f_2\) are called p,c. Neither c nor
the paper's scalar parent \(F_2\) is a bimetric parameter or the CD
curvature coefficient. Differentiate in x at fixed u before substitution.
The useful identities are

\[
p^2=f/2,\qquad
2p-cx+2F_4x^2=1,\qquad \frac{p_x}{p}=-\frac1{2w}.
\tag{3}
\]

Direct substitution into the independently transcribed source
Eqs. (4.5), (4.9)--(4.12) gives the five target relations

\[
\begin{split}
f_{\rm out}&=\frac{w}{2h},\qquad \alpha_1=0,\\
\alpha_3&=\frac1{hx},\\
\alpha_4&=\frac{11x-4h+4}{4hxw},\\
\alpha_5&=-\frac1{hxw},
\end{split}
\qquad \alpha_2=-\alpha_1=0.
\tag{4}
\]

Independently insert
\(F_{2,\rm r}=-1/2+(1-X)/(2h)\), \(A_1=0\), \(A_3=1/(hX)\)
into the frozen Ia completion and then apply (1). It gives exactly (4).
The tests also compare this independent target directly with the original
CD specification and Ia API, not the primary affine dictionary.
No field equation, background value \(X=1\), or Hubble denominator is
used in these identities. This check of normal-form coefficients alone
does not establish the full distortion Hessian's invertibility.

## 3. The printed Q2 issue cannot be ignored

Let a subscript u denote the clock derivative in these units. With the
denominator in (3), the lower-order source expressions take the form

\[
P=F_{2,\rm parent}+3x(p_u-F_3x)^2,\qquad
Q_1=-2f_u+4p(p_u-F_3x).
\]

The coefficient required by the literal connection calculation is

\[
Q_{2,\rm corr}
=\frac{2f_u}{x}-\frac{4(p-3xp_x)(p_u-F_3x)}x.
\tag{5}
\]

The PDF instead prints \(p-3p_x\) in this bracket. Our lower-order theorem
uses (5), so it requires the separate literal action audit; it is not a
consequence of the printed formula without correction. The independent
[connection engine](../src/p8_affine/connection.py) supplies that separate
verification for the lift. The wrapper must check the actual coefficient
bridge, not only trust a source label.

A small independent control already distinguishes the two formulas.
For the pure \(p(\phi,x)R_\Gamma\) action, the conformal Palatini solution
has metric \(\bar g=p g\), up to the projective gauge. Its curvature
contraction gives

\[
p\,g^{\mu\nu}R_{\mu\nu}[\bar g]
=pR[g]-3\Box p+\frac{3}{2p}(\partial p)^2.
\]

Since
\(\partial_\mu p=p_u v_\mu+2p_x\phi_{\mu\nu}v^\nu\), the non-boundary
completion contains \(6p_up_x\,v^\mu\phi_{\mu\nu}v^\nu/p\).
The unspecialized corrected source formula, whose second term also has
the denominator \(2p-cx+2F_4x^2=2p\) in this pure-pR case, gives this Q2;
the printed bracket gives an extra factor \(1/x\).
The direct Fraction contraction test retains both terms in
\(\partial p\). This is a restricted calibration, not a substitute for
the full projective quotient calculation.

## 4. Exact lower-order match modulo a stated divergence

Set \(q=Q_1\), with

\[
q_x+\left(\frac1{2x}-\frac{3p_x}{2p}\right)q
=3\frac{p_x}{p}f_u,\qquad q(u,-1)=0.
\tag{6}
\]

Define the remaining parent functions by

\[
F_3=-\frac{q+f_u}{4px},\qquad
F_{2,\rm parent}
=F_{\rm repo}(u,-x)+xq_u
-\frac{3x(q+2f_u)^2}{16p^2}.
\tag{7}
\]

The relation \(f_u=4pp_u\) makes the first source expression exactly
\(Q_1=q\). Equation (5) gives
\(Q_{2,\rm corr}-2q_x=-2\mathcal E_q\), where \(\mathcal E_q\) is the
left side of (6) minus its right side. The other expression becomes
\(P=F_{\rm repo}(u,-x)+xq_u\).

For the Levi-Civita derivative in the paper convention,

\[
\nabla_\mu(q\,v^\mu)
=q\Box\phi+xq_u+2q_xv^\mu\phi_{\mu\nu}v^\nu.
\tag{8}
\]

Hence, on (6), the complete lower-order density
\(P+Q_1\Box\phi+Q_2v^\mu\phi_{\mu\nu}v^\nu\) equals
\(F_{\rm repo}(u,-x)\) plus (8). This is a compact-variation boundary,
not use of the clock or metric equations. After (1), it gives the actual
repository scalar F and K=0 on the open tube. The proof works for any
smooth F_repo, in particular the frozen rational CD function; it does
not need to expand that large polynomial again.

The incorrect printed bracket leaves the exact residual

\[
Q_{2,\rm print}-Q_{2,\rm corr}
=\frac{3(1-x)p_x}{px}(q+2f_u),
\]

generically nonzero off the trajectory. At x=-1, q and f_u both vanish,
so a background-only check would miss this defect. Also q_u(u,-1)=0
because the boundary condition holds identically in u. None of these
trajectory identities is substituted before the x differentiation.

## 5. Global tube and analytic ODE existence

The entire closed CD tube corresponds to

\[
u\in\mathbb R,\qquad -11/10\le x\le-9/10.
\tag{9}
\]

Since h>=1 and \(w/h=1-(1+x)/h\), continuously on (9),

\[
\frac9{10}\le\frac wh\le\frac{11}{10},\quad
\frac9{20}\le f\le\frac{11}{20},\quad
\frac9{40}\le p^2\le\frac{11}{40},\quad
\frac9{20}<p<\frac{11}{20}.
\tag{10}
\]

The first bounds are checked coefficientwise by writing their cleared
differences as sums of h-1 and distance to the x endpoints. All x, p and
w denominators are bounded away from zero; the denominator in (3) is one.
The additional factor seen in the separate connection quotient satisfies

\[
8p^2-1\ge4/5>0.
\tag{11}
\]

The weak inequality is necessary: equality occurs at u=0,x=-9/10.
It is not strictly greater than 4/5 on the closed tube. On the larger
open neighborhood \(-6/5<x<-4/5\), for every real u, the corresponding
factor is strictly greater than 3/5. This gives a genuine smooth open
coefficient domain, not merely an along-trajectory reconstruction.
Bounding one factor does not independently prove the connection rank;
that remains the separate literal quotient audit's obligation.

Both coefficient and forcing in the linear equation (6) are real analytic
on this open neighborhood. Its positive integrating factor is

\[
\mu(u,x)=\sqrt{-x}\,p(u,x)^{-3/2},\qquad
q(u,x)=\mu(u,x)^{-1}
\int_{-1}^{x}\mu(u,s)\,3\frac{p_s(u,s)}{p(u,s)}f_u(u,s)\,ds.
\tag{12}
\]

Differentiation proves (6) and its initial condition. The integral is over
a finite x segment wholly in the domain; its coefficients are analytic
jointly in u and s. Parameter-dependent integration therefore gives a
jointly real analytic q. The usual integrating-factor uniqueness proves
that these local analytic functions agree everywhere on the strip.
Every finite u is included. There is no time evolution through an H=0
division, no finite-u blowup of these coefficient functions, and no
additional integration constant to be selected at a crossing.

There is also a useful uniform bound. The exact polynomial identity

\[
9h^2-h_u^2=9(1+u^2)^4(u^2-1)^2\ge0
\]

gives \(|h_u|\le3h\), so on (9)

\[
|f_u|\le\frac3{20h},\quad
\left|\frac{p_x}{p}\right|\le\frac5{9h},\quad
|\mathcal A_{\rm ODE}|\le\frac{25}{18},\quad
|\mathcal B_{\rm ODE}|\le\frac1{4h^2}.
\]

Using \(|x+1|\le1/10\), the integral solution or Gronwall yields

\[
|q|\le\frac1{40h^2}e^{5/36}
<\frac1{40h^2}\frac1{1-5/36}
=\frac9{310h^2}<\frac1{32}.
\tag{13}
\]

The geometric-series exponential majorant is valid because \(5/36<1\).
Analyticity supplies the q_u in (7); (13) alone is not promoted to a
bound on every derivative or on quantum/strong-coupling effects.

## 6. Evidence and scope

The module audit.py imports neither the primary affine dictionary nor the
connection solver. Its source substitutions, original-CD completion,
dense Fraction signature/product-rule fixtures and continuous domain
identities are separate checks. The scientific tests explicitly reject
a wrong A3 signature, the uncorrected printed Q2, an unsupported strict
endpoint margin, floating-point/bool inputs and invalid exact domains.
Finite fixtures supplement, rather than replace, the continuous proof.

This audit establishes the coefficient/IBP/domain part of the proposed
lift, conditional on the separately verified literal connection
elimination and its regular quotient. It is not a new proof of scalar
kinetic positivity, a propagating healthy heavy sector, a loop-stable
completion, a physical cutoff, an X=0 vacuum or the adopted V/G tests.
Connection-dependent matter and extra curvature/derivative operators
would require their own source and rank checks. Original P8 is not
closed by this classical change of variables and auxiliary-field audit.
