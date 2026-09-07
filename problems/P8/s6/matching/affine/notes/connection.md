# Unrestricted connection: direct action, quotient, and Schur audit

## 1. Claim and source conventions

This component proves pointwise algebraic elimination of all 64 components of
the affine distortion for the stated CD lift, modulo exactly four projective
directions on its regular domain. It also reconstructs every quadratic and
linear scalar-Hessian coefficient from the stationary action. It does not
prove a propagating heavy-field interpretation, a kinetic sign, a UV
completion, or the separate vacuum and gravitational-error requirements of
the adopted S6 contract.

The primary action and index definitions are
[Aoki--Shimada, arXiv:1806.02589v2](https://arxiv.org/pdf/1806.02589v2),
equations (2.1)--(2.3), (2.16), (3.16)--(3.17), and (4.1)--(4.4).
Their convention is used literally here: signature `(-+++)`,
`X = g^{ab} phi_a phi_b` with no factor of one half, and the derivative
index of `Gamma^a_{bc}` is **last**. In particular,

\[
 \nabla^\Gamma_c\phi_b=\partial_c\phi_b-\Gamma^a{}_{bc}\phi_a,
 \qquad \kappa^a{}_{bc}=\Gamma^a{}_{bc}-\{^a{}_{bc}\}_g.
\]

No lower-index symmetry, torsion-free condition, or metric compatibility is
imposed. The double-dual definition of the affine Einstein tensor in (4.2)
must not be replaced by a symmetrized affine Ricci tensor.

Write the connection-dependent action as

\[
 pR_\Gamma+cG^{ab}_\Gamma\phi_a\phi_b
       +F_3L_3^\Gamma+F_4L_4^\Gamma.
 \tag{1}
\]

An independent scalar function `F2` and minimally coupled first-derivative
free-scalar matter do not depend on the connection and can be added back
unchanged. In particular, this audit introduces neither a matter metric
change nor a hidden connection coupling for the frozen free chi field.
Matter that does couple to the connection is outside this calculation.

At a chosen point take Levi-Civita normal coordinates and the orthonormal
timelike scalar rest frame

\[
 g_{ab}=\operatorname{diag}(-1,1,1,1),\quad
 v_a=\phi_a=(s,0,0,0),\quad s>0,\quad X=-s^2,
 \quad H_{ab}=\nabla_a\nabla_b\phi.
 \tag{2}
\]

All ten independent entries of the symmetric **Levi-Civita** Hessian `H`
remain arbitrary. Thus this is not a homogeneous-background truncation.
At fixed timelike gradient every such point can be put in this frame, so
equality of the resulting tensorial scalar polynomials for every `H`
establishes the local covariant identity on the timelike domain.

The code uses `p_phi`, `p_X`, `c_phi`, `c_X`, and `F3` as independent
forcing data even after substituting the values of `c` and `F4` for CD.
The actual lift's derivative relations are imposed separately in the
coefficient dictionary; none is silently used to make a source equation
vanish.

## 2. Literal quadratic polynomial and all boundary terms

With the conventions above,

\[
 R^a{}_{bcd}(\Gamma)=R^a{}_{bcd}(g)
 +\nabla_c\kappa^a{}_{bd}-\nabla_d\kappa^a{}_{bc}
 +\kappa^a{}_{ec}\kappa^e{}_{bd}
 -\kappa^a{}_{ed}\kappa^e{}_{bc}.
 \tag{3}
\]

The first pair of connection derivatives is integrated by parts before
specializing the scalar gradient. The coefficient derivative for the
`p R_Gamma` term is

\[
 p_{,a}=s(p_\phi\delta_{a0}-2p_XH_{a0}).
 \tag{4}
\]

Consequently its quadratic and linear terms are, with
`eta_a = (-1,1,1,1)_a`,

\[
 \begin{split}
 L_{p,2}&=p\sum_{abe}\eta_b
 (\kappa^a{}_{ea}\kappa^e{}_{bb}
             -\kappa^a{}_{eb}\kappa^e{}_{ba}),\\
 L_{p,1}&=\sum_{ab}\eta_b
 (-p_{,a}\kappa^a{}_{bb}+p_{,b}\kappa^a{}_{ba}).
 \end{split}\tag{5}
\]

Let `epsilon^{0123}=1`; the sign choice of this epsilon cancels between
the two factors. The double-dual curvature coefficient is

\[
 C^{abcd}=\frac14\sum_{\gamma\alpha\beta}
 \eta_a\eta_\gamma
 \epsilon^{\gamma\alpha ab}\epsilon^{\gamma\beta cd}
 c\,v_\alpha v_\beta.
 \tag{6}
\]

It gives the literal source
`-partial_c C^{abcd} kappa^a_{bd}
+partial_d C^{abcd} kappa^a_{bc}`. The derivative used is

\[
 \partial_e(cv_\alpha v_\beta)=
 s(c_\phi\delta_{e0}-2c_XH_{e0})s^2\delta_{\alpha0}\delta_{\beta0}
 +cs(H_{e\alpha}\delta_{\beta0}+\delta_{\alpha0}H_{e\beta}).
 \tag{7}
\]

Only after this operation can its quadratic piece be simplified to

\[
 L_{c,2}=\frac{cs^2}{2}\sum_{ije}
 (\kappa^i{}_{ei}\kappa^e{}_{jj}
               -\kappa^i{}_{ej}\kappa^e{}_{ji}),\qquad i,j=1,2,3.
 \tag{8}
\]

Putting the gradient in the rest frame before differentiating its
coefficient would lose the last two terms of (7).

For the two literal epsilon Galileons, define the generally nonsymmetric
spatial matrix

\[
 \mathsf H_{ij}=H_{ij}-s\kappa^0{}_{ji}.
\]

Direct epsilon contraction gives

\[
 L_3^\Gamma=2s^2\operatorname{tr}\mathsf H,
 \qquad
 L_4^\Gamma=s^2[(\operatorname{tr}\mathsf H)^2
                         -\operatorname{tr}(\mathsf H^2)].
 \tag{9}
\]

The second trace is not a Euclidean Frobenius square. An independent
nonsymmetric matrix test distinguishes these expressions. Raising a scalar
gradient before taking its affine derivative is a different operation in
this nonmetric theory and is not used.

There remains a Levi-Civita curvature term. The identity

\[
 R_{ab}v^av^b=(\Box\phi)^2-H_{ab}H^{ab}
 -\nabla_a(v^a\Box\phi-H^{ab}v_b)
 \tag{10}
\]

turns (1) into curvature coefficient `f=p-cX/2` plus the bare density

\[
 \begin{split}
 L_0={}&2F_3s^2\operatorname{tr}H_{\rm sp}
 +F_4s^2[(\operatorname{tr}H_{\rm sp})^2
                              -\operatorname{tr}(H_{\rm sp}^2)]\\
 &+c[(\Box\phi)^2-H_{ab}H^{ab}]
 +c_\phi[X\Box\phi-vHv]
 +2c_X[(vHv)\Box\phi-vH^2v].
 \end{split}\tag{11}
\]

Here `vHv = v^a H_ab v^b` and
`vH2v = v^a H_ac g^{cd} H_db v^b`. The derivatives of `c` in (11)
are required integration-by-parts terms, not optional scheme choices.

`literal()` constructs (5)--(9) directly, returns the 64-vector source
`J`, and uses the Hessian normalization

\[
 L_\kappa=L_0+J^T\kappa+\tfrac12\kappa^TM\kappa.
 \tag{12}
\]

The factor of one half is fixed by differentiating the literal quadratic
polynomial. The independent tests rebuild its curvature and epsilon
contractions on dense, unrestricted distortions, rather than importing a
published on-shell connection.

## 3. Four projective directions and the complete quotient

For each `c=0,...,3`, the projective column is

\[
 (G_c)^a{}_{bd}=\delta^a_b\delta_{dc}.
\]

Direct calculation gives `M G=0` and `G^T J=0` already for generic
`p,c,F4` and generic derivative forcing. Impose the trace gauge

\[
 \sum_a\kappa^a{}_{ac}=0,
 \quad\kappa^0{}_{0c}=-\sum_{i=1}^3\kappa^i{}_{ic}.
 \tag{13}
\]

Keep the other 60 components in increasing `(a,b,c)` index order
`16a+4b+c`. Let `E` be this 64-by-60 embedding. It has full rank,
`G^T E=0`, and `[E G]` has rank 64. Each projective orbit has exactly
one representative in (13), because its trace changes by `4 U_c`.

For the CD pointwise values use

\[
 c=\frac{4p^2-2p}{s^2},\qquad
 F_4=\frac{1/2-2p^2}{s^4},\qquad f=2p^2,
 \qquad\Delta:=2p-cX+2F_4X^2=1.
 \tag{14}
\]

The quotient Hessian `M_q=E^T M E` is independent of `s`. Its nonzero-entry
graph separates into nine blocks. In the component order generated by
`quotient()`, their sizes and determinants are:

| Size | Determinant |
| --- | --- |
| 9 | `8192 p^14` |
| 9 | `-32 p^6` |
| 6 | `p^4 (8p^2-1)` |
| 6 | `p^4 (8p^2-1)` |
| 9 | `8192 p^14` |
| 6 | `p^4 (8p^2-1)` |
| 9 | `8192 p^14` |
| 3 | `-16 p^6` |
| 3 | `-16 p^6` |

Off-block entries vanish exactly, and multiplying these determinants gives

\[
 \det M_q=-2^{52}p^{72}(8p^2-1)^3.
 \tag{15}
\]

Thus `p>0`, `s>0`, `p^2 != 1/8` imply a unique gauge-fixed solution
for every source in (12). No restriction to an ansatz of sourced tensor
structures was used to infer this rank.

At `p=sqrt(2)/4`, the three size-six blocks each acquire one additional
null vector; all other blocks remain invertible. The quotient rank is 57,
and the unrestricted rank is also 57 (seven null directions including
the four projective ones). The actual source is orthogonal to all three
additional null vectors, and the regular sourced particular solution
continues finitely there. This does **not** provide a unique elimination:
one may add any of the three extra null vectors. In particular, the
denominators in the particular connection formula of the source's
Appendix B are not by themselves a full 60-direction rank argument.
This is a checked exceptional-domain control, not a broader verdict
about that paper or a dynamical instability.

## 4. Uniform closed-tube inverse bound

The CD lift satisfies, including the closed tube endpoints,

\[
 \frac9{40}\le p^2\le\frac{11}{40},\qquad
 \frac9{20}<p<\frac{11}{20},\qquad
 8p^2-1\ge\frac45.
 \tag{16}
\]

The last inequality is not made strict on the closed tube. The bound
`inverse_bound()` is an independent rational-entry proof, not a sampled
matrix norm. Each block inverse is explicitly computed and its product
with the block checked to be the identity. After cancellation, the
denominators of its entries have only a numerical factor and factors
`p` and `8p^2-1`. For a polynomial numerator `sum_j n_j p^j`, use

\[
 \left|\sum_jn_jp^j\right|\le
       \sum_j|n_j|(11/20)^j.
 \tag{17}
\]

For a denominator `c p^a (8p^2-1)^b`, use
`|c| (9/20)^a (4/5)^b` as its strictly positive lower bound. Unknown
denominator factors raise an error instead of receiving an unproved
bound. Sum the resulting entry bounds in every row of every block.
The largest of the 60 exact row sums is

\[
 \boxed{\ \|M_q^{-1}\|_\infty
          \le\frac{880109}{36000}<25\ }.
 \tag{18}
\]

Here the norm is the induced component maximum norm (absolute row-sum
norm) on the displayed trace-gauge coordinates and their dual source
coordinates in the orthonormal scalar rest frame, with the action
normalization of (12). Consequently
`||kappa_q||_infinity <= (880109/36000) ||E^T J||_infinity`.
The tensor reconstructed by `E` obeys the additional elementary bound
`||E kappa_q||_infinity <= 3 ||kappa_q||_infinity`.
The norm is not uniform under arbitrarily large Lorentz boosts or
arbitrary changes of component coordinates. Neither it nor an
indefinite algebraic Hessian is a physical propagator, kinetic residue,
or mass/gap statement.

## 5. Literal Schur elimination and coefficient reconstruction

On the regular domain,

\[
 \kappa_*=-E M_q^{-1}E^T J,\qquad
 L_{\rm red}=L_0+\tfrac12J^T\kappa_*.
 \tag{19}
\]

The code checks all 64 equations `M kappa_*+J=0`, not just the
60 gauge-projected equations. This also follows from the complement
in section 3 and `G^T J=0`. Substitution into the original polynomial
agrees exactly with (19).

Use the five independent invariant basis elements

\[
 H_{ab}H^{ab},\quad(\Box\phi)^2,\quad
 (vHv)\Box\phi,\quad vH^2v,\quad(vHv)^2.
 \tag{20}
\]

`coefficients()` extracts their coefficients and those of `Box phi`,
`vHv`, and the constant term from the polynomial in all ten `H_ab`.
It then reconstructs the complete density, whose difference is zero.
This last check excludes an untested Hessian monomial or an accidental
homogeneous-only agreement.

Set `t=p_phi+F3 s^2` and

\[
 A=-c_Xs^4+4p^2+4pp_Xs^2-2p-2p_Xs^2.
\]

The result (with the connection-independent parent scalar term omitted)
is

\[
 \begin{split}
 P_{\rm corr}&=-3s^2t^2,\\
 Q_1&=4pt-2p_\phi-c_\phi s^2,\\
 Q_2&=\frac{4(p+3s^2p_X)t-2p_\phi-c_\phi s^2}{s^2},\\
 \alpha_1&=\alpha_2=0,\qquad
 \alpha_3=-\frac{2A}{s^4},\qquad
 \alpha_5=-\frac{4p_XA}{ps^4},\\
 \alpha_4&=\frac{2}{ps^4}
 \big[-c_Xps^4+2c_Xp_Xs^6+4p^3-4p^2p_Xs^2-2p^2\\
 &\hspace{32mm}-2pp_X^2s^4+2pp_Xs^2+4p_X^2s^4\big].
 \end{split}\tag{21}
\]

These are independently compared with all eight scalar coefficients
in the source's (4.6)--(4.12), under (14). All comparisons vanish after
the single `Q2` correction justified next. Neither these published
coefficients nor the source's on-shell connection formula is used to
build, invert, or reconstruct (12).

## 6. Direct lower-order control of the printed Q2 term

In general symbols, the directly confirmed expression is

\[
 Q_2=\frac{2f_\phi}{X}
 -\frac{4(p-3Xp_X)(p_\phi-F_3X)}{X\Delta}.
 \tag{22}
\]

Equation (4.8) of the specified v2 PDF instead prints `p-3 p_X` in
the second numerator. On the CD values the direct result minus that
printed expression is

\[
 \frac{12p_X(s^2+1)(F_3s^2+p_\phi)}{s^2}.
 \tag{23}
\]

It is generically nonzero; setting `p_X=0` would hide the discrepancy.
This correction is not selected by agreement with the desired CD
coefficient. There is a separate direct Palatini check: set
`c=F3=F4=c_phi=c_X=0`, leaving arbitrary `p(phi,X)>0`. Define
`p_a=partial_a p`. The trace-gauge connection

\[
 \kappa^a{}_{bc}
   =\frac{\delta^a_c p_b-g_{bc}p^a}{2p}
 \tag{24}
\]

is obtained from the conformal Levi-Civita connection by a projective
shift. Direct substitution satisfies all 64 Euler equations of the
literal Palatini action and yields

\[
 L_{\rm red}-pR(g)=\frac{3}{2p}(\partial p)^2,
 \quad Q_1=0,\quad Q_2=\frac{6p_\phi p_X}{p}.
 \tag{25}
\]

The printed expression instead gives `6 p_phi p_X/(p X)` in this
control. Equation (22) gives (25). This fixes precisely the missing
factor of `X`; it is not a claim that other source results fail.

## 7. Replay and boundary of the result

The public `checks()` returns 30 named exact residuals, including full
matrix residuals. The ordinary tests independently rebuild literal
curvature/epsilon contractions, check projective complementarity and
the entire block decomposition, test the exceptional kernel and
source compatibility, and verify the Palatini normalization. Closed
parameter endpoints and inexact, nonreal, zero, and singular inputs
have explicit controls. The validator permits regular points below
`p^2=1/8` without interpreting them as healthy; it rejects the
exceptional value itself when invertibility is requested.

For a smooth CD lift in its timelike tube, the unique gauge-fixed
solution and reduced action depend smoothly on the local metric and
scalar jets because (15) has no zero there. This is exact local
auxiliary elimination, not a derivative expansion or an assertion
that a massive connection has been integrated out below a cutoff.
The separate principal/signature/lower-order dictionary still must
be applied to match the frozen physical-g/free-chi CD action. The
algebraic inverse alone settles none of the additional physical
V/G requirements, nor a stationary-vacuum extension outside this
timelike tube.
