# Literal-action domain proof

All source equation numbers below refer to the v1 HTML. The signature is
`(-,+,+,+)`. A subscript on a scalar is an ordinary derivative in flat
coordinates; its gradient is raised only where displayed. Work on an open
set with strictly timelike gradient when dividing by `X`.

## 1. Complete leading higher-derivative action

Freeze `a=a1(phi0)`, `g=g1(phi0)` and `f=1/2-g != 0`. Substituting source
(19) into (4) gives

\[
 A_3=\frac{r}{X}+O(1),\quad
 A_4=-\frac{r}{X}+O(1),\quad
 A_5=-\frac{r^2}{2fX}+O(1),\qquad r=2g-a.
\]

These limits are evaluated from the complete rational coefficients, not by
retaining a selected operator. Both printed and earlier-Ia versions of A4
have these residues; their exact difference is regular as specified in
the formulation. The denominator factors other than `X` have nonzero
limits at this point.

For `phi=phi0+epsilon psi`, the respective elementary invariants scale as
`L1,L2=O(epsilon^2)`, `L3,L4=O(epsilon^4)` and `L5=O(epsilon^6)`, while
`X=epsilon^2 Xpsi`. Hence the singular part has the directional limit

\[
 \epsilon^{-2}\sum_{i=3}^5 A_iL_i\longrightarrow
 r\frac{L_3[\psi]-L_4[\psi]}{X_\psi}.
\]

This holds uniformly with the necessary derivatives on compact sets of
smooth directions whose `Xpsi` is bounded strictly below zero, at regular
coefficient points. Coefficient variation in phi is higher order. A5
does not enter this limit. The constant-coefficient contribution
`A1(L1-L2)` is a flat divergence. `F2 R` vanishes on the fixed-flat-metric
restriction. At smooth coefficient points, the rest of `F` and a fixed
smooth chi restriction contribute only ordinary lower-derivative Taylor
terms, with a bilinear quadratic part.

If the full action had a smooth vacuum extension, its restriction to fixed
metric and scalar directions would have a bilinear second differential,
up to divergences. Subtracting these ordinary quadratic contributions
would make the displayed `r L` quadratic too. Sections 2 and 3 contradict
this for `r != 0`. Metric constraints or mixing cannot turn a nonbilinear
restriction of a purported smooth full action into a bilinear one. This
does not claim anything about a different reduced theory with a restricted
field domain excluding the constant-clock configurations.

At `phi0=0`, the source coefficients need not have an ordinary smooth
Taylor expansion; we do not use one. Section 5 gives an exact smooth
restriction `chi=B(phi)` on which all their Q-dependent reconstructed
terms cancel. It establishes the same full-action obstruction at
`chi0=B(0)`. At other chi values no existing second variation of the
literal nonsmooth potential is presumed. The no-open-neighborhood result
in section 4 needs only the generic smooth field points away from zero.

## 2. Full symmetric-Hessian Euler check

In a `2+1` slice embedded in four dimensions, let `v=(p,q,s)` denote the
**covariant** gradient, `eta=diag(-1,1,1)` and `H` the full symmetric
Hessian with six independent entries. Then

\[
 X=v^T\eta v,\quad
 \mathcal L(v,H)=
 \frac{[(\eta v)^TH(\eta v)]\operatorname{tr}(\eta H)
       -(\eta v)^TH\eta H(\eta v)}{X}.
\]

Differentiate with respect to every independent symmetric H entry before
setting `H=diag(d0,d1,d2)`. For a locally quadratic profile the third and
fourth jets vanish. The exact Euler expression reduces to

\[
 E(\mathcal L)=
 -\sum_i d_i\mathcal L_{v_iv_i}
 +\sum_i d_i^2\mathcal L_{H_{ii}v_iv_i}
 +\sum_{i<j}d_id_j\mathcal L_{H_{ij}v_iv_j}
 =-\frac{2d_0d_1d_2}{X}.
\]

There is no extra factor of two in the last sum: varying one independent
off-diagonal H entry already varies both matrix slots. Prematurely setting
off-diagonal entries to zero would lose actual variational terms.

At the origin choose `psi1=t+(t^2+x^2+y^2)/2` and
`psi2=2t+(t^2+x^2+y^2)/2`. The two outputs and their sum-direction output
are respectively `2,1/2,16/9`; the additivity defect is `-13/18`.
All directions, and their sum, have timelike gradient on a common small
neighborhood. They can be embedded in smooth local variations. An Euler
operator obtained from a quadratic form would be linear and cannot have
this defect. A boundary term has identically zero Euler expression, so
integration by parts cannot remove it. By contrast the rank-one homogeneous
profile `H=diag(1,0,0)` gives zero: a coefficient pole or a homogeneous-only
test would not have proved this result.

## 3. Independent action-level bulk/IBP control

Put

\[
 J^\mu=\psi^\mu\Box\psi-\psi^{\mu\nu}\psi_\nu,
 \qquad D=(\Box\psi)^2-\psi_{\mu\nu}\psi^{\mu\nu}.
\]

Full third-jet differentiation verifies
`div J=D` and `(partial_mu X)J^mu=2(L3-L4)`. Consequently

\[
 \mathcal L=\tfrac12\operatorname{div}(\log|X|J)
             -\tfrac12\log|X|D.
\]

Take `psi=t+l T(t)(cos x+cos y)`, with real nonzero smooth compactly
supported T, periodic x,y and an optional unit z period. Small enough l
keeps the gradient uniformly timelike. The direct rational-action expansion,
with normalized spatial average, gives

\[
 \langle\mathcal L\rangle=
 l^2(TT''+T'^2)
 +l^4(\tfrac12T^3T''+T^2T'^2)+O(l^5).
\]

The cubic average is zero. Compact time support makes the quadratic
coefficient a boundary term and the quartic integrated coefficient

\[
 -\tfrac12\int T^2T'^2\,dt<0.
\]

Nonzero compactly supported smooth T cannot have `T T'=0` everywhere,
which proves strictness. A separate Fraction Laurent-polynomial constant
term calculation uses `g=cos x+cos y` and gives

\[
 \langle g^4\rangle=9/4,\quad
 \langle g^2|\nabla g|^2\rangle=3/4,\quad
 \langle|\nabla g|^4\rangle=5/4,\quad
 \langle g^2[(\Delta g)^2-|D^2g|^2]\rangle=1,\quad
 \langle|\nabla g|^2[(\Delta g)^2-|D^2g|^2]\rangle=0.
\]

Applying these to the log/IBP representation gives the independent raw
quartic average
`3T^3T''/4+7T^2T'^2/4+9TT'^2T''/4+3T'^4/4`.
The identities `int T^3T''=-3 int T^2T'^2` and
`int TT'^2T''=-int T'^4/3` reproduce exactly the same nonzero bulk term.

For the five profiles `t+kl T g`, with `k=2,-2,1,-1,0`, take weights
`1,1,-4,-4,6`. Every quadratic functional has zero such finite difference;
here its leading value is `24` times the nonzero quartic coefficient.
All profiles have the same boundary jets, so boundary-function contributions
also cancel. This is a second, action-level obstruction, not just the
pointwise density failing a parallelogram law.

The torus is an exact Fourier calculation device. Section 2 already proves
the local Minkowski result. If a compact-spatial Minkowski bulk control is
desired, multiply the periodic perturbation by a smooth cutoff on a large
integer-period box with fixed-width transition. All derivatives are bounded
independently of box size, the interior bulk is of volume order, and the
transition contribution is of surface order. For sufficiently large boxes
the bulk remains nonzero. No global torus vacuum or S-matrix is assumed.

## 4. Why isolated zeros cannot furnish an open vacuum chart

For the fixed benchmark write `z=phi/tau`,

\[
 a(z)=\sigma(z)(1+z^2)^{1/6}+\sigma(-z)(1+z^2)^{1/20},\quad h=a'/a,
\]
\[
 N=z^2\tanh(z+1/10)+\tanh z,\quad D=(1+z^2)h,\quad Z=9D-2N.
\]

The positive scale factor, N,D,Z and
`g1=(2/3)sech^2(z+1/10)` are real analytic at every finite real z.
On the regular coefficient chart,
`a1/g1=(2N/D-1)/4` and `r/g1=Z/(4D)`.

Exact local differentiation gives `D'(0)=13/60`, `N'(0)=1` and
`Z'(0)=-1/20`. Thus D and Z are not identically zero. Their real zero
sets have no finite accumulation point unless the corresponding analytic
function vanishes identically. At z=0 the a1 quotient is removable, with
`a1/g1=107/52` and `r/g1=-3/52`. Other D zeros need not be declared
removable for the argument: deleting their isolated set still leaves
generic regular points in every open interval. Similarly f is analytic,
nonconstant, and has no open zero set.

It follows that every open finite-phi interval contains points with all
displayed coefficients regular, `phi != 0`, `f != 0` and `r != 0`.
Sections 1–3 forbid a smooth local full-action extension at those points
if it agrees with the literal action for all sufficiently small nonzero
X in the same neighborhood (agreement up to harmless boundary terms does
not help). Hence no open smooth vacuum neighborhood of that unchanged
kind exists. At `r=0` the checked leading obstruction vanishes; this gate
does not classify the next order there. At `f=0` the proof's rational
nondegenerate expansion is unavailable and the usual positive Einstein
vacuum premise already fails. No existence or absence of isolated formal
stationary points is asserted.

For clarity the often-written equations
`f0-W2 chi0^2=0`, `f0'-W2' chi0^2=0`, `W2 chi0=0` are only formal
zeroth-derivative necessities for a hypothetical regular extension retaining
those values. They are not the actual Euler equations of the undefined
literal X=0 action. Solving them alone would not supply a vacuum chart.

## 5. Fractional tilt and full source cancellation

Choose `b=9/500` for this audit of the printed Q formula, giving
`alpha=2(1+b)=509/250`. The source leaves b adjustable; this is an explicit
audit choice, not a claimed printed numerical benchmark. Source (21)
gives `n_s=1-2b*10/9=24/25` for this choice. This algebraic dictionary is
not an observational-fit certification.
At z=0,

\[
 Q=\tfrac12-\tfrac14|z|^\alpha+O(|z|^{\alpha+1}).
\]

Thus Q is C2, with second derivative locally Holder of exponent `9/250`,
but not C3. Its leading third derivative in phi is
`-alpha(alpha-1)(alpha-2) sgn(z)|z|^(alpha-3)/(4 tau^3)`.
There is no such cusp conclusion for the smooth `b=0` control.

Let `B(phi)=bar chi(phi)>0`. Source (23c) gives

\[
 W_2=-\frac{QB''+(3HQ+Q')B'}{B}.
\]

For the source example `lambda=7/50`, `chi_bias=29/1000`, `p=c=10`,
`B(0)=29/1000+(7/50)sech(10)>0` and
`tau B'(0)=(7/500)sech(10)tanh(10)>0`.
The leading nonsmooth W2 contribution is
`alpha B'(0) sgn(z)|z|^(alpha-1)/(4 tau B(0))`.
It is C1 but not C2. W2 itself does not diverge. Constant B is an explicit
exclusion control, giving W2=0.

It would be incorrect to infer a physical singularity from W2 or from the
raw reconstructed f0 alone. Source (24) contains
`f0=f0_regular+W2 B^2`, `f1=f1_regular-Q B'^2`. Combining them with **both**
chi source terms gives the exact Q-dependent part

\[
 L_Q=Q[X B'^2-(\partial\chi)^2]+W_2(B^2-\chi^2).
\]

This vanishes identically on `chi=B(phi)`, for arbitrary off-shell phi and
metric. More generally set `xi=chi-B(phi)`, a smooth invertible field
change. A single explicit IBP gives

\[
 L_Q\simeq-Q(\partial\xi)^2-W_2\xi^2
 +2\xi\{(Q'B'+QB'')(X+1)+QB'(\Box\phi+3H)\}.
\]

The linear term vanishes on the reconstructed clock background
`X=-1, Box phi=-3H`. The code checks the complete identity without
discarding source terms. The stated regularity concerns this literal
coefficient chart; integration by parts can change which coefficient
contains a derivative of Q. It is not a reduced kinetic, instability,
high-frequency, or field-redefinition-invariant health verdict. A C-infinity
coefficient repair would be a separately named modified model.

## 6. Scope of the matching conclusion

The theorem is invariant under adding local boundary terms and cannot be
evaded by inspecting homogeneous clock perturbations alone. It is not
invariant under changing the action off the rolling tube, which the frozen
S6.1 context explicitly permits as a new candidate. A smooth splice can
leave the source's `X=-1` tube unchanged and alter a neighborhood of X=0;
this gate neither constructs nor rules out a healthy matched parent of
that modified theory. It does not impose real analyticity on that parent.

The chi sector here is interacting and phi-dependent, unlike frozen M1.
The source's nonlinear cosmological checks and cubic/bispectrum scale
estimates are different questions from a common stationary-vacuum EFT,
quartic scattering, finite-gravity subtraction errors, and UV completion.
No result on those questions, or on all DHOST bounces, is inferred.
