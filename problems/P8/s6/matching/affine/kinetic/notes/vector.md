# Selective quotient vector: exact elimination and rolling linearization

## 1. Claim, normalization, and frozen input

This component considers one **new** action: add

\[
 -\frac{\lambda}{4}F_{\mu\nu}(V)F^{\mu\nu}(V),\qquad
 V_\mu=\kappa^a{}_{\mu a}-\frac14\kappa^a{}_{a\mu},\qquad
 F_{\mu\nu}(V)=\partial_\mu V_\nu-\partial_\nu V_\mu
 \tag{1}
\]

to the frozen S6.37 affine lift. Here initially `M=tau=1`, and the source
signature is `(-+++)`. The derivative index of the affine connection is
last: `kappa^a_bc=Gamma^a_bc-LC(g)^a_bc`. There is no torsion restriction.
The original physical metric and free canonical chi remain unchanged.
The distortion traces in (1) are tensors, and their projective changes
cancel: `kappa^a_bmu -> kappa^a_bmu+delta^a_b U_mu` adds `U_mu-U_mu=0`
to `V_mu`. Thus (1) depends on the 60-dimensional projective quotient,
not the projective Maxwell spectator studied separately.

The exact inputs are the frozen `p8_affine.connection`, `dictionary`,
and `lower` modules, especially the full unrestricted curvature/epsilon
action and integration-by-parts derivation in `../notes/connection.md`
(path relative to this child root). That source note audits the literal
Aoki--Shimada action and conventions. No published on-shell connection
formula or torsion-free ansatz is substituted here. The new vector curl
and its properties below are independently calculated additions, not a
kinetic theorem attributed to that paper.

The conclusions are exact elimination of the other 56 quotient
components, a continuous inverse/sign bound on the closed CD tube,
preservation of the original rolling background, and its complete
linearized vector source. A sign of an algebraic Hessian or an isolated
Proca mass is **not** a conclusion about the full coupled constraints,
kinetic health, heavy gap, cutoff, or UV completion.

## 2. Direct quotient Schur action and its 56-dimensional complement

In a Levi-Civita normal frame with scalar gradient
`v_a=(s,0,0,0)`, `s>0`, the frozen literal action is

\[
 L=L_0+J^T\kappa+\tfrac12\kappa^TM\kappa.
 \tag{2}
\]

All 64 components and all ten independent symmetric Levi-Civita
Hessian entries are retained. Let `E` be the frozen 64-by-60 trace-gauge
embedding, `M_q=E^T M E`, `J_q=E^T J`, and let `T` denote the explicit
4-by-64 contraction in (1). Set

\[
 N=TE,\qquad W=M_q^{-1}N^T,\qquad D=NW.
 \tag{3}
\]

`vector_map()` constructs `T` by its indices and checks its four
projective null directions. `schur()` solves every block of the actual
frozen `M_q` against `N^T`. It does not assume an invariant-sector ansatz
or a desired value of `D`. The resulting four-vector matrix is

\[
 \begin{split}
 D&=\operatorname{diag}(d_0,d_s,d_s,d_s),\\
 d_0&=\frac{3(p^3+2p^2+p-1)}{2p},\qquad
 d_s=\frac{(4p-5)(4p-1)}{32p^2}.
 \end{split}\tag{4}
\]

Write `L=W D^{-1}` and `Pi=I-LN`. The replay checks, entry by entry,

\[
 \begin{gathered}
 M_qW=N^T,\quad NL=I_4,\quad M_qL=N^TD^{-1},\\
 N\Pi=0,\quad\Pi^2=\Pi,\quad\operatorname{tr}\Pi=56,\\
 L^TM_qL=D^{-1},\quad\Pi^TM_qL=0.
 \end{gathered}\tag{5}
\]

Since an idempotent has rank equal to its trace, its image is exactly
the 56-dimensional `ker N`. The restriction of `M_q` to that kernel is
nondegenerate whenever `M_q` and `D` are invertible. Indeed, if
`z in ker N` is `M_q`-orthogonal to `ker N`, then `M_q z=N^T ell` for
some four-vector `ell`. Thus `z=W ell` and `0=Nz=D ell` imply `z=0`.
This proves invertibility of the fixed-vector complement without
silently selecting or dropping any of its equations.

Let `k_*=-M_q^{-1}J_q`, `V_*=Nk_*`. At fixed `V`, the unique
stationary representative is

\[
 k(V)=k_*+L(V-V_*).
 \tag{6}
\]

It has `Nk(V)=V`, and its residual `M_q k(V)+J_q=N^TD^{-1}(V-V_*)`
annihilates every fixed-`V` variation. Completing the square in (2)
gives the exact reduced action

\[
 L_{\mathrm{CD/M1}}+
 \tfrac12(V-V_*)^TD^{-1}(V-V_*)
 -\frac\lambda4F(V)^2.
 \tag{7}
\]

The last term contains no other quotient coordinate, so the same
56 algebraic equations remain exact after adding it. These are local
changes of connection coordinates at fixed metric and clock; they do
not replace the physical matter metric or alter chi's source coupling.

### Continuous closed-tube control

On the actual closed tube, `9/40 <= p^2 <= 11/40` and `p>0`. In
particular

\[
 \frac{47}{100}<p<\frac{11}{20},\qquad
 \frac9{40}-\left(\frac{47}{100}\right)^2=\frac{41}{10000},\quad
 \left(\frac{11}{20}\right)^2-\frac{11}{40}=\frac{11}{400}.
 \tag{8}
\]

The polynomial `f(p)=p^3+2p^2+p-1` has positive derivative
`3p^2+4p+1` for positive `p`, and `f(47/100)=15623/10^6`. Equations
(4) and (8) therefore prove the strict, continuous bounds

\[
 d_0>\frac{46869}{1100000}>0,\qquad
 -d_s>\frac{(5-4(11/20))(4(47/100)-1)}{32(11/20)^2}
 =\frac{14}{55}>0,
 \qquad \|D^{-1}\|_\infty<\frac{1100000}{46869}<24.
 \tag{9}
\]

The norm is the maximum absolute row-sum in this declared orthonormal
scalar-rest frame, not a Lorentz-frame-independent norm or a physical
gap. The parent already proves `M_q` invertible on this tube. Its
coarser `p>=9/20` bound must **not** be used to prove the sign of (4):
`f(9/20)<0`. This is an explicit omission/control test.

At `p=1/2`, which holds on the entire original clock trajectory,

\[
 D=\frac38\operatorname{diag}(1,-1,-1,-1),\qquad
 \tfrac12\Delta V^TD^{-1}\Delta V
 =\frac43(\Delta V_0^2-|\Delta\mathbf V|^2).
 \tag{10}
\]

The sign in (10) is useful for the subsequent constraint calculation;
it is not that calculation.

## 3. Covariant stationary vector from every connection source

Put `x=v^2<0`, `H_{mu nu}=nabla_mu nabla_nu phi`, and use lower vector
components. Covariance and the complete rest-frame calculation give

\[
 V_{*\mu}=a v_\mu+b H_{\mu\nu}v^\nu
       +c v_\mu\Box\phi+d v_\mu(v^\alpha H_{\alpha\beta}v^\beta).
 \tag{11}
\]

The letters `a,b,c,d` in (11) are **vector coefficients**, not the
scale factor or the parent's affine coupling named `c`. The module
constructs `V_*=T kappa_*` from the frozen full 64-vector solution.
It extracts (11) from the independent constant, `H11`, `H00`, and
`H01` coefficients, then reconstructs all four components with every
Hessian entry arbitrary. Before applying derivative relations, it finds

\[
 \begin{split}
 a={}&-\frac3{4p}
 (4F_3p^2x+4F_3px-4p^2p_\phi-4pp_\phi+p_\phi),\\
 b={}&-\frac{c_Xx^2-2p^2+2pp_Xx+p-4p_Xx}{px},\\
 c={}&\frac{(2p-1)(2p+3)}{2x},\\
 d={}&\frac{2c_Xx^2-4p^3+12p^2p_Xx-8p^2
                      +16pp_Xx+5p-11p_Xx}{2px^2}.
 \end{split}\tag{12}
\]

Here `c_X` on the right is the derivative of the **parent** affine
coupling. Its independent `c_phi` forcing was retained in the full
source; its cancellation in this projection is an output, not an input.
Rest-frame covariance is sufficient: any timelike scalar gradient can
be put in this frame, and the complete independent Hessian determines
the tensor identity. No homogeneity assumption was used in (11)--(12).

For the actual CD lift define `h=(1+u^2)^3`, and retain its smooth lower
coefficient `q` from the frozen integrating-factor ODE. The substitutions
are

\[
 \begin{gathered}
 p^2=\frac{h-1-x}{4h},\quad p_X=-\frac1{8hp},\quad
 f_\phi=\frac{h_\phi(1+x)}{2h^2},\quad p_\phi=\frac{f_\phi}{4p},\\
 c_X=\frac{2(1-4p)p_X}{x}-\frac{2(p-2p^2)}{x^2},\qquad
 F_3=-\frac{q+f_\phi}{4px}.
 \end{gathered}\tag{13}
\]

The denominator of `F3` is the product `4*p*x`, not `4*p_X`.
Substitution directly into (12) gives

\[
 \begin{split}
 a&=\frac{3(p+1)q}{4p}
       +\frac{3(8p^2+8p-1)f_\phi}{16p^2},\\
 b&=-\frac{2p-1}{x}-\frac{3p+1}{4hp^2},\qquad
 c=\frac{(2p-1)(2p+3)}{2x},\\
 d&=-\frac{4p^2-1}{2x^2}-\frac{12p^2-7}{16hp^2x}.
 \end{split}\tag{14}
\]

The imposed lower boundary function is `q(u,-1)=0`. Its ODE forcing
vanishes at `x=-1`, hence `q_X(u,-1)=0`; also `q_phi(u,-1)=0` by
differentiating that identically zero boundary function. Thus

\[
 (a,b,c,d)\big|_{x=-1}=(0,-5/(2h),0,-1/h),\qquad
 a_X\big|_{x=-1}=\frac{15h_\phi}{8h^2},\quad
 c_X^{\mathrm{vector}}\big|_{x=-1}=\frac1h.
 \tag{15}
\]

These are **total** coefficient derivatives including `p_X`. Merely
holding (15)'s point values fixed would lose the rolling lapse terms.

## 4. Exact original background and full ADM first variation

The original solution has `phi=u`, physical scale factor
`a_FLRW=(1+u^2)^2`, and `Hubble=4u/(1+u^2)`. On its entire trajectory
`x=-1`, `p=1/2`, and the parent's `c`, `F3`, `F4`, `p_phi`, `c_phi`
vanish. In a scalar-rest normal frame `H0mu=0`, while its spatial
Hessian need not vanish. Direct substitution into all 64 original
source entries gives `J=0` and `kappa_*=0`, not just `T kappa_*=0`.
Consequently `V_*=V=0` and `F(V)=0` preserve every background field
equation in the new action: both added terms and their first variations
vanish. This uses the actual rolling background, not an off-shell flat
metric with frozen coefficients.

For clarity now denote `Hubble` by `Hcal`. In unitary gauge take the
full first-order ADM metric

\[
 N=1+\epsilon n,\quad N_i=\epsilon B_i,\quad
 \gamma_{ij}=a_{\rm FLRW}^2(1+2\epsilon\zeta)\delta_{ij}.
 \tag{16}
\]

For the scalar shift, `B_i=partial_i beta`; keeping arbitrary `B_i`
also checks its cancellation. The inverse has
`g00=-1+2 eps n`, `g0i=eps B_i/a_FLRW^2`. Differentiating this metric
before inserting its background gives

\[
 \begin{split}
 H_{00}&=-\epsilon\dot n,\qquad
 H_{0i}=-\epsilon(\partial_i n+\mathcal H B_i),\\
 H_{ij}&=-a_{\rm FLRW}^2\mathcal H\delta_{ij}
  +\epsilon\{\partial_{(i}B_{j)}
  -a_{\rm FLRW}^2(2\mathcal H\zeta+\dot\zeta-2\mathcal H n)\delta_{ij}\}.
 \end{split}\tag{17}
\]

The raised scalar gradient also varies. Its `delta v^j=B_j/a_FLRW^2`
contribution cancels the shift in `H_{i nu}v^nu`, giving

\[
 \begin{gathered}
 H_{i\nu}v^\nu=\epsilon\partial_i n,\qquad
 vHv=-\epsilon\dot n,\qquad x=-1+2\epsilon n,\\
 \Box\phi=-3\mathcal H+\epsilon\{
 \dot n-3\dot\zeta+6\mathcal H n
            +a_{\rm FLRW}^{-2}\partial_iB_i\}.
 \end{gathered}\tag{18}
\]

In particular `delta c_vector=2n/h` multiplies the **nonzero background**
`Box(phi)=-3 Hcal`. Combining all terms in (14)--(18), before using the
clock relation, yields

\[
 \delta V_{*0}=-\frac3{2h}\dot n
       +\left(\frac{15h_\phi}{4h^2}-\frac{6\mathcal H}{h}\right)n,
 \qquad \delta V_{*i}=-\frac5{2h}\partial_i n.
 \tag{19}
\]

For the actual profile `h_phi/h=3 Hcal/2`, this simplifies to

\[
 \boxed{\delta V_{*0}=-\frac3{2h}\dot n
                        -\frac{3\mathcal H}{8h}n,\qquad
        \delta V_{*i}=-\frac5{2h}\partial_i n.}
 \tag{20}
\]

There is no `beta` or `zeta` dependence anywhere along the trajectory.
Their absence is the result of the full ADM and coefficient variations,
not a gauge truncation setting them to zero. Since the background vector
is zero, perturbing the frame does not add a hidden background-vector
term to (20). At the bounce (20) is `(-3 n_dot/2,-5 grad(n)/2)`.

## 5. Exact-gradient change and the constraint boundary

Taking the curl of (20) gives

\[
 \delta F(V_*)_{0i}=-h^{-1}\partial_i\dot n
                     +\frac{33\mathcal H}{8h}\partial_i n.
 \tag{21}
\]

Define the invertible linear, derivative-dependent change

\[
 W=V+d\left[\frac{3n}{2h}\right].
 \tag{22}
\]

It is a gradient of the **whole product**, not `V+(3/(2h))dn`.
Commuting partial derivatives proves `F(W)=F(V)` exactly at this
perturbative order, including time-dependent `h`. Its mass differences
are

\[
 (V-V_*)_0=W_0+\frac{21\mathcal H}{8h}n,\qquad
 (V-V_*)_i=W_i+h^{-1}\partial_i n.
 \tag{23}
\]

Using (10) in the background coordinate frame, the complete quadratic
added density, including its volume factor, is therefore

\[
 a_{\rm FLRW}^3\left[-\frac\lambda4F(W)^2
  +\frac43\left\{\left(W_0+\frac{21\mathcal H}{8h}n\right)^2
       -a_{\rm FLRW}^{-2}(W_i+h^{-1}\partial_i n)^2\right\}\right].
 \tag{24}
\]

The metric contractions in `F(W)^2` use the background metric in this
quadratic expression. No lapse velocity remains in (24). In particular,
a raw `n_dot^2` term before (22) does not alone prove that the DHOST
constraint is lost. Equation (22) is not a point transformation and is
not asserted to be a global nonlinear reduction to a spectator. The
full lapse/shift constraints, reduced kinetic matrix, and spatial
principal coefficients remain a separate coupled calculation. Equally,
(24) is not a proof of health merely because an isolated vector has the
usual center Proca sign. None of these statements changes the original
physical `g` or `chi`.

## 6. Physical units and independent controls

With physical `M^2>0`, `tau>0`, and dimensionless positive physical curl
coupling `zeta_phys`, use `u=t/tau` and normalized connection components
`V_norm=tau V_phys`. Factoring `M^2 tau^2` from the action gives

\[
 \lambda=\frac{\zeta_{\rm phys}}{M^2\tau^2},\qquad
 D_{\rm phys}=M^{-2}D_{\rm norm},\qquad
 m_{\rm isolated,phys}^2=\frac{8M^2}{3\zeta_{\rm phys}},\quad
 m_{\rm isolated,norm}^2=\tau^2m_{\rm isolated,phys}^2.
 \tag{25}
\]

The isolated mass is just the coefficient comparison of (10) with a
canonical Maxwell term after freezing the metric and other scalar
variables. It is not a rolling coupled gap. For example `M^2=3`,
`tau=2`, `zeta_phys=5` give `lambda=5/12`, physical mass squared `8/5`,
and normalized mass squared `32/5`. Thus restoring units is not the
replacement `lambda=zeta_phys` at arbitrary `M,tau`.

The ordinary tests separately invert the complete 60-by-60 quotient
at the center and both exact tube endpoints, test a dense generic
source and all 56 projected Euler directions, and derive the ADM
Hessian/shift cancellations before the vector simplification. Controls
omit coefficient variations, use the inadequate lower-`p` bracket,
or replace (22) by the wrong time-dependent partial-product shift.
Every such omission is distinguished by a nonzero exact expression.
Exact real finite input guards reject strings, booleans, binary floats,
nonpositive units, and unproved/outside tube values. The 30 named checks
also retain the complete matrix residuals, not just their determinants.
