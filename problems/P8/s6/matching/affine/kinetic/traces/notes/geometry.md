# Two projectively invariant traces and the rolling quadratic action

## 1. Scope and conventions

The frozen input is the unrestricted 64-component S6.37 affine lift of
the original physical CD/M1 action, not a torsion-free specialization.
The derivative index of `kappa^a_bc=Gamma^a_bc-LC(g)^a_bc` is last. Use
source signature `(-+++)`, `x=g^{mu nu}phi_mu phi_nu`, and the timelike
scalar rest frame. At first set `M=tau=1`.

This component treats the constant real combination

\[
 \begin{split}
 V_\mu&=\kappa^a{}_{\mu a}-\tfrac14\kappa^a{}_{a\mu},\\
 U_\mu&=g_{\mu a}g^{bc}\kappa^a{}_{bc}-\tfrac14\kappa^a{}_{a\mu},\\
 T_\mu&=A V_\mu+B U_\mu,\qquad (A,B)\ne(0,0),
 \end{split}\tag{1}
\]

and the new term `S_K[T]=-(lambda/4) integral sqrt(-g) F(T)^2`, where
`lambda` is any nonzero real constant. It proves the exact trace maps,
their on-clock Schur action, the complete rolling source, and the
zero-Schur **quadratic** multiplier identity. Coupled dynamical verdicts
are separate: neither an algebraic Hessian sign nor a raw lapse-velocity
square is used here as an unproved kinetic or UV conclusion.

The literal action, all boundary terms, and exact projective quotient
are sourced by the frozen `p8_affine.connection` module and its proof
`../../notes/connection.md`. The actual lower-coefficient ODE and the
full first-order ADM Hessian are the frozen S6.37/S6.38 inputs, with
`../notes/vector.md` proving the latter. These paths are relative to
this child root. In particular the new two-trace result is a direct
extension of a replayed action calculation, not a published propagator
imported from a restricted metric-affine model.

Under the unrestricted projective change
`kappa^a_bmu -> kappa^a_bmu+delta^a_b Z_mu`, each first trace in (1)
changes by `Z_mu` and its subtracted trace changes by the same amount.
Both are covectors. `traces()` constructs their 8-by-64 contraction
matrix directly and checks all four projective columns. An independent
test also transforms every distortion index by an exact rational
Lorentz boost; this checks the metric contractions and lower-index
placement of `U`, rather than fitting its temporal component alone.

## 2. Literal quotient inversion and the eight-dimensional Schur matrix

Let `E` be the frozen trace-gauge embedding into the 64 components. It
retains all 60 quotient coordinates, with no lower-index symmetry. The
connection action, after its documented boundaries, is

\[
 L=L_0+J_q^T k+\tfrac12 k^T M_q k,\qquad
 M_q=E^TME,\quad J_q=E^TJ.
 \tag{2}
\]

On the entire closed CD tube `9/40 <= p^2 <= 11/40`, `p>0`, the frozen
action has exactly four projective null directions and invertible
`M_q`. Its literal block matrices, not a selected connection ansatz,
are used below. Stack the two full maps as `T_8`, and write

\[
 N_8=T_8E,\quad W_8=M_q^{-1}N_8^T,\quad D_8=N_8W_8.
 \tag{3}
\]

`schur()` separately solves every block against **all eight** columns
and checks `M_q W_8=N_8^T`. Before specializing the scalar coefficients,
the nonzero diagonal spacetime blocks are

\[
 \begin{array}{ll}
 D_{VV,00}=\dfrac{3(p^3+2p^2+p-1)}{2p},&
 D_{UU,00}=\dfrac{3(p-1)^2}{2},\\[4pt]
 D_{VU,00}=\dfrac{3(2p^3-2p-1)}{4p},&
 D_{VV,ii}=D_{UU,ii}=\dfrac{(4p-5)(4p-1)}{32p^2},\\[4pt]
 D_{VU,ii}=-\dfrac{16p^2-40p-5}{32p^2}.&
 \end{array}\tag{4}
\]

The original rolling trajectory has `p=1/2` at every time, not merely
at the bounce. With trace order `(V,U)` the result there is

\[
 D_8=\frac38
 \begin{pmatrix}1&-7\\-7&1\end{pmatrix}
 \mathbin{\otimes}\eta_+,\qquad
 \eta_+=\operatorname{diag}(1,-1,-1,-1)=-g_{\rm rest}.
 \tag{5}
\]

The explicit minus sign relative to the source metric is essential.
This matrix has rank eight and determinant `531441/256`. It therefore
cannot be read as two independent equal-sign Proca forms: the cross
entry is `-7` times the individual temporal entry.

For `C=(A I_4, B I_4)`, define `N=CN_8`, `W=W_8 C^T`, on clock. Then

\[
 M_qW=N^T,\qquad NW=W^TM_qW=\gamma\eta_+,
 \qquad \gamma=\frac38(A^2-14AB+B^2).
 \tag{6}
\]

Independently of `gamma`, `N` has rank four whenever `(A,B)!=0`.
For `A!=0`, use the four retained columns
`kappa^{a(mu)}_{mu a(mu)}`, `a(mu)=(mu+1) mod 4`. Their minor is `A I`.
If `A=0`, use `kappa^mu_{b(mu)b(mu)}` with the same cyclic choice; its
minor is `B diag(g_mu_mu g_b_b)`, invertible for `B!=0`.
These columns are not among the four trace-gauge omitted entries.
The code checks both minors exactly, so no Schur invertibility is
silently used to prove rank at the exceptional combinations.

## 3. A multiplier identity retaining every connection equation

Put `k_*=-M_q^{-1}J_q` and `y=k-k_*`. The stationary traces are
`T_*=N k_*`. Introduce an independent `T` and a four-component
multiplier `ell` before performing any division by `gamma`:

\[
 L=L_{\rm CD/M1}+\tfrac12y^TM_qy+L_K[T]
                 +\ell^T(T-T_*-Ny).
 \tag{7}
\]

The connection equation is `M_q y-N^T ell=0`, hence uniquely
`y=W ell`. Substitution gives

\[
 L=L_{\rm CD/M1}+L_K[T]+\ell^T(T-T_*)
                      -\tfrac12\ell^T D\ell,
 \qquad D=\gamma\eta_+.
 \tag{8}
\]

`action_identity()` proves both the action substitution and the
component identities

\[
 M_qW=N^T,\qquad M E W=T_{\rm full}^T.
 \tag{9}
\]

The second is the full 64-component equation. Its four projective
equations have zero source; none of the other 60 is dropped. Writing
`E_K[T]` for the variational derivative of the actual kinetic functional,
normalized by the same volume element as (7), the other equations of
(8) are

\[
 T-T_*-D\ell=0,\qquad E_K[T]+\ell=0.
 \tag{10}
\]

The notation `E_K` retains its full differential operator, with all
time-dependent metric coefficients and boundaries of the action. No
fixed-frequency or slowly-varying approximation is made: varying
`T=N(u)y+T_*` simply gives `N(u)^T E_K`, after taking the functional
derivative with respect to `T`. Thus (9)--(10) also hold in a rolling
frame, rather than only for a matrix frozen in time.

For `gamma!=0`, eliminating `ell=D^{-1}(T-T_*)` in (8) gives

\[
 L=L_{\rm CD/M1}+L_K[T]
                +\frac1{2\gamma}(T-T_*)^T\eta_+(T-T_*).
 \tag{11}
\]

The sign is plus. Equation (11) is the constrained Schur action;
its vector kinetic sign and the lapse/shift constraints must still be
analyzed together to decide a physical spectrum.

### Null Schur combinations

If `gamma=0` and the combination is nonzero, necessarily `B!=0` and
`A/B=7+4sqrt(3)` or `7-4sqrt(3)`. In particular `A-B!=0`.
There is no fixed-`T` 56-dimensional complement inverse in this case.
Indeed `ker N` still has dimension 56, but its `M_q`-radical is exactly
the four-dimensional image of `W`:

\[
 NW=0,\qquad W^TM_qW=0.
 \tag{12}
\]

`W` has rank four by (9) and the rank of `N`. Conversely if a vector
`z in ker N` is `M_q`-orthogonal to `ker N`, then `M_q z=N^T ell`,
so `z=W ell`. This establishes the radical and the nondegenerate
52-dimensional quotient, rather than mistaking the four null vectors
for omitted healthy fields.

Equations (8)--(10) remain valid and now imply exactly

\[
 T=T_*,\quad \ell=-E_K[T_*],\quad y=-W E_K[T_*].
 \tag{13}
\]

All connection equations are satisfied, and `y^T M_q y=0`. Hence the
stationary quadratic action is precisely

\[
 L_{\rm CD/M1}^{(2)}+L_K[T_*^{(1)}].
 \tag{14}
\]

The implementation checks (12)--(14) as polynomial identities modulo
`A^2-14AB+B^2`, and tests both exact real roots separately. It never
substitutes a divergent `1/gamma` expression.

This is explicitly a **quadratic rolling-background** assertion.
The original trajectory has `kappa_*=0`; the quadratic action contains
`M_q` and the trace maps at that background, and its first-order source
`J_q^(1)`. Their field-dependent coefficient changes in `M_q k^2`
are cubic, so no such term is lost from (7) at this order. However the
general expressions (4) do not obey `D=0` on an open field tube when
the constants select a null on-clock combination. An exact test at
`p=51/100` distinguishes this case. Neither nonlinear constant-rank
elimination nor a smooth inverse of the exceptional complement is
claimed.

## 4. Full stationary source and rolling lapse variation

The complete `U_*` is obtained by contracting the frozen 64-component
solution, keeping all ten entries of the scalar Hessian arbitrary.
It is reconstructed in the covariant basis

\[
 U_{*\mu}=a_U\phi_\mu+b_U H_{\mu\nu}\phi^\nu
       +c_U\phi_\mu\Box\phi+d_U\phi_\mu(\phi^\alpha H_{\alpha\beta}\phi^\beta).
 \tag{15}
\]

These coefficients are distinct from the scale factor and the parent
affine coefficient named `c`. Let `h=(1+u^2)^3`, and let `q(u,x)` be
the frozen smooth lower-order coefficient with `q(u,-1)=0`. With the
actual CD substitutions

\[
 p^2=\frac{h-1-x}{4h},\quad p_X=-\frac1{8hp},\quad
 f_\phi=\frac{h_\phi(1+x)}{2h^2},\quad
 F_3=-\frac{q+f_\phi}{4px},
 \tag{16}
\]

the all-source calculation gives

\[
 \begin{split}
 a_U&=\frac{3(p-1)q}{4p}
        +\frac{3(8p^2-8p-1)f_\phi}{16p^2},\\
 b_U&=\frac{2p-1}{x}+\frac{3p+1}{4hp^2},\qquad
 c_U=\frac{(2p-1)^2}{2x},\\
 d_U&=-\frac{4p^2-1}{2x^2}-\frac{12p^2+1}{16hp^2x}.
 \end{split}\tag{17}
\]

The code first retains generic derivative forcing in the literal
source and only then makes (16) and the parent `c_X` substitution.
The reconstruction residual has four components and all independent
Hessian coefficients, so (17) is not fitted only to the bounce lapse.

On clock, the lower ODE gives `q_X(u,-1)=0`. The values and total
coefficient derivatives needed for a lapse perturbation are

\[
 (a_U,b_U,c_U,d_U)=(0,5/(2h),0,1/h),\quad
 a_{U,X}=-\frac{9h_\phi}{8h^2},\quad c_{U,X}=0.
 \tag{18}
\]

Use the full S6.38 ADM inverse/Hessian with unitary `phi=u`,
`N=1+n`, `N_i=B_i`, and spatial metric
`a_FLRW^2(1+2zeta)delta_ij`, where `a_FLRW=(1+u^2)^2` and
`Hcal=4u/(1+u^2)`. In particular `delta x=2n`,
`H_{i nu}phi^nu=partial_i n`, and `phi.H.phi=-n_dot`. The cancellation
of the scalar shift in `H.phi` requires varying the raised gradient
as well as the coordinate Hessian. The frozen ADM calculation retains
arbitrary shift, its derivatives, and spatial curvature jets before
this simplification.

Equation (18) therefore gives, before using `h_phi/h=3Hcal/2`,
`delta Ustar_0=3n_dot/(2h)-9h_phi*n/(4h^2)`. The complete result is

\[
 \begin{split}
 \delta V_{*0}&=-\frac3{2h}\dot n-\frac{3\mathcal H}{8h}n,
 &\delta V_{*i}&=-\frac5{2h}\partial_i n,\\
 \delta U_{*0}&=+\frac3{2h}\dot n-\frac{27\mathcal H}{8h}n,
 &\delta U_{*i}&=+\frac5{2h}\partial_i n.
 \end{split}\tag{19}
\]

No shift or `zeta` term survives, but the background coefficient
variations do survive in the temporal lapse coefficient. In particular
`Ustar=-Vstar` is false off center: their sum is
`(-15Hcal*n/(4h),0,0,0)` at first order.

For the general combination define the true gradient shift

\[
 \mathcal W=T+d\left[\frac{3(A-B)}{2h}n\right].
 \tag{20}
\]

It is the derivative of the **whole product**, including `h(u)`;
therefore `F(mathcal W)=F(T)`. With lower coordinate spatial components,

\[
 (T-T_*)_0=\mathcal W_0+d_Tn,\quad
 (T-T_*)_i=\mathcal W_i+e_T\partial_i n,\qquad
 d_T=\frac{3(7A+3B)\mathcal H}{8h},\quad e_T=\frac{A-B}{h}.
 \tag{21}
\]

For `gamma!=0`, (11) thus adds the quadratic density divided by
`a_FLRW^3`

\[
 -\frac\lambda4F(\mathcal W)^2+
 \frac1{2\gamma}\left[(\mathcal W_0+d_Tn)^2
       -a_{\rm FLRW}^{-2}(\mathcal W_i+e_T\partial_i n)^2\right].
 \tag{22}
\]

For `gamma=0`, substitution (14) instead uses the actual source curl

\[
 F(T_*^{(1)})_{0i}=-e_T\partial_i\dot n
            +\frac{3(11A-B)\mathcal H}{8h}\partial_i n,
 \qquad F(T_*^{(1)})_{ij}=0.
 \tag{23}
\]

For a real normalized scalar Fourier amplitude and spatial momentum
squared `q_spatial=k_comoving^2/a_FLRW^2>0`, its density divided by
`a_FLRW^3` is

\[
 \frac{\lambda q_{\rm spatial}}2
 \left[-e_T\dot n+\frac{3(11A-B)\mathcal H}{8h}n\right]^2.
 \tag{24}
\]

This is the exact quadratic action to which the metric/matter scalar
constraints must be applied, not a conclusion based only on seeing
`n_dot^2`. The code names `q_spatial` separately from the coefficient
function `q(u,x)` in (16). It does not freeze the time dependence of
the physical momentum when deriving the subsequent Euler equations.

All 64 parent sources vanish on the original rolling solution and
`kappa_*=0` there. Thus both traces and their curls vanish; each new
constant kinetic addition has zero first variation on that same
background. The entire physical metric, clock and free chi are kept.

## 5. Units, domains, and evidence boundary

For positive physical `M^2,tau`, with constant dimensionless `A,B`,

\[
 T_{\rm norm}=\tau T_{\rm phys},\qquad
 D_{\rm phys}=M^{-2}D_{\rm norm},\qquad
 \lambda=\frac{\zeta_{\rm phys}}{M^2\tau^2}.
 \tag{25}
\]

Here `zeta_phys` is a nonzero real dimensionless four-dimensional
curl coefficient. `M^2 tau^2` is factored from the action. The geometry
APIs allow both signs of `zeta_phys`, both null real combinations, and
all other nonzero real `(A,B)`; they reject the trivial zero map, zero
kinetic coupling, inexact floats, booleans, strings, nonfinite/nonreal
values and nonpositive units. No division by `gamma` occurs in the
shared public identities or parameter guard.

The independent tests reconstruct dense unrestricted traces, perform
a nontrivial Lorentz transformation, directly invert the complete
center quotient, and verify both minors, both null roots, every
multiplier Euler equation and the wrong nonlinear continuation control.
The scientific scope is the original source-preserving rolling action
and its stated kinetic promotions. A full two-by-two kinetic-form
classification, a vacuum extension, frequency window, heavy gap or
UV completion does not follow from these geometry identities alone.
