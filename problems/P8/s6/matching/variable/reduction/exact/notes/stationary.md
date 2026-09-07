# Exact isotropic own-f branch on the fixed physical probe

This S6.33 calculation keeps the unchanged VARIABLE action, its positive
canonical clock composed with `theta=T`, the physical `g`, and the free
`chi` coupling. It solves only the **own-f** equations on the off-shell
probe `g=eta`. It does not impose the `g` or clock equations and is not
another full-parent/CD background. An exact homogeneous metric branch and
a functional inverse for its tensor perturbations are different objects.

The convention is P8(b)'s `(+---)`, `R_B=-6(Hdot+2H^2)` and
`EH=-M^2 R_B/2`. The physical-g proper clock is `T=tau*u`, with `M,tau>0`.
The actual uniform theorem uses

\[
 0<\delta\le1/100,\qquad c=2+\delta,\qquad |u|\le1/100.             \tag{1}
\]

The face `delta=0` in the desingularized proof is not an action at `u=0`.
The physical matter metric and its proper-time dictionary remain fixed.

## 1. Literal equations and their undivided identity

For `g=eta`, let `f=diag(N^2,-b^2,-b^2,-b^2)`, with `b,N>0`. Primes in
this section mean physical-g proper-time derivatives. The positive root
has `e1=N+3b`, `e4=Nb^3`. The exact own-f density is therefore

\[
 L_f=-3M^2bb'^2/N-2\beta_1N-6\beta_1b-2\beta_4Nb^3.                \tag{2}
\]

The covariant Einstein density differs from its first term by
`D_T(3M^2 b^2 b'/N)`. The API retains this boundary. The physical Einstein
term, `beta0`, and both matter terms have no own-f variation; the lapse
has not been gauge-fixed before variation. Literal lapse and common-scale
variations give

\[
 \begin{split}
 C&=3M^2bb'^2/N^2-2\beta_1-2\beta_4b^3,\\
 E&=6M^2bb''/N+3M^2b'^2/N-6M^2bb'N'/N^2-6\beta_1-6\beta_4Nb^2.
 \end{split}                                                       \tag{3}
\]

`E` sums the three equal spatial variations; the other own-f components
vanish by isotropy. Direct differentiation proves

\[
 C'-(b'/N)E=-2[\beta_1'+\beta_4'b^3-3\beta_1b'/N].                 \tag{4}
\]

The prescribed clock is held fixed in the variation, but its coefficient
derivatives are retained in (4). No Hubble rate or metric velocity has
been divided out.

## 2. Center data without assuming even metrics

The coefficient profiles are even. Their positive algebraic root is
`r^3=-beta1/beta4`, with `r(0)=2`. Assume only a regular solution with `b`
of class `C2`, `N` of class `C1`, `b(0)=r(0)` and `N>0`. Neither metric
function is initially assumed even. Evaluating `C=0` first gives
`3M^2 r0 b'(0)^2/N0^2=0`; hence the equation itself forces `b'(0)=0`.
Writing `Z=N0/r0`, the spatial equation yields

\[
 b''_0=\beta_{1,0}Z(1-Z)/M^2.                                     \tag{5}
\]

Differentiating `beta1+beta4*r^3=0`, or using its second-order Taylor
expansion, shows that the coefficient of `T^2` in the lapse constraint is

\[
 C_{[2]}=3M^2r_0b''_0{}^2/N_0^2
                   +3\beta_{1,0}(b''_0-r''_0)/r_0.
\]

This step needs no third derivative of `b` or second derivative/parity of
`N`. Substituting (5) gives the exact factor

\[
 M^2r_0 C_{[2]}/3
 =\beta_{1,0}[\beta_{1,0}(1-Z)-M^2r''_0].                          \tag{6}
\]

For the stipulated nonzero link,

\[
 Z=1-M^2r''_0/\beta_{1,0},\qquad b''_0=Zr''_0.                    \tag{7}
\]

There has been no division by `b''`. The actual profiles have
`beta1(0)=(M^2/tau^2)*32/[c(c-2)]` and
`r''(0)=-8(c+2)/(c*tau^2)`, giving

\[
 N_0=c^2/2,\qquad b''_0=-2c(c+2)/\tau^2.                           \tag{8}
\]

An odd lapse can satisfy these center equations alone. A test supplies
such a trial and detects its nonzero order-`u^2` source-balance residual.
Evenness below follows from the full equations, not the center data alone.

## 3. Fixed-c existence, uniqueness and automatic parity

Put `v=u^2` and `beta_i=(M^2/tau^2) B_i(v)`. The metric can still be an
arbitrary function `b(u)`. Define

\[
 F(v,b)=2(B_{1,v}+B_{4,v}b^3)/(3B_1).                              \tag{9}
\]

Equation (4) gives `b_u/N=u F(v,b)`, since `B1` is nonzero. The lapse
constraint then implies the algebraic equation

\[
 \mathcal A(v,b)=3bvF(v,b)^2-2(B_1+B_4b^3)=0.                      \tag{10}
\]

At `(v,b)=(0,2)`, for every fixed admissible `c`, its analytic coefficients
satisfy

\[
 \mathcal A_b=96/[c(c-2)]\ne0,\qquad F_0=-4(c+2)/c\ne0.           \tag{11}
\]

The analytic implicit-function theorem gives a unique analytic `b=b(v)`
near the center. **Every** regular positive solution through `b(0)=2`
must satisfy (10), and therefore coincides with this branch locally.
In particular `b(u)` is even as a consequence. On a punctured interval,
consistency fixes

\[
 N(v)=2b_v(v)/F(v,b(v)).                                            \tag{12}
\]

Its nonzero denominator extends (12) analytically through the center,
with the positive value (8). Continuity identifies the original lapse
there too. Thus `N` is also even, without a parity premise.

Conversely, (10)-(12) make `C=0` and the bracket in (4) zero. For `u!=0`,
`b_u/N=uF` is nonzero, so (4) gives `E=0`. Smoothness extends the spatial
equation to `u=0`. No division by `u` or `b_u` is used to assert an equation
at the center. This proves actual local own-f existence and uniqueness,
not only necessary center jets.

## 4. Joint desingularization and the uniform box

The naive rational equation contains `c-y=delta+8v+...`. Clearing it while
forgetting the branch can destroy the IFT Jacobian. Define instead

\[
 \begin{split}
 d&=1+v,&y&=2/d^4,&D&=c-y,\\
 Q&=32(1-v)/(cd^{14}),&J&=cd^4(1-7v)+12v,\\
 R&=d^8J/[8c(1-v)],&L&=R_v/R,&H&=D(J_v/J-6/d)-D_v.
 \end{split}                                                       \tag{13}
\]

`Q` is a rational numerator, not the tensor below; `R` is not curvature.
The unchanged profiles obey `B1=Q/D`, `B4=-QR/D`. Extract the factor

\[
 1-Rb^3=vD\zeta,\quad b=[(1-vD\zeta)/R]^{1/3},\quad
 F=2(v\zeta H-L)/3,\quad \mathcal T=3bF^2/(2Q).                    \tag{14}
\]

Use the positive cube root. Expressions (9) and (14) agree exactly where
the original coefficients are defined, and the lapse equation is
`C=2vQ(T-zeta)`. The regular branch satisfies `zeta=T`, including its
analytic center value, rather than discarding that condition at `v=0`.
At the joint center,

\[
 (b,Q,R,L,H,F,\zeta)=(2,16,1/8,12,-8,-8,12),\qquad
 \partial_\zeta(3bF^2-2Q\zeta)=-32.                               \tag{15}
\]

The independent exact-Fraction proof in [intervals.md](intervals.md)
establishes, on the entire box `0<=v<=1/10000`, `0<=delta<=1/100`,
`11<=zeta<=13`, that `11<T<13`, `|T_zeta|<1/100`, all denominator/root
conditions hold, `F<0`, and `19/10<b,N<21/10`. Banach contraction gives
a unique fixed point for each parameter pair; local analytic-IFT branches
patch by that uniqueness. The actual derivative and lapse are

\[
 \zeta_v=\mathcal T_v/(1-\mathcal T_\zeta),\qquad
 N=2(b_v+b_\zeta\zeta_v)/F,                                       \tag{16}
\]

where `b_v` on the right is a partial derivative. Dropping the implicit
term is an actual nonzero error, not a permissible clock choice.

For (1), `D>=delta>0`, so the original equations are defined throughout
the window. This construction is the fixed-c germ of section 3. A regular
solution through `b(0)=2` continued throughout the window cannot first
leave its `zeta` box: `zeta=(1-Rb^3)/(vD)` is continuous for `u!=0`, has
the interior center limit, and either first boundary would contradict
the strict self-map. No second regular continuation of the germ is added.

At `delta=v=0`, however, `D=0` and `Q=16`: the literal `B1` diverges.
Only the desingularized metric and its coefficient jets extend there.
Neither an action at `c=2` nor a bounded potential/Green inverse follows.

## 5. Exact center jets and the first retained metric

Direct implicit differentiation and a separate Fraction series solve give

\[
 \begin{split}
 b_{uu}(0)&=-2c(c+2),\\
 b_{uuuu}(0)&=-12(9c^5+14c^4-82c^3+44c^2+168c-256)/c,\\
 N_{uu}(0)&=(27c^4-12c^3-236c^2+592c-512)/4.
 \end{split}                                                       \tag{17}
\]

The joint limiting checks include `zeta_v=204`, `zeta_c=0`, `b_v=-8`,
`b_vv=-56`, `N_uu=16`. A physical time derivative of order `n` adds
`tau^-n`; `b,N` themselves are dimensionless. In particular

\[
 f_{00}(0)=N_0^2=c^4/4,                                           \tag{18}
\]

not the original rolling-g solution's `c^2`. On this unit-volume probe
the frozen first own-f correction instead gives

\[
 (f_0+f_2)_{00}(0)=2c^2-4,\qquad
 f_{00}(0)-(f_0+f_2)_{00}(0)=(c^2-4)^2/4.                           \tag{19}
\]

Indeed its central Schouten component is `P00=-2r''/r`, hence
`f2_00=-2M^2 r^2 r''/beta1`. Equation (7) makes the exact defect
`r^2(Z-1)^2`. This compares one metric component, not all coefficient
derivatives, CD DHOST functions, or a convergent derivative expansion.

## 6. Exact tensor action and two distinct inner operators

For homogeneous norm-two coordinates take

\[
 g=\operatorname{diag}(1,-e^q,-e^{-q},-1),\qquad
 f=\operatorname{diag}(N^2,-b^2e^Q,-b^2e^{-Q},-b^2).
\]

Both volumes are independent of their tensor amplitudes. The exact root
trace is `N+b[1+2cosh((Q-q)/2)]`. Including the physical-g Einstein term,
the literal quadratic action is

\[
 L_T=M^2q'^2/4+K_fQ'^2/4-\beta_1b(Q-q)^2/2,
 \qquad K_f=M^2b^3/N.                                             \tag{20}
\]

The `beta4` determinant and homogeneous canonical clock/free-chi
densities have no tensor term in this parametrization. The Einstein
boundary depends only on isotropic volume and has no missing tensor
boundary. Own-f isotropic corrections beginning at tensor order two do
not alter (20), since their linear action coefficient vanishes on the
exact own-f background. The physical-g/clock background is still off
shell; varying (20) does not solve its remaining isotropic equations.

Set `nu=2beta1*b`. Prescribing the physical perturbation `q` gives

\[
 (K_fQ')'+\nu(Q-q)=0.                                              \tag{21}
\]

Varying `q` too gives `M^2 q''+nu(q-Q)=0`. At the center,

\[
 K_f=16M^2/c^2,\quad \tau^2\nu/K_f=8c/(c-2),\quad
 \tau^2\nu(1/M^2+1/K_f)=8(c^2+16)/[c(c-2)].                        \tag{22}
\]

For fixed compact `x`, put `u=sqrt(delta)*x`, taking delta small enough
that `delta*x^2<=1/10000`. Analyticity on the closed joint box, locally
through its boundary, supplies bounded `v,c` derivatives and

\[
 b=2+O(\delta),\quad N=2+O(\delta),\quad K_f=4M^2+O(\delta),\quad
 D/\delta=1+8x^2+O(\delta).                                       \tag{23}
\]

These constants may depend on the fixed compact interval. Since
`K_f,u=2u K_f,v`, the inner drift `sqrt(delta)*K_f,u/K_f=O(delta)`.
Equation (21) therefore tends to

\[
 Q_{xx}+16(Q-q)/(1+8x^2)=0.                                        \tag{24}
\]

Canonical normalization `Y=sqrt(K_f/2)*Q` adds the u-time pump
`-(sqrt(K_f))_uu/sqrt(K_f)`. Its inner contribution is multiplied by
delta and vanishes as `O(delta)`, since both first and second v derivatives
are bounded. It cannot change 16 in (24).

Varying both tensor amplitudes instead gives

\[
 q_{xx}+64(q-Q)/(1+8x^2)=0,\qquad
 Q_{xx}+16(Q-q)/(1+8x^2)=0.                                        \tag{25}
\]

Thus `q+4Q` has zero second derivative, while `Q-q` has coefficient
`80/(1+8x^2)`. The prescribed-q own-f inverse and the coupled-relative
problem are distinct. Replacing 16 by 80 in (24) changes the question.
The shared limiting number also does not identify this finite-c off-shell
action with the original rolling-g tensor system.

Coefficient convergence gives ordinary ODE solution convergence on each
fixed compact inner interval when the inner-coordinate initial data and,
in the prescribed-q problem, the prescribed inputs converge there. Merely
finite but delta-dependent data need not converge. This is not a fixed
physical-window bound.
Homogeneous tensor data remain arbitrary until specified. An exact
isotropic branch selects no retarded, advanced, Feynman or other Green
function, nor a source-preserving inverse norm or heavy-data suppression.

## 7. Exceptional branches and scope controls

No `b''` division occurs. Dividing by nonzero `beta1` is explicit. If
`beta1(0)=0`, the positive turning-point lapse equation also requires
`beta4(0)=0`; the root and IFT premises then fail. This does not exclude
that different degenerate system. If `r''=0`, (7) gives `Z=1,b''=0`, but
(12) may need higher-order desingularization. A genuine own-f control is
constant positive `r`, arbitrary nonzero `beta1(T)`, `beta4=-beta1/r^3`,
and `b=N=r`: all own-f equations vanish. No full-parent vacuum is asserted.

The actual profile has nonzero `r''` and `F0`. A separate fixed-local
`c=1` algebra control yields `N0=1/2`, `b_uu=-6`; it is outside the
uniform positive-delta box and implies no tensor health. Negative lapses,
zero scales, other roots, off-domain continuations and added operators
are not silently included.

No full CD solution, scalar/vector stability, complete-parent matching,
controlled curvature expansion, Lorentz-invariant vacuum, loop bound,
physical cutoff or UV completion follows. Original P8 remains open.
This branch and its operator are exact calibration targets for future
approximate elimination and explicit state choices.
