# Proof: sourced variable coefficients, an actual local bounce, and its inner tensor operator

## 1. Conventions and source-aware undivided identities

We use the actual P8(b) matter frame, signature `+---`, curvature
`R_B=-6(DH+2H^2)`, and Einstein action `-G sqrt|g| R_B/2`. Consequently
`3G H^2=rho` and `-2G DH=rho+p`. These are not the opposite curvature/Einstein
tensor conventions used in some photon-track sources. All `G`, `M`, and
`tau` below are constant and positive. Interaction coefficients have mass
dimension four; the interaction normalization is `-2 beta`, as in S6.17.

First allow one canonical clock on either endpoint and let
`beta_n=beta_n(phi_g,phi_f)`. Write

\[
y=a_f/a_g>0,\qquad c=N_f/N_g>0,\qquad
q_g=D_g\phi_g,\quad q_f=D_f\phi_f,
\]
\[
U=\beta_0+3\beta_1y+3\beta_2y^2+\beta_3y^3,
\quad V=\beta_1+3\beta_2y+3\beta_3y^2+\beta_4y^3,
\quad P=2(\beta_1+2\beta_2y+\beta_3y^2).
\]

The literal homogeneous interaction density before fixing either lapse is

\[
-2\bigl[N_ga_g^3\beta_0+
\beta_1(N_fa_g^3+3N_ga_g^2a_f)
+3\beta_2(N_fa_g^2a_f+N_ga_ga_f^2)
+\beta_3(3N_fa_ga_f^2+N_ga_f^3)+\beta_4N_fa_f^3\bigr].
\]

Differentiating the four coframes, with the clocks held fixed, gives

\[
\rho_{Ig}=2U,\quad \rho_{If}=2V/y^3,\qquad
\nu_{Ig}=(y-c)P,\quad \nu_{If}=\frac{(c-y)P}{cy^3}.
\tag{1}
\]

Here `nu=rho+p`. With a separately written clock potential `W_i`, its equation
is `D_i q_i+3H_i q_i+W_{i,phi}=F_i`, where literal clock variation gives

\[
F_g=-2(U_{,\phi_g}+cV_{,\phi_g}),\qquad
F_f=-\frac{2(U_{,\phi_f}+cV_{,\phi_f})}{cy^3}.
\tag{2}
\]

`D_g y=y(cH_f-H_g)`. Add the clock's balance `q_i F_i` to the interaction
balance, rather than pretending that the clock is separately conserved:

\[
\mathcal C_g=D_g\rho_{Ig}+3H_g\nu_{Ig}+q_gF_g=c\mathcal B,
\qquad
\mathcal C_f=D_f\rho_{If}+3H_f\nu_{If}+q_fF_f=-\frac{\mathcal B}{cy^3},
\]
\[
\boxed{\mathcal B=3P(yH_f-H_g)
+2(U_{,\phi_f}q_f-V_{,\phi_g}q_g).}
\tag{3}
\]

Thus `N_g^2 a_g^3 C_g+N_f^2 a_f^3 C_f=0`. The bimetric on-shell constraint
is the undivided equation `B=0`, not `P(yH_f-H_g)=0`. On a finite tree whose
edge coefficients depend only on its endpoint clocks, the same incidence
argument can set the total edge flux to zero, but it still sets (3), not the
old constant-coefficient velocity factor, to zero. No division by `P`, `H`,
or a clock speed is used in this identity.

For context only, an Einstein-frame canonical clock with other matter on
`h_i=A_i(phi_i)^2 g_i` also has a trace exchange. In four dimensions, with
`alpha=d log A/d phi`, `rho_E=A^4 rho_h`, `p_E=A^4 p_h`, the ordinary-matter
balance is `alpha q (rho_E-3p_E)` and the clock force includes its negative.
They cancel locally when both are retained. The physical Hubble rate is then
`H_h=A^-1(H_g+alpha q)`, not `H_g`. Our specified fixture has `A=1` and keeps
the literal physical metric `g`; it does not silently transfer a conformal
frame theorem. Source audit and prior-art attribution are in `sources.md`.

## 2. The actual canonical-clock/free-matter solution

Take the action in the formulation with no separate `W_phi`, no `f` source,
and beta independent of the free `chi`. Set `u=T/tau`, `d=1+u^2` and

\[
a=d^2,\quad b=2d^{-2},\quad N_g=1,\quad N_f=c,\quad
y=2d^{-4},\quad h=4u/d,\quad h'=4(1-u^2)/d^2.
\]

The two physical endpoint Hubble rates and derivative clocks are

\[
H_g=h/\tau,\quad H_f=-h/(c\tau),\qquad
D_f H_f=-h'/(c^2\tau^2).
\tag{4}
\]

Define the profiles in the formulation. In `M=tau=1` notation put
`Ubar=b0+3b1 y`, `Vbar=b1+b4 y^3`, `Pbar=2b1`. Direct substitution gives
the four independent isotropic Einstein equations:

\[
3h^2=\bar n/2+2\bar U,\quad
3h^2/c^2=2\bar V/y^3,\quad
-2h'=\bar n+(y-c)\bar P,\quad
2h'/c^2=(c-y)\bar P/(cy^3).
\tag{5}
\]

The fourth equation determines `b1`; adding the weighted two null equations
determines `nbar=2(y^3/c-1)h'`. The two lapse equations determine `b0,b4`.
The scalar identities, checked independently rather than inferred only from
the metric equations, are

\[
(a^3\chi_{,T})_{,T}=0,
\]
\[
\frac12\bar k'+3h\bar k
+2\bigl[b_0'+(c+3y)b_1'+cy^3b_4'\bigr]=0.
\tag{6}
\]

Equation (6) is the clock equation multiplied by its strictly positive clock
speed, with the common physical units restored in every term. In the beta
derivatives in (6), the coefficient profile is differentiated through its
clock dependence; the explicit `y` in a partial derivative of `U` or `V` is
held fixed. The modified Bianchi identity is exactly

\[
3\bar P(-yh/c-h)-2(b_1'+y^3b_4')=0.
\tag{7}
\]

At generic nonzero `u`, either term in (7) is nonzero. Omitting the clock
exchange is a tested negative control. It is precisely why S6.17 cannot be
applied to this action. Homogeneity, isotropy, and the diagonal root make
the remaining vector/off-diagonal/shear background equations vanish. The
separately authored audit checks the lapse and spatial equations explicitly.
An optional constant free `f` scalar obeys its unsourced equation identically.

On `|u|<=1/10`, `d<=101/100` and the exact rational margins imply

\[
d^{12}<8/7,\quad y^3>7,\quad h'>3,\quad y>1,\quad
0<\chi_{,u}^2/M^2\le 1/100.
\]

For `c=1`, `nbar>36`, hence `kbar>3599/100`. For `2<c<=4`,
`nbar>2(7/4-1)3=9/2`, hence `kbar>449/100`. Every scale/lapse is positive;
for `c=1`, `c-y<0` stays separated from zero, and for fixed `c>2`,
`c-y>=c-2>0`. All rational profiles are analytic on a neighborhood of the
closed interval for each fixed admitted `c`. Their positive analytic square
root gives a strictly monotone analytic clock; the analytic inverse function
theorem supplies `u(phi/M)` on its image, including local extensions at the
endpoints. This constructs actual coefficient functions in the action.
No uniform boundedness of those functions as `c` tends to `2` is asserted.

At `u=0`, the exact `(b0,b1,b4,kbar)` values are

\[
c=1:\ (178,-32,4,5599/100),\qquad
c=4:\ (-26,4,-1/2,799/100).
\]

All physical equations restore with `beta=M^2 b/tau^2`,
`phi_{,T}^2=M^2 kbar/tau^2`, `chi_{,T}=M/(10 tau a^3)`.
Although the target scale factor and free rolling matter are realized locally,
the interacting canonical clock is not the old unitary DHOST clock and no
old C/D physical operator or field map has been constructed.

## 3. Literal transverse shifts and two-TT action

For one TT polarization normalized by `e_ij e_ij=1`, let `gamma_g,gamma_f`
denote the metric amplitudes. Literal diagonal exponential metric variables
use roots `c,y exp(Delta/2),y exp(-Delta/2),y`, with norm-two polarization
`Delta=gamma_f-gamma_g`. Expanding all five elementary polynomials gives
the norm-two coefficient `-N_g a^3 mu Delta^2/4`, where

\[
\mu=2y[\beta_1+\beta_2(c+y)+\beta_3cy].
\tag{8}
\]

For the present beta1 family, `mu=yP`. Direct extrinsic-curvature evaluation
gives the two time-kinetic terms; a literal spatial curvature calculation for
`diag(exp(gamma(z)),exp(-gamma(z)),1)` gives `R_3=-(gamma_{,z})^2/2`.
The resulting unit-normalization physical-time quadratic action is

\[
S_T^{(2)}=\frac18\int dT\,a^3\left\{
M^2[\dot\gamma_g^2-q^2\gamma_g^2]
+\frac{M^2y^3}{c}[\dot\gamma_f^2-(c/y)^2q^2\gamma_f^2]
-\mu(\gamma_f-\gamma_g)^2\right\},\quad q=k_{\rm com}/a.
\tag{9}
\]

Both homogeneous canonical matter fields have exactly unit determinant under
these traceless exponential spatial perturbations. They add no tensor
pressure Hessian to (8). No external tensor source or state is being imposed.
This assertion concerns the specified scalar action, not a generic material
with a linear TT anisotropic stress.

For a transverse relative ADM shift `s`, a direct two-by-two square root has
trace `sqrt((c+y)^2-a^2 y^2 s^2/N_g^2)`. Expanding its elementary polynomials
gives the literal quadratic coefficient

\[
C_s=\frac{a^5y^2P}{2N_g(c+y)}.
\tag{10}
\]

The EH vector squares have positive coefficients. In a normalization where
their action is `A k^2(E'-s_g)^2+B k^2(F'-s_f)^2+C_s(s_f-s_g)^2`, elimination
of the two shifts gives

\[
K_{E-F}=\frac{AB C_s k^2}{ABk^2+(A+B)C_s},\qquad A,B>0.
\tag{11}
\]

For `C_s<0`, it is negative once the denominator is positive at sufficiently
large `k`; for `C_s>0`, it is positive for every nonzero `k`. This statement
alone says nothing about the remaining gradient, scalar constraints, a
below-cutoff instability rate, or the applicability of arbitrarily high `k`.
At the center the literal coefficients are `-128/3` for `c=1` and `8/3` for
`c=4` in `M=tau=N_g=a=1`. The signs persist throughout the admitted local
intervals because `P` has the sign of `c-y` there.

The full two-TT principal squared speeds relative to actual `g` are `1` and
`c^2/y^2`. For `c>2`, the latter exceeds one everywhere on this local interval.
The leading locked ansatz has squared speed

\[
c_L^2=\frac{1+cy}{1+y^3/c},\qquad
c_L^2-1=\frac{y(c^2-y^2)}{c+y^3}>0.
\tag{12}
\]

Equation (12) is an exact coefficient for a locked configuration, not a
controlled integrated-out light mode. Such a claim would require the full
time-dependent operator, initial data, physical source and error budget.

## 4. Exact moving-weight canonical map

Use bars implicitly for `M=tau=1`, and prime for `d/du`. Define

\[
K_1=a^3/8,\quad K_2=b^3/(8c),\quad K_\Sigma=K_1+K_2,
\quad K_R=K_1K_2/K_\Sigma,\quad w_i=K_i/K_\Sigma,
\]
\[
v=w_1\gamma_g+w_2\gamma_f,\quad \Delta=\gamma_f-\gamma_g,
\quad l=f_\Sigma v,\quad Q=f_R\Delta,
\quad f_\Sigma=\sqrt{2K_\Sigma},\quad f_R=\sqrt{2K_R}.
\tag{13}
\]

Restoring physical fields multiplies both `f` by `M`; the action in `u` then
has overall factor `1/tau`. Equivalently use physical time, with connections
divided by `tau` and squared frequencies divided by `tau^2`. This preserves
the actual physical amplitudes and does not identify `u` momentum with
physical-time canonical momentum without its factor of `tau`.

Set `theta_i=f_i'/f_i`, `N_i=theta_i'+theta_i^2`,
`omega=w_2'/(2 sqrt(w_1w_2))`, and

\[
D=\sqrt{w_1w_2}(c^2/y^2-1),\quad
c_H^2=w_2+w_1c^2/y^2,\quad \bar q^2=\kappa^2/a^2,
\quad \kappa=\tau k_{\rm com},
\]
\[
\mathfrak m^2=\frac{\tau^2\mu}{M^2}(1+c/y^3)
=\frac{2y h'(y^3+c)}{c(c-y)}.
\tag{14}
\]

Here `kappa` is momentum (the canonical module's `kbar` symbol), distinct
from the earlier scalar kinetic profile `kbar`. The last quantity is an
algebraic relative mass coefficient, not an already
proven rolling gap. From the literal kinetic action the exact transformed
Lagrangian, before a time boundary is removed, is

\[
\mathcal L=\tfrac12\bigl[(l'-\theta_\Sigma l-2\omega Q)^2
+(Q'-\theta_R Q)^2
-\bar q^2c_L^2l^2-2\bar q^2DlQ
-(\mathfrak m^2+\bar q^2c_H^2)Q^2\bigr].
\tag{15}
\]

The normalization boundary is `-(theta_Sigma l^2+theta_R Q^2)'/2`. After it
is retained or accounted for, the two equations are

\[
l''+V_{LL}l+\mathcal B^*Q=0,\qquad
Q''+V_{HH}Q+\mathcal B l=0,
\tag{16}
\]
\[
V_{LL}=\bar q^2c_L^2-N_\Sigma,\quad
V_{HH}=\mathfrak m^2+\bar q^2c_H^2-N_R-4\omega^2,
\]
\[
\mathcal B=2\omega\partial_u-2\omega\theta_\Sigma+\bar q^2D,
\quad
\mathcal B^*=-2\omega\partial_u-2\omega'-2\omega\theta_\Sigma+\bar q^2D.
\tag{17}
\]

The formal-adjoint surface term is `(2 omega l Q)'`. In particular, even if
`omega=0` at a center, its derivative cannot be discarded. A future retarded
elimination would be `Q=Q_hom-G_R B l`, retaining homogeneous heavy data and
the chosen retarded inverse. It is not obtained by interpreting a causal
inverse as a symmetric single-copy action.

For this family, exactly

\[
\theta_\Sigma=\frac{6u(cd^{12}-8)}{d(cd^{12}+8)},\quad
\theta_R=-\theta_\Sigma,\quad
\omega=-\frac{24\sqrt{2c}\,u d^5}{cd^{12}+8}.
\]

At `u=0`,

\[
\mathfrak m_0^2=\frac{16(c+8)}{c(c-2)},\quad
N_{\Sigma0}=\frac{6(c-8)}{c+8},\quad
N_{R0}=\frac{6(8-c)}{c+8},\quad
\omega_0=0,\quad \omega'_0=-\frac{24\sqrt{2c}}{c+8}.
\tag{18}
\]

For example at `c=4`, the zero-momentum heavy diagonal is `24-2=22`, while
the locked squared speed is `3`. Neither center number is a full rolling
gap or a physical light-only matching theorem.

The inverse physical map is explicit:

\[
\gamma_g=l/f_\Sigma-w_2 Q/f_R,\qquad
\gamma_f=l/f_\Sigma+w_1 Q/f_R.
\tag{19}
\]

It must be retained when preparing physical data or a physical source. Both
canonical fields are tensor amplitudes; the scalar clock perturbations are
not being normalized or solved by (13).

## 5. Exact nonadiabatic ratios and shrinking inner limit

Put `delta=c-2=epsilon^2`, `u=epsilon x`, with `epsilon>0`. For a fixed finite
`kappa`, the explicit algebraic coefficient can be written

\[
\mathfrak m^2(u,c)=
\frac{16(1-u^2)(8+cd^{12})}{c d^{14}(cd^4-2)}.
\tag{20}
\]

Its center grows as `80/delta`, but `c-y=delta+8u^2+O(u^4)` varies on the
shrinking physical time scale `sqrt(delta) tau`. Literal derivatives give

\[
\frac{(\mathfrak m^2)''_0}{(\mathfrak m_0^2)^2}
=-\frac{c(7c^2+146c-240)}{8(c+8)^2}\ \longrightarrow\ -\frac15.
\tag{21}
\]

The coefficient is even in `u`. For its positive square root at `c>2`,
`omega_m''/omega_m^3` is half (21) at the center and tends to `-1/10`.
These are dimensionless physical adiabatic ratios: the powers of `tau`
cancel. `omega_m` here is a frequency diagnostic, not the rotation `omega`.

The actual canonical correction `V_HH-mathfrak m^2` is analytic at
`(u,c)=(0,2)`. Its value and second derivative there are exactly

\[
\kappa^2-18/5,\qquad -24\kappa^2/5-252/5.
\tag{22}
\]

Thus the finite correction cannot change the double-pole/squared-simple-pole
ratio (21). The actual canonical diagonal has the same limits `-1/5` and
`-1/10`, for its real positive square root at sufficiently small `delta`.
The independent audit also differentiates the full canonical diagonal
directly at zero momentum. No invariant spectral-gap theorem follows from
the diagonal or its sign; the mixed operator remains part of the dynamics.

Taking all `u` derivatives in (17) before changing the time variable, then
multiplying (16) by `epsilon^2`, gives coefficient limits

\[
\epsilon^2 V_{HH}\to\frac{80}{1+8x^2},\quad
\epsilon^2 V_{LL}\to0,\quad
\pm2\epsilon\omega\to0,
\]
\[
\epsilon^2(\bar q^2D-2\omega\theta_\Sigma)\to0,\quad
\epsilon^2(\bar q^2D-2\omega'-2\omega\theta_\Sigma)\to0.
\tag{23}
\]

Hence the complete limiting two-field equations, not just a frozen mass
estimate, are

\[
l_{,xx}=0,\qquad Q_{,xx}+\frac{80}{1+8x^2}Q=0.
\tag{24}
\]

The normalization itself is nonsingular on this inner chart:
`w1->1/5`, `w2->4/5`, `f_Sigma->M sqrt(5)/2`, `f_R->M/sqrt(5)`.
The physical projection approaches `gamma_g=2(l-2Q)/(M sqrt(5))`, so the
relative mode is not projected away. This is an exact field-map statement,
not a specified source response or a comparison of equal physical data.
At `epsilon=0` the coefficient functions diverge; there is no regular
finite-coefficient parent action at that parameter value.

## 6. Uniform convergence on each fixed inner interval

Fix `X<infinity`, any bounded finite momentum set, and an initial `x_*` in
`[-X,X]`. Restrict `0<epsilon<=min(1,1/[10 max(1,X)])`; then the entire
inner interval is within the certified clock domain and `2<c<=3<=4`.
After substitution, the only apparently singular denominator from the
mass is removable:

\[
\frac{(2+\epsilon^2)(1+\epsilon^2x^2)^4-2}{\epsilon^2}
=1+8x^2+O(\epsilon^2).
\tag{25}
\]

The exact quotient is a polynomial in `epsilon^2,x^2` with nonnegative
coefficients and constant-in-`epsilon` term `1+8x^2>=1`. The coefficientwise
Fraction engine verifies both division and leading term. Other denominators
are positive powers of `d`, `c`, or `cd^12+8`, and the radical `sqrt(2c)` is
on its positive analytic branch. Every scaled coefficient in (23) therefore
extends smoothly through `epsilon=0` on the fixed compact interval, uniformly
with any fixed finite number of `x` derivatives. This proves uniform, not
merely pointwise, coefficient convergence there.

For `Z=(l,Q,l_x,Q_x)`, (16) is a first-order system `Z_x=A_epsilon(x)Z`.
The preceding continuation gives `sup||A_epsilon-A_0||->0` and a uniform
finite coefficient bound on this compact set. Variation of constants and
Gronwall then imply uniform convergence of its Cauchy propagator to that of
(24), and hence of solutions for convergent initial data. For example the
propagator difference is bounded by
`2X exp(4X C_X) sup||A_epsilon-A_0||` when both coefficient norms are at most
`C_X`; this is a qualitative compact-set estimate, not a computed numeric
epsilon threshold. The same argument allows a uniformly convergent bounded
forcing specified in these normalized variables.

Nothing here makes `X` increase as `epsilon` decreases. In physical units
the window shrinks to `|T|<=epsilon tau X`. A fixed physical CD window, a
long-time or momentum-growing regime, and an actual physical-source error
would require additional uniform bounds and source/data dictionaries.

## 7. A derivative-omission control, not a matching exclusion

For the limiting operator with constant normalized forcing,

\[
R_{,xx}+\frac{80}{1+8x^2}R=1,
\]

the exact particular solution is `(1+8x^2)/96`. Algebraic inversion instead
gives `R_alg=(1+8x^2)/80`, which leaves residual `1/5`; their ratio is `5/6`.
The residual and ratio are exact polynomial identities. The two particular
solutions have not been assigned equal retarded initial data, and this
constant forcing has not been identified with a conserved physical TT
stress. Accordingly this is an omission control for a nonvanishing
derivative term, not an observable matching exclusion.

Likewise the locked center excess satisfies
`(c_L^2(0)-1)/(c-2)->4/5`. Algebraic derivative effects need not be negligible
at that same small parameter. A genuine light-only S6 calculation still has
to prepare heavy data, retain the physical map and source, and bound the
finite-time causal inverse. The present result identifies that obligation;
it neither asserts its success nor proves that every nonadiabatic reduction
must fail.

## 8. Evidence and theorem boundary

The primary engine differentiates the literal coframe density, square root,
curvature and canonical maps. A separate Fraction/Taylor engine recomputes
actual background equations and center derivatives, and a sparse polynomial
engine expands elementary polynomials and the removable inner denominator.
The verifier compares their outputs, not just pre-entered verdicts. A
separately authored audit supplies additional independent covariant/action
and full canonical controls. Exact arithmetic plus the written compactness
and inverse-function arguments are the evidence boundary; no formal proof
assistant or exhaustive numerical scan is claimed.

The result demonstrates an actual local positive-canonical-kinetic bounce
outside the frozen constant-beta/separate-stress premise, then computes a
specific tensor obstruction to inferring adiabaticity from mass height.
It does not certify full scalar constraints or health, complete-bounce
matching, cutoff, loop closure, UV completion, or the original DHOST
operators. It does not change the completed scoped photon objective or
frozen linear classification. Original P8 remains open.
