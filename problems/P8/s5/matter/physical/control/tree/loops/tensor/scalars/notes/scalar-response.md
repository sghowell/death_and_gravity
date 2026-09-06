# Coupled scalar Weyl response: geometry, phase reduction and finite-band proof

## 1. Source and convention boundary

The calculation concerns a new, constant-coefficient local candidate in the
frozen physical matter frame, not the complete renormalized matter loop.
The CD/M1 action and S5.10 lineage remain immutable. Signature is `+---`;
changing the sign convention of the Riemann tensor does not change `C^2`.
The matter perturbation is `delta chi=M*s`, so the baseline scalar action
has an overall `M^2`, while the added Weyl functional has coefficient `c_C`.
Consequently the equation/phase expansion parameter has units of length
squared: `beta=c_C/M^2`.

As a primary-source cross-check, Salvio, *Inflationary Perturbations in No-Scale
Theories*, [arXiv:1703.08012](https://arxiv.org/pdf/1703.08012), eqs. (2.18),
(3.1) and (4.1), uses the same metric signature and gives the scalar Weyl
action proportional to the squared spatial Laplacian of `Phi+Psi`. His
coefficient is `c_C=-1/(2 f_2^2)` modulo the Euler density; translating it
gives `4 c_C/3` below. His Einstein-frame background and extra-mode quantum
interpretation are not transplanted to CD/M1.

For perturbative reduction, see Parker and Simon,
[gr-qc/9211002](https://arxiv.org/pdf/gr-qc/9211002), and the field-redefinition
discussion in Solomon and Trodden,
[1709.09695](https://arxiv.org/pdf/1709.09695), section II.C. Glavan,
[1710.01562](https://arxiv.org/pdf/1710.01562), explains relevant distinctions
between reduced actions/equations and their exact resummations. The identities
here are derived directly rather than inferred from those sources' models.

## 2. Full scalar Weyl quadratic functional

Divide the physical metric by its background `a^2`. Its linear scalar part is

\[
h_{00}=2A,\quad h_{0i}=-\partial_i B_c,\quad
h_{ij}=-2\zeta_{\rm phys}\delta_{ij}-2\partial_i\partial_j E.
\]

Under a conformal-time displacement `T`, `A` changes by `-T'-calH*T`,
`zeta_phys` by `-calH*T`, and `B_c-E'` by `T`. Therefore

\[
\mathcal L=A-\zeta_{\rm phys}+(B_c-E')'
\]

is gauge invariant. A direct four-index contraction with arbitrary scalar
Hessian jets (take the wavevector along `z` by rotation) gives the pointwise
identity

\[
C_{\rm lin}^2=\frac43
(A_{,zz}-\zeta_{{\rm phys},zz}+B_{c,\eta zz}-E_{,\eta\eta zz})^2.
\]

Equivalently the scalar magnetic Weyl tensor vanishes and its electric part
is one half the traceless spatial Hessian of `Lensing`, up to the irrelevant
overall curvature sign. Both routes include the shift and scalar shear.
For a single spatial Fourier mode of squared norm one, conformal invariance
of `sqrt(|g|)*C^2` in four dimensions yields

\[
\Delta S_S^{(2)}=\frac{4c_C}{3}\int d\eta\,k_{\rm com}^4\mathcal L^2
=\frac{4c_C}{3}\int dt\,a^3q^2\mathcal L^2.
\]

No background or scalar constraint equation was used in this geometry.
`C=0` on FLRW removes the background and first variation, not this second
variation. The positive conformal-scalar control `A=zeta`, zero shift/shear
has zero Weyl, whereas a static inhomogeneous lapse has a nonzero result.

In the old clock and scalar spatial gauge, `E=0`. Since physical
`b=a^2*beta_shift`, `B_c=b/a` and `B_c'=bdot-Hb`. The auxiliary spatial map
is `zeta_phys=v_old+delta*n`, so

\[
\boxed{\mathcal L=(1-\delta)n-v_{\rm old}+\dot b-Hb.}
\]

Omitting `delta*n` changes the physical operator; it is not a convention.

## 3. Regular auxiliary reduction before any velocity-chart inversion

Suppress the overall `M^2*a^3`. After the pinned matter boundary, the old
quadratic Lagrangian is

\[
\begin{split}
L_0={}&-3\dot v^2+\Sigma n^2+6\theta n\dot v
+\tfrac12\dot s^2+wn\dot s-3l\dot v\,s\\
&+q[v^2+2\Lambda nv-\tfrac12s^2+2\theta nb-2\dot v b-lbs],\\
w={}&l(3\delta-1),\qquad J=\Sigma+3\theta^2-w^2/2.
\end{split}
\]

Legendre-transform only `vdot,sdot`, whose Hessian `diag(-6,1)` is
invertible independently of the crossing:

\[
p=-6\dot v+6\theta n-3ls-2qb,\qquad P_s=\dot s+wn.
\]

The resulting Hamiltonian has exact auxiliary Hessian

\[
H_{0,(n,b)(n,b)}=\operatorname{diag}(-2J,-2q^2/3),\qquad
\det=4Jq^2/3.
\]

The stationary solution, with no `theta` or `Lambda` denominator, is

\[
\begin{split}
R&=-\theta p/2+\Lambda qv+wP_s/2-3l\theta s/2,\\
n_0&=-R/J,\quad b_0=-p/(2q),\\
h_0&=-qv^2+P_s^2/2+(q/2-3l^2/4)s^2-lps/2+R^2/J.
\end{split}
\]

The code checks this full Hamiltonian against the independently pinned S5.6
unitary and gamma targets. In gamma variables
`p=-2q*Q_g`, `v=P_g/(2q)`, with the indispensable density/time generator
`-H*Q_g*P_g`. These are phase charts. The positivity restrictions on the
separate finite-q velocity charts do not invalidate regular stationary phase
elimination at `theta=0`.

At first beta order, substitution of the old stationary `n0,b0` in the new
functional is sufficient: the linear change of the old action in auxiliary
directions vanishes identically. This is a perturbative stationary-action
identity, not exact elimination of the new `bdot` dynamics or a choice of
its extra initial data.

## 4. Derivative phase map and rank-one coupled Hamiltonian

After old auxiliary substitution **off shell**, retain `pdot`:

\[
\begin{split}
\mathcal L_{\rm off}&=(1-\delta)n_0-v-\frac{\dot p+Hp}{2q},\\
E_p&=\dot p+3Hp-2q(v+\Lambda n_0),\\
\mathcal L_0&=(1-\delta-\Lambda)n_0-2(v+Hb_0),\\
\mathcal L_{\rm off}&=\mathcal L_0-E_p/(2q).
\end{split}
\]

Here `qdot=-2Hq`. On CD, `Lambda=1-3delta`, hence
`Lensing0=2delta*n0-2(v+H*b0)`. The old variational derivative in `v` is
`-M^2*a^3*E_p`, including the expanding momentum density. The new action is

\[
\beta a^3\left[\tfrac43q^2\mathcal L_0^2
-\tfrac43q\mathcal L_0E_p+\tfrac13E_p^2\right]
\]

after dividing by `M^2`. Its last two terms are removed through first order
by the full derivative field map

\[
v_{\rm old}=v+\beta(E_p-4q\mathcal L_0)/3,
\]

holding `p,s,P_s` unchanged at this stage. The old canonical one-form produces
the boundary `beta*a^3*p*(E_p-4q*Lensing0)/3`; it is not silently a canonical
map in the original variables. Omitting `E_p/3` leaves the exact action
residual `beta*a^3*E_p^2/3` off shell. On the perturbative branch
`E_p=O(beta)`, that map reduces to `v_old=v-4beta*q*Lensing0/3` through first
order, giving the **specified representative**

\[
\boxed{h_{\rm red}=h_0+\beta h_1,\qquad
h_1=-\tfrac43q^2\mathcal L_0^2.}
\]

It acts on the full four-dimensional two-scalar phase vector. For example,
at the bounce `delta=1/2`, `theta=H=0`, `lambda=-1/2`, `l=1/10` in local units,

\[
\mathcal L_0=(-2+q/(2J))v-P_s/(40J).
\]

Dropping the matter channel loses a nonzero exact term even at the crossing.
If `Lensing0=r^T X`, the correction Hessian is
`-(8/3)q^2*r*r^T`, rank one (or zero at a degenerate row), and its canonical
generator is Hamiltonian. The phase ODE is first order with four data; it
does not assert that the resummed original candidate has only four data.

## 5. Physical reconstruction, not a changed matter frame

For physical lapse/shift reconstruction, unlike the first-order effective
action, the auxiliary correction is required. Before the derivative phase
map the Weyl Euler derivatives divided by `a^3` are

\[
\delta_nL_1=\tfrac83q^2(1-\delta)\mathcal L,
\qquad \mathrm{EL}_bL_1=-\tfrac83q^2\dot{\mathcal L}.
\]

The second identity uses `(a^3*q^2)'=-H*a^3*q^2`; neglecting this weight
produces an incorrect extra Hubble term. Inverting the old diagonal Hessian
gives

\[
n_{1,\rm before}=\frac{4q^2(\delta-1)}{3J}\mathcal L_0,
\qquad b_1=4D_t\mathcal L_0.
\]

The lapse must also be pulled through `v_old=v+beta*v1`. Since
`partial_v n0=-Lambda*q/J`, the total first-order lapse is

\[
n_1=\frac{4q^2(\Lambda+\delta-1)}{3J}\mathcal L_0
=-\frac{8\delta q^2}{3J}\mathcal L_0\quad\text{on CD}.
\]

Thus `zeta_phys=v+delta*n0+beta*(v1+delta*n1)`. The matter field `s` itself
has no map at this step. Its **physical** normalized momentum density is
`pi_chi/(M*a^3)=l+P_s+3l*v_old`, so its first-order map includes `3l*v1`.
The pre-mixed auxiliary metric momentum is `p+3l*s`; these are the old
S5.6 boundary coordinates, not an invented matter decoupling.

There is a direct velocity check, beyond the auxiliary stationarity test:
the reduced `sdot`, with the reconstructed lapse, obeys the original matter
relation `P_s=sdot+w*n_phys` through first order. Differentiating the branch
`v_old` map, including `qdot`, similarly reproduces the original metric
relation `p=-6*vdot_old+6theta*n_phys-3l*s-2q*b_phys`. Both are exact symbolic
first-order identities; omitting reconstruction is not innocuous.

`D_t` in `b1` uses the old full coupled phase equations and all coefficient
derivatives. `old_flow_derivative` supplies this operator with explicit
Hubble, theta, J, delta and Lambda jets, `ldot=-3Hl`, `qdot=-2Hq`, and both
momentum-density drifts. Compact substitutions must satisfy
`Lambda_dot=-3*delta_dot`. In moving local units `qbar=ell^2*q`, the correct
drift is `ell*d(qbar)/dt=-6x*qbar`, not zero. These time weights are checked
against the pinned S5.7 coefficient module.

The gamma **canonical** coordinate remains `-p/(2q)`. It equals only `b0`,
not the corrected physical shift `b0+4beta*D_t Lensing0`. No physical metric
redefinition is used to discard the matching obligation.

## 6. Full coupled normalization and explicit representative ODE

Use fixed centre time and spatial units `ell0`, and define
`s_ell=ell/ell0`. In these units the unitary gravitational coordinate is `v`
and the gamma coordinate is `b0/ell0`; their conjugate densities scale as
`ell0*p` and `ell0^2*P_b`, respectively. The matter pair is
`s,ell0*P_s`. The expansion coefficient is `epsilon0=c_C/(M*ell0)^2`.

Take the pinned inverse canonical map, not a velocity replacement:

\[
Q_{\rm chart}=a^{-3/2}T^{-1}Y,\qquad
p_{\rm chart}=a^{-3/2}T^T(P-SY).
\]

`T=Tbar*diag(s_ell^-d,1)`, with chart weight `d=0` in unitary and `d=1`
in gamma; `S=Sbar/s_ell`. The exact old Lagrangian and Hamiltonian are

\[
L_0=\tfrac12|\dot Y+\Omega Y|^2-\tfrac12Y^TWY,
\quad H_0=\tfrac12(P^TP+Y^TWY)-P^T\Omega Y.
\]

The dot here is fixed dimensionless time. `Omega` is antisymmetric. If
`Lensing0=cQ^T*Qchart+cp^T*pchart`, define

\[
r_P=T c_p,\quad r_Y=T^{-T}c_Q-S^T r_P,\quad
\mathcal L_0=a^{-3/2}(r_Y^TY+r_P^TP).
\]

Volume cancels in `a^3*Lensing0^2`. With `r=(rY,rP)` and the standard
symplectic matrix `Jcan`, the actual reduced equation is

\[
\frac{dX}{d\sigma}=\left[
\begin{pmatrix}-\Omega&I\\-W&-\Omega\end{pmatrix}
-\frac{8\epsilon_0}{3}q^2\mathbb J\,r r^T\right]X.
\]

`canonical_generator` implements this factored expression. No principal-limit
`K=G`, eigenmode diagonalization, omission of `S`, or frozen-frequency
propagator enters it. In particular `P` is not `Ydot`.

## 7. Certified finite-band coefficient and generator bounds

Take `abs(t-t0)<=ell0/100`, normalize `a(t0)=1`, and impose centre momentum
`L<=k_com*ell0<=4U`, `L=10^11`, `U=10^12`. The old covering proof gives a
single valid chart over each full window: gamma for `abs(x0)<=9/50`, unitary
otherwise. It also gives `1/2<=a,s_ell<=2`, and hence

\[
q_{\rm fixed}\ge L^2/4>1,\quad
\bar q\ge L^2/16>10^{20},\quad
q_{\rm fixed},\bar q\le(16U)^2<10^{27}=:Q_*.
\]

The window avoids velocity-chart poles. This covers all finite centre times;
it does not follow one fixed comoving momentum uniformly to both infinities.
The old energy metric is `G=diag(W,I)` and
`qI/2 <= W <= 3qI/2` in fixed centre units. The exact free propagator has
norm at most `sqrt(11/7)<4/3` between times in the same window, measured in
the corresponding time-dependent old energy norms.

The compact CD coefficient bounds and `s_ell` give
`|H|,|theta|<=8`, `|l|<=1/5`, `|Lambda|<=1`, `delta<=1/2`, `J>=1/40`.
The exact CD rows are

\[
\begin{array}{c|cc|cc}
&c_{Q1}&c_{Q2}&c_{p1}&c_{p2}\\ \hline
u&-2-2\delta\Lambda q/J&3\delta l\theta/J&\delta\theta/J+H/q&\delta l\Lambda/J\\
g&-2\delta\theta q/J-2H&3\delta l\theta/J&-1/q-\delta\Lambda/J&\delta l\Lambda/J
\end{array}
\]

Positive Laurent-monomial arithmetic gives unitary row norms at most
`138q,172` and gamma row norms at most `432q,25`. Use the common bounds
`432q,172`. The pinned full whitening gives `||T||<46`, `||T^-1||<=20`, and
the exact old symmetric-boundary majorant is

\[
\|S\|\le2q C_S,\qquad C_S=42828424431/2000000.
\]

It is replayed from `coefficient_bound(z*Sbar)` in both charts; it is not
inferred from a few sampled matrices. Thus

\[
\|r_P\|\le46\cdot172<8000,\quad
\|r_Y\|\le[20\cdot432+2C_S\cdot46\cdot172]q<4\cdot10^8q,
\]

and the dual energy-row norm satisfies

\[
D_r:=r_Y^TW^{-1}r_Y+r_P^Tr_P<4\cdot10^{17}q.
\]

For the rank-one generator without `epsilon0`, its **exact** energy norm is

\[
\|A_1\|_G=\tfrac83q^2
\sqrt{(r_P^TWr_P+r_Y^Tr_Y)D_r}.
\]

The first factor squared is at most `lambda_max(W)*D_r`. Therefore

\[
\boxed{\|A_1\|_G<3\cdot10^{18}q^{7/2}<10^{114}.}
\]

The substantial momentum enhancement is retained. No claim that these
constants are sharp is needed. Exact rational margin checks certify all
numerical enlargements. The new normalized momentum Hessian is
`I-(8epsilon0/3)q^2*rP*rP^T`; the stated coefficient gate bounds its correction
norm below `1/2`. This establishes positivity of this chosen representative's
velocity Hessian on the band, not absence of heavy modes of an exact
higher-derivative system.

## 8. Finite-window evolution and first Born remainder

Let `U0` be the exact full coupled free propagator. In its interaction picture,
the perturbation norm is at most `(11/7)*10^114*|epsilon0|`. The full window
has dimensionless length `1/50`; thus the integral is bounded conservatively
by `mu=10^114*epsilon`, with `epsilon=|c_C|/(M*tau)^2` since `ell0>=tau`.
Require the independent input `epsilon<=10^-116`, so `mu<=1/100`.

The ordered-integral expansion, with no commutativity assumption, gives

\[
\|U_{\rm red}-U_0\|\le\tfrac43(e^\mu-1)
\le\tfrac43\frac\mu{1-\mu}<10^{115}\epsilon,
\]

and the difference after its first time-ordered/Born term is at most

\[
\tfrac43(e^\mu-1-\mu)
\le\tfrac43\frac{\mu^2}{2(1-\mu)}<10^{229}\epsilon^2.
\]

Norms are operator norms from the centre's old energy space to the observation
time's old energy space. These comparisons use equal **reduced canonical**
centre data. They are exact bounds for the selected linear reduced ODE;
the second estimate is not a bound on unknown beta-squared action terms or
on an exact branch of the unreduced equations.

## 9. Branch phase map and its inverse

There is also a useful controlled bookkeeping map. The branch map changes
only old `v`, by `-4epsilon0*q*Lensing0/3`, before the gamma and whitening maps.
In the unitary chart its normalized image column is
`[T*e1; S*T*e1]`, whose energy norm is at most `92q*(1+C_S)`. In gamma it
changes `P_g` by `-8epsilon0*q^2*Lensing0/3`, with image column
`[0; T^-T*e1]` of norm at most `20`. Since `sqrt(D_r)<10^9*sqrt(q)`, both give

\[
X_{\rm old}=(I+\epsilon_0F_{\rm map})X_{\rm red}+O(\epsilon_0^2),\quad
\|F_{\rm map}\|_G<3\cdot10^{15}q^{5/2}<10^{84}.
\]

Here `X_old` is the old normalization of `(v_old,s,p,P_s)`, not a new
normalization chosen to diagonalize the corrected operator. Set
`eta=10^84*epsilon<=10^-32`. The **chosen map** `I+epsilon0*Fmap` has inverse
norm at most `1/(1-eta)`. It is legitimate to use this inverse to impose equal
pulled-old-phase centre data for this chosen representative. Its propagator
is `A(t)*Ured(t,t0)*A(t0)^-1`, not `Ured` with unchanged data. Its deviation
from `U0` is bounded by

\[
\tfrac43\frac{(e^\mu-1)+\eta(e^\mu+1)}{1-\eta}
<6\mu+8\eta<10^{116}\epsilon.
\]

Treating this first-order map exactly is a declared higher-order
representative choice. It neither removes the formal `O(epsilon0^2)`
ambiguity nor turns equal old auxiliary phase data into equal physical
`zeta,n,b` data. The separate reconstruction in section 5 is indispensable.

## 10. Interpretation and verification boundary

For independent `|c_C|<=1`, `M*tau=10^324`, one has
`epsilon<=10^-648`, `mu<=10^-534`, reduced error at most `10^-533`, first
Born remainder at most `10^-1067`, branch map at most `10^-564`, and chosen
pulled-phase error at most `10^-532`. The scalar gate is intentionally much
stronger than S5.10's tensor-only gate; blindly using the latter fails an
explicit negative control.

The direct geometry, partial Legendre Hessian, stationary substitutions,
derivative map, auxiliary reconstruction, old Hamiltonian bridge and generic
normalization are exact symbolic identities. The coefficient bounds are
positive rational polynomial/Laurent majorants with replayed denominators;
they are not numerical sampling claims. Independent tests supply a second
four-index/gauge/action/rank-one route. The stationary-action and finite-time
energy/Dyson implications are written mathematical arguments, not formalized
proofs of PDE or UV statements.

These results advance the same candidate's two physical scalar channels at
first order. They do not establish full quantum superluminality, a UV front
velocity, a physical EFT ghost, complete candidate nonlinear/tree control,
vector/extra-mode health, a corrected vacuum or backreaction solution.
Finite matching coefficients and omitted nonlocal, anomaly, state and mixed
loop terms remain independently necessary. A constant coefficient is not a
spacetime-dependent renormalization scale, and isolated running cannot supply
the matching assumption used in the conditional estimates.
