# S6.14 proof: actual matter cone and the surviving symmetric extension

## 1. Action, conventions and hypotheses

Retain the full action and coframe gradient convention certified by S6.13:

\[
S=-\frac G2\int\det e\,R[g]-\frac F2\int\det v\,R[f]
-2\int\det u\,[B+\operatorname{tr}(u^{-1}(p_g e+p_f v))]
-2\int(b_g\det e+b_f\det v)+\epsilon S_m[h,\psi].
\]

Here g=e^T eta e, f=v^T eta v and the **actual physical metric** is
h=u^T eta u. G,F have mass dimension two; B,p_g,p_f,b_g,b_f and rho,p
have dimension four. Frames and epsilon are dimensionless. All parameters
are constant. The endpoint terms are an explicitly different action from
the restricted source model screened in S6.13, and do not enter the u equation.

The curvature convention is the original P8(b) convention:
R^rho_(sigma mu nu)=partial_mu Gamma^rho_(nu sigma)-partial_nu
Gamma^rho_(mu sigma)+..., R_B=-6(Hdot+2H²), and G_B00=+3H².
EH=-G R_B/2 has a positive tensor kinetic term. This differs from the
R_FK=+6(...), G_FK00=-3H² convention in the P8(a) FK calculations; the
opposite curvature sign cannot be copied while retaining this EH coefficient.
The old P8 formulation and S6.13 are hash-pinned by the verifier.

Use a regular common spatially flat homogeneous chart

\[
e=\mathrm{diag}(n_e,a_e,a_e,a_e),\quad
v=\mathrm{diag}(n_v,a_v,a_v,a_v),\quad
u=\mathrm{diag}(n_u,a_u,a_u,a_u),
\]

with all six entries positive. The diagonal positive chart satisfies the
parent Lorentz/root conditions. The claims apply on actual smooth solutions,
not on arbitrary points satisfying only the auxiliary equation.

For canonical homogeneous scalars, allow a positive field-space metric
K_IJ(psi) and potential V(psi). Then

\[
\rho+p=K_{IJ}D_T\psi^I D_T\psi^J,\qquad
dT=n_u dt,\qquad n_h=\epsilon(\rho+p)\geq0.
\]

Strict rolling means n_h>0, including any nonzero frozen free-chi velocity.
An epsilon prefactor can be absorbed into canonical matter normalization;
it cannot be taken to zero at fixed physical chi action and velocity to
erase n_h. Negative quantum null stress, fluids with an independent TT
response, derivative matter couplings and extra u kinetic terms are not
covered by the tensor-action corollary. The background cone identity itself
only uses the stated homogeneous stress components.

## 2. Full u equation gives a weighted cone identity

The covariant Euler map is varied before any FLRW restriction. Its matter
gradient on the homogeneous chart is

\[
J_u=\mathrm{diag}(-a_u^3\rho,n_u a_u^2 p,n_u a_u^2 p,n_u a_u^2 p).
\]

Set P_g=p_g a_e/a_u, P_f=p_f a_v/a_u. Normalize E_u,00 by -2a_u³
and any E_u,ii by -2n_u a_u²; these nonzero factors give respectively

\[
B+3(P_g+P_f)=-\epsilon\rho/2,\qquad
B+P_g(c_e+2)+P_f(c_f+2)=\epsilon p/2,
\tag{1}
\]

where c_e=a_u n_e/(n_u a_e), c_f=a_u n_v/(n_u a_v).
These are the two Einstein tensor propagation speeds measured in physical
h time and spatial length. Subtract (1), without a Hubble quotient, lapse
gauge change, matter expansion or division by B:

\[
P_g(c_e-1)+P_f(c_f-1)=n_h/2.
\tag{2}
\]

For p_g,p_f>0, both weights are positive. Consequently

\[
\max(c_e,c_f)\ \geq\ 1+\frac{n_h}{2(P_g+P_f)}>1
\tag{3}
\]

when matter rolls. Only at least one speed is constrained: it is incorrect
to infer that both exceed one. Neither the Einstein coefficients, B nor
endpoint cosmological terms enter (2). The result holds at a bounce without
using any quotient of Hubble rates.

This is a full-parent principal statement. A large relative algebraic mass
can change which channels remain in a finite-band light theory. Equation
(3) alone does not decide an asymmetric light-only reduction. For example,
the purely quadratic two-field model K=diag(100,1), gradient=diag(25,4)
has principal squared speeds 1/4 and 4, but a positive very heavy relative
spring locks the light squared speed to 29/101. This countercontrol is not
asserted to be a solution of our parent; it diagnoses the invalid inference.

## 3. Literal tensor elimination and actual physical source

For one real unit polarization E, tr E=0 and tr E²=1, perturb the spatial
coframes by exp(gamma_i E/2). At quadratic order, tr(u^-1 e) contributes
(a_e/a_u)(gamma_e-gamma_u)²/8. The three determinants are fixed along
these exponential tensor paths. The homogeneous canonical scalar kinetic
term depends on u's time component, not on these spatial tensors. It has
no independent linear TT anisotropic perturbation. Thus neither its action
nor the endpoint terms add TT derivatives or a tensor potential here.

Divide the algebraic TT density by the background physical volume. With
the explicitly normalized external source S_j=int sqrt|h| j gamma_u,

\[
L_{\rm alg}/\sqrt{|h|}
=-\frac14[P_g(\gamma_e-\gamma_u)^2+P_f(\gamma_v-\gamma_u)^2]
+j\gamma_u.
\tag{4}
\]

The exact linear source-dependent auxiliary response is

\[
\gamma_u=\bar\gamma+\frac{2j}{P_g+P_f},\qquad
\bar\gamma=\frac{P_g\gamma_e+P_f\gamma_v}{P_g+P_f}.
\tag{5}
\]

For positive links the tensor denominator is nonzero, independent of any
unproved scalar auxiliary inverse. Completing the square yields

\[
L_{\rm alg}/\sqrt{|h|}
=-\frac{P_gP_f}{4(P_g+P_f)}(\gamma_e-\gamma_v)^2
+j\bar\gamma+\frac{j^2}{P_g+P_f}
-\frac{P_g+P_f}{4}\left(\gamma_u-\bar\gamma-\frac{2j}{P_g+P_f}\right)^2.
\tag{6}
\]

This identity retains the external probe's algebraic contact. The derivative
of the reduced source functional with respect to j gives (5), not merely
bar-gamma. Its local contact does not remove the common propagating response.

For precise stress normalization, with the canonical stress convention
delta S_probe=-(1/2)int sqrt|h| T^mu nu delta h_mu nu and
delta h_ij=-a_u² gamma_u E_ij, define the physical orthonormal TT amplitude
Pi=T^(hat i hat j)E_ij. Then **j=Pi/2**. With instead a plus sign in the
definition of the variational source stress, j=-Pi/2. The literal j action
in (4) fixes all replay coefficients, irrespective of this naming convention.
A spatial Fourier TT probe with zero density/flux, zero trace and transverse
spatial divergence is separately conserved on the FLRW background for any
smooth time profile: the temporal connection term is proportional to its
trace. It is not the internal scalar source, and its amplitude need not obey
a positivity condition. The response is linear in the external probe and
second order in the perturbation/source generating functional.

## 4. Exchange symmetry places the faster cone in the physical common mode

Now impose G=F, p_g=p_f=q>0, b_g=b_f and the exact invariant background
branch e=v=r. Write r=diag(n,a,a,a), u=diag(N,A,A,A)r with N,A>0.
Then P_g=P_f=q/A and c_e=c_f=A/N. Equation (2) becomes

\[
\frac1N-\frac1A=\frac{n_h}{4q},\qquad
c_T=\frac AN=1+\frac{A n_h}{4q}>1.
\tag{7}
\]

No vacuum branch expansion or small-source approximation is used. For
gamma=(gamma_e+gamma_v)/2, delta=(gamma_e-gamma_v)/2, (5) gives
gamma_u=gamma+A j/q. The common tensor is therefore an actual physical
metric perturbation, with unsuppressed source coefficient one. Odd tensors
have zero linear physical projection; exchange symmetry forbids common/odd
mixing even on a time-dependent background. Unequal Einstein coefficients
would instead introduce a velocity cross term and are excluded from this
stronger corollary, though not from (2).

The full source-free quadratic tensor action, with coordinate-time dots, is

\[
S_T^{(2)}=\int dt\,d^3x\,\frac{a^3}{8n}
\left\{2G[\dot\gamma^2+\dot\delta^2]
-\frac{2Gn^2k^2}{a^2}(\gamma^2+\delta^2)
-4q n^2NA^2\delta^2\right\}.
\tag{8}
\]

The Einstein-Hilbert derivative terms follow directly from the ADM tensor
shear and spatial curvature; exponential traceless tensors do not change
the background volume. Isotropy separates these modes from all linear
scalar/vector constraints. The positive-root coframe constraints hold to
the required tensor order. No independent matter TT response is dropped.

Changing to actual physical variables dT=Nn dt, a_h=Aa gives for the common
mode, and identically for the derivative part of the relative mode,

\[
S_\gamma^{(2)}=\frac18\int dT\,d^3x\,a_h^3
\left[G_{T,h}(D_T\gamma)^2-F_{T,h}\frac{k^2}{a_h^2}\gamma^2\right],
\quad G_{T,h}=\frac{2GN}{A^3},\quad F_{T,h}=\frac{2G}{NA}.
\tag{9}
\]

Both coefficients are positive; their ratio is c_T²=A²/N². The source
terms are int dT a_h³[j gamma+A j²/(2q)]. Thus eliminating the relative
tensor cannot repair (7) on this exactly symmetric branch. This is stronger
than an observation about a full-parent heavy cone: the common physical
channel already has it.

Do not identify a cosmological channel with a literal stationary mass pole.
Let f_c=sqrt(a_h³ G_T,h)/2, Y=f_c gamma and theta=D_T log f_c. The exact
canonical action includes the boundary -(theta Y²)'/2 and squared frequency
c_T² k²/a_h²-f_c''/f_c. Time dependence of this normalization does not change
the k² coefficient. The relative algebraic mass is 2qA²/(GN)>0 but its
canonical frequency also contains the normalization term. No Green-function,
adiabaticity or heavy-gap theorem is inferred from that algebraic mass.

## 5. A full local rolling solution makes the corollary non-vacuous

Use the **same** separately named flat-vacuum extension as S6.13:
B=-6q, b_g=b_f=-q, q,G,epsilon>0, one free canonical scalar V=0.
Choose r proper time t (n=1). On A>1 the auxiliary equations give

\[
\rho=p=\frac{12q(A-1)}{\epsilon A},\quad
N=\frac{A}{6A-5},\quad
H_r^2=\frac{2q(A^3-1)}{3G}.
\tag{10}
\]

Select the positive square root H_r. For any A0>1 and a0>0 solve

\[
A'=-\frac{6H_r A(A-1)}{6A-5},\qquad
a'=H_r a,\qquad
\psi'=N\sqrt{2\rho}.
\tag{11}
\]

The vector field is analytic on this open domain, including the positive
square roots, so the ordinary local existence/uniqueness theorem applies.
All scales/lapses and the scalar kinetic term remain positive for a small
interval by continuity. This is a local actual solution, not a prescription
of unrelated instantaneous acceleration and density.

Differentiating the constraint, or substituting (11) directly, gives

\[
H_r'=-\frac{6qA^3(A-1)}{G(6A-5)},\qquad
3G H_r^2=2qA^3-2q,\qquad
G(2H_r'+3H_r^2)=2qNA^2-2q.
\tag{12}
\]

These are both full e/v lapse and spatial equations in the B curvature
convention. The u equations are (10). Finally the actual scalar current
(Aa)^3 D_T psi is conserved: its logarithmic t derivative is
3(H_r+A'/A)+rho'/(2rho)=0. All off-diagonal coframe equations vanish by
isotropy. The root audit independently substitutes the full covariant
Einstein and matter tensors into every component of all three Euler maps.

At A0=2 one has N0=2/7, rho0=6q/epsilon,
H_r0²=14q/(3G), H_r0'=-48q/(7G) and c_T0=7.
For q=G=1, G_T,h=1/14 and F_T,h=7/2. Arbitrarily near-vacuum choices
A0>1 also exist and have c_T>1 tending to one as A0 tends to one.
The physical Hubble rate is H_h=H_r/A>0; **this is not a bounce**.
It does not prove a global interval, a subcutoff signal advance, scalar/vector
health or the existence of an original CD/M1 matching trajectory.

## 6. Exact quantitative gaps and conditional correction budgets

At any symmetric rolling point define delta=A n_h/(4q)>0. Then

\[
c_T-1=\delta,\quad c_T^2-1=2\delta+\delta^2,\quad
F_{T,h}-G_{T,h}=G_{T,h}\delta(2+\delta)>0.
\tag{13}
\]

Suppose a separately specified physical-frame two-derivative correction has
|delta G_T|<=E_G, |delta F_T|<=E_F, and G_T,h-E_G>0. Then

\[
F_{\rm new}-G_{\rm new}\geq
G_{T,h}\delta(2+\delta)-E_G-E_F.
\tag{14}
\]

A positive right side preserves the superluminal cone. Hence E_G+E_F at
least the uncorrected positive gap is necessary, not sufficient, to obtain
a positive-kinetic subluminal correction within that representative. The
API rejects floats, nonfinite parameters, wrong signs and missing exact
domains. It does not manufacture a bound on actual omitted operators,
background displacement, loops, nonlocal response or dispersive terms.
Those require a physical field/source dictionary and a separate calculation.

For a conditional connection to the frozen CD/M1 free scalar, preserve its
actual canonical action and velocity D_T chi=M_target/(10tau) at the bounce.
If all additional matter is positive canonical, n_h>=M_target²/(100tau²).
If one additionally identifies the extension's vacuum common Planck
coefficient 2G with M_target² and its vacuum relative mass m0²=2q/G, then

\[
c_T-1\geq\frac{A}{100(m_0\tau)^2}.
\tag{15}
\]

The one-scalar example is not thereby identified with the old clock-plus-chi
theory. These are explicit hypothetical matching assumptions, not a supplied
matching map or an inferred cutoff. A large finite mass can make the excess
small, but cannot change its strict sign at a rolling point in this exact
positive-link action. Higher-derivative/EFT or quantum completions may change
the coefficients and require their own quantitative test.

## 7. Controls and the precise conclusion

Negative links are not silently excluded by a vacuum-health assertion.
The **actual** mixed-sign extended action p_g=-2,p_f=1,B=3,b_g=2,b_f=-1
has e=v=u=I, zero matter and all three full flat Euler matrices zero.
Its exact vacuum auxiliary solution is u=2e-v. The eliminated potential
has relative TT spring q_eff=2p_g p_f/(p_g+p_f)=4>0 and positive Einstein
kinetics yield m_FP²=4(1/G+1/F)>0. This is only a flat quadratic control,
not a scalar/rolling/full-health claim. It suffices to show why positivity
of both links is an explicit premise, not deduced from that flat mass.

A separate mixed-sign **auxiliary-only point** has a_e=a_v=a_u=n_u=1,
n_e=1/2,n_v=1,p_g=-2,p_f=1,B=5/2,epsilon=1,rho=p=1.
It solves both u equations with c_e=1/2,c_f=1 and n_h=2. It is not claimed
to solve e/v equations; the root audit verifies that those fail for the
stated static frames. This diagnoses the positive-weight premise in (3).

Nonrolling n_h=0 has no strict gap; the positive-link flat vacuum is a
nonempty control. Zero link weights can destroy the particular auxiliary
inverse and are not admitted by (7). Negative lapse/scale charts, anisotropic
matter perturbations, unequal/exchange-breaking branches, added u EH or
other operators, and quantum stress are outside the corresponding theorem.
The scalar/vector spectrum is not decided. The unchanged original CD/M1
action has the physical-cone requirement; this exact symmetric canonical
extension cannot itself supply that requirement at rolling points. The
result is not a proof of a closed causal curve, a universal UV inconsistency,
or exclusion of all allowed S6 light-only matching constructions. No
stationary positivity amplitude or cutoff is inferred. S6 and P8 remain open.
