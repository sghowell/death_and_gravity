# Exact tensor/vector derivation and the bounded obstruction

## 1. Variables, units and physical frame

Use the pinned +--- action and curvature conventions. Write g with lapse
Ng and scale a, f with lapse Nf and scale b, and define

\[
y=b/a>0,\quad c=N_f/N_g>0,\quad
r=\alpha+\beta y,\quad s=\alpha+\beta c.
\]

The physical metric is g_eff=g(alpha I+beta sqrt(g^-1 f))². Its lapse and
scale are Ne=Ng*s and Ae=a*r. Physical cosmic time satisfies dT=Ne*dt.
G and F have mass dimension two, m4, rho, p, mu and Xi have dimension four.
Here beta without a subscript is a positive matter-coupling constant, not
one of the five HR potential coefficients beta_n. The matter variables are
rho=(D_eff chi)²/2+V, p=(D_eff chi)²/2-V, n=rho+p.

These formulas use the full background equations, not independently
conserved fictitious induced g/f sources. A canonical scalar has no
independent linear transverse vector or TT matter fluctuation. Its
homogeneous stress still contributes to both quadratic metric sectors.

## 2. Literal TT action including that stress

For D=diag(1,-1,0), tr(D²)=2, take spatial metrics
a² exp(hD), b² exp(gamma D), and set d=gamma-h. Their determinants are
exactly independent of h,gamma. The four square-root eigenvalues are
c, y exp(d/2), y exp(-d/2), y. Therefore their elementary generating
polynomial is

\[
(1+cz)(1+yz)[1+2y\cosh(d/2)z+y^2z^2].
\]

The interaction density is minus Ng*a³*m4 times the corresponding beta_n
sum. The effective spatial volume divided by a³ is

\[
r[\alpha^2+\beta^2y^2+2\alpha\beta y\cosh(d/2)].
\]

The effective lapse remains Ng*s for this ansatz. Holding chi and its
coordinate-time derivative fixed, the matter density is that volume times
Ng*s*p. Expanding both terms gives L2=-Ng*a³*mu*d²/4, where

\[
\boxed{\mu=y\{m4[\beta_1+\beta_2(y+c)+\beta_3cy]
                    -\alpha\beta r s p\}.}
\]

The literal ADM extrinsic-curvature combination for g is
tr(K²)-tr(K)²=-6Hg²+hdot²/(2Ng²), and similarly for f. For the positive
spatial metric diag(exp(h(z)),exp(-h(z)),1), a direct Christoffel contraction
gives R3=-(partial_z h)²/2. These yield the kinetic and gradient terms,
respectively. Dividing the norm-two polarization result by two gives, per
TT tensor contraction,

\[
S_T^{(2)}={1\over8}\int dt\,d^3k\left[
{Ga^3\over N_g}|\dot h|^2-GN_ga k^2|h|^2
+{Fb^3\over N_f}|\dot\gamma|^2-FN_fb k^2|\gamma|^2
-N_ga^3\mu|h-\gamma|^2\right].
\]

The Einstein coefficients and principal gradient terms are positive in the
stipulated domain. Mu is a relative algebraic stiffness. No rolling mass
eigenvalue, canonical gap or instability timescale is identified with mu/G.

## 3. Finite transverse shift from a literal square root

Use a common coordinate transformation to set the g coordinate shift to
zero, and write the f relative coordinate shift along x as sigma. This
computes the coefficient of the gauge-invariant relative shift. Set
w=b² sigma²/Ng². After a diagonal basis conjugation, the t-x block of g^-1 f
has the representative

\[
M=\begin{pmatrix}c^2-v^2&-vy\\vy&y^2\end{pmatrix},\quad v^2=w,
\quad\det M=c^2y^2.
\]

Its square root on the branch continuous from the positive FLRW root is

\[
S_2={M+cyI\over L},\quad L=\sqrt{(c+y)^2-w}.
\]

The identity (M+cyI)²=[(c+y)²-w]M is checked coefficientwise. For small
shifts, L is positive and the root is analytic; at zero shift,
L=c+y and dL/dw=-1/[2(c+y)]. There is no assumption about the sign of mu.
The two other square-root eigenvalues equal y, giving the generating
polynomial (1+yz)²(1+Lz+cy z²).

Define D_A=alpha²+alpha*beta*L+beta²*c*y and
W=alpha²+beta²*y²+2alpha*beta*y(c+y)/L. Direct determinant and inverse-metric
algebra gives

\[
\sqrt{|g_{\rm eff}|}=N_ga^3r^2D_A,\qquad
g_{\rm eff}^{00}={W\over N_g^2D_A^2}.
\]

In taking the shift variation, chi_dot is held fixed at its background
value chi_dot²=Ng²*s²*n, rather than holding the perturbed matter density
fixed. The matter density divided by Ng*a³ is consequently

\[
r^2\left[{s^2n W\over2D_A}-{(\rho-p)D_A\over2}\right].
\]

Its derivative, added to the interaction derivative, gives the coefficient
Z_v/[2(c+y)] of w, with

\[
\boxed{Z_v=P+\alpha\beta r^2\rho+
 {\alpha\beta r s y\over c+y}(\rho+p),\quad
\Xi={2y^2Z_v\over c+y},\quad
P=m4(\beta_1+2\beta_2y+\beta_3y^2).}
\]

This is the definition used at mu=0. Away from such a zero it agrees with
the source shorthand Xi=mu/cV². In particular P>0, rho>0 and n>0 suffice
for Xi>0, but no such sign theorem is asserted for arbitrary sources or
HR coefficients.

## 4. Complete vector constraint elimination, including mu=0

Use spatial vectors h_ij=i(k_i E_j+k_j E_i)/2 and the corresponding S for
f, with each transverse to k. Let B=a*sigma_g/Ng and b_v=b*sigma_f/Nf be
the source shift variables. Define u_s=2Ng*B/a and v_s=2Nf*b_v/b. They are
twice the coordinate shifts, not the physical lapse ratio c. At any k>0,
the action for one transverse polarization is

\[
{Ga^3\over8N_g}\left[
A_0(\dot E-u_s)^2+B_0(\dot S-v_s)^2+C_0(u_s-v_s)^2
-{N_g^2\mu\over2G}k^2(E-S)^2\right],
\]

\[
A_0={k^2\over2},\quad B_0={Fy^3k^2\over2Gc},\quad
C_0={a^2\Xi\over2G}.
\]

The shift term follows literally from section 3. The EH shear kinetic term
has the factor k²/2 in this spatial-vector convention. The relative
spatial-potential term follows from the TT relative quadratic contraction:
|(i/2)(k_i(E-S)_j+k_j(E-S)_i)|²=k²|E-S|²/2. The common-vector gauge terms
cancel on the background equations. This is also precisely the translated
unreduced source equation 4.8.

For general positive A0,B0,C0 the shift matrix is
[[A0+C0,-C0],[-C0,B0+C0]], with determinant
D0=A0*B0+A0*C0+B0*C0>0. Completing its full positive square leaves

\[
{A_0B_0C_0\over D_0}(\dot E-\dot S)^2
=(A_0^{-1}+B_0^{-1}+C_0^{-1})^{-1}(\dot E-\dot S)^2.
\]

This identity retains both shifts, is not a large-k approximation, and
does not divide by mu. Set q=k(E-S), holding its nonzero comoving k fixed.
The exact reduced quadratic action is

\[
\boxed{L_V=K\dot q^2-Uq^2,\quad
K={Ga^3\mathcal C\over16N_g},\quad U={N_ga^3\mu\over16},\quad
\mathcal C=\left[1+{Gc\over Fy^3}+{Gk^2\over a^2\Xi}\right]^{-1}.}
\]

It is regular with K>0 even at mu=0 as long as Xi>0. The homogeneous k=0
mode is outside this q chart and outside the vector-principal conclusion.

## 5. Physical clock, principal speed and the gap boundary

Changing first to physical time gives K_T=Ne*K=G*a³*s*C/16 and
U_T=U/Ne=a³*mu/(16s). The exact physical-time equation is

\[
q_{TT}+{(K_T)_T\over K_T}q_T+{U_T\over K_T}q=0,
\quad {U_T\over K_T}={\mu\over Gs^2}
\left[1+{Gc\over Fy^3}+{Gk^2\over a^2\Xi}\right].
\]

Comparing its large-k coefficient with k_phys²=k²/(a*r)² gives

\[
\boxed{c_{V,{\rm eff}}^2=\left({r\over s}\right)^2{\mu\over\Xi}.}
\]

The first-derivative coefficient has a finite large-k limit on a smooth
Xi>0 domain and cannot alter this leading k² coefficient. No EFT-validity
claim is attached to sending k to infinity. This is a formal principal
coefficient of the stated classical quadratic system.

In either clock, for an action K*qdot²-U*q², define Q=sqrt(2K)*q and
f=(sqrt K)dot/sqrt K. Removing only the explicit total derivative gives

\[
L={1\over2}[\dot Q^2-\Omega^2 Q^2]
-{d\over dt}\left({fQ^2\over2}\right),\qquad
\boxed{\Omega^2={U\over K}-{(\sqrt K)^{\ddot{}}\over\sqrt K}}.
\]

Use K_T,U_T and T derivatives for a physical oscillator. The K derivatives
contain volume, lapse, ratio and momentum dependence; they are not dropped
or borrowed from the old CD/M1 oscillator. For example K_T=exp(2lambda*T),
U_T=0 gives Omega²=-lambda², not zero. This is only an identity/control,
not a background solution asserted for this theory. Thus mu=0 does not
alone prove a vanishing canonical gap or strong coupling. A rolling gap,
growth rate or finite-energy estimate would require additional analysis.

## 6. Analytic jets and the genuinely independent scale

For the S6.6 parameter family, take G=F=M², m4=M²*m², alpha=beta=1,
beta_n=(0,0,1,0,0), and set u=mT. Bars divide densities by M²*m² and
Hubbles by m. The frozen gate chose m=1/tau, but this equality must not
be imposed when assessing other parents. On its pressure branch,

\[
\bar p={2y\over(1+y)^2},\quad
{\mu\over m4}={y(y-1)(y-c)\over1+y}.
\]

Use the pinned exact reconstruction ODE with h=H_eff/m and initial
y=1, rho_bar=1/2, h=0. Let A=h'(0), allowing any analytic prescribed h.
The initial X=Hg/m=sqrt(7/3), Y=Hf/m=-sqrt(7/3), Ng=Nf=1/2 and n_bar=1
are independent of A. Differentiating the actual ODE and lapse quotient
without any Hubble division gives

\[
X'=Y'=-2,\quad y'=-\sqrt{7/3},\quad\bar\rho'=0,\quad
c'=-{5+12A\over3\sqrt{7/3}}.
\]

At the initial point both factors y-1 and y-c vanish. Consequently unknown
second derivatives of y and c drop out of the leading stiffness term:

\[
\boxed{{\mu\over m4}=(1/3-2A)u^2+O(u^3),\qquad
{\Xi(0)\over m4}=6,\qquad
c_{V,{\rm eff}}^2={1-6A\over18}u^2+O(u^3).}
\]

The free solution's actual-clock jet is A=7/36, not a chosen time parameter.
Thus mu/m4=-u²/18+O(u³) and cV,eff²=-u²/108+O(u³). The frozen reconstructed
CD solution has A=4, giving mu/m4=-23u²/3+O(u³) and
cV,eff²=-23u²/18+O(u³).

The inherited solutions and the relevant root/quotient functions are
analytic on an open interval about zero, with positive lapses, density,
null source and Xi. If f(u)=a*u²+O(u³), a<0, continuity of f(u)/u² extended
by a at zero gives a two-sided punctured neighborhood with f<0. Applying
this to either coefficient proves the stated local principal-gradient
obstruction while K remains positive. This is an existence proof of a
negative interval, not a computed lower bound on its width, its accumulated
growth, or the wavenumbers lying below an actual background-dependent EFT
cutoff.

For independently scaled CD,

\[
h(u)={4u\over(m\tau)^2+u^2},\quad A={4\over(m\tau)^2},\quad
{\mu\over m4}={(m\tau)^2-24\over3(m\tau)^2}u^2+O(u^3).
\]

Hence m*tau>sqrt(24) reverses this leading coefficient. The local analytic
reconstruction remains a different, potential-reconstructed candidate;
it is not the frozen action or a free-M1 matching solution. The equality
case is undecided by these jets. Even strict positivity of the leading
term certifies neither the scalar sector nor finite-window vector health,
canonical gaps, interaction suppression or a cutoff hierarchy.

## 7. Verification and omission boundaries

The main engine differentiates the literal tensor/source action, computes
the matrix-root identities, completes the two-shift square and
differentiates the pinned actual-time reconstruction. A separate sparse
Fraction polynomial engine checks the TT series, the matter inverse-metric
quotient rule, both shift equations and their completed square, and the
bounce jet in the exact root²=7/3 quotient. These are coefficientwise
identities, not interpolation from rational time samples.

Independent tests also retain the original source action and physical
clock. Controls detect omission of the pressure term, the shift matter
term, inverse-effective-metric variation, composite clock/ruler factors,
the unprovided extra source mass factor, and the rolling normalization
connection. Positive Xi and k>0 are explicit constraint-domain hypotheses;
a negative-Xi fixture is deliberately outside them. Certificate replay
compares every source/doc/test hash and all exact outputs with the stored
report, recursively replaying the pinned prior lineage without writes.

This gate does not derive the coupled scalar quadratic action, its regular
constraints or gradient matrix. It does not compute the composite
coupling's BD-mode cutoff on these rolling backgrounds. Those are research
obligations, not conclusions to infer from the vector signs or from the
positivity of a source paper's scalar kinetic matrix alone.
