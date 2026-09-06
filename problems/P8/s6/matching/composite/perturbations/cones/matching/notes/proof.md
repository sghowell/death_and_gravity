# Exact rolling two-TT reduction and limited hierarchy screens

## 1. Domain, physical time and normalization

Adopt the fixed action and physical composite metric in
[FORMULATION.md](../FORMULATION.md). Work on a regular common-flat
positive-root interval with positive finite lapses and scales. All
coefficients needed below are sufficiently differentiable; the centre
examples are analytic. A prime in sections 1–5 denotes \(d/dT\), where
\(dT=N_e\,dt=N_gs\,dt\), and \(q=k/A_e\) for fixed comoving \(k\).
The two TT polarizations obey identical independent copies of the same
quadratic action; the following normalization is per real contraction.

The literal source-aware TT action, already independently varied and
pinned in S6.7, is
\[
 L=K_1h'^2+K_2\gamma'^2
 -q^2(K_1c_g^2h^2+K_2c_f^2\gamma^2)-U(\gamma-h)^2,
\]
\[
 K_1=\frac{Ga^3s}{8},\quad
 K_2=\frac{Fa^3y^3s}{8c},\quad
 U=\frac{a^3\mu}{8s},\quad
 c_g^2=(r/s)^2,\quad c_f^2=(cr/(ys))^2,
\]
\[
 \mu=y\{m_4[\beta_1+\beta_2(y+c)+\beta_3cy]-\alpha\beta rs\,p\}.
\]
The pressure term is retained. The transformation to physical time
multiplies the coordinate-time kinetic coefficients by \(N_e\) and
divides the stiffness by \(N_e\). No assumption about a late-time
proportional background is used.

Define
\[
 K_\Sigma=K_1+K_2,\quad K_R=\frac{K_1K_2}{K_\Sigma},\quad
 w_i=\frac{K_i}{K_\Sigma},\quad
 f_\Sigma=\sqrt{2K_\Sigma},\quad f_R=\sqrt{2K_R}.
\]
All of these kinetic factors are positive. The weighted physical field
and relative field, followed by their canonical normalization, are
\[
 u=w_1h+w_2\gamma,\quad\delta=\gamma-h,\quad
 l=f_\Sigma u,\quad H=f_R\delta .
\]
The inverse map is
\[
 h=\frac l{f_\Sigma}-w_2\frac H{f_R},\qquad
 \gamma=\frac l{f_\Sigma}+w_1\frac H{f_R}.
\]
In particular, this inverse map, and not just a locked field, is needed
when comparing physical tensor data.

Let
\[
 \theta_\Sigma=f_\Sigma'/f_\Sigma,\quad \theta_R=f_R'/f_R,\quad
 N_\Sigma=f_\Sigma''/f_\Sigma,\quad N_R=f_R''/f_R,
\]
\[
 \omega=\sqrt{w_1w_2}\,\partial_T\log\sqrt{K_2/K_1}
 =\frac{f_\Sigma w_2'}{2f_R}.
\]
This follows from
\(w_2'=w_1w_2\,\partial_T\log(K_2/K_1)\). Differentiating the
weighted field rather than freezing its weights gives the exact identity
\[
 K_1h'^2+K_2\gamma'^2
 =K_\Sigma(u'-w_2'\delta)^2+K_R\delta'^2.
\]
It is checked by two separate symbolic/coefficientwise implementations.

The rotated spatial and algebraic data are
\[
 c_L^2=w_1c_g^2+w_2c_f^2
 =\left(\frac r s\right)^2\frac{G+Fcy}{G+Fy^3/c},
\]
\[
 c_H^2=w_2c_g^2+w_1c_f^2,\qquad
 D=\sqrt{w_1w_2}(c_f^2-c_g^2),
\]
\[
 m_{\rm alg}^2=\frac{U}{K_R}
 =\frac{\mu}{s^2}\left(\frac1G+\frac{c}{Fy^3}\right).
\]
The last expression is the relative algebraic stiffness divided by its
kinetic coefficient. It is not by definition a rolling physical gap.
The exact identities below allow \(m_{\rm alg}^2\) to be any real value.

## 2. Full canonical action, boundary and equations

Substitution gives, before the self-normalization integrations by parts,
\[
 L_{\rm pre}=\tfrac12(l'-\theta_\Sigma l-2\omega H)^2
 +\tfrac12(H'-\theta_RH)^2
 -\tfrac12[q^2c_L^2l^2+2q^2D\,lH+
 (m_{\rm alg}^2+q^2c_H^2)H^2].
\]
Define
\[
 V_{LL}=q^2c_L^2-N_\Sigma,\qquad
 V_{HH}=m_{\rm alg}^2+q^2c_H^2-N_R-4\omega^2.
\]
The exact boundary relation is
\[
 L_{\rm pre}=L_{\rm post}
 -\partial_T\frac{\theta_\Sigma l^2+\theta_RH^2}{2},
\]
\[
 L_{\rm post}=\tfrac12(l'^2+H'^2-V_{LL}l^2-V_{HH}H^2)
 -2\omega H l'-(q^2D-2\omega\theta_\Sigma)lH .
\]
Neither \(\omega'\) nor the self-normalization curvatures may be
discarded. With
\[
 B=2\omega(\partial_T-\theta_\Sigma)+q^2D,\qquad
 B^*=-2\omega\partial_T-2\omega'-2\omega\theta_\Sigma+q^2D,
\]
the literal surface identity is
\[
 fBg-(B^*f)g=\partial_T(2\omega fg).
\]
The two Euler equations are therefore
\[
 D_Ll+B^*H=0,\qquad D_HH+Bl=0,\qquad
 D_L=\partial_T^2+V_{LL},\quad D_H=\partial_T^2+V_{HH}.
\]
The post-boundary momenta and Hamiltonian are
\[
 p_l=l'-2\omega H,\qquad p_H=H',
\]
\[
 {\cal H}=\tfrac12[p_l^2+p_H^2+V_{LL}l^2+
 (m_{\rm alg}^2+q^2c_H^2-N_R)H^2]
 +2\omega p_lH+(q^2D-2\omega\theta_\Sigma)lH .
\]
A canonical diagonal potential is not alone the spectrum of this coupled
time-dependent system. In particular \(V_{HH}<0\) is not used below as a
proof of physical growth or the absence of every possible heavy response.

## 3. Causal equation Schur formula and controlled-approximation inputs

Fix an initial time and let \(G_R\) be the retarded inverse of \(D_H\)
with zero initial response. For the actual heavy initial data,
\[
 H=H_{\rm hom}-G_RBl,\qquad
 (D_L-B^*G_RB)l+B^*H_{\rm hom}=0.
\]
This is an exact equation identity. It follows by solving the second
Euler equation and substituting in the first. The homogeneous term cannot
be omitted unless the corresponding physical/canonical initial data are
specified and matched. A causal inverse is not symmetric under
interchange of its time arguments. Substituting it into an ordinary
single-copy quadratic action would symmetrize its kernel upon variation;
that is not the above initial-value equation.

On an interval where \(m_{\rm alg}^2\ne0\), consider the formal first
inverse-mass representative
\[
 H_0=-m_{\rm alg}^{-2}Bl,\qquad
 D_Ll-B^*(m_{\rm alg}^{-2}Bl)=0 .
\]
The light equation residual of \((l,H_0)\) is exactly zero and the heavy
equation residual is exactly
\[
 R_H=-(\partial_T^2+V_{HH}-m_{\rm alg}^2)
       [m_{\rm alg}^{-2}Bl].
\]
Thus derivatives of the mass, moving eigenvector, normalization, physical
wave number and rotated gradient coefficients remain present. The
corresponding formal local first-order correction to the light action is
\[
 \Delta L_{\rm local}=\frac{(Bl)^2}{2m_{\rm alg}^2}.
\]
Its Euler contribution, with the sign convention above, is
\(-B^*(m_{\rm alg}^{-2}Bl)\). This local expansion is not an exact
retarded effective action.

The zeroth locked action is
\(K_\Sigma u'^2-q^2K_\Sigma c_L^2u^2\), or its canonically normalized
version with \(N_\Sigma\). Its leading locked speed is \(c_L\), not the
larger of the full two-TT principal speeds. Interpreting this expression
as a controlled light EFT needs quantitative information, for example:

- A uniform positive algebraic mass \(m_{\rm alg}^2\geq M_*^2\), if the
  proposed construction is specifically a mass-led derivative expansion.
- A chosen finite momentum band and small \(q_{\max}/M_*\), bounded
  gradient weights, small mixing/normalization rates relative to \(M_*\),
  and appropriate derivatives of those rates, the mass and \(q\).
- Small physical background curvature relative to \(M_*^2\) when claiming
  a covariant slow-background expansion.
- An actual heavy or full coupled propagator bound on the chosen window,
  a light solution/derivative class controlling \(R_H\), and matched or
  suppressed heavy homogeneous data.

These are construction inputs, not automatic consequences of a positive
number \(m_{\rm alg}^2\) at one time. For illustration only, a consistent
counting \(\partial_Tl=O(\epsilon M_*l)\),
\(\omega,\theta=O(\epsilon M_*)\), \(q=O(\epsilon M_*)\), with analogous
coefficient derivative bounds, gives \(Bl=O(\epsilon^2M_*^2l)\),
\(H_0=O(\epsilon^2l)\), and \(R_H=O(\epsilon^4M_*^2l)\).
No constants making this counting uniform are established here.

More precisely, if the full canonical propagator obeys a supplied norm
bound \(C_I\) on the time window, the exact variation-of-constants formula
for the full solution and the phase representative of \((l,H_0)\) gives
\[
 \|z-z_0\|\leq C_I\left(\|z_i-z_{0i}\|+
       \int_I |R_H|\,dT\right)
\]
in a compatible phase/source norm. Additional inverse-map bounds then
control the physical \(h,\gamma\). This conditional identity explains
what a residual proof would need; no advantageous \(C_I\) or small
physical-data error is certified in this gate.

## 4. Exact zero-momentum charge and the necessary countercontrol

At \(k=0\) the physical common field is cyclic:
\[
 L=K_\Sigma(u'-w_2'\delta)^2+K_R\delta'^2-U\delta^2,\qquad
 J=2K_\Sigma(u'-w_2'\delta)
  =f_\Sigma(p_l-\theta_\Sigma l),\qquad J'=0.
\]
At fixed \(J\), varying the relative field, then substituting \(H=f_R\delta\),
gives
\[
 H''+(m_{\rm alg}^2-N_R)H=-\frac{w_2'}{f_R}J.
\]
Thus the fixed-charge \(k=0\) Routh frequency is
\(\Omega_R^2=m_{\rm alg}^2-N_R\). It can be positive when \(\mu=0\),
purely from the rolling normalization. This exact charge reduction is a
useful control against both zero-stiffness and negative-diagonal
overclaims. It is not the finite-\(k\) coupled spectrum or a positive
retarded Green bound.

The same caution applies to a large mixing coefficient: at \(k=0\),
\(Bl=2\omega f_\Sigma u'\). It vanishes for a constant common physical
field. Large canonical coefficients do not establish heavy excitation for
every choice of physical light state.

## 5. Dimensions

With dimensionless scale factors, \(G,F\) have mass dimension two and
\(\mu,U\) have dimension four. \(K_i,K_\Sigma,K_R\) have dimension two;
\(f_\Sigma,f_R\) have dimension one. \(q,\theta_i,\omega\) have
dimension one; \(N_i,m_{\rm alg}^2,\Omega_R^2,\dot H_e\) have dimension
two. The ratios compared below are dimensionless in the same physical
clock. Fourier-volume conventions do not change these ratios.

## 6. Full local CD-centre reconstruction

For this section, primes mean \(d/du\), with \(u=mT\). Fix the beta2 model
and exact states in the formulation. Bars on pressure and density are
suppressed below: they are divided by \(M^2m^2\).
On the pressure branch, the pinned reconstruction functions are
\[
 p=\frac{2y}{(1+y)^2},\quad
 X=\sqrt{y^2+(1+y)^3\rho/3},\quad
 Y=-\sqrt{y^{-2}+(1+y)^3\rho/(3y^3)},
\]
\[
 y'=\frac{y(1+y)[XY-h(X+Y)]}{X-yY},\qquad
 \rho'=-3h(\rho+p),
\]
\[
 N_g=\frac{(1+y)h-yY}{X-yY},\qquad
 N_f=\frac{X-(1+y)h}{X-yY},\qquad c=N_f/N_g.
\]
The lapse gauge is \(N_g+N_f=1\), so \(T\) is already physical time up to
the explicit \(u=mT\) scaling. At either specified initial state,
\(X>0>Y\), \(X-yY>0\), \(N_g,N_f>0\), and \(\rho+p>0\).
The radicands and denominators remain positive in some neighbourhood.
The analytic ODE theorem supplies a local solution for every finite
\(A\) with the prescribed analytic CD \(h\). A canonical scalar and local
potential follow from \(\chi'=M\sqrt{\rho+p}\) and
\(V=M^2m^2(\rho-p)/2\), using the locally invertible scalar clock.
The potential generally depends on the prescribed history; no one
globally free scalar action for all \(A\) is asserted.

To obtain the normalization jets, let \(z=(y,\rho,h)\). At the centre
\[
 z'=(y',0,A),\quad
 z''=(y'',-3A(\tfrac12+p),0),\quad
 y''=(\partial_z y')z'.
\]
For example,
\[
 c'=(\partial_z c)z',\qquad
 c''=(\partial_z c)z''+
       \sum_{ij}(\partial_i\partial_jc)z_i'z_j'.
\]
The \(c''\), \(\rho''\), \(y''\), and prescribed \(h''=0\) terms are
retained in the implementation and independent audit.
Up to constants which disappear under logarithmic differentiation,
\[
 f_\Sigma^2\propto
 A_e^3\frac{(1+c)(c+y^3)}{c(1+y)^3},\qquad
 f_R^2\propto
 A_e^3\frac{(1+c)y^3}{(1+y)^3(c+y^3)}.
\]
The common volume contributes \(3A/2\) to each \(\theta_i'\) at the
centre. All reported first-order rates are divided by \(m\), and
second-order rates/masses by \(m^2\).

## 7. Symmetric state: zero algebraic mass, but nonzero Routh frequency

For \(y_0=1,\rho_0=1/2\), set \(\sigma=\sqrt{7/3}\). Direct ODE
differentiation gives
\[
 y'=-\sigma,\quad c'=-\frac{5+12A}{3\sigma},\quad
 y''=y'^2,\quad c''=c'^2,\quad
 \theta_\Sigma=\theta_R=\omega'=0,
\]
\[
 \omega=\frac{3A-4}{3\sigma},\quad m_{\rm alg}^2=0,\quad
 c_L^2=c_H^2=1,\quad D=0,
\]
\[
 N_\Sigma=N=\frac{144A^2-6A+67}{84},\qquad
 N_R=N-4\omega^2=\frac{9A}{2}-\frac94.
\]
These equations are in the dimensionless units of section 6. Both
canonical Euler diagonals at the centre are \(q^2-N\). However,
\[
 \Omega_R^2=\frac94-\frac{9A}{2},\qquad
 N=\frac{12}{7}(A-\tfrac1{48})^2+\frac{51}{64}>0,
\]
\[
 3N-\Omega_R^2=\frac{36A^2+30A+1}{7}>0\quad(A\geq0).
\]
For the independently scaled candidate \(m\tau=5\), \(A=4/25\),
\[
 N=\frac{43579}{52500},\qquad
 \Omega_R^2=\frac{153}{100}>0 .
\]
Therefore \(\mu=0\) does not eliminate all canonical frequencies. The
positive Routh squared frequency is nevertheless never parametrically larger
than the locked inertia curvature \(N\) in this split. The rigorous
mass-led obstruction is the zero algebraic mass; no full finite-\(k\)
gap or dynamical instability inference is drawn from the negative
diagonal or the Routh comparison.

## 8. Asymmetric state: exact all-\(A\) hierarchy bounds

For \(y_0=2,\rho_0=1/2\), put \(d=\sqrt{442}\). Direct reconstruction
and canonical normalization give
\[
 c=d/13,\qquad
 \frac{m_{\rm alg}^2}{m^2}=\frac{107-5d}{42},\qquad
 c_L^2=\frac{46631-1928d}{6517},
\]
\[
 \frac{\omega^2}{m^2}
 =\frac{[(376+16d)A-1037]^2}{16(60996+3047d)}.
\]
The exact rational interval
\[
 \frac{10511}{500}<d<\frac{106}{5}
\]
is proved by squaring positive endpoints. Applying this interval
coefficientwise to affine expressions in \(d\) proves
\[
 0<\frac{m_{\rm alg}^2}{m^2}<\frac9{200},\qquad 0<c_L^2<1.
\]
Thus this state passes a leading locked-cone test. That fact alone
does not identify a controlled low-energy theory.

If \(0\leq A\leq1/4\), the looser \(d<22\) yields
\[
 1037-(376+16d)A>855,\qquad
 16(60996+3047d)<2048480,
\]
\[
 \omega^2/m^2>
 \frac{731025}{2048480}>\frac13 .
\]
If \(A\geq1/4\), instead \(\dot H_e/m^2=A\geq1/4\).
Combining either case with \(m_{\rm alg}^2/m^2<9/200\) proves
\[
 \frac{\max\{\dot H_e,\omega^2\}}{m_{\rm alg}^2}
 >\frac{1/4}{9/200}=\frac{50}{9}\qquad(A\geq0).
\]
Hence neither \(M\) nor changing \(m\tau\) can make both curvature and
mixing parametrically small relative to the algebraic mass at these
specified centre data.

For completeness, the full second-jet polynomials are
\[
 \frac{N_\Sigma}{m^2}=
 \frac{67835+1568d}{79781}A^2+
 \frac{-3796510+198215d}{1839656}A+
 \frac{-300037201409+14441671901d}{6490306368},
\]
\[
 \frac{\Omega_R^2}{m^2}=
 \frac{-1299+24d}{722}A^2-
 \frac{228939+795d}{70756}A+
 \frac{257325553-2446624d}{142644096}.
\]
The same exact root interval proves that all three coefficients of
\(N_\Sigma\) are positive, with constant term greater than \(1/2\).
It proves that the linear and quadratic coefficients of \(\Omega_R^2\)
are negative, with constant term less than \(3/2\). Therefore
\[
 N_\Sigma/m^2>\tfrac12,\quad
 \Omega_R^2/m^2<\tfrac32,\quad \Omega_R^2<3N_\Sigma
\]
for \(A\geq0\). All coefficient enclosures and strict margins are
recomputed in the certificate. This last comparison again concerns the
specified fixed-charge \(k=0\) split, not the finite-\(k\) spectrum.

## 9. Verification and exact limit of the result

The verifier replays the immutable parent report, checks the full action
and derivative residuals, evaluates the full ODE jets, and proves the
strict rational enclosures. Separate Fraction polynomial identities and
an independent physical-clock/action/jet audit guard normalization,
source, time-boundary and algebra mistakes. Controls show that omitting
the rotating weight, heavy homogeneous data or mass derivatives changes
the result. The positive Routh-frequency control blocks the incorrect
converse that zero stiffness implies no canonical frequency.

This is a falsifiable screen of a particular parametric construction:
the specified states cannot realize the stated mass-led slow-background
hierarchy just by varying the interaction-to-duration ratio. It is not a
theorem that every nonadiabatic light-only finite-window reduction fails,
nor that all composite parents, other initial states, other HR potentials,
or ultraviolet completions fail. S6.8's full principal-contract theorem
is not silently imposed on the light-only problem. Finding an actual
bounded finite-band heavy response with suitable state matching, or
establishing a more general obstruction to it, remains research work.
The omitted scalar/vector/cutoff and original-frame matching obligations
are not closed here.
