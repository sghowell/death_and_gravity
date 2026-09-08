# Coupled scalar constraint and kinetic audit

This note concerns exactly the S6.37 auxiliary action with the one added
term −zeta_physical F(V)²/4. It is not an audit of a general connection
kinetic completion. The algebraic quotient reduction supplied by the
companion vector calculation is independently interfaced here by
contracting the full frozen 64-component stationary connection. The
original physical metric, clock and free M1 scalar are retained.

The result is a failure of positive scalar kinetic energy for this literal
finite-derivative extension at sufficiently high spatial momentum on each
punctured side of the bounce. It is not a determination that these modes
lie in an independently justified low-energy EFT band.

## 1. Units, gauges, and the original action

Set M=tau=1 initially and let u be the original clock/proper time on the
background. Write

    a=(1+u²)², h=(1+u²)³, H=4u/(1+u²),
    theta=u(4h−1)/(1+u²)^4, Lambda=1−3/(2h),
    delta=1/(2h), ell=1/(10a³), w=ell(3delta−1).

J is exactly the rational function in the original CD/M1 certificate.
Its numerator is an even polynomial with nonnegative coefficients and
strictly positive constant; its denominator is 800(1+u²)^18. Hence J>0
at every finite real u. Also theta=0 precisely at u=0, without using H as
a denominator.

Use the regular clock gauge delta_phi=0, which is legitimate because the
actual clock derivative is one. Fix the scalar spatial shear to zero.
The physical spatial curvature perturbation is zeta_metric=v+delta*n,
where N=1+n, and b=a² psi is the scalar shift. Let s=delta_chi and
q=k_comoving²/a²>0. This q is spatial momentum squared, not the auxiliary
coefficient function used in the S6.37 lower-order dictionary.

After the explicitly verified original DHOST time boundaries and the
free-matter boundary 3a³ell*v*s, the original action density divided by a³
is

    L_CD = −3 v_dot² + S n² + 6theta n v_dot
           + s_dot²/2 + w n s_dot − 3ell v_dot s
           + q[v²+2Lambda n v−s²/2
               +2theta n b−2v_dot b−ell b s],
    S=J+w²/2−3theta².

This is the full quadratic scalar/constraint action, not just its
principal part. Its derivation is frozen in
[the original action proof](../../../../../notes/s1-derived-action.md)
and the metric-derived coupled-matter module. In particular the actual
matter current a³ell=1/10, its lapse coupling w and its shift coupling
−q ell b s are not omitted. The coordinates v and the shift chart used
below are established constrained variables, not a change of matter frame.

## 2. Independent rolling vector reconstruction

In the source signature (−+++), contract the literal full solution as

    Vstar_mu = sum_a kappa_star^a_(mu a)
               −(1/4)sum_a kappa_star^a_(a mu).

All lower connection indices are independent. In the orthonormal frame
aligned with the timelike clock, the coefficient of each spatial H_ii in
Vstar_0 is −(2p−1)(2p+3)/(2s_clock). It vanishes at p=1/2, s_clock=1.
Thus arbitrary first-order spatial metric velocities and shift terms in
those H_ii do not survive. The clock-normal mixed Hessian is
H_normal,i=−partial_i n at first order. Retaining just the coordinate
H_0i and forgetting its normal projection would leave a spurious H*b.

The full coefficient variations must still multiply the nonzero
background H_ii=−H. With x=−N^−2,

    p_N=−1/(2h), (p_phi)_N=h_dot/(2h²), (J3)_N=h_dot/(2h²).

The last identity uses the exact boundary data of the lower coefficient
ODE, q_coefficient(u,−1)=q_coefficient,x(u,−1)=0, not its erroneous
extension to all x. Contracting first and then differentiating gives

    delta Vstar_0 = −3 n_dot/(2h)
                   +[15h_dot/(4h²)−6H/h]n
                 = −3 n_dot/(2h)−3Hn/(8h),
    delta Vstar_i = −5 partial_i n/(2h),

where h_dot/h=3H/2. These are derived in scalar.vector_reconstruction()
from the frozen connection solution, independently of the companion
vector module's rolling formula. The background itself has Vstar=0.

Define alpha=3/(2h) and the true exact gradient shift

    W_mu=V_mu+partial_mu(alpha*n).

The exterior field strength is unchanged. Using alpha*d n alone would
miss alpha_dot*n in W0 and would not preserve F away from the center.
The quadratic algebraic term becomes

    (4/3)[(W0+d n)²−a^−2(Wi+e partial_i n)²],
    d=21H/(8h), e=1/h.

There is no lapse velocity in this representation, and no additional
scalar shift term. This local, invertible derivative transformation is
used at the level of the full constrained action; it does not justify
discarding a constraint or claiming the raw lapse-velocity square is a
new independent degree of freedom. Since both Vstar and F vanish on the
actual background, the added term has zero first variation there.

## 3. Temporal vector and lapse/shift constraints

For the longitudinal vector set Wi=partial_i sigma. At quadratic order
the additional action divided by a³ is exactly

    L_extra=(4/3)(W0+d n)²−(4/3)q(sigma+e n)²
             +(zeta*q/2)(sigma_dot−W0)².

Here and below zeta is the normalized positive coupling defined in
section 7. For every q,zeta>0, the W0 equation is regular and gives

    W0=(3zeta*q*sigma_dot−8d n)/(8+3zeta*q),
    L_extra=C(sigma_dot+d n)²−(4/3)q(sigma+e n)²,
    C=4zeta*q/(8+3zeta*q)>0.

The primary lapse-velocity degeneracy survives at this quadratic order.
For example, before the gradient shift at the center, the apparent
kinetic addition is proportional to
(sigma_original,dot+3n_dot/2)² after W0 elimination; its two-velocity
Hessian has rank one, not two. The null direction is (2,−3).

Most importantly, the scalar shift constraint remains exactly

    v_dot=theta*n−ell*s/2.

No field equation for the free scalar or lapse has been dropped. For
theta nonzero this equation determines n; the lapse equation determines
b. Jointly solving n,b,W0 gives the same reduced action as this order of
elimination. The ordinary tests independently perform that joint solve
at exact positive- and negative-time points.

## 4. Direct physical scalar kinetic verdict on the punctured chart

For any fixed finite u nonzero, substitute
n=(v_dot+ell*s/2)/theta into the full action. In the configuration fields
(v,s,sigma), the complete velocity quadratic form is

    L_kin = (1/2)(s_dot+(w/theta)v_dot)²
            +C(sigma_dot+(d/theta)v_dot)²
            +[J−(4/3)q e²]v_dot²/theta².

This is an exact congruence by an invertible triangular velocity map.
Equivalently the vector and matter pivots are C and 1/2, and the final
Schur pivot is

    K_clock,Schur=[J−(4/3)q e²]/theta².

Therefore for every finite u nonzero, zeta>0 and

    q > 3J(u)h(u)²/4,

the reduced three-scalar velocity form has precisely two positive and
one negative eigenvalue. Equality is a kinetic degeneracy, not a strict
negative-pivot case. This conclusion includes all lapse, scalar shift,
clock and free-chi constraints of the stated action. It is not inferred
from an isolated vector mass or a fixed-metric clock truncation.

The criterion is a direct configuration-velocity test in a regular
chart at every point where it is used. In particular it fails the
positive scalar kinetic requirement for the literal finite-derivative
parent. At the bounce J=1199/800, e=1, so the threshold tends to
3597/3200. Given any fixed q at the center strictly above that value,
continuity and q(u)=k_comoving²/a(u)² give punctured intervals on both
sides with the strict negative pivot. The singular unitary chart is
never evaluated at u=0 in this proof.

## 5. Regular center Hamiltonian and the time derivatives not discarded

An independent Legendre calculation retains n before lapse elimination.
Use density-normalized momenta P=s_dot+w n and
Pi=2C(sigma_dot+d n). The shift constraint implies the original curvature
momentum is −2a³q b. Define

    J_eff=J−(4/3)q e²,
    R_eff=q(theta*b+Lambda*v)+wP/2−3ell*s*theta/2
           +dPi/2−(4/3)q e sigma,
    C0=P²/2+q ell b s+q s²/2−q v²−3ell²s²/4
        +Pi²/(4C)+(4/3)q sigma².

The Hamiltonian before eliminating n is C0−J_eff n²−2R_eff n,
so n=−R_eff/J_eff and H=C0+R_eff²/J_eff where J_eff is nonzero.
The canonical shift momentum is P_b=2a³q v. Its time-dependent swap
contributes −H*b*P_b to the canonical Hamiltonian; this boundary term
must not be omitted.

At u=0 the velocity matrix in the regular (b,s,sigma) chart is

    Kbb=−2q(8q−9)/(19q−18),
    Kbs=3q/[20(19q−18)],
    Kss=(3800q−3597)/[400(19q−18)], Ksigma,sigma=C,

with zero remaining cross entries. Its two-scalar determinant is
−q(3200q−3597)/[400(19q−18)]. The lapse pole q=3597/3200 and the
velocity-chart pole q=18/19 are distinct. For q at least 8 this block
has one negative eigenvalue, whereas the original CD/M1 shift-chart
matrix is positive at q>6. This is corroboration, not the proof in
section 4: a canonical coordinate/momentum exchange can turn a
potential sign into an instantaneous velocity sign. The original
baseline itself illustrates why small-q chart inertia is insufficient.

The center jets used in scalar.center_pencil() are

    H_dot=4, theta_dot=3, d_dot=21/2,
    e_dot=alpha_dot=J_dot=w_dot=Lambda_dot=ell_dot=a_dot=q_dot=0.

The code builds the full canonical Hamiltonian before this substitution,
including its moving boundary. With H0 its six-dimensional phase Hessian,
H1 is computed as (4 partial_H+3 partial_theta+(21/2)partial_d)H0 before
evaluation. If the Hamiltonian blocks are (coordinates first)

    F0 = [[L,A],[Cblock,−L^T]],

the actual center configuration equation is Q_ddot=Mq Q+Mv Q_dot with

    Mv=(L A+A_dot−A L^T)A^−1,
    Mq=L_dot+L²+A Cblock−Mv L.

Its identity against F1+F0² is checked exactly. Replacing F1 by zero
does not pass this check. No frozen characteristic polynomial from this
center Cauchy matrix is asserted to be a uniform physical cone or an
eigenfrequency cutoff. Those stronger claims are unnecessary for the
direct punctured kinetic failure.

## 6. Controls and exclusions from the verdict

For a fixed metric/lapse, the isolated longitudinal Proca pivot C is
positive for zeta>0. That check alone misses the negative term generated
by the coupled lapse-gradient constraint. Conversely zeta<0 already
has a wrong transverse Maxwell velocity sign, so it is not a positive
kinetic control for this construction.

At zeta=0, W0=−d n and sigma=−e n are both algebraic; substitution sets
the entire extra action to zero. The vector pair disappears, and the
original CD action returns. This is a change of constraint rank, not a
continuous three-propagating-scalar positivity argument. The same
zero-coupling test rules out concluding that an arbitrarily weak positive
kinetic promotion is the unchanged auxiliary theory.

The exact center Gamma-chart pole values are not a UV scale. In
particular small positive zeta can place a problematic eigenfrequency
outside a separately prescribed EFT domain. No such frequency, cutoff,
strong-coupling, loop, or remainder contract has been demonstrated by
this calculation. The conclusion rejects the stated action as an
everywhere-positive scalar-kinetic finite-derivative parent. It does not
exclude all controlled light-only EFT completions or change the original
CD/M1 target, whose frozen result remains intact.

## 7. Restoring physical scales

Let physical time be tau*u, phi=tau*u, chi=M*chi_bar and
V_normalized=tau*V_physical. Factoring M² tau² out of the action gives

    zeta_normalized=zeta_physical/(M² tau²),
    q_normalized=(tau*k_physical)².

All zeta and q in this note and the symbolic scalar APIs are these
normalized quantities. The physical Maxwell coefficient is
dimensionless. The isolated center Proca mass, which is not the coupled
scalar health criterion, is

    m_Proca,physical²=8M²/(3zeta_physical),
    m_Proca,normalized²=8/(3zeta_normalized).

For example M²=3,tau=2,zeta_physical=5 gives normalized zeta=5/12,
physical mass squared 8/5 and normalized mass squared 32/5. The common
positive action prefactor cannot change any of the kinetic inertias.

The public checks expose the literal contraction, all algebraic
constraints, kinetic square completion, both Hamiltonian routes and
center coefficient jets. Exact-domain guards precede cached symbolic
results and reject inexact, nonreal, zero-mode and wrong-chart inputs.
No certificate is generated by the three files owned by this audit.
