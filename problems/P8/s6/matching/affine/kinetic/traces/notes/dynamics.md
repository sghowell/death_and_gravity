# Complete rank-one two-trace curl dynamics

This independently authored argument covers every nonzero real constant
rank-one curl term built from the two specified projective-invariant
traces. It keeps the original physical metric, clock and free M1 scalar.
It concerns the actual rolling quadratic action of the S6.37 lift; it
does not assert a nonlinear open-tube inverse for the new kinetic theory.
The separate rank-two matrix problem is not used as a premise here.

## 1. Action, source and normalization

Use the frozen source convention (−+++), with the last connection index
the derivative index. Define

    V_mu=kappa^a_(mu a)−kappa^a_(a mu)/4,
    U_mu=g_(mu a) g^(bc) kappa^a_(bc)−kappa^a_(a mu)/4,
    T_mu=A V_mu+B U_mu,

where A and B are real constants. The only new term is −zeta F(T)²/4.
All coefficients in the following formulas are dimensionless after
factoring M² tau² from the action:

    zeta=zeta_physical/(M² tau²), q=(tau*k_physical)²>0.

The fixed parameters A,B are dimensionless. When (A,B)=(0,0), or when
zeta=0, the term is identically absent and the old auxiliary result is
unchanged. Those are positive controls, not propagating branches.

The module constructs both four-row trace maps directly in all 64
components and contracts every block of the frozen 60-component inverse.
On the entire rolling trajectory p=1/2 the resulting two-trace response is

    D_VV=D_UU=(3/8)eta, D_VU=D_UV=−(21/8)eta,
    eta=diag(1,−1,−1,−1).

Consequently D_T=gamma*eta, with

    gamma=3(A²−14AB+B²)/8.

The signature of eta here is the displayed response-matrix signature;
the spacetime source convention remains (−+++). These formulas are not
obtained by substituting an isolated Proca mass or by restricting torsion.
The generic eight-row trace map has rank eight and kills all four
projective columns. The actual all-component geometry is also checked by
the separately authored companion calculation and root interfaces.

## 2. Full rolling source before any center specialization

Let a=(1+u²)², h=(1+u²)³ and H=4u/(1+u²). Contracting the full
stationary connection and differentiating its coefficients gives

    Vstar_0=−3 n_dot/(2h)−3H n/(8h),
    Vstar_i=−5 partial_i n/(2h),
    Ustar_0=+3 n_dot/(2h)−27H n/(8h),
    Ustar_i=+5 partial_i n/(2h).

In particular the complete Ustar variation is not simply −Vstar.
Both traces have zero coefficients of the spatial-Hessian perturbation
at p=1/2, so no spatial metric velocity or scalar shift was suppressed
to obtain these identities. The lapse coefficient still receives the
variation of that coefficient times the nonzero background Hessian.
The U trace coefficient is −(2p−1)²/(2s_clock); its first p derivative
also vanishes there, unlike the V trace coefficient. The coefficient
source then supplies the displayed −27H/(8h) term. Both derivatives
p_phi,N and J3_N equal h_dot/(2h²), as derived from the complete lower
dictionary with its exact clock data, not by assuming the lower
coefficient is identically zero off the clock.

Use the exact gradient change W=T+d(alpha*n), with

    alpha=3(A−B)/(2h),
    d=3(7A+3B)H/(8h), e=(A−B)/h.

The exterior curl is unchanged, including the alpha_dot*n part of the
temporal component. The algebraic term for gamma nonzero is therefore

    [(W0+d n)²−a^−2(Wi+e partial_i n)²]/(2gamma).

All trace sources and their curls vanish on the exact background, so
the new action has zero first variation there. In the null case the
actual constrained trace has curl

    F(Tstar)_(0i)=partial_i[−e n_dot+f n],
    f=d−e_dot=3(11A−B)H/(8h).

This includes time-dependent coefficients before evaluation. It is not
a frozen-center substitution.

## 3. Nonzero-Schur constraint chart

Retain the original spatial gauge, physical clock gauge delta_phi=0,
and free matter perturbation s=delta_chi. Write the physical spatial
curvature as zeta_metric=v+delta*n, with delta=1/(2h), and b=a² psi.
The original action divided by a³ is

    L_CD=−3v_dot²+S n²+6Theta n v_dot+s_dot²/2
          +w n s_dot−3ell v_dot s
          +q[v²+2Lambda n v−s²/2
              +2Theta n b−2v_dot b−ell b s],
    S=J+w²/2−3Theta², ell=1/(10a³), w=ell(3delta−1).

This is the full previously derived metric/free-matter constraint
action, not a spectator-scalar approximation. J>0 for every finite u;
Theta vanishes only at u=0. Its all-real polynomial proof and exact
background are frozen inputs, checked by the parent interfaces.

For Wi=partial_i sigma, the added density is

    L_extra=(W0+d n)²/(2gamma)−q(sigma+e n)²/(2gamma)
              +zeta*q*(sigma_dot−W0)²/2.

At gamma nonzero and 1+gamma*zeta*q nonzero, its W0 equation gives

    W0=(gamma*zeta*q*sigma_dot−d n)/(1+gamma*zeta*q),
    L_extra=C(sigma_dot+d n)²−q(sigma+e n)²/(2gamma),
    C=zeta*q/[2(1+gamma*zeta*q)].

No lapse velocity remains. The scalar shift equation is unchanged:
v_dot=Theta*n−ell*s/2. On Theta nonzero it fixes n; the lapse equation
then reconstructs b. Jointly eliminating n,b,W0 from the literal
unreduced action independently reproduces the same kinetic form:

    L_kin=(1/2)[s_dot+(w/Theta)v_dot]²
           +C[sigma_dot+(d/Theta)v_dot]²
           +[J−q e²/(2gamma)]v_dot²/Theta².

Thus the three scalar pivots are C, 1/2, and
[J−q e²/(2gamma)]/Theta². The transformation completing these squares
is triangular and invertible on the stated chart. This is a direct
configuration-velocity verdict after all constraints, not an
instantaneous sign of a canonically exchanged momentum Hessian.

For zeta>0 the two nonzero-Schur cases are exhaustive:

| Schur sign | Exact sufficient momentum condition | Negative pivot |
| --- | --- | --- |
| gamma>0 | q>2gamma J/e² | clock Schur |
| gamma<0 | q>−1/(gamma zeta) | C |

In the first case e cannot vanish: A=B would give
gamma=−9A²/2, inconsistent with gamma>0. The other two pivots are
positive. In the second case the temporal denominator is strictly
negative and nonzero, C<0, while J−q e²/(2gamma)>0 and the matter
pivot is positive. Each branch therefore has exactly one negative
scalar kinetic direction above its displayed finite threshold.

For zeta<0 and gamma nonzero, the transverse trace components instead
give the direct negative velocity coefficient zeta/2. Their stationary
sources vanish in the transverse sector, their kinetic velocities are
nondegenerate, and they do not couple to a temporal vector multiplier
or to the metric-vector constraint on this rolling background. This
argument is deliberately not used for gamma=0, where the trace itself
is constrained.

## 4. The zero-Schur branch is not a limit of that chart

For nonzero (A,B), gamma=0 means A/B=7±4sqrt(3); B is then nonzero
and A−B is nonzero. No division by gamma is permitted.

Let M be the invertible frozen quotient matrix, N the four-row T map,
and kappa_star the original stationary solution. At quadratic order
the full connection equation after adding the kinetic action is

    M(kappa−kappa_star)+N^T E_kin[T]=0.

Project it with N M^−1. Since N M^−1 N^T=0 in this branch, it forces
T=Tstar, irrespective of derivatives in E_kin. The remaining unique
gauge-fixed connection is

    kappa=kappa_star−M^−1 N^T E_kin[Tstar].

This satisfies every connection equation. The auxiliary correction
evaluated there is (1/2)E_kin^T N M^−1 N^T E_kin=0. Hence the actual
reduced quadratic action is exactly CD/M1 plus

    a³*zeta*q*[−e n_dot+f n]²/2.

There is no independent Proca homogeneous solution to append, nor a
retarded inversion. The value zeta*q/2 obtained by formally putting
gamma=0 in the preceding C formula misses the new constraint and is
not a kinetic reduction of this branch. An independent two-variable
null-Schur model checks this distinction; the actual four-trace
projection and all-component equation are supplied by the full
geometry, not inferred from that toy.

The projected response vanishes on the rolling quadratic coefficients.
It is not asserted to vanish on a nonlinear open field tube.

## 5. Complete null-branch highest-derivative Legendre map

The added null-branch term is independent of b. For any finite u
nonzero its equation again fixes

    n=(Q+ell*s/2)/Theta, Q=v_dot,
    n_dot=(v_ddot+ell*s_dot/2+ell_dot*s/2)/Theta
            −Theta_dot*(Q+ell*s/2)/Theta².

The lapse Euler equation now contains n_ddot, but its coefficient of b
is still exactly 2a³qTheta. It reconstructs b as a function of the
remaining jets rather than imposing another constraint on them. The
clock equation is supplied by the diffeomorphism identity after the
metric, connection and free-matter equations, since the actual clock
derivative is nonzero. Neither metric auxiliary equation has been
discarded by substituting the shift constraint.

Set

    a_e=−e/Theta, b_e=a_e*ell/2,
    r=−e*ell_dot*s/(2Theta)
       +[e*Theta_dot/Theta²+f/Theta](Q+ell*s/2).

The complete reduced Lagrangian can be written

    L=a³[L0(v,Q,s)+L1(v,Q,s)s_dot+s_dot²/2
          +zeta*q*(a_e*v_ddot+b_e*s_dot+r)²/2],
    L1=w*n.

L0 is the remaining full CD scalar density, not just a potential on
the background. The joint highest-derivative Hessian in (v_ddot,s_dot)
is

    a³ [[zeta*q*e²/Theta², zeta*q*e²*ell/(2Theta²)],
        [zeta*q*e²*ell/(2Theta²), 1+zeta*q*e²*ell²/(4Theta²)]],

with determinant a^6*zeta*q*e²/Theta². It is nonzero for either allowed
sign of zeta. The free canonical matter term is essential to this
joint rank calculation; dropping it makes the highest block rank one
and is an explicit failed control. No frozen-matter approximation is
used to infer nondegeneracy.

Let P1=partial L/partial v_ddot and Ps=partial L/partial s_dot, and
introduce the independent Ostrogradsky momentum
P0=partial L/partial Q−d_u P1. The exact inverse is

    s_dot=Ps/a³−ell*P1/(2a³)−L1,
    v_ddot=[P1/(a³*zeta*q*a_e)−b_e*s_dot−r]/a_e.

Therefore, in the physical constrained phase variables
(v,Q,s;P0,P1,Ps), the Hamiltonian is

    H=P0*Q+[Ps−ell*P1/2−a³L1]²/(2a³)
        +P1²/(2a³*zeta*q*a_e²)−(r/a_e)P1−a³L0.

The invertible highest-derivative map leaves no primary relation on
P0. Both metric auxiliary equations have already been reconstructed.
For fixed Q nonzero the independent P0 term makes H unbounded in both
directions, for either sign of zeta. This is the nondegenerate
higher-derivative/Ostrogradsky verdict in the actual coupled system;
it does not follow from the sign of an isolated mass or from a
frozen-center Cauchy polynomial. Arbitrarily close punctured times
suffice, and no Theta=0 substitution is made.

## 6. Exhaustiveness, controls and physical scope

Every nonzero pair (A,B) has gamma positive, negative, or zero. The
two signs of a nonzero real zeta are covered above in each case. Thus
every nontrivial rank-one term in this constant two-trace curl family
fails the stated positive-kinetic/no-Ostrogradsky parent requirement
on the original rolling solution. A=B=0 or zeta=0 instead returns the
unchanged auxiliary theory and changes the propagating rank. The
exceptional temporal divisor 1+gamma*zeta*q=0 is excluded rather than
crossed; the strict momentum conditions avoid it.

Rescaling (A,B) to t(A,B), with nonzero real t, and zeta to zeta/t²
leaves the added action unchanged. Gamma scales by t², e and d by t,
and both q e²/(2gamma) and gamma*zeta are invariant. In the null case
zeta*e² is invariant. The exact checks and proportional-parameter
controls retain these normalizations.

For M²,tau>0 and physical zeta, the isolated nonnull mass control would
be M²/(gamma*zeta_physical), not a mass inferred from an auxiliary
Hessian alone. For A=1,B=0,M²=3,tau=2,zeta_physical=5 it equals 8/5,
with normalized coupling 5/12 and normalized mass squared 32/5. At
gamma=0 no such mass is assigned. These values are normalization
controls, not proofs that an unhealthy frequency lies below a cutoff.

The theorem concerns the literal finite-derivative action and its
actual constrained quadratic modes. It does not establish an EFT
frequency range, prove an all-orders instability, exclude other
kinetic contractions or fields, or provide a nonlinear open-tube
inverse. No V/G/B or original P8 closure follows from this family
screen. The target physical metric/free-M1 action and every frozen
ancestor remain unchanged.
