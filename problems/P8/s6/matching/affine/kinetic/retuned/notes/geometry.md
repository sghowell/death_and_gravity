# Source-centered matching, all-64 rank and the primary constraint

## 1. Literal new action, not a reinterpretation of S6.39

Let kappa_star[g,phi] be the complete frozen S6.37 stationary distortion
and T=V+U. The specified source Tstar is T[kappa_star], a local covariant
function, not an external rolling profile. Add

    Delta L_mass = -mu*g^(ab)(T-Tstar)_a(T-Tstar)_b/2, mu=11/9,

after extracting the original M². The full physical volume and original
free chi are retained. Both the linear term and the Tstar² counterterm
are part of this action. Dropping either changes the matching problem.

At zero curl coefficient, the added action and its first variation
vanish identically at kappa=kappa_star. If the new connection Hessian
is invertible on the same projective quotient, its exact elimination
therefore returns the entire original CD/M1 action, including its
off-trajectory scalar coefficients and boundaries. This is a new mass
matrix, outside the unchanged-mass curl family excluded by S6.39.

The primary-source review of
[Ikeda, arXiv:2311.11104, sections IV and V](https://arxiv.org/pdf/2311.11104)
is a methodological warning, not a construction imported here: projective
symmetry alone does not establish vector-tensor degeneracy. Its
vector-tensor families and kinematic restrictions are not substituted
for the original scalar, metric, matter and unrestricted connection.

## 2. Exact mass update and uniform rank

Let E be the original 64-by-60 trace-gauge embedding, M the quotient
Hessian, and N the four-row T map. Direct contraction of every original
block gives, before specializing to the clock,

    D=N M^-1 N^T=diag(Dt,Ds,Ds,Ds),
    Dt=3(2p³-1)/p, Ds=(8p+5)/(8p²).

With eta=diag(1,-1,-1,-1), the literal new Hessian is
M_new=M+mu*N^T*eta*N. The exact response is

    D_new=diag(gamma_t,-gamma_s,-gamma_s,-gamma_s),
    gamma_t=(1/Dt+mu)^-1, gamma_s=(mu-1/Ds)^-1.

All 64 stationary Euler equations and all 64 lifted source equations
are checked, not only their four trace contractions. The determinant
ratio is (1+mu*Dt)(1-mu*Ds)^3. At p=1/2 the old response is
-(9/2)*eta, whereas the new response is exactly eta.

The original closed p² tube is inside 47/100<p<53/100. Monotonicity
and exact rational endpoint inequalities give -6<Dt<-3 and 4<Ds<6.
Consequently

    18/19<gamma_t<9/8, 18/19<gamma_s<36/35.

Neither the old nor new quotient loses rank on that tube. The Woodbury
inverse is M^-1-M^-1 N^T(eta/mu+D)^-1 N M^-1. In the same specified
quotient infinity norm as S6.37, ||M^-1||<25, ||N||=12, ||N^T||=2,
and the middle inverse has norm below 11/24. Thus

    ||M_new^-1|| < 25+25²*12*2*(11/24) = 6900.

This is a fixed-frame algebraic inverse bound, not a propagator or a
Lorentz-invariant operator norm. Independent dense 60-dimensional
inversions and determinant comparisons test three exact p values.

The lift W=M_new^-1 N^T and right inverse W*D_new^-1 split off a
56-dimensional complement. Its projector, cross Euler equations,
right-inverse identity and retained quadratic action D_new^-1 are
verified exactly. Nondegeneracy follows from the full quotient and
four-trace ranks. With a curl term retained, this eliminates only the
complement, leaving a dynamical T field, not an auxiliary T inverse.

## 3. Complete source, including off-clock Hessian dependence

In the timelike frame phi_a=(s,0,0,0), x=-s² and H_ab=nabla_a nabla_b phi,
the four-component exact source is (Tstar_0,0,0,0), where

    Tstar_0 = 6p*J3*s³ +3(4p²-1)*p_phi*s/(2p)
               -3(4p²-1)*p_X*s*H00/p
               -(4p²-1)*(H11+H22+H33)/s.

Its independent covariant reconstruction is

    Tstar_mu=phi_mu*(alpha+c*Box(phi)+d*vHv),
    alpha=6p*J3*s²+3(4p²-1)*p_phi/(2p),
    c=-(4p²-1)/s²,
    d=-(4p²-1)/s^4-3(4p²-1)*p_X/(p*s²).

In particular the H.v term cancels between V and U. The coefficients
of the remaining Hessians do not vanish off the clock. Their effect
must be tested against the original degeneracy, not discarded.

## 4. Primary null and the nonlinear spatial chart

Keep H00 and the full spatial trace T_H independent. Including the
actual fR time boundary, f=2p², f_X=-1/(2h), gives the complete
original trace-velocity square

    L_trace=-(-3H00*s²+8*T_H*h*p²)²/(48*h²*p²*s²).

Its null is (delta H00,delta T_H)=(1,3s²/(8hp²)). The derivative-bearing
part of Tstar_0 is precisely

    -(4p²-1)*(-3H00*s²+8*T_H*h*p²)/(8hp²s).

It annihilates that same null. The centered mass term therefore
preserves the primary lapse degeneracy exactly, not just at p=1/2.
Omitting the fR time boundary leaves another null in the old trace
block but misaligns it with this source; the combined action then
fails this degeneracy test. That stronger negative control is replayed.

Use the unchanged nonlinear chart
h_physical=(4p²)^(-1/2)*h_hat, omega=-log(4p²)/4. Since
omega_X=1/(16hp²) and omega_phi=-p_phi/(2p), the source becomes simply

    Tstar_normal=(4p²-1)*K_hat+6p*J3*s³.

Every lapse velocity cancels, including the p_phi terms. The physical
metric in the matter and Maxwell action is still h_physical; h_hat
is only a local change of variables. Coordinate vector components can
be retained, so the exterior curl introduces no metric or lapse time
derivative through a moving orthonormal coframe.

The added relative trace-velocity coefficient is
3(4p²-1)²/(8 gamma_t p²)<19/1080<1/50 on the whole tube. The six metric,
three vector and one matter velocities have a nonsingular ten-by-ten
block for positive curl: five metric shear directions, the three
vector velocities and matter are positive; the unreduced gravitational
trace remains negative. This is not a physical ghost test. The lapse,
shift and vector temporal variables still have no time derivatives.

Primary degeneracy and this velocity rank do not, alone, prove the
complete nonlinear secondary-constraint rank, physical degree count,
hyperbolicity or nonlinear health. Those questions remain separate.
