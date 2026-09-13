# Whole constrained Euler system and local unforced existence

Use configuration coordinates x=(log(a_hat),chi,hbar) and velocities
v=(H,mc,mh). At fixed lapse the normalized full homogeneous velocity
Hessian is K0=diag(-6D,Z,Z). Its lapse/velocity cross vector is
q=(6Theta,w_M1,w_H), and L0_NN=2Cnn. These statements are checked by
direct differentiation of the complete homogeneous action, not from a
selected lapse or heavy equation.

At the constraint data of the preceding note, Cnn>0. The implicit-
function theorem therefore solves E_N=L0_N=0 uniquely for
N=nstar(u,x,v) on a sufficiently small neighborhood of those data.
The derivative of the reduced Lagrangian with respect to x or v equals
the original derivative evaluated at nstar because its lapse derivative
is zero. Its entire velocity Hessian, apart from the positive common
kappa a^3 weight, is the Schur complement

    K_reduced=K0-q q^T/(2Cnn).

With J=Cnn+3Theta^2/D-(w.w)/(2Z), exact whole-matrix identities give

    det(K_reduced)=-6D Z^2 J/Cnn,
    K_reduced^-1=K0^-1+K0^-1 q q^T K0^-1/(2J).

The initial strict J,Cnn margins thus make this THREE-dimensional
velocity system regular. The negative homogeneous gravitational entry
does not prevent ordinary local ODE existence; it is not by itself a
ghost or high-frequency stability claim. All three accelerations are
solved simultaneously from the reduced Euler equations.

To make the local-existence argument explicit, choose a small closed
coordinate/velocity box around the admitted datum, inside the implicit-
lapse domain and the strict coefficient margins. The reduced Euler
right-hand side and its first derivatives are bounded there because
the full coefficient functions have the required differentiability and
the displayed inverse has no zero denominator. Extend to a slightly
larger such box, and choose a time length small enough that the bounded
right-hand side stays in the inner box and its Lipschitz constant times
that length is below1. The integral-map Picard iteration is a contraction.
It gives a unique local solution with the stated initial data. Repeated
differentiation gives the regularity needed by the S254 principal proof.
The constants may depend badly on the fixed heavy mass and the data;
no numerical or uniform-mass lifetime is claimed.

Reconstruct N using its implicit function. Differentiating E_N=0 gives
the same lapse derivative as the simultaneous full Euler/constraint
system. The reduced Euler equations imply the original three homogeneous
Euler equations because the omitted nstar derivatives multiply E_N=0.
The lapse constraint therefore holds throughout, not only initially.
The independent finite diagnostic integrates these two formulations
separately, checking all fields, the lapse and conserved M1 charge; it
is not a numerical simulation of the giant-mass physical theory.

This homogeneous solution lifts to a solution of the FULL LOCAL classical
field equations. In the full covariant action, isotropy and homogeneity
make spatial-vector, tracefree-metric and spatial-gradient Euler components
zero; lapse and volume equations give the remaining metric components.
Both matter Euler equations hold. The exact aligned vector and affine
complement reconstruction from the source note gives their Euler equations.
Finally, vary the full covariant action under a compactly supported
diffeomorphism. Integrating the metric/vector/connection Lie-derivative
terms by parts yields the Noether identity. Once their Euler expressions
and the two matter expressions vanish, it reduces to
E_clock partial_mu(phi)=0. Since phi=u and partial_u(phi)=1, the clock
Euler expression vanishes too. Thus gauge fixing has not discarded an
independent background equation.

This reasoning concerns the local classical action with its fixed
coefficient functions. The nonlocal quantum effective action and its
normal-vector/light/heavy onepoints are different equations. No quantum
stationarity is smuggled into the classical Noether argument.
