# Freeze physical Euler coefficients, not the background action

Use the physical fields f=(n,zeta) and the literal two-current
operator derived in S6.73 after setting the independent
prepared matter charge to zero. The reconstructed matter
rate is -3ell*v-w*n, not its negative.

At the bounce the exact background gives
H=Theta=0, Theta'=3, ell=1/10, w=1/20,
J_e=1199/800+4*10^-6,
delta=1/2, delta'=0 and delta''=-3.
The implementation differentiates the actual frozen rational
background to obtain these jets.

The physical operator is obtained with v=zeta-delta*n and
the full time-dependent Euler chain rule before its
coefficients are evaluated at the bounce. Applying only
then the constant T congruence gives

E0=[[749377/250000,-3/200],[-3/200,-9/100]],
E1=0, E2=diag(0,6).

In particular the Theta' and delta'' terms must not be
discarded separately. Freezing a differently chosen
time-dependent coordinate system first is another diagnostic,
not a legitimate proof that these terms never mattered.

## Local Gaussian operator

The S6.67 report contains the finite physical Euler-current
coefficient rows, not a newly chosen action or subtraction.
It is hash-pinned here and fully rebuilt by the parent chain.
Sum all adiabatic orders, evaluate the coefficients at u=0,
and apply the constant T congruence.

With the symbol m still unevaluated, the result is

C0_NN=(697897*m^4+855392*m^2+10020608)/52488,
C0_NZ=(6867*m^4+39776*m^2-84864)/324,
C0_ZN=(763*m^4+1760*m^2+8960)/36,
C0_ZZ=3*(15*m^4+160*m^2-384)/2;

C2_NN=(35917*m^2-7740)/6561,
C2_NZ=2*(1243*m^2-8676)/243,
C2_ZN=2*(1243*m^2-972)/243,
C2_ZZ=4*(5*m^2-36);

C1=C3=0,
C4=Freg=[[-7357/6561,-562/243],[-562/243,-4]].

Every row is independently replayed against these compact
polynomials. The original differential operator is weighted
self-adjoint; its point-frozen lower coefficient matrices
need not individually be symmetric. They are not symmetrized.

At m=1000, the largest absolute C0 entry is 22500239999424,
below 10^15, and the largest absolute C2 entry is 19999856,
below 10^9. The exact C4 is included inside Breg below, not
also counted in the lower-order error.

## Fixed tadpole, normalized only once

Write rho,p for the fixed selected vector background values.
The fixed scalar tadpole cancels that vector gradient, and
its physical Euler Hessian is
[[-rho-p,3rho],[3rho,-9p]].
Its constant T congruence is

[[2rho-13p/4,3rho-9p/2],[3rho-9p/2,-9p]].

These are Euler coefficients, not the physical stress
Jacobian. Each entry has absolute value at most 9(E+P),
where E,P are S6.68's physical profile bounds. Those bounds
already include the powers of L and 64*pi^2 from the
physical stress normalization. Multiplying by 64*pi^2*L^2
defines Ctad in the normalized loop bracket, exactly once.

Using pi^2<10 gives

abs(Ctad_ij) <= 5760*L^2*(E+P)
             = 84242272974349549/1296 < 10^15.

The bound is uniform in the stated constant profile box and
therefore includes the actual values at the bounce.
No derivative of a reselected stress profile is inserted.
The nonlinear metric-chart Hessian contact is not added:
we use physical Euler coefficients and a constant congruence,
and the total vector-plus-fixed-tadpole background gradient
already vanishes.
