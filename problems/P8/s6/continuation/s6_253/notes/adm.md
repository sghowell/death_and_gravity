# Entire homogeneous-background scalar ADM block

Use additive lapse perturbation n, logarithmic spatial perturbation v,
M1 perturbation sigma and contravariant shift divergence b. Background
N,a_hat,H=dot(a_hat)/a_hat and m=dot(chi) are independent coordinate-time
variables; H and m are not divided by N. Put q=P^2/a_hat^2. The exact
CD+M1 ADM density per kappa*N*sqrt(hhat) is

    M/2*(Kij Kij-K^2)+B*K+Fhat+C3*R3
      +U/(2N^2)*(chi_dot-N^i partial_i chi)^2
      -Cchi/2*hhat^ij partial_i chi partial_j chi.

For hhat_ij=a_hat^2 exp(2v) delta_ij,
R3=a_hat^-2 exp(-2v)[-4 Delta v-2(partial v)^2]. The expansion retains
shift transport in both K and the matter velocity. Its integrations by
parts use integral(N^i partial_i v)=-integral(v b), and the analogous
matter identity. The scalar shift shear-square difference is a spatial
boundary at quadratic order on a homogeneous background. No individual
transport term is silently set to zero.

Define background functions and lapse derivatives before restriction:

    D=M/N, Z=U/N, C=N*C3, Y=N*Cchi, c=Z*m,
    L0=-3D H^2+3B H+N Fhat+Z m^2/2,
    Theta=-H D_N+B_N/2, w=m Z_N,
    Cnn=(-3H^2 D_NN+3H B_NN+(N Fhat)_NN+Z_NN m^2/2)/2,
    Qn=-3H^2 D_N+3H B_N+(N Fhat)_N+Z_N m^2/2,
    J=Cnn+3Theta^2/D-w^2/(2Z),
    cv=-18D H+9B,
    Vvv=9L0/2-(d_t+3H)cv/2,
    Vvs=-3(d_t+3H)c,
    Lnv=3Qn+4(N C3)_N q.

The WHOLE normalized quadratic base Lagrangian is

    -3D vdot^2+Z sigmadot^2/2+6Theta n vdot+w n sigmadot
    +(J-3Theta^2/D+w^2/(2Z))n^2
    -3c vdot sigma+Lnv n v+Vvv v^2+Vvs v sigma
    +b(2D vdot-2Theta n+c sigma)+2C q v^2-Y q sigma^2/2.

The raw expression instead contains cv*v*vdot+3c*v*sigmadot and9L0*v^2/2.
The difference is precisely the total time derivative

    d_t[kappa*a_hat^3*(cv*v^2/2+3c*v*sigma)].

This boundary must be transported with the canonical state and the other
inherited chart boundaries. It is not a license to reset endpoint data.

At the actual current reference, D=Z=1, m=ell,
w=-ell E, C=1/2, Y=1, Qn=Tcorr, Vvv=9A/2, Vvs=0 and J=Jc.
parent.reference_data verifies every coefficient, including the full
fixed-profile lapse pivot. In particular Vvs is zero only on that
reference; its physical derivatives are not zero. The independent raw
ADM expansion checks the entire expression with independent family jets.
Separate direct parent-path tests differentiate all18 resulting coefficient
and scale rows through second order at three rational time fixtures.
The independent profile functions in those diagnostics test the algebra;
they are not replacement physical stress profiles or a new reference state.

This base block alone is not the full off-reference Gaussian. The complete
source-square/longitudinal-vector block is included next. Off-shell
background terms and their time derivatives are retained, not replaced
by reference equations during variation.
