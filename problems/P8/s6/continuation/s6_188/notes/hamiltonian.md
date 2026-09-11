# Exact constrained anisotropic Proca Hamiltonian

Let E=exp gamma,Q=E^-1 and det E=1. The positive spatial
metric is a^2 E and the volume is a^3. The complete parent
on this clock history has X=R=1,trace K=3H,uHu=0 and
full affine source S=0, as in S6.187. Therefore its
canonical vector action is precisely ordinary Proca
with mass1000, not a scalar or longitudinal-only surrogate.

With coordinate A_i and canonical pi^i,

    pi=a Q(A'-grad A0),
    div pi+a^3 m^2 A0=0.

In the Hamiltonian the auxiliary terms are
-A0 div pi-a^3 m^2 A0^2/2. Solving the constraint gives
+(div pi)^2/(2a^3 m^2). Thus at each comoving internal
momentum k, with C_k A=k cross A,

    H_gamma=1/2[pi^T K_gamma pi+A^T V_gamma A],
    K_gamma=a^-1 E+kk^T/(a^3 m^2),
    V_gamma=a m^2 Q+a^-1 C_k^T E C_k.

For complex Fourier amplitudes use Hermitian quadratic
forms with the usual real-field pairing. Equivalently
each real cosine/sine oscillator sector has the displayed
three-mode Hamiltonian; there is no fourth oscillator.

To check the magnetic term independently, for any
symmetric Q use the polynomial identity

    C_k^T adj(Q) C_k=(k^T Q k)Q-(Qk)(Qk)^T.

Here adj(Q)=E since det Q=1. Direct multiplication also
gives K_gamma V_gamma=(m^2+a^-2 k^T Q k)I. This is an
instantaneous algebraic identity. Time-dependent canonical
normalization adds its generator; this identity alone
is not a dynamical frequency, gap or stability theorem.

At gamma=0, K0 has transverse eigenvalue1/a and
longitudinal eigenvalue omega^2/(a m^2), while V0 has
eigenvalues a omega^2 and a m^2 respectively. Their
canonical normalizers are exactly the old
gT^2=a and gL^2=a m^2/omega^2.

All mass, magnetic, electric and constrained terms are
positive. For ||gamma||op<=delta, their separate
exponential bounds yield, in the quadratic-form order,

    exp(-delta)M0 <= M_gamma <= exp(delta)M0,
    M=diag(V,K).

This uses the full exponentials, not a quadratic
approximation to the Hamiltonian.
