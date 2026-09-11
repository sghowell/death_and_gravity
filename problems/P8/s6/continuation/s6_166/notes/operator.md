# First-quantized operator-norm limits and their errors

On the one-particle Hilbert space L2(R^3,C^4), let

    H(t)=alpha.p+beta M(t),
    H_asym=alpha.p+beta M_asym.

The mass is bounded and real. The Hamiltonians are selfadjoint
on the common first-order Sobolev domain. Fourier decomposition
and the unitary two-level evolution define the propagator U.
The difference H(t)-H_asym is a bounded multiplication operator,
with norm |M(t)-M_asym|, uniformly in momentum.

For either tail use

    W(T)=U(0,T) exp(-i H_asym T).

Its strong derivative on the common domain is

    W'(T)=i U(0,T)(H(T)-H_asym) exp(-i H_asym T).

Duhamel's formula extends the integral estimate to all of L2.
The real-profile inequality from S6.165 gives, for |T|>=tau,

    integral_(|T|)^infinity |M(t)-M_asym|dt
       <= Delta tau^8/(56|T|^7).

The orientation of the negative tail changes no norm.
Consequently W(T) is operator-norm Cauchy. Its adjoint is
norm Cauchy too, so the limit W_asym and its adjoint satisfy
both inverse identities. W_asym is a unitary ONE-PARTICLE
operator, not merely a modewise limiting matrix.

If P_asym is either asymptotic spectral projector,

    c(T)=U(0,T)P_asym U(T,0)
        =W(T)P_asym W(T)^*.

The two-factor difference identity and ||P_asym||=1 give

    ||c_asym-c(T)|| <= Delta tau^8/(28|T|^7).

The limits are complementary orthogonal projectors, so they
give the same pure quasifree CAR states as the earlier
mode construction and the external theorem's projector limits.

At the actual Delta<=3e197, tau=1e-100, one has
Delta tau^8<=3e-603. Thus at |T|=1,

    ||W_asym-W(T)|| <=3e-603/56 <1e-604,
    ||c_asym-c(T)|| <=3e-603/28 <2e-604.

This bound is independent of momentum. It is not a bound
on derivatives of a distributional kernel or on stress.
A unitary on the one-particle space need not be implementable
as a unitary in the chosen global infinite-volume Fock
representation. That separate claim is not made.

Independent finite-interval two-level calculations compare
transported asymptotic projectors at several endpoints,
with step refinement and norm preservation controls.
They test a toy parameter regime, not the actual extreme
frequencies, and are not validated integration.
