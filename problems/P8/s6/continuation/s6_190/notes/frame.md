# Full canonical matrix frame and pure-state covariance

Retain the actual S6.188 Hamiltonian with M=diag(V,K),
K=a^-1 exp(gamma)+kk^T/(a^3 m^2) and
V=a m^2 exp(-gamma)+a^-1 Ck^T exp(gamma) Ck.
Its proved identity K V=omega^2 I, with
omega^2=m^2+a^-2 k^T exp(-gamma) k, concerns instantaneous
positive matrices. It is not a dynamical polarization gap.

Let B be the unique real positive square root of K.
The time-dependent canonical transformation q=B^-1 A,
p=B pi preserves the symplectic form. Direct differentiation,
including the derivative of the frame, gives

    q'=p-Lq,  p'=-omega^2 q+L^T p,  L=B^-1 B'.

Set Q=sqrt(omega)q, P=p/sqrt(omega), rho=omega'/(2omega).
The real generator has blocks
[rho I-L, omega I; -omega I, L^T-rho I].
Consequently b=(Q+iP)/sqrt(2) and its conjugate satisfy

    b'=(-i omega I+R)b+S bdag,
    R=(L^T-L)/2, S=rho I-(L+L^T)/2.

R is real skew, S real symmetric. All three physical modes
are retained. General time-dependent matrices need not
commute, and no derivative of a chosen polarization
eigenbasis or separation among polarization frequencies
is used.

Write U,V for the coefficients of the same initial
annihilators in b,bdag. The preserved column CCR gives

    Udag U-Vdag V=I, U^T V-V^T U=0.

U is invertible. The graph r=V U^-1 is complex symmetric,
I-rdag r is positive, and ||r||op<1. Direct differentiation
yields the exact nonlinear equation

    r'=2i omega r+[R,r]+S-r S r.

The same CCR gives U Udag=(I-rdag r)^-1. With
F(r)=[I+r; -i(I-r)]/sqrt(2), the complete real
symmetrized six-quadrature covariance is

    Sigma(r)=Re[F(r)(I-rdag r)^-1 F(r)dag].

This includes cross-covariances. Irrelevant common
positive-frequency phases and unitary column rotations
cancel, not physical polarization mixing. This is the
actual propagated pure Gaussian state; a finite reference
introduced later is not substituted for its initial data.

For ||r||,||s||<=rho0=1/10, let C(r)=(I-rdag r)^-1.
Then ||C||<=1/(1-rho0^2)<2,
||F||^2=||I+rdag r||<=1+rho0^2<4, and
||F(r)-F(s)||=||r-s||. The inverse difference identity gives

    ||C(r)-C(s)|| <= 2rho0/(1-rho0^2)^2 ||r-s||
                  < ||r-s||.

Expanding the two F factors and the inverse separately
bounds the covariance difference by
(4+4+4)||r-s||<16||r-s||. Taking real parts cannot
increase its operator norm on real quadratures.
