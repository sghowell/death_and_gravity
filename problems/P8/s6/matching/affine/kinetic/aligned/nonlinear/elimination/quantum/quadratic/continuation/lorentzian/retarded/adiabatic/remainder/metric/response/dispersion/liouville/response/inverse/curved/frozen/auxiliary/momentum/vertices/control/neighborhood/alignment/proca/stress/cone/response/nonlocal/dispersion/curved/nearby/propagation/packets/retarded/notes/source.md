# Actual local relational source and the regular low-frequency chart

Use the unchanged S6.88 actual classical solution with central lapse
N(0)=1+10^-6 and W=Pi_W=0. Put T=10^-7 and work on the inner
real interval I=[-T/2,T/2]. No old quantum state or fixed profile
is transferred to this solution.

The background-dependent scalar

    O=chi-chi_background(phi)

vanishes on the background. Its linear perturbation is
delta chi-(chi_background'/phi_background')delta phi; it is
gauge invariant at this order and equals chi in unitary gauge.
The local probe action is integral sqrt(-g) J O. With J first
order, its quadratic unitary-gauge density is exactly

    R^3 N e^3 J chi.

There is no linear lapse, shift or metric source from a nonzero
background O. The clock and matter probe forces have vanishing
background Noether combination, as checked in S6.90. Thus using
the original linear constraints remains legitimate. J is an
external classical scalar probe, not a new UV-complete apparatus.

## Original canonical equations

Before the q-dependent clock momentum chart, X=(v,chi,p,P) has

    H2=m(P-3ell v)^2+(2/3)a ell p chi
       +(gq-aell^2)chi^2+2rqv^2-F1^2/(2h),
    F1=alpha p+beta(P-3ell v)+4r_N qv.

The canonical one-form has measure R^3. Consequently

    X'=M_X X+N e^3 e_4 J,
    M_X=Omega Hessian(H2)-diag(0,0,3Hhat,3Hhat).

source.py recovers this literal Hamiltonian from S6.89 and checks
that M_X is a polynomial of degree at most two in q. Its only
lapse denominator is the nonzero actual h. In particular, neither
a finite-q Legendre determinant nor Hhat is inverted.

For q=k^2/R^2>0 the full S6.90 transformation is

    Y=R^(3/2)(-kp/(2q), k chi, 2qv+Bp/2, P),
    B=alpha/r_N.

The derivative of R, q and B is included. Direct multiplication
verifies that (D'+D M_X)D^-1 is exactly the S6.90 Laurent
generator, with every lower-order term. It also verifies

    Y_source'=R^(3/2) N e^3 e_4 J,
    O_hat=R^(-3/2) Y_2/k.

These are absolute source/output normalizations, not a relative
packet estimate near a zero of a sine. The smooth polynomial
extension of M_X to k=0 is used only to bound a compact Fourier
ball. No separate claim about the strictly homogeneous scalar
constraint sector is needed at that measure-zero point.

## Retarded solution

For each k>0 solve the original first-order ODE with zero
canonical data before the source. It exists uniquely on I.
The high-frequency estimates below are polynomial bounds on
the fundamental transfer; on compact frequency balls the
ordinary ODE estimate is finite. Smooth compact sources have
rapidly decreasing spatial Fourier transforms. All finite
space and time derivatives of the sourced solution are therefore
well defined; time differentiation only adds polynomial
coefficients in k. This defines the classical retarded response

    O(t,x)=integral_{s<t} ds d^3y G(t,s;x-y) J(s,y).

The factors of physical spacetime volume are already in G's
source normalization. A detector integrated with physical
volume differs from a coordinate-volume test by a positive
smooth factor, which leaves compact support and the existence
of a nonzero pairing unchanged.
