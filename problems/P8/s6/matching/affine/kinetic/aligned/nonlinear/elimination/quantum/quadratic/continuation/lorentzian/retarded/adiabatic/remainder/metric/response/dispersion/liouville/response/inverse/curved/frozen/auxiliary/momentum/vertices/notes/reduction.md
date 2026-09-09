# Actual stationary and spatial physical reduction

The S6.75 Hamiltonian per hat spatial volume is
H(N,Y)=A(N)+L(N,Y)+Q(N,Y). There are nine independent
invariants, nine linear monomials and four quadratic
monomials, in addition to the background A.

Here L and Q mean degree one and two in those independent
invariants, not physical field degree. Every invariant
has zero physical background. Shear, electric/magnetic
squares, vector mass and matter-gradient squares start
at physical degree two. Curvature, trace and matter
momentum deviations, and vector momentum divergence
can have linear terms. Thus an omitted invariant term
of degree at least five cannot contribute to H4 after
full physical substitution.

## Original lapse jets, not replacement clock jets

For I(u,1)=0, differentiate its original defining primitive:

I_N(u,N)=-N^-2 I_s(u,N^-1).

The fixed lower endpoint implies
I_phi(u,N)=partial_u I(u,N) at fixed N.
All its derivatives through N order four are retained.
In particular I_phi does not vanish at the bounce merely
because I(u=0,N) vanishes. The closed S6.75 expression
is independently replayed coefficient by coefficient.

The original lower vector source satisfies

Q_x=forcing(u,x)-coefficient(u,x)*Q,
x=-N^-2, x_N=2N^-3, Q(u,-1)=0.

A finite ordinary Taylor recursion integrates this ODE
through N order four. Symmetric nilpotent labels provide
exact derivatives: the coefficient of any k distinct
labels in f(1+sum epsilon_i) is f^(k)(1).
No interpolation or floating evaluation is used.

The complete I derivative must be simplified before
converting to QQ(u): individual algebraic time prefactors
may be irrational while their product is rational.
The implementation differentiates that complete expression.
No algebraic branch is changed to force a cancellation.

The resulting 14 coefficient rows each have five lapse
derivatives. Direct closed-bounce differentiation checks
all 70 values. A1=0 and A2=-2J_e are replayed at symbolic
time, with J_e=J+4*10^-6/h^2 and J_e>1/40 on I.

At all symbolic times, all 14 on-clock monomials equal
the complete S6.76 on-clock Hamiltonian. Its four possible
linear physical lapse-force coefficients also agree,
as does A2. Together with the identical spatial geometry,
momentum solution and moving boundaries, this proves
the quadratic bridge for every mode and direction before
any phase-pair sampling. Additional full physical samples
at -1/2,0,1/2 test that bridge independently.

## Stationary series and full physical geometry

Let A_k,L_k,Q_k denote actual N derivatives at N=1.
Then, with invariant weights retained,

n1=-L1/A2,
F2=A3*n1^2/2+L2*n1+Q1,
n2=-F2/A2,
n3=-(A3*n1*n2+A4*n1^3/6+L2*n2
      +L3*n1^2/2+Q2*n1)/A2.

The reduced Hamiltonian has invariant pieces

H0=A0,
H1=L0,
H2=Q0-L1^2/(2A2),
H3=A3*n1^3/6+L2*n1^2/2+Q1*n1,
H4=A4*n1^4/24+L3*n1^3/6+Q2*n1^2/2-F2^2/(2A2).

These generic formulas were independently derived before
the present model. Here their actual physical substitution
also checks the lapse force through degree n-1 and compares
direct stationary substitution with H through n=2,3,4.
Only the linear physical lapse is needed for H2.

The metric inverse, determinant volume, Christoffel symbols
and scalar curvature are all expanded from
g=(1+2v)I+gamma_TT. No curvature shortcut is used.
S6.76 supplies the complete matter-vector York recursion.

With M=pi*g/sqrt(g), the actual invariants are

dp=2 tr(M)/3+2H,
dc=(ell+P_chi)/sqrt(g)-ell,
j=div(Pi)/sqrt(g),
Sigma=tr(M^2)-tr(M)^2/3,
e=g_ij Pi^i Pi^j/(zeta*det(g)),
b=zeta*g^ik*g^jl*F_ij*F_kl,
V=g^ij W_i W_j,
Z=g^ij partial_i chi partial_j chi,
R=R3[g].

Here zeta=10^-6. The background matter momentum density
and all density factors are retained. Multiply H by
sqrt(g), then retain the original moving canonical terms

-2H*pi:g
-(H'+3H^2)*(6v+3v^2-gamma:gamma/2)
-ell*P_chi.

Time derivatives are taken before evaluation at a chosen
time. The current canonical momenta are those of S6.76,
before comparison with the old mixed-boundary variables.
No further matter momentum shift is made in these vertices.

A possible fourth-order longitudinal metric momentum
contracts only the isotropic background first variation
in H4 and hence vanishes by flat trace. This is why the
third-order spatial solution is sufficient; the zero
total wavevector is never inverted.

## Independent vector quartic and Fourier reality

At the bounce, take flat spatial metric, no matter
fluctuation, W=0 and transverse electric momenta.
The matter/vector momentum current vanishes, j=0,
and e=Pi^2/zeta. The first lapse force is then e/4,
of physical degree two. Consequently

H4=(Pi^2)^2/(64*J_e*zeta^2).

For four labelled unit Pi_z legs in a nonexceptional
planar quadrilateral, the integrated coefficient is

24/(64*J_e*zeta^2)
=187500000000000000/749377.

The general reduction exactly reproduces this nonzero
coefficient. Omitting the lapse response misses it.

Mixed fixtures need not have a real individual Fourier
coefficient: odd spatial derivatives produce imaginary
terms. The real functional instead obeys
H(k1,...,kn)^*=H(-k1,...,-kn).
Both complete momentum reversal and permutation of all
labels/fields are checked for mixed cubic and quartic
fixtures. The general property follows from real source
coefficients, real spatial derivatives, the real even
York inverse and commutative multilinear extraction.
