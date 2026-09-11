# Complete specified local tensor Hessian

Set R0=6(H'+2H^2) and A=5m^2/12-R0/36.
The constant m^4 term has no unimodular tensor
variation. Since the first TT variation of R is
zero, the R^2 quadratic term is twice R0 times
the second scalar-curvature term. The exact local
quadratic action is consequently

    S_loc^(2)=1/(64pi^2) integral a^3 {
        A tr[gamma_t^2-a^-2(grad gamma)^2]
        -(1/60)tr[(D gamma)^2]}.

This is the whole tensor Hessian of the specified
finite local density, not the whole Gaussian
determinant. In the physical unit-Frobenius canonical
chart gamma=2h/sqrt(kappa), its prefactor becomes
1/(16pi^2 kappa).

With the actual a^3 dt dx pairing, the formal adjoint is

    Ddag=partial_t^2+5H partial_t+2H'+6H^2-a^-2 Delta.

It is not D. Define L=partial_t^2+3H partial_t-a^-2 Delta
and L_A=a^-3 partial_t(a^3 A partial_t)-A a^-2 Delta.
Then L_A=A L+A'partial_t. Independent Euler variation
of the literal quadratic action gives the local
canonical correction

    Qloc=(L_A+Ddag D/60)/(8pi^2 kappa).

The full Euler sign is -Lh-Qloc h+source=0 before
adding the remaining determinant response. It is
fixed by action variation, not the noise convention.

The complete fourth-order composition is

    Ddag D=partial_t^4+6H partial_t^3
       +(4H'+11H^2)partial_t^2
       +(H''+7HH'+6H^3)partial_t
       -2a^-2 Delta partial_t^2-2H a^-2 Delta partial_t
       +a^-4 Delta^2.

The pure Delta zeroth-time coefficient cancels after
differentiating a^-2. All derivatives and contacts
in this local operator are retained. On a Fourier
component the same expression follows from two
conformal wave operators; independent compact
action variations verify its pairing.

No order reduction, extra-pole deletion or four-data
branch selection is performed. Such an operator
identity and norm do not supply a full inverse.
