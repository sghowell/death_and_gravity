# Physical stress reconstruction, including every moving factor

For a sector and physical output, the frozen per-comoving-mode
readout is

    Q=(A |p_u|^2+B Omega^2 |v_u|^2)/(2a^3).

In S6.68 the coefficient variations delta A,delta B already
include the variation of the physical output normalization:
-3zeta for energy, -(n+3zeta) for pressure, and the corresponding
moving Hamiltonian weights. In particular all second mass
vertices remain. The a in the formula below is the background
a, so this normalization is not varied a second time.

Let D=k^2+U, Omega^2=f^2 D, and put
Y=(S,C,P)=(|v_sigma|^2,Re(p_sigma v_sigma*),|p_sigma|^2).
Then

    Q=f[A P+B D S]/(2a^3).

At fixed physical u,

    delta |v_u|^2=(delta S|_sigma+xi S'-r S)/f,
    delta |p_u|^2=f(delta P|_sigma+xi P'+r P),
    delta Omega^2=2r Omega^2+f^2 delta U|_u.

Substitution gives the exact physical variation

    delta Q=f/(2a^3) {
       A delta P|_sigma+B D delta S|_sigma
       +(delta A+r A)P+(delta B+r B)D S+B delta U|_u S
       +xi(A P'+B D S') }.

The history term here differentiates the covariance only.
The other coefficients are evaluated at fixed u and have their
own explicit variations. Replacing this history by xi Q' without
also changing those coefficient variations would double count
background derivatives.

The code constructs the response row, fixed-output contact row
and history row multiplying Y. All four sector/output rows
independently match the physical-time expression using the
original frequency variation and readout contacts. Multiplicity
is still two transverse plus one longitudinal, the comoving
measure is unchanged, and the common normalization by L^-2
is applied exactly as in S6.68.

## An independent action-vertex check

For the longitudinal coordinate before canonical normalization,
its acoustic pump is C_k. Express the varied physical action
in the unvaried acoustic coordinate. With chi=delta log A,
r=delta log f and delta U=U(2zeta+alpha n), its insertion is

    delta L_L =[(2chi-delta U/D-r)p_sigma^2
                       -(2chi+r)D v_sigma^2]/2.

Inserting chi=zeta+(alpha+beta)n/4 reproduces the original
current weights:

    p^2: (-1+alpha z)n+(1+2z)zeta,
    D v^2: -(1+beta)n-zeta, z=k^2/(k^2+U).

The transverse action similarly gives -r_T p^2/2 and
-[r_T D+delta U_T]v^2/2. These independently recover its two
Hamiltonian vertices. This first insertion check does not
replace the second-source contacts in the full readout above.

The covariance equation and this algebra do not by themselves
supply a bound for the ultraviolet momentum integral. That is
why the original physical subtraction is retained rather than
replacing it with a new auxiliary scalar prescription.
