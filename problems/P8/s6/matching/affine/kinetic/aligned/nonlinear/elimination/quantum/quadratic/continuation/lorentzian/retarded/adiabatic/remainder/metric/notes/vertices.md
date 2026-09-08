# Independent physical vertices

Work in the original metric with lapse N, scale a_hat and comoving
momentum held fixed under variation. Write q=|k|^2/a_hat^2.
The physical oscillator Hamiltonian is

    H_mode=(P^2/g^2+g^2 Omega^2 Q^2)/2.

Before choosing moving canonical variables,

    g_T^2=a_hat^(D-2)/N,
    Omega_T^2=N^2(q+m^2 b_m),

    g_L^2=a_hat^D m^2 a_m q/[N(q+m^2 a_m)],
    Omega_L^2=N^2 b_m(q/a_m+m^2).

The longitudinal formulas are used at positive q and continued by
the regular isotropic canonical limit. No independent zero-momentum
longitudinal direction is selected.

Let D_Z=a_hat partial_a_hat-2q partial_q. The coordinate-density
currents are J_N=-partial_N H and J_Z=-D_Z H. They obey

    J_N=-a_hat^D rho, J_Z=D N a_hat^D p.

For p_mode=P/g=v'-d v, their readout is
(A p_mode^2+B Omega^2 v^2)/2, with

    A_X=partial_X log(g^2),
    B_X=-partial_X log(g^2 Omega^2).

These identities are independently differentiated directly from H,
before restricting to N=1. They include the contact terms.
The code compares every baseline and varied readout weight against
its separately written closed expression.

## Moving normalization

At the clock let z=q/(q+m^2), alpha=a_m,N, beta=b_m,N.
For sources n and zeta,

    delta log Omega_T=[1+beta(1-z)/2]n-z zeta,
    delta log Omega_L=[1+(beta-alpha z)/2]n-z zeta,

    delta log g_T=[(D-2)zeta-n]/2,
    delta log g_L=[(D-2+2z)zeta+(alpha z-1)n]/2.

The rate and frequency-rate variations are the full time derivatives
of these expressions. In particular, derivatives of the fixed
alpha and beta profiles are not discarded. With d=g'/g and
lambda=Omega'/Omega,

    delta U=(delta d)'+2d delta d, delta lambda=(delta log Omega)'.

Differentiating the Riccati expansion determines delta P2,delta P4
and all readout coefficients. The current coefficient at half-order
j also receives (1-2j) delta log Omega times its baseline value.
Transverse multiplicity remains D-1.

The existing second mass derivatives, for h=(1+u^2)^3, are

    a_m,NN=-4/(3h)-8/(9h^2),
    b_m,NN=-28/(27h)+152/(729h^2).

They enter the differentiated current weights even though they do
not enter the first frequency variation. Omitting them loses a
finite lapse contact. Its zero-order retained-minus-omitted
unit-source fixture at the bounce is -826/243 in the common
mass-fourth normalization.
