# Physical constraint, contact and static normalization

Use ds^2=-du^2+a^2 dx^2 and physical spatial coordinate components W_i.
For independent mass coefficients a_m,b_m the Lagrangian density is

    L=a(W_i'-partial_i W_0)^2/2-F_ij F_ij/(4a)
      +a^3 m^2 a_m W_0^2/2-a m^2 b_m W_i W_i/2.

The canonical momentum is pi^i=a(W_i'-partial_i W_0). Integrating
the spatial derivative by parts in the Hamiltonian gives
-W_0 partial_i pi^i-a^3 m^2 a_m W_0^2/2. Hence

    W_0=-partial_i pi^i/(a^3 m^2 a_m),
    H=pi_i pi_i/(2a)+F_ij F_ij/(4a)+a m^2 b_m W_i W_i/2
      +(partial_i pi^i)^2/(2a^3 m^2 a_m).

This elimination is pointwise in the external mass source, even
when that source has arbitrary spatial dependence.

For a homogeneous source and physical q=|p_com|^2/a^2, the transverse
coordinate is W_T and the longitudinal coordinate is sigma,
W_i=i(p_com)_i sigma. Their Hamiltonians are

    H_T=Pi_T^2/(2a)+a(q+m^2 b_m)W_T^2/2,
    H_L=Pi_sigma^2/(2a^3 q)+Pi_sigma^2/(2a^3 m^2 a_m)
        +a^3 m^2 b_m q sigma^2/2.

The q=0 coordinate sigma is not used as an independent zero-mode
chart; the reconstructed field two-point function has a regular
isotropic limit. With the fixed original-clock canonical map
v=g w, p=Pi/g and omega^2=m^2+q,

    g_T^2=a, g_L^2=a^3 m^2 q/omega^2,
    J_T=-partial_n H_T=-m^2 beta v_T^2/2,
    J_L=-partial_n H_L=(alpha z p_L^2-beta omega^2 v_L^2)/2,
    z=q/omega^2.

Do not vary the moving canonical map when differentiating these
fixed-clock Hamiltonian source vertices. The independent frozen
physical-Hamiltonian check gives the same source generator.

Although the original mass insertion is linear in n, its eliminated
Hamiltonian is not. The local response contact is

    partial_n J|_0=-partial_n^2 H|_0
      =-alpha^2(partial_i pi^i)^2/(a^3 m^2),
    contact per longitudinal momentum=-alpha^2 z |p_L|^2.

There is no transverse mass-only contact. The displayed coincidence
expectation is bare and divergent, not a pre-renormalized observable.

## Flat independent check

For each flat oscillator, omega_T(n)^2=q+m^2(1+beta n) and
omega_L(n)^2=(1+beta n)[q/(1+alpha n)+m^2].
The static susceptibility is -partial_n^2(omega_sector/2)|_0.
The Abel-integrated bubble yields m^4 beta^2/(8 omega^3) for T
and omega(beta+alpha z)^2/8 for L. Adding the longitudinal contact
-alpha^2 z omega/2 reproduces the exact energy derivative in both
sectors. Without it the longitudinal answer is too large by
alpha^2 z omega/2. This is an algebraic nonzero control, not a
bound on the time-dependent selected state.
