"""Physical temporal constraint, positive reference energy and causal estimates."""

from functools import cache

import sympy as s


@cache
def data():
    a, m, K, z = s.symbols("a m kappa zeta", positive=True)
    A0, j0, divp = s.symbols("A0 J0 div_pi", real=True)
    A, pi, j, grad0, F2 = s.symbols("A_i pi_i J_i partial_i_A0 F_ij_squared", real=True)
    velocity = s.Symbol("A_i_dot", real=True)
    # One electric component plus the full magnetic square notation.
    L = (
        a * (velocity - grad0) ** 2 / 2
        - F2 / (4 * a)
        + a**3 * m * m * A0 * A0 / 2
        - a * m * m * A * A / 2
        - a**3 * A0 * j0
        + a * A * j
    )
    pi_actual = s.diff(L, velocity)
    H = s.expand((pi * velocity - L).subs(velocity, pi / a + grad0))
    H_integrated = H.subs(pi * grad0, -A0 * divp)
    A0star = (j0 - divp / a**3) / (m * m)
    Epart = pi * pi / (2 * a) + F2 / (4 * a) + a * m * m * A * A / 2
    reduced = Epart + (divp - a**3 * j0) ** 2 / (2 * a**3 * m * m) - a * A * j
    S0, Ssp2 = s.symbols("S0 spatial_S_squared", real=True)
    contact = -K * a**3 * (S0 * S0 - Ssp2 / (a * a)) / 2
    source_contact = (a**3 * j0 * j0 / (2 * m * m) + contact).subs(
        {m * m: 1 / z, j0: s.sqrt(K / z) * S0}, simultaneous=True
    )
    Ee, EB, Em, ED, Hubble = s.symbols(
        "E_electric E_magnetic E_mass E_div Hubble", real=True
    )
    Jnorm, DJnorm, G0norm = s.symbols(
        "spatial_J_norm div_J_norm grad_J0_norm", nonnegative=True
    )
    forcing = s.sqrt(a * Jnorm**2 + DJnorm**2 / (a * m * m) + a * G0norm**2 / (m * m))
    B = s.Rational(5, 4) ** 6
    return {
        "canonical_parameters": "A=sqrt(kappa*zeta)W, J=sqrt(kappa/zeta)S, zeta=1e-6 and mass1000",
        "uneliminated_physical_Lagrangian": L,
        "temporal_constraint_solution": A0star,
        "reduced_driven_Hamiltonian_without_light_contact": reduced,
        "separate_retained_light_contact_Hamiltonian": contact,
        "temporal_source_constant_plus_full_contact": s.simplify(source_contact),
        "positive_reference_energy": "E0=integral[pi^2/(2a)+F_ij^2/(4a)+a*m^2*A_i^2/2+(div pi)^2/(2a^3*m^2)] d^3x",
        "exact_explicit_time_energy_derivative": Hubble * (-Ee - EB + Em - 3 * ED),
        "exact_source_power": "integral[pi.J_sp+a A_sp.grad J0+(div pi)(div J_sp)/(a^2 m^2)] d^3x; the antisymmetric curl contracted with grad J0 integrates to zero",
        "forcing_norm": forcing,
        "differential_bound": "v=sqrt(2E0): v'<=3|H|v/2+forcing_norm; use sqrt(2E0+epsilon^2) before epsilon->0, not division by zero energy",
        "exact_CD_propagation_majorant": B,
        "retarded_energy_bound": "v(t)<=integral_(t0)^t exp[(3/2) integral_s^t |H|] forcing_norm(s) ds <= (5/4)^6 integral forcing_norm",
        "constraint_boundary": "No conserved-source or transverse-current assumption. J0 and div J_sp remain in the temporal equation and forcing norm.",
        "checks": {
            "canonical_electric_momentum": s.expand(pi_actual - a * (velocity - grad0)),
            "literal_temporal_Euler_constraint": s.factor(
                s.diff(H_integrated, A0) + divp + a**3 * m * m * A0 - a**3 * j0
            ),
            "literal_temporal_elimination": s.factor(
                H_integrated.subs(A0, A0star) - reduced
            ),
            "positive_energy_source_decomposition": s.expand(
                reduced
                - Epart
                - divp**2 / (2 * a**3 * m * m)
                + j0 * divp / (m * m)
                - a**3 * j0 * j0 / (2 * m * m)
                + a * A * j
            ),
            "full_light_contact_temporal_constant_cancellation": s.simplify(
                source_contact - K * a * Ssp2 / 2
            ),
            "longitudinal_energy_weight": s.diff(a**-3, a) * a + 3 * a**-3,
            "CD_absolute_H_integral": s.integrate(
                4 * s.Symbol("t", real=True) / (1 + s.Symbol("t", real=True) ** 2),
                (s.Symbol("t", real=True), 0, s.Rational(1, 2)),
            )
            * 2
            - 4 * s.log(s.Rational(5, 4)),
            "exact_CD_propagation_exponential": s.exp(6 * s.log(s.Rational(5, 4))) - B,
            "maximum_spatial_volume": s.Rational(25, 16) ** 3 - B,
        },
    }


@cache
def modes():
    a, m, k = s.symbols("a m k", positive=True)
    q, p, j0, j, H = s.symbols("A pi J0 J Hubble", real=True)
    E_L = p * p / (2 * a) + a * m * m * q * q / 2 + k * k * p * p / (2 * a**3 * m * m)
    dq = p / a + k * k * p / (a**3 * m * m) - k * j0 / (m * m)
    dp = -a * m * m * q + a * j
    explicit = s.diff(E_L, a) * a * H
    power = p * j + k * k * p * j / (a * a * m * m) - a * k * q * j0
    actual = s.diff(E_L, q) * dq + s.diff(E_L, p) * dp + explicit
    E_T = p * p / (2 * a) + (a * m * m + k * k / a) * q * q / 2
    dqT = p / a
    dpT = -(a * m * m + k * k / a) * q + a * j
    actualT = s.diff(E_T, q) * dqT + s.diff(E_T, p) * dpT + s.diff(E_T, a) * a * H
    time, omega, nu = s.symbols("time omega nu", positive=True)
    particular = (s.sin(nu * time) - (nu / omega) * s.sin(omega * time)) / (
        omega * omega - nu * nu
    )
    resonant = (s.sin(omega * time) - omega * time * s.cos(omega * time)) / (
        2 * omega**2
    )
    return {
        "longitudinal_energy": E_L,
        "longitudinal_dq": dq,
        "longitudinal_dp": dp,
        "longitudinal_source_power": power,
        "transverse_energy": E_T,
        "transverse_dq": dqT,
        "transverse_dp": dpT,
        "finite_time_resonant_response": resonant,
        "scope": "Independent real sin/cos normalized Fourier modes retain the temporal source; arbitrary all-momentum energy proof is not replaced by this grid.",
        "checks": {
            "longitudinal_full_energy_identity": s.factor(actual - explicit - power),
            "transverse_full_energy_identity": s.factor(
                actualT - s.diff(E_T, a) * a * H - p * j
            ),
            "finite_time_resonant_limit": s.simplify(
                s.limit(particular, nu, omega) - resonant
            ),
            "resonant_forced_equation": s.simplify(
                s.diff(resonant, time, 2)
                + omega * omega * resonant
                - s.sin(omega * time)
            ),
            "resonant_zero_initial_position": resonant.subs(time, 0),
            "resonant_zero_initial_velocity": s.diff(resonant, time).subs(time, 0),
        },
    }
