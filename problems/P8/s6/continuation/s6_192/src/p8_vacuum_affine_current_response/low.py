"""Momentum-independent covariance growth and exact finite-band response."""

from functools import cache

import sympy as s
from p8_vacuum_affine_matrix_adiabatic import initial
from p8_vacuum_affine_matrix_response_tail import mixed

from . import vertices

MASS = initial.jets.MASS
PARTITION = initial.PARTITION
SIGMA = (4 * s.Integer(10) ** 12, 4 * s.Integer(10) ** 13, 3 * s.Integer(10) ** 14)
CURRENT = (s.Integer(10) ** 14, s.Integer(10) ** 15, s.Integer(10) ** 16)


def generator(omega, rotation, squeeze):
    unit = s.eye(rotation.rows)
    return (
        (rotation + squeeze)
        .row_join(omega * unit)
        .col_join((-omega * unit).row_join(rotation - squeeze))
    )


def covariance_vector_field(A, C):
    return A * C + C * A.T


@cache
def constants():
    c = mixed.constants()
    a1 = 3 + (c["R"][0, 1] + c["S"][0, 1]) / MASS
    a2 = 6 + (c["R"][0, 2] + c["S"][0, 2]) / MASS
    g0, g1, g2 = vertices.G
    s0, s1, s2 = SIGMA
    current = (
        9 * g0 * s0,
        9 * ((g0 + g1) * s0 / MASS + g0 * s1),
        9
        * ((2 * g0 + 2 * g1 + g2) * s0 / MASS**2 + 2 * (g0 + g1) * s1 / MASS + g0 * s2),
    )
    integral = tuple(
        s.Rational(4, 18) * CURRENT[a] * PARTITION ** (4 + a) / (4 + a)
        for a in range(3)
    )
    return {
        "generator_parameter_bounds_over_nu": (a1, a2),
        "covariance_initial_trace": s.Integer(108),
        "covariance_zero_growth_display_before_rounding": 108 * 3**22,
        "covariance_first_second_before_rounding": (
            SIGMA[0] * 8,
            SIGMA[0] * (56 / MASS + 64),
        ),
        "current_integrand_before_rounding": current,
        "finite_band_complete_current_integrals": integral,
    }


@cache
def data():
    R = s.Matrix([[0, 1, -2], [-1, 0, 3], [2, -3, 0]])
    S = s.Matrix([[2, 1, 0], [1, -1, 2], [0, 2, 3]])
    w = s.Symbol("omega", positive=True)
    A = generator(w, R, S)
    J = s.zeros(3).row_join(s.eye(3)).col_join((-s.eye(3)).row_join(s.zeros(3)))
    eps = s.Symbol("epsilon", real=True)
    A1 = generator(s.Integer(2), R / 2, S / 3)
    A2 = generator(s.Integer(3), R / 3, S / 5)
    C0 = s.diag(1, 2, 3, 4, 5, 6)
    C1 = s.diag(S, S)
    C2 = s.diag(S * S, S * S)
    Ap = A + eps * A1 + eps**2 * A2 / 2
    Cp = C0 + eps * C1 + eps**2 * C2 / 2
    rhs = covariance_vector_field(Ap, Cp).applyfunc(s.expand)
    first = covariance_vector_field(A, C1) + covariance_vector_field(A1, C0)
    second = (
        covariance_vector_field(A, C2)
        + 2 * covariance_vector_field(A1, C1)
        + covariance_vector_field(A2, C0)
    )
    t, u, v, rate, b1, b2 = s.symbols("t u v rate b1 b2", positive=True)
    first_time = s.integrate(2 * b1, (u, 0, t))
    second_time = s.integrate(2 * b2 + 8 * b1**2 * u, (u, 0, t))
    n, k = s.symbols("nu K", positive=True)
    checks = {
        "fast_momentum_and_rotation_are_skew": A + A.T - 2 * s.diag(S, -S),
        "full_balanced_generator_is_symplectic": A * J + J * A.T,
        "complete_first_covariance_parameter_equation": rhs.applyfunc(
            lambda x: x.coeff(eps, 1)
        )
        - first,
        "complete_second_covariance_parameter_equation": rhs.applyfunc(
            lambda x: 2 * x.coeff(eps, 2)
        )
        - second,
        "first_Duhamel_common_exponent_integral": first_time - 2 * b1 * t,
        "second_Duhamel_common_exponent_integral": second_time
        - 2 * b2 * t
        - 4 * b1**2 * t**2,
        "no_spurious_composed_full_slab_exponent": s.exp(rate * (t - u))
        * s.exp(rate * (u - v))
        * s.exp(rate * v)
        - s.exp(rate * t),
    }
    for a in range(3):
        checks[f"finite_band_radial_integral_{a}"] = s.integrate(
            n ** (3 + a), (n, 0, k)
        ) - k ** (4 + a) / (4 + a)
    c = constants()
    return {
        "state": "The same S55 all-order state has initial balanced covariance trace<=108; the initial covariance and its first two amplitude derivatives are unchanged because the whole initial shear-jet neighborhood is common.",
        "uniform_propagator": "The exact balanced generator has skew fast part and squeeze norm<=11, so ||U(t,s)||<=exp(11|t-s|) for every momentum; no exponential in momentum or amplitude is used.",
        "full_actual_covariance_displays": SIGMA,
        "covariance_powers": "||Sigma^(a)|| < SIGMA[a]*nu_minus^a for a0,1,2, uniformly on the unit CD slab. Common-exponent Duhamel estimates retain the second-order iterated first response.",
        "actual_current_integrand_displays": CURRENT,
        "current_powers": "|partial_epsilon^a J_D|<CURRENT[a]*nu_minus^(1+a); this is unrenormalized at each finite momentum, not an integrable ultraviolet claim.",
        "constants": c,
        "checks": checks,
        "gates": {
            "generator_first_second_displays": c["generator_parameter_bounds_over_nu"][
                0
            ]
            < 4
            and c["generator_parameter_bounds_over_nu"][1] < 28,
            "actual_initial_trace_growth": c[
                "covariance_zero_growth_display_before_rounding"
            ]
            < SIGMA[0],
            "covariance_first_second_displays": all(
                c["covariance_first_second_before_rounding"][a] < SIGMA[a + 1]
                for a in range(2)
            ),
            "complete_physical_current_integrand_displays": all(
                c["current_integrand_before_rounding"][a] < CURRENT[a] for a in range(3)
            ),
            "same_fixed_positive_mass_and_analysis_partition": MASS == 1000
            and PARTITION == 10**16,
        },
    }
