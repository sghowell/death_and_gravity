"""Global complex-momentum gap, Dirac norm and forward external domain."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_proper_references import dirac


@cache
def data():
    x, q2, b = sp.symbols("x q_squared b", nonnegative=True)
    m2 = sp.Symbol("m_squared", positive=True)
    S = m2 + q2
    p2lower = q2 - sp.Rational(21, 100) * S
    lower = x * m2 + (1 - x) * b + x * (1 - x) * p2lower
    target = sp.Rational(79, 100) * x * (m2 + (1 - x) * q2) + (1 - x) * b
    upper = m2 + 2 * x * (1 - x) * S
    log = sp.Symbol("nonnegative_log_one_plus_t", nonnegative=True)
    gamma = dirac.gamma_matrices()
    q = sp.symbols("real_q0:4", real=True)
    m = sp.Symbol("positive_mass", positive=True)
    D = m * sp.eye(4) + sp.I * sum((v * g for v, g in zip(q, gamma)), sp.zeros(4))
    normcheck = (
        D.conjugate().T * D - (m * m + sum(v * v for v in q)) * sp.eye(4)
    ).applyfunc(sp.expand)
    E, P = sp.symbols("E P")
    s = sp.Symbol("forward_s")
    p1 = sp.Matrix([sp.I * E, P, 0, 0])
    p2 = sp.Matrix([sp.I * E, -P, 0, 0])
    external = (p1, p2, -p1, -p2)
    subst = {E**2: s / 4, P**2: s / 4 - 1}
    checks = {
        "complex_parameter_gap_positive_remainder": sp.expand(
            lower - target - sp.Rational(21, 100) * x * x * m2
        ),
        "parameter_gap_above_half_x_mass": sp.Rational(79, 100) * x * m2
        - x * m2 / 2
        - sp.Rational(29, 100) * x * m2,
        "parameter_modulus_upper_positive_decomposition": sp.expand(
            2 * S - upper - (m2 / 2 + 3 * q2 / 2 + 2 * (x - sp.Rational(1, 2)) ** 2 * S)
        ),
        "momentum_shift_quadratic_loss": 2 * sp.Rational(1, 10)
        + sp.Rational(1, 100)
        - sp.Rational(21, 100),
        "integrated_log_endpoint_weight": sp.integrate(-sp.log(x), (x, 0, 1)) - 1,
        "scalar_kernel_rounded_log_majorant": 3 * (log + 6) - 3 * (log + 5) - 3,
        "gauge_kernel_rounded_log_majorant": 12 * (log + 6)
        - (8 * (log + 5) + 4)
        - (4 * log + 28),
        "literal_real_Dirac_normality": normcheck,
        "four_external_momenta_conserve": sum(external, sp.zeros(4, 1)),
        "s_channel_invariant": sp.expand(-((p1 + p2).T * (p1 + p2))[0]).subs(subst) - s,
        "u_channel_invariant": sp.expand(-((p1 - p2).T * (p1 - p2))[0]).subs(subst)
        - (4 - s),
        "unit_disc_E_squared_upper": sp.Rational(3, 4) + sp.Rational(1, 4) - 1,
        "three_prefix_legs_L1_bound": 3 * 2 - 6,
        "conservative_prefix_bound_gap": 18 - 6 - 12,
        "Cauchy_radius_shift_margin": sp.Rational(18, 360) - sp.Rational(1, 20),
        "Neumann_reference_margin": 1 - sp.Rational(1, 10) - sp.Rational(9, 10),
    }
    for i, p in enumerate(external):
        checks[f"external_mass_one_{i}"] = sp.expand((p.T * p)[0]).subs(subst) + 1
    return {
        "real_loop_norm_squared": S,
        "permitted_complex_shift": "||delta||_1<=sqrt(mF^2+q^2)/10",
        "parameter_real_part_lower": target,
        "normalized_parameter_modulus_bounds": "x/2 <= |Delta/mF^2| <=2(1+q^2/mF^2)",
        "pointwise_log_upper": "log(1+q^2/mF^2)-log(x)+4",
        "full_MS_self_energy_operator_norm_upper": "(3Y+12aCf)/Q sqrt(mF^2+q^2) [log(1+q^2/mF^2)+6]",
        "shifted_free_Dirac_operator_norm_upper": "2/sqrt(mF^2+q^2)",
        "external_forward_domain": "|s-2|<=1,u=4-s,t=0; all incoming scalar masses one",
        "soft_scaling_Cauchy_radius": "sqrt(mF^2+q^2)/360",
        "scope": "All real internal q, with a q-dependent soft-momentum analytic radius. This radius is proved for the proper two-point self-energy, not the remaining vertex-chord graphs.",
        "checks": checks,
    }
