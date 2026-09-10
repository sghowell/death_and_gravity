"""Full tensor contraction and fixed-negative-transfer one-loop pole subtraction."""

from functools import cache

import sympy as sp

from . import triangles


@cache
def data():
    E, p, Q = sp.symbols("P_energy P_spatial transfer_spatial", real=True)
    F, D, v, t = sp.symbols("spin_two_F stress_D crossing_v transfer_t", real=True)
    Pl = sp.Symbol("positive_Planck_mass", positive=True)
    eta = sp.diag(1, -1, -1, -1)
    P1 = sp.Matrix([E, p, 0, 0])
    P2 = sp.Matrix([E, -p, 0, 0])
    q = sp.Matrix([0, 0, 0, Q])
    T1 = 2 * F * (P1 * P1.T) + D * (q * q.T - eta * t) / 2
    T2 = 2 * F * (P2 * P2.T) + D * (q * q.T - eta * t) / 2
    full = sp.expand(
        sp.trace(T1 * eta * T2 * eta) - sp.trace(T1 * eta) * sp.trace(T2 * eta) / 2
    )
    aa = 1 - t / 4
    reduced = sp.expand(
        full.subs({E * E: (aa + v / 2) / 2, p * p: (v / 2 - aa) / 2, Q * Q: -t})
    )
    target = (
        F * F * (v * v - 2 * aa * aa)
        + F * D * aa * t
        - sp.Rational(3, 8) * D * D * t * t
    )
    eP, eQ, e2, local = sp.symbols(
        "e_dot_P e_dot_q e_squared local_metric_vertex", real=True
    )
    projected = 2 * F * eP**2 + D * (eQ**2 - e2 * t) / 2 + local * e2
    hbar = sp.Symbol("loop_order_marker", real=True)
    f1, slope, remainder = sp.symbols(
        "one_loop_subtracted_form_factor slope_at_zero transfer_remainder", real=True
    )
    amplitude = -((1 + hbar * f1) ** 2) / (Pl * Pl * t)
    one_loop = sp.diff(amplitude, hbar).subs(hbar, 0)
    expansion = one_loop.subs(f1, t * slope + t * t * remainder)
    d = triangles.data()
    g = sp.Symbol("equal_mass_cubic_squared", positive=True)
    single_identical_slope = (
        g * d["equal_mass_two_distinct_field_integral"] / (192 * sp.pi**2)
    )
    return {
        "literal_stress_projector_numerator": reduced,
        "fixed_t_forward_coefficient": -F * F / (Pl * Pl * t),
        "one_loop_pole_subtracted_coefficient": one_loop,
        "finite_negative_t_vertex_correction": -2 * slope / (Pl * Pl),
        "finite_t_convergence_error": -2 * t * remainder / (Pl * Pl),
        "single_identical_field_equal_mass_diagnostic_c2": -4
        * single_identical_slope
        / (Pl * Pl),
        "scope": "Only the one-loop matter-vertex dressing of the t-channel graviton pole. Tree pole is removed at fixed negative t before its limit. The f1-squared term is two-loop and not resummed. Other exchange channels, graviton loops and high-energy contour terms remain separate.",
        "checks": {
            "literal_all_component_stress_propagator_contraction": sp.factor(
                reduced - target
            ),
            "actual_crossing_coefficient_of_full_tensor_numerator": sp.diff(
                reduced, v, 2
            )
            / 2
            - F * F,
            "null_projection_removes_D_and_metric_local_terms": sp.expand(
                projected.subs({eQ: 0, e2: 0}) - 2 * F * eP**2
            ),
            "one_loop_not_squared_one_loop_resummation": sp.factor(
                one_loop + 2 * f1 / (Pl * Pl * t)
            ),
            "fixed_negative_t_subtraction_before_limit": sp.factor(
                expansion + 2 * slope / (Pl * Pl) + 2 * t * remainder / (Pl * Pl)
            ),
            "known_single_identical_scalar_equal_mass_normalization_control": sp.simplify(
                -4 * single_identical_slope / (Pl * Pl)
                + (45 - 8 * sp.pi * sp.sqrt(3)) * g / (1296 * sp.pi**2 * Pl * Pl)
            ),
        },
    }
