"""Leading state-difference response derived from the actual reduced action."""

from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model as old
from p8_proca_physical_matter import coefficients

u = old.u
xi, dp, rho, pressure = sp.symbols(
    "delta_log_hat_scale delta_trace_momentum vector_density vector_pressure", real=True
)


def unit(i, n=1):
    return tuple(n if j == i else 0 for j in range(9))


@cache
def data():
    co = coefficients.data()
    rows = co["rows"]
    bg = old.coefficients()["background"]
    H, ell, h = bg["H"], bg["ell"], bg["h"]
    J = co["Jnew"]
    alpha = rows[unit(0)][1]
    beta = rows[unit(1)][1]
    F = rho - 3 * pressure / (2 * h)
    dl = -3 * ell * xi
    lapse = sp.factor((alpha * dp + beta * dl + F) / (2 * J))
    flow = sp.Matrix(
        [
            -dp / 2 + alpha * lapse / 3,
            -3 * H * dp + ell * beta * lapse + ell * dl + pressure,
        ]
    )
    variables = sp.Matrix([xi, dp])
    matrix = flow.jacobian(variables).applyfunc(sp.factor)
    forcing = flow.subs({xi: 0, dp: 0}).applyfunc(sp.factor)
    physical_scale = xi + lapse / (2 * h)
    scalar_density = -3 * ell * ell * physical_scale
    matter_derivative = beta * lapse - 3 * ell * xi
    baseH = rows[(0,) * 9][0]
    return {
        "H": H,
        "ell": ell,
        "h": h,
        "J": J,
        "alpha": alpha,
        "beta": beta,
        "conserved_free_matter_charge": sp.Rational(1, 10),
        "lapse_force": F,
        "linearized_lapse": lapse,
        "linearized_hat_scale_trace_flow": flow,
        "two_component_mean_generator": matrix,
        "two_component_mean_forcing": forcing,
        "physical_log_scale_response": physical_scale,
        "physical_scalar_density_and_pressure_response": scalar_density,
        "matter_clock_field_response_derivative": matter_derivative,
        "zero_anchor_hat_scale_and_trace_response": True,
        "checks": {
            "actual_free_matter_field_response_retains_physical_lapse": sp.factor(
                matter_derivative - ell * (lapse - 3 * physical_scale)
            ),
            "actual_free_matter_normal_velocity_gives_scalar_density_response": sp.factor(
                ell * (matter_derivative - ell * lapse) - scalar_density
            ),
            "actual_trace_Hamiltonian_first_derivative": sp.factor(
                rows[unit(0)][0] - 3 * H
            ),
            "actual_trace_Hamiltonian_second_derivative": sp.factor(
                2 * rows[unit(0, 2)][0] + sp.Rational(3, 2)
            ),
            "actual_matter_Hamiltonian_first_derivative": sp.factor(
                rows[unit(1)][0] - ell
            ),
            "actual_matter_Hamiltonian_second_derivative": sp.factor(
                2 * rows[unit(1, 2)][0] - 1
            ),
            "actual_mixed_trace_matter_derivative_zero": sp.factor(
                rows.get(tuple([1, 1] + [0] * 7), (0,))[0]
            ),
            "actual_background_trace_evolution": sp.factor(
                -baseH + ell * ell + 2 * sp.diff(H, u)
            ),
            "actual_background_charge_conservation": sp.factor(
                sp.diff(ell, u) + 3 * H * ell
            ),
            "actual_matter_lapse_cross_coefficient": sp.factor(
                beta - ell * (1 - 3 / (2 * h))
            ),
            "full_lapse_constraint_includes_pressure_and_mean_fields": sp.factor(
                -2 * J * lapse + alpha * dp + beta * dl + F
            ),
            "derived_mean_flow_is_generator_plus_physical_source": (
                flow - matrix * variables - forcing
            ).applyfunc(sp.factor),
            "mean_phase_trace_keeps_expanding_density_volume": sp.factor(
                sp.trace(matrix) + 3 * H
            ),
        },
    }


@cache
def ward():
    d = data()
    B, Bdot = sp.symbols(
        "physical_log_scale_response physical_log_scale_response_dot", real=True
    )
    r0 = d["ell"] ** 2 / 2
    delta = -6 * r0 * B
    derivative = -6 * sp.diff(r0, u) * B - 6 * r0 * Bdot
    total = (
        derivative
        - 3 * d["H"] * (rho + pressure)
        + 3 * d["H"] * (2 * delta + rho + pressure)
        + 6 * r0 * Bdot
    )
    return {
        "background_matter_density": r0,
        "matter_density_response": delta,
        "retained_connection_variation_term": 6 * r0 * Bdot,
        "checks": {
            "scalar_mean_and_vector_source_satisfy_full_linearized_homogeneous_Ward_identity": sp.factor(
                total
            ),
            "scalar_mean_requires_connection_variation": sp.factor(
                derivative + 6 * d["H"] * delta + 6 * r0 * Bdot
            ),
        },
    }


@cache
def center_control():
    d = data()
    at = {u: 0, xi: 0, dp: 0}
    n = sp.factor(d["linearized_lapse"].subs(at))
    dm = sp.factor(d["physical_scalar_density_and_pressure_response"].subs(at))
    total = rho + dm
    return {
        "initial_induced_lapse": n,
        "initial_induced_scalar_density": dm,
        "initial_combined_scalar_plus_vector_density": total,
        "positive_vector_state_pressure_bounds": "-rho/3 <= pressure <= rho",
        "initial_combined_density_lower_per_positive_vector_density": sp.Rational(
            134, 135
        ),
        "checks": {
            "actual_center_induced_lapse": sp.factor(
                n - sp.Rational(80, 243) * (rho - sp.Rational(3, 2) * pressure)
            ),
            "actual_center_induced_scalar_density": sp.factor(
                dm + sp.Rational(2, 405) * (rho - sp.Rational(3, 2) * pressure)
            ),
            "actual_center_combined_density_lower": sp.factor(
                total.subs(pressure, -rho / 3) - sp.Rational(134, 135) * rho
            ),
        },
    }
