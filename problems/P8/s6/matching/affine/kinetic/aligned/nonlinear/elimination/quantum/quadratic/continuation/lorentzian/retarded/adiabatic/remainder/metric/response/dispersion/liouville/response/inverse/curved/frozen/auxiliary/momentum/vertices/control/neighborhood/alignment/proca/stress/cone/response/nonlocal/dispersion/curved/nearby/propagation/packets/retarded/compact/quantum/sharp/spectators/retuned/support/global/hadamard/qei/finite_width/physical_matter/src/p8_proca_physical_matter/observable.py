"""Physical normal-frame density and pressure, before and after lapse reduction."""

from functools import cache

import sympy as sp


@cache
def data():
    mass_scale, time_scale = sp.symbols("positive_M positive_tau", positive=True)
    action_factor = mass_scale**2 * time_scale**2
    N, h = sp.symbols("N h", positive=True)
    ell, dc, grad = sp.symbols("ell dc spatial_gradient_squared", real=True)
    ratio = (h - 1 + N**-2) / h
    conformal = ratio ** (-sp.Rational(1, 4))
    rho = (ell + dc) ** 2 / (2 * conformal**6) + grad / (2 * conformal**2)
    pressure = (ell + dc) ** 2 / (2 * conformal**6) - grad / (6 * conformal**2)
    eps = sp.Symbol("epsilon", real=True)
    n1, n2, c1, c2, g2 = sp.symbols("n1 n2 dc1 dc2 gradient2", real=True)

    def coefficients(expr):
        # First differentiate in independent invariant/lapse variables.
        value = sum(
            sp.diff(expr, N, j).subs(N, 1) * (N - 1) ** j / sp.factorial(j)
            for j in range(3)
        )
        series = sp.Poly(
            sp.expand(
                value.subs(
                    {
                        N: 1 + eps * n1 + eps * eps * n2,
                        dc: eps * c1 + eps * eps * c2,
                        grad: eps * eps * g2,
                    },
                    simultaneous=True,
                )
            ),
            eps,
        )
        return tuple(sp.factor(series.nth(j)) for j in range(3))

    rr, pp = coefficients(rho), coefficients(pressure)
    kinetic2 = (
        ell * c2
        + c1 * c1 / 2
        - 3 * ell * c1 * n1 / h
        - 3 * ell * ell * n2 / (2 * h)
        + sp.Rational(3, 4) * ell * ell * (1 + 3 * h) * n1 * n1 / h**2
    )
    A2, A3, F1, F2, L2, Q1 = sp.symbols("A2 A3 F1 F2 L2 Q1", real=True, nonzero=True)
    lapse1 = -F1 / A2
    lapse2 = -(F2 + A3 * lapse1 * lapse1 / 2 + L2 * lapse1 + Q1) / A2
    force = (
        A2 * (eps * lapse1 + eps**2 * lapse2)
        + A3 * (eps * lapse1) ** 2 / 2
        + eps * F1
        + eps**2 * F2
        + eps * L2 * (eps * lapse1)
        + eps**2 * Q1
    )
    return {
        "N": N,
        "h": h,
        "ell": ell,
        "dc": dc,
        "gradient": grad,
        "inherited_overall_action_factor": action_factor,
        "physical_stress_multiplier": mass_scale**2 / time_scale**2,
        "physical_density": rho,
        "physical_isotropic_pressure": pressure,
        "rho_orders": rr,
        "pressure_orders": pp,
        "second_order_lapse_coefficient_in_density": sp.diff(rr[2], n2),
        "universal_first_lapse": lapse1,
        "universal_second_lapse": lapse2,
        "checks": {
            "physical_stress_restoration_retains_overall_action_and_length_scale": sp.factor(
                action_factor / time_scale**4 - mass_scale**2 / time_scale**2
            ),
            "physical_conformal_volume_is_cube": sp.factor(
                conformal**3 - ratio ** (-sp.Rational(3, 4))
            ),
            "physical_normal_velocity_is_canonical_density_divided_by_actual_volume": sp.factor(
                rho - (ell + dc) ** 2 / (2 * conformal**6) - grad / (2 * conformal**2)
            ),
            "full_density_background": rr[0] - ell * ell / 2,
            "full_pressure_background": pp[0] - ell * ell / 2,
            "actual_density_linear_lapse_coefficient": sp.factor(
                rr[1] - ell * c1 + 3 * ell * ell * n1 / (2 * h)
            ),
            "actual_full_quadratic_density": sp.factor(rr[2] - kinetic2 - g2 / 2),
            "actual_full_quadratic_pressure": sp.factor(pp[2] - kinetic2 + g2 / 6),
            "quadratic_density_pressure_difference_is_actual_gradient": sp.factor(
                rr[2] - pp[2] - 2 * g2 / 3
            ),
            "nonzero_second_order_lapse_density_coefficient": sp.factor(
                sp.diff(rr[2], n2) + 3 * ell * ell / (2 * h)
            ),
            "full_second_order_stationary_lapse_force": sp.expand(force),
        },
    }
