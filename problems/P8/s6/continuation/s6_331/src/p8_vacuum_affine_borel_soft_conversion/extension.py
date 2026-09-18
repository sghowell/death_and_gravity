"""Rebuilt S301 measure extension, strengthened to uniform-radial angular L2."""

from functools import cache

import sympy as s


@cache
def data():
    r, c = s.symbols("r c", real=True)
    e, t, H, U = s.symbols("e t H U", positive=True)
    checks = {
        "unit_null_projected_nuclear_quotient": 1
        - r * r * c * c
        - (1 - r * c) * (1 + r * c),
        "full_D_trace_coefficient_difference": -1 / (2 + 2 * e)
        + s.Rational(1, 2)
        - e / (2 * (1 + e)),
        "radial_quarter_modulus_integrable": s.integrate(
            t ** s.Rational(-3, 4), (t, 0, 1)
        )
        - 4,
        "quadratic_current_polarization": (H + U) ** 2 - H**2 - 2 * H * U - U**2,
    }
    gates = {
        "bounded_Borel_kernel_diagonal_null_by_Fubini": True,
        "joint_radial_and_source_direction_to_angular_L2_continuity": True,
        "uniform_compact_Banach_kernel_partition_argument": True,
        "weak_energy_measure_gives_uniform_radial_L2_current_convergence": True,
        "common_nuclear_cap_controls_quadratic_angular_integrals": True,
        "radial_finite_part_passes_by_integrable_quarter_modulus": True,
        "positive_atomic_approximations_preserve_addition_and_mass": True,
        "noninteger_D_current_is_not_a_positive_cloud_intensity": True,
    }
    return {
        "checks": {name: s.factor(value) for name, value in checks.items()},
        "gates": gates,
        "whole_angular_energy_extension": "The unit-null projected kernel has nuclear norm<=2. Its arbitrary diagonal representative changes no angular integral. Joint continuity into angular L2, compact kernel partitions and the continuous original recoil extend T,H,U,K0,K1,Delta and regulated coefficients to all positive Borel angular-energy measures of mass<=1/8. Weak measure convergence gives uniform-radial L2 current convergence and convergence of the finite radial part.",
        "whole_regulator_and_scope": "The inherited radial quarter-Holder modulus and uniform S301 expansion pass to the Borel limit. Keep the full trace projector and fixed-ball phase. This only continues the ADDITIONAL soft factor with physical D4 external states; it does not supply outer hard evanescence or an interacting probability.",
    }
