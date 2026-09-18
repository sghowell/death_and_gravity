"""Total-energy rather than particle-count coefficient and difference bounds."""

from functools import cache

import sympy as s

from . import kernel, source


def coefficient_upper():
    return s.Integer(320) / source.KAPPA ** s.Rational(3, 2)


def coefficient_change_upper(total_energy):
    return 25000 * kernel.energy(total_energy) / source.KAPPA ** s.Rational(3, 2)


@cache
def data():
    R = s.Rational(1, 8)
    gates = {
        "mixed_pair_cap": 3 * 2 * 4 + 2 + 2 * 4 + 2 == 36,
        "massive_real_cap": s.Rational(6 * 56 * 102, 16 * 9) == 238,
        "mixed_real_cap": s.Rational(4 * 36 * 4, 4 * 9) == 16,
        "phase_real_cap": 49 * (24 + s.Rational(3, 2) * R) / 36 < 33,
        "full_real_cap": 238 + 16 * R + 5 * R * R / 36 + 33 < 275,
        "massive_imaginary_cap": s.Rational(2 * 2 * 102, 8 * 3) == 17,
        "mixed_imaginary_cap": s.Rational(2 * 2 * 36, 8 * 3) == 6,
        "imaginary_phase_cap": s.Rational(49 * 4, 4 * 3) < 17,
        "full_imaginary_cap": 17 + 6 * R + R * R / 6 + 17 < 35,
        "whole_coefficient_cap": 275 + 35 < 320,
        "ratio_difference_cap": 4 * 4 + 4 * 4 * 16 == 272,
        "B_difference_cap": 12 * 16 + 3 * 272 + 12 == 1020,
        "real_massive_difference_cap": s.Rational(
            6 * (161 * 26 * 102 + 56 * 2040), 16 * 9
        )
        < 23000,
        "new_null_real_difference_cap": 16 + 5 * R / 36 < 17,
        "phase_difference_cap": (866 * 25 + 48 * s.Rational(43, 2)) / 36 < 631,
        "real_difference_cap": 23000 + 17 + 631 < 24000,
        "imaginary_massive_difference_cap": s.Rational(4 * 102 + 2 * 2040, 8 * 3)
        == 187,
        "new_null_imaginary_difference_cap": 6 + R / 6 < 7,
        "imaginary_phase_difference_cap": s.Rational(866 * 4, 4 * 3) < 289,
        "imaginary_difference_cap": 187 + 7 + 289 < 500,
        "whole_difference_cap": 24000 + 500 < 25000,
        "all_finite_multiplicities_and_angular_approaches": True,
        "bound_not_a_uniform_full_loop_remainder": True,
    }
    E, W, Q2 = s.symbols("E W Q2", real=True)
    checks = {
        "exact_outgoing_pair_gamma": ((2 * E - W) ** 2 - (W * W - Q2) - 2) / 2
        - (2 * E * E - 2 * E * W + Q2 / 2 - 1),
        "outgoing_pair_gamma_lower_endpoint": 2 * s.Rational(5, 4) ** 2
        - 2 * s.Rational(5, 4) * R
        - 1
        - s.Rational(29, 16),
        "delta_projected_current": s.Integer(2 * 3 * 144 + 2 - 866),
        "real_phase_derivative_budget": s.Rational(2 * 4 * 5, 2)
        + s.Rational(3, 2)
        - s.Rational(43, 2),
        "symmetric_pair_difference": s.Integer(2 * 1020 - 2040),
    }
    return {
        "checks": {k: s.factor(v) for k, v in checks.items()},
        "gates": {k: bool(v) for k, v in gates.items()},
        "whole_uniform_dimensionless_bounds": (320, 25000),
        "whole_kinematic_bounds": (51, 102, 36, 866, 272, 1020, 2040, 26),
        "whole_uniform_statement": "Every finite positive null-radiation state with R<=1/8 on the same original recoil has |C_sigma|<320 and |C_sigma-C_Born|<25000R at fixed E/u, uniformly in directions off exactly soft-collinear atoms and in all approaches. The physical log coefficient has an additional kappa^(-3/2).",
        "whole_real_difference_budget": (23000, 17, 631, 24000),
        "whole_imaginary_difference_budget": (187, 7, 289, 500),
        "whole_proof_boundary": "Analytic compact-domain Doppler, TT contraction and derivative bounds supply uniformity. Exact arithmetic gates record their slack; calibrations are not the uniform proof. No finite hard-loop remainder, simultaneous multiple-soft theorem or all-N hard measure is implied.",
    }
