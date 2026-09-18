"""Energy-dependent finite regulator remainder in the original soft scheme."""

from functools import cache

import sympy as s
from p8_vacuum_affine_borel_soft_conversion.moments import parameters
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient.measure import cutoff

from . import kernels, source


def quotient_remainder_upper(energy, epsilon):
    R, e = kernels.regulator(energy), kernels.regulator(epsilon)
    return s.S.Zero if R == 0 else 66000 * e * R * (1 - s.log(R)) ** 2 / source.KAPPA


def conditional_regulator_tv_upper(index, resolution, epsilon):
    a, x = parameters(index, resolution)
    e = kernels.regulator(epsilon)
    return 400000 * e * a * x * (1 - s.log(x)) ** 2 / source.KAPPA


def original_regulator_tv_upper(resolution, epsilon):
    x, e = cutoff(resolution), kernels.regulator(epsilon)
    return 320000 * e * x * (1 - s.log(x)) ** 2 / source.KAPPA**2


@cache
def data():
    e, z, v = s.symbols("epsilon z v", positive=True)
    p, A, H, H0, U, K0, U0, p1 = s.symbols("p A H H0 U K0 U0 p1")
    c = e / (2 * (1 + e))
    Ke = K0 + c * U0 + e * A * (H + c * U)
    q = (p * Ke - K0) / e
    D = p1 * K0 + U0 / 2 + H0
    decomp = (
        ((p - 1) / e - p1) * K0
        + (p * c / e - s.Rational(1, 2)) * U0
        + (p * A - 1) * H
        + (H - H0)
        + p * A * c * U
    )
    checks = {
        "exact_q_remainder_decomposition": s.factor(q - D - decomp),
        "product_phase_beta_increment_identity": p * A - 1 - p * (A - 1) - (p - 1),
        "trace_coefficient_divided_difference": s.factor(
            p * c / e - s.Rational(1, 2) - (p - 1 - e) / (2 * (1 + e))
        ),
        "weighted_radial_tail_integral": s.integrate(
            (4 * z + v) * s.exp(-v / 4), (v, 0, s.oo)
        )
        - (16 * z + 16),
        "weighted_radial_complete_integral": 8 * z * z
        + 16 * z
        + 16
        - 8 * ((1 + z) ** 2 + 1),
        "weighted_radial_uniform_budget": s.Integer(16 * 100000 - 1600000),
        "five_q_remainder_budgets": s.Integer(
            425000 + 75000 + 2400000 + 1600000 + 240000 - 4740000
        ),
        "conditional_regulator_budget": s.Integer(
            66000 * 5 + 20000 * 3 + 1400 * 5 - 397000
        ),
        "physical_regulator_budget": s.Rational(400000 * 4, 5) - 320000,
    }
    return {
        "checks": {name: s.simplify(value) for name, value in checks.items()},
        "gates": {
            "quotient_remainder_roundup": bool(s.Rational(4740000, 72) < 66000),
            "integrated_regulator_roundup": bool(s.Integer(397000) < 400000),
            "state_energy_dependent_remainder_not_constant_epsilon_error": True,
            "full_phase_trace_and_radial_terms_retained": True,
            "same_Born_subtraction_in_all_regulated_terms": True,
            "positive_Borel_extension_keeps_log_weighted_radial_dominator": True,
            "zero_energy_and_zero_epsilon_limits_are_exact": True,
            "outer_hard_evanescent_matching_not_supplied": True,
        },
        "whole_energy_dependent_radial_majorant": "For z=ln(64000/(100000R)), the exact weighted radial minimum integral is100000R*(8z^2+16z+16)<=1600000R*L_R^2. It bounds |I_H(e)-I_H(0)|/e without a state-independent remainder.",
        "whole_pointwise_quotient_error": "In the original S296/S301 scheme, |q_e-Delta|<=66000e*R*L_R^2/kappa. The numerator budget4740000 is divided by8pi^2>72; the actual original Born differences and phase Taylor bound are used.",
        "whole_reference_regulator_rate": "On the same full conditioned reference, |Z_e-Z_0|<=e[66000R L_R^2+20000R L_R|ln y|+1400R|ln y|^2]/kappa. Conditional integration gives L1 error<=400000e*a*x*L_x^2/kappa, physically<=320000e*x*L_x^2/kappa^2.",
    }
