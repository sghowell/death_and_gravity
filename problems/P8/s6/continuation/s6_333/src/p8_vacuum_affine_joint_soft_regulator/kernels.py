"""Uniform regulated remaining-energy kernels, including epsilon zero."""

from functools import cache

import sympy as s
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient.measure import cutoff
from p8_vacuum_affine_radiative_soft_state_transfer.continuity import (
    require_radiated_energy,
)


def regulator(value):
    return require_radiated_energy(value)


def remaining_kernels(epsilon, remaining):
    e, y = regulator(epsilon), cutoff(remaining)
    if e == 0:
        return s.S.One, s.log(y)
    g = y ** (2 * e)
    return g, (g - 1) / (2 * e)


@cache
def data():
    e, y, q, da, z = s.symbols("epsilon y q da z", positive=True)
    g = y ** (2 * e)
    h = (g - 1) / (2 * e)
    checks = {
        "same_state_regulated_mark_decomposition": ((da + 2 * e * q) * g - da) / (2 * e)
        - q * g
        - da * h,
        "power_kernel_derivative": s.diff(g, y) - 2 * e * y ** (2 * e - 1),
        "log_kernel_derivative": s.diff(h, y) - y ** (2 * e - 1),
        "log_kernel_zero_regulator": s.limit(h, e, 0) - s.log(y),
        "power_kernel_zero_regulator": s.limit(g, e, 0) - 1,
        "exp_remainder_value_at_zero": (s.exp(-z) - 1 + z).subs(z, 0),
        "exp_remainder_first_derivative_at_zero": s.diff(s.exp(-z) - 1 + z, z).subs(
            z, 0
        ),
        "exp_remainder_second_derivative": s.diff(s.exp(-z) - 1 + z, z, 2) - s.exp(-z),
    }
    return {
        "checks": {name: s.simplify(value) for name, value in checks.items()},
        "gates": {
            "power_between_zero_and_one_in_physical_remaining_domain": True,
            "regulated_log_dominated_by_absolute_log": True,
            "log_kernel_increment_bounded_by_logratio": True,
            "power_kernel_increment_bounded_by_two_epsilon_logratio": True,
            "log_kernel_remainder_bounded_by_epsilon_log_squared": True,
            "epsilon_zero_separately_defined_without_division": True,
        },
        "whole_exact_decomposition": "Z_e=q_e*g_e+deltaa*h_e, q_e=(deltaa_e-deltaa)/(2e), g_e=y^(2e), h_e=(g_e-1)/(2e). At e=0 use q_0=deltaDelta,g_0=1,h_0=ln y.",
        "whole_kernel_bounds": "For0<y<=y+t<=x<=1/8 and0<=e<=1/8: 0<=g_e<=1,|h_e|<=|ln y|,|h_e(y+t)-h_e(y)|<=ln(1+t/y),|g_e(y+t)-g_e(y)|<=2e*ln(1+t/y). Also|g_e-1|<=2e*|ln y| and|h_e-ln y|<=e*|ln y|^2.",
    }
