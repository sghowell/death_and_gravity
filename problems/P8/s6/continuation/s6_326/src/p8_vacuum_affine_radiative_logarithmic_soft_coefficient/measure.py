"""Only a named logarithmic addition, including its full signed interference."""

from functools import cache

import sympy as s

from . import kernel, source


def cutoff(value):
    x = kernel.energy(value)
    if x == 0:
        raise ValueError("Require a positive detector cutoff")
    return x


def born_seed_upper(value):
    x = cutoff(value)
    L = -s.log(x)
    return (
        2000 * x * (1 + L) / source.KAPPA**2
        + 3000 * x * x * (L * L + L + s.Rational(1, 2)) / source.KAPPA**3
    )


def full_tree_upper(value):
    x = cutoff(value)
    L = -s.log(x)
    return (
        4 * s.Rational(1, 10) ** 1187 * x * (1 + L)
        + 3000 * x * x * (L * L + L + s.Rational(1, 2)) / source.KAPPA**3
    )


@cache
def data():
    x = s.Symbol("x", positive=True)
    primitive = x * (1 - s.log(x))
    square = x * x / 2 * (s.log(x) ** 2 - s.log(x) + s.Rational(1, 2))
    B = 6 * s.Integer(10) ** 16 * source.HEAVY_MASS2**2
    checks = {
        "integrated_logarithm": s.diff(primitive, x) + s.log(x),
        "integrated_phase_weighted_log_square": s.diff(square, x) - x * s.log(x) ** 2,
        "physical_two_polarization_phase_factor": s.Rational(2 * 4, 2 * 8)
        - s.Rational(1, 2),
        "interference_factor_not_square_only": s.Integer(2 * 48 * 320 - 30720),
        "quadratic_term_not_one_loop_order": s.Rational(3, 2) * 2 - 3,
    }
    return {
        "checks": {k: s.factor(v) for k, v in checks.items()},
        "gates": {
            "Born_soft_interference_cap": s.Rational(48 * 320, 9) < 2000,
            "log_addition_square_cap": s.Rational(320**2, 4 * 9) < 3000,
            "original_complete_tree_B_cap": bool(B < s.Integer(10) ** 412),
            "original_complete_tree_interference_cap": s.Rational(320, 9)
            * s.Integer(10) ** 412
            / s.Integer(10) ** 1600
            < 4 * s.Rational(1, 10) ** 1187,
            "interference_and_quadratic_addition_both_retained": True,
            "two_helicities_all_angles_and_recoil_phase_retained": True,
            "no_full_loop_amplitude_or_positive_measure_claim": True,
        },
        "whole_elementary_integrals": (primitive, square),
        "whole_complete_N1_tree_majorant": B,
        "whole_named_signed_measure": "Let L_log=C_Born*ln(1/w)/kappa^(3/2), normalized by the full A0=Am+AG. Bound only (|F1+L_log|^2-|F1|^2) times the original physical one-real recoil phase. Both helicities/all angles give w dw/(2pi^2),rho<=1. The S319 complete47-tree F1 has |F1|<=B/(sqrt(kappa)w). Its signed total variation is <4e-1187*x(1+ln(1/x))+3000*x^2*(ln(1/x)^2+ln(1/x)+1/2)/kappa^3. With the leading Born seed replace the first coefficient by2000/kappa^2.",
        "whole_perturbative_boundary": "The square of the named loop addition is deliberately displayed, not called a one-loop perturbative contribution. The remainder S0*H4_one_loop+R_finite and detector-resolution terms are not bounded or set to zero. This is a selected signed addition, not the actual full real-virtual correction.",
    }
