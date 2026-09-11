"""Explicit full-spacetime fourth-Sobolev norm of the fixed local response."""

from functools import cache

import sympy as s
from p8_vacuum_affine_proca_gaussian.bridge import KAPPA, MASS

from .operator import R0, A, H, t


@cache
def data():
    h = s.Rational(8, 5)
    hp = s.Integer(4)
    hpp = s.Integer(12)
    Amax = s.Integer(500000)
    Aprime = s.Integer(4)
    laplacian_count = s.Integer(3)
    bilaplacian_count = s.Integer(15)
    fourth = (
        1
        + 6 * h
        + (4 * hp + 11 * h * h)
        + (hpp + 7 * h * hp + 6 * h**3)
        + 4
        + 4 * h
        + 4
    )
    wave = (1 + 3 * h + 2) * Amax + Aprime
    total = (wave + s.Rational(160, 60)) / 72
    checks = {
        "actual_H_absolute_bound_factor": s.factor(
            s.Rational(8, 5) - H - 4 * (1 - 2 * t) * (2 - t) / (5 * (1 + t * t))
        ),
        "actual_Hprime_absolute_formula": s.factor(
            s.diff(H, t) - 4 * (1 - t * t) / (1 + t * t) ** 2
        ),
        "actual_Hsecond_absolute_formula": s.factor(
            s.diff(H, t, 2) + 8 * t * (3 - t * t) / (1 + t * t) ** 3
        ),
        "actual_Rprime_absolute_formula": s.factor(
            s.diff(R0, t) - 48 * t * (5 - 7 * t * t) / (1 + t * t) ** 3
        ),
        "actual_Aprime_control": s.factor(s.diff(A, t) + s.diff(R0, t) / 36),
        "full_fourth_order_coefficient_sum": fourth - s.Rational(18817, 125),
        "full_A_wave_coefficient_sum": wave - 3900004,
        "complete_laplacian_multiindex_count": sum(1**2 for _ in range(3))
        - laplacian_count,
        "complete_bilaplacian_multiindex_count": 3 + 3 * 2**2 - bilaplacian_count,
        "canonical_local_operator_display": s.Rational(60000, KAPPA)
        - 6 * s.Rational(1, 10**796),
    }
    return {
        "domain": "Actual CD slab [-1/2,1/2] and arbitrary compact spatial/time TT test; fixed mass1000 and kappa1e800",
        "norm": "||h||H4^2=sum_(j+|alpha|<=4)||partial_t^j partial_x^alpha h||L2(dt dx;F)^2, with every spatial multiindex counted once",
        "background_bounds": {
            "abs_H": h,
            "abs_Hprime": hp,
            "abs_Hsecond": hpp,
            "abs_R0_upper": s.Integer(72),
            "abs_Rprime_upper": s.Integer(120),
            "abs_A": Amax,
            "abs_Aprime": Aprime,
        },
        "laplacian_norm": "||Delta partial_t^j h||<=sqrt(3)||h||H4<2||h||H4 for j<=2; ||Delta^2 h||<=sqrt(15)||h||H4<4||h||H4",
        "full_DdagD_actual_coefficient_sum": fourth,
        "full_DdagD_display": s.Integer(160),
        "full_A_wave_coefficient_sum": wave,
        "local_response_upper_before_inverse_kappa": total,
        "complete_canonical_local_response_display": "||Qloc h||L2<6e-796 ||h||H4",
        "scope": "A bound for the complete specified finite-local tensor Hessian with full spatial derivatives. It is not the nonlocal/state-dependent response, a finite-coupling feedback bound, a cutoff or a pole-removal prescription.",
        "checks": checks,
        "gates": {
            "positive_A_from_actual_mass_and_curvature": s.Rational(5, 12) * MASS**2
            - s.Rational(72, 36)
            > 0,
            "actual_A_below_display": s.Rational(5, 12) * MASS**2 < Amax,
            "actual_R0_upper_enclosure": 24 * (1 + 7 * s.Rational(1, 4)) < 72,
            "actual_Rprime_upper_enclosure": 48 * s.Rational(1, 2) * 5 <= 120,
            "actual_Aprime_below_display": s.Rational(120, 36) < Aprime,
            "laplacian_count_below_two_squared": laplacian_count < 2**2,
            "bilaplacian_count_below_four_squared": bilaplacian_count < 4**2,
            "full_fourth_order_coefficient_below_display": fourth < 160,
            "pi_squared_lower": s.pi**2 > 9,
            "full_operator_coefficient_below_sixty_thousand": 0 < total < 60000,
            "same_fixed_hierarchy": KAPPA == 10**800,
        },
    }
