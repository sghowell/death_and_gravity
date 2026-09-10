"""Exact massive polynomial parameter and Euclidean kernel calibrations."""

from fractions import Fraction

import sympy as sp
from p8_exceptional_vacuum import analytic


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational")
    return sp.Rational(value)


def point(lam=analytic.VACUUM_LAMBDA_BAR, gamma=analytic.FIXED_GAMMA):
    lam, gamma = rational(lam), rational(gamma)
    if not lam > gamma > 0:
        raise ValueError(
            "Require positive couplings and heavy mass strictly above two-light threshold"
        )
    D = 2 * lam / gamma
    mu2 = D + 2
    G2 = 2 * lam * D**3
    quartic = G2 * (3 * D - 2) / D**2
    delta = quartic - 3 * G2 / mu2
    return {
        "lambda": lam,
        "gamma": gamma,
        "D": D,
        "heavy_mass_squared": mu2,
        "cubic_squared": G2,
        "quartic": quartic,
        "stable_margin": delta,
        "kernel_upper": (quartic - G2 / mu2) / 2,
        "tree_forward_b2": 4 * lam,
    }


def kernel_point(momentum_squared=0, field=1):
    y, phi = rational(momentum_squared), rational(field)
    if y < 0:
        raise ValueError("Require nonnegative Euclidean momentum squared")
    d = point()
    F = d["kernel_upper"] - d["cubic_squared"] / (y + d["heavy_mass_squared"])
    return {
        "momentum_squared": y,
        "field": phi,
        "F": F,
        "Schur_kernel": y + 1 + F * phi * phi,
        "log_increment": F * phi * phi / (y + 1),
        "rational_one_loop_remainder_upper": d["kernel_upper"] ** 3 * phi**6 / 1728,
    }


def bad_cases():
    invalid = (True, False, 1.0, sp.Float(1), "1", sp.I, sp.oo, -sp.oo, sp.zoo, sp.nan)
    rows = []
    for i, v in enumerate(invalid):
        rows += [
            ("lambda_" + str(i), point, (v,)),
            ("gamma_" + str(i), point, (1, v)),
            ("momentum_" + str(i), kernel_point, (v,)),
            ("field_" + str(i), kernel_point, (0, v)),
        ]
    rows += [
        ("parameter_domain_" + str(i), point, args)
        for i, args in enumerate(((0, 1), (-1, 1), (1, 0), (1, -1), (1, 1), (1, 2)))
    ]
    rows += [("negative_Euclidean_momentum", kernel_point, (-1,))]
    return rows
