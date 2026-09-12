"""Exact finite radial moments and absolute finite-polynomial coefficient sum."""

from functools import cache

import sympy as s
from p8_vacuum_affine_reference_spatial_analyticity import domain
from p8_vector_state.comparison import AMAX

from . import degrees

FINITE = s.Integer(10) ** 50


def require_order(order):
    if (
        isinstance(order, bool)
        or not isinstance(order, (int, s.Integer))
        or order not in (4, 5, 6, 7)
    ):
        raise ValueError("Only exact integer integrable radial orders4..7")
    return int(order)


def moment(order):
    return _moment(require_order(order))


@cache
def _moment(order):
    m = domain.MASS
    A = AMAX
    return {
        4: A**3 / (8 * s.pi * m),
        5: A**3 / (6 * s.pi**2 * m * m),
        6: A**3 / (32 * s.pi * m**3),
        7: A**3 / (15 * s.pi**2 * m**4),
    }[order]


def upper(order):
    order = require_order(order)
    m = domain.MASS
    A = AMAX
    return {
        4: A**3 / (24 * m),
        5: A**3 / (54 * m * m),
        6: A**3 / (96 * m**3),
        7: A**3 / (135 * m**4),
    }[order]


def finite_cell_bound(endpoint, spatial):
    j, n = degrees.require_cell(endpoint, spatial)
    if (j, n) not in degrees.finite_cells():
        raise ValueError("This cell still belongs to the candidate UV sector")
    return degrees.row_majorant(j, n) * upper(j + n - 1)


@cache
def constants():
    return {
        "complete_finite_polynomial_upper": sum(
            finite_cell_bound(j, n) for j, n in degrees.finite_cells()
        ),
        "all_ten_finite_cell_bounds": {
            f"{j}_{n}": finite_cell_bound(j, n) for j, n in degrees.finite_cells()
        },
    }


@cache
def data():
    checks = {}
    m = domain.MASS
    A = AMAX
    for order in range(4, 8):
        beta = (
            A**3
            * m ** (3 - order)
            * s.gamma(s.Rational(order - 3, 2))
            / (8 * s.pi ** s.Rational(3, 2) * s.gamma(s.Rational(order, 2)))
        )
        checks[f"complete_massive_radial_integral_order_{order}"] = s.simplify(
            beta - moment(order)
        )
    checks["all_finite_coefficients_restored"] = (
        len(constants()["all_ten_finite_cell_bounds"]) - 10
    )
    return {
        "moment": "The full comparison majorant integral J_s=integral d^3k/(2pi)^3 nu^-s is Amax^3 m^(3-s) Gamma((s-3)/2)/(8pi^(3/2)Gamma(s/2)). All four exact moments s4..7 and rational upper bounds are retained.",
        "high_band": "The actual finite coefficients are integrated over the fixed |k|>=m band. Extending a positive majorant to all k gives the displayed J_s upper bounds; it does not replace the original physical/reference definition or remove a low-mode term.",
        "finite_polynomial": "After dominated cutoff removal, every retained coefficient is integrated on a P-independent fixed band. The ten-cell sum is therefore a genuine finite spatial polynomial of degree<=4, with source time jets<=4 and background-dependent coefficient tensors. This statement applies only to the absolutely convergent cells.",
        "bound": "The full sum is below1e50||D||L2 X46[Gamma]. It must be restored to the known response because it was subtracted in S203; it is not deleted or identified with a freely adjustable covariant counterterm.",
        "constants": constants(),
        "checks": checks,
        "gates": {
            "complete_finite_polynomial_display": constants()[
                "complete_finite_polynomial_upper"
            ]
            < FINITE,
            "all_moments_positive": all(moment(order) > 0 for order in range(4, 8)),
            "candidate_UV_not_sent_to_finite_integrals": all(
                j + n - 1 > 3 for j, n in degrees.finite_cells()
            ),
            "no_endpoint_zero_cell_declared_finite": all(
                j > 0 for j, n in degrees.finite_cells()
            ),
        },
    }
