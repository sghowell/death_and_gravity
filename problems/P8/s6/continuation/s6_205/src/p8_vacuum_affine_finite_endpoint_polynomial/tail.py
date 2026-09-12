"""Original both-created-mode tail for all ten finite Taylor cells."""

from functools import cache

import sympy as s
from p8_vacuum_affine_reference_spatial_analyticity import domain
from p8_vector_state.comparison import AMAX

from . import degrees, radials

TAIL = s.Integer(10) ** 55


def low_transfer_cell(endpoint, spatial):
    j, n = degrees.require_cell(endpoint, spatial)
    if (j, n) not in degrees.finite_cells():
        raise ValueError("Only absolutely integrable Taylor cells have this tail")
    order = j + n - 1
    return (
        degrees.row_majorant(j, n)
        * AMAX**order
        * 2 ** (order - 3)
        * domain.MASS ** (4 - order)
        / (18 * (order - 3))
    )


@cache
def constants():
    low = sum(low_transfer_cell(j, n) for j, n in degrees.finite_cells())
    high = 2 * radials.constants()["complete_finite_polynomial_upper"]
    return {
        "complete_low_external_transfer_tail_upper": low,
        "complete_high_external_transfer_tail_upper": high,
        "full_finite_polynomial_removed_union_tail_upper": low + high,
    }


@cache
def data():
    k, K = s.symbols("k K", positive=True)
    checks = {}
    for order in range(4, 8):
        checks[f"complete_half_band_radial_tail_order_{order}"] = s.integrate(
            k ** (2 - order), (k, K / 2, s.oo)
        ) - 2 ** (order - 3) * K ** (3 - order) / (order - 3)
    checks["largest_extra_external_tail_degree"] = (
        max(n + 1 for j, n in degrees.finite_cells()) - 5
    )
    return {
        "same_regulator": "Each finite Taylor integrand retains chi_K(P,k)=1_(|k|<=K,|-k+P|<=K) and the high band|k|>=m. The finite-K integrated result is not assumed polynomial or local.",
        "low_transfer": "For |P|<=K/2, a removed pair forces |k|>K/2. For s=j+n-1, integrate nu^-s<=Amax^s |k|^-s. At K>=m, the tail is at most Amax^s 2^(s-3)m^(4-s)/(18(s-3)K), using pi^2>9.",
        "high_transfer": "For |P|>K/2, the full moment bound is multiplied by2|P|/K. The resulting spatial degree n+1 is at most5. The same X46 norm covers both regimes, with all source time jets retained.",
        "full_tail": "Adding every finite cell and both transfer regimes gives |P_fin-P_fin,K|<1e55||D||L2 X46[Gamma]/K. No sharp domain derivative, one-leg replacement or physical cutoff is introduced.",
        "continuum_locality_boundary": "Only after absolute dominated removal can one integrate the finite coefficient tensors on the fixed P-independent high band and call the result a finite spatial polynomial. This argument does not extend to the fifteen candidate UV cells.",
        "constants": constants(),
        "checks": checks,
        "gates": {
            "complete_low_and_high_tail_display": constants()[
                "full_finite_polynomial_removed_union_tail_upper"
            ]
            < TAIL,
            "every_low_transfer_cell_positive": all(
                low_transfer_cell(j, n) > 0 for j, n in degrees.finite_cells()
            ),
            "same_six_spatial_derivatives_suffice": max(
                n + 1 for j, n in degrees.finite_cells()
            )
            <= 6,
            "all_four_radial_tail_orders_retained": {
                j + n - 1 for j, n in degrees.finite_cells()
            }
            == {4, 5, 6, 7},
        },
    }
