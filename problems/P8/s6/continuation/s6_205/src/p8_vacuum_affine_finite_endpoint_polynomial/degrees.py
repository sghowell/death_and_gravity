"""Individual complete endpoint/time-jet rows and finite Taylor-degree cells."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference
from p8_vacuum_affine_reference_spatial_analyticity import domain
from p8_vacuum_affine_retarded_reference_boundary import majorants


def require_cell(endpoint, spatial):
    for value in (endpoint, spatial):
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, s.Integer))
            or not 0 <= value <= 4
        ):
            raise ValueError("Require exact integer endpoint and spatial orders0..4")
    return int(endpoint), int(spatial)


def finite_cells():
    return tuple((j, n) for j in range(5) for n in range(5) if j + n >= 5)


def ultraviolet_cells():
    return tuple((j, n) for j in range(5) for n in range(5) if j + n <= 4)


def row_majorant(endpoint, spatial):
    j, n = require_cell(endpoint, spatial)
    return (
        9
        * reference.PAIR_REF**2
        * sum(majorants.recurrence()[j])
        / domain.INNER**j
        / domain.DELTA**n
    )


def jet_majorant(endpoint, spatial, time_jet):
    j, n = require_cell(endpoint, spatial)
    if (
        isinstance(time_jet, bool)
        or not isinstance(time_jet, (int, s.Integer))
        or not 0 <= time_jet <= j
    ):
        raise ValueError("Require an exact source time-jet order0..endpoint")
    return (
        9
        * reference.PAIR_REF**2
        * majorants.recurrence()[j][int(time_jet)]
        / domain.INNER ** (j - time_jet)
        / domain.DELTA**n
    )


@cache
def data():
    finite = finite_cells()
    uv = ultraviolet_cells()
    checks = {
        "complete_twenty_five_cell_partition": len(finite) + len(uv) - 25,
        "no_overlap_between_finite_and_candidate_UV": len(set(finite) & set(uv)),
        "finite_cell_row_counts": s.Matrix(
            [sum(j == row for j, n in finite) for row in range(5)]
        )
        - s.Matrix([0, 1, 2, 3, 4]),
        "candidate_UV_cell_row_counts": s.Matrix(
            [sum(j == row for j, n in uv) for row in range(5)]
        )
        - s.Matrix([5, 4, 3, 2, 1]),
        "all_forty_finite_time_jet_entries": sum(j + 1 for j, n in finite) - 40,
        "all_thirty_five_candidate_UV_time_jet_entries": sum(j + 1 for j, n in uv) - 35,
        "all_seventy_five_endpoint_Taylor_time_jet_entries": sum(
            j + 1 for j in range(5) for n in range(5)
        )
        - 75,
        "finite_terms_have_no_zero_spatial_degree": sum(n == 0 for j, n in finite),
    }
    return {
        "individual_row": "For endpoint j0..4 and source time jet r<=j, the complete joint-complex pre-current coefficient is bounded by9 Cref^2 c_jr Rtime^(-(j-r)) nu^(1-j). The full nine-pair product rule and all time jets remain.",
        "spatial_Taylor": "The nth homogeneous spatial Taylor term, n0..4, is bounded by9 Cref^2 c_jr Rtime^(-(j-r)) delta^-n |P|^n nu^(1-j-n). Cauchy acts on pre-current coefficients, before the imaginary part, not on tests, Borel preparation or a sharp regulator domain.",
        "complete_row_bound": "Summing source time jets and using Rtime<1 gives D_jn=9 Cref^2 sum_r c_jr Rtime^-j delta^-n. No test derivative or polarization is omitted.",
        "finite_cells": finite,
        "candidate_UV_cells": uv,
        "finite_degree_test": "The radial measure gives large-k power k^(3-j-n). Exactly the ten cells j+n>=5 have an absolutely integrable majorant, with radial exponent s=j+n-1 in4..7. The fifteen others are only candidate UV cells; cancellations may still improve them, and no actual divergence coefficient is asserted.",
        "complete_finite_row_majorants": {
            f"{j}_{n}": row_majorant(j, n) for j, n in finite
        },
        "zero_transfer": "Every finite cell has n>=1, so its continuum polynomial vanishes at P=0. This restoration does not change the already established homogeneous zero-transfer result.",
        "checks": checks,
        "gates": {
            "all_ten_finite_cells": len(finite) == 10,
            "all_fifteen_candidate_cells": len(uv) == 15,
            "finite_radial_exponents": {j + n - 1 for j, n in finite} == {4, 5, 6, 7},
            "complete_jet_rows_bounded": all(
                sum(jet_majorant(j, n, r) for r in range(j + 1)) <= row_majorant(j, n)
                for j in range(5)
                for n in range(5)
            ),
            "same_joint_spatial_radius": domain.DELTA == s.Rational(1, 10**6),
            "same_complete_time_radius": domain.INNER == s.Rational(1, 20000),
            "finite_restoration_vanishes_at_zero_transfer": all(
                n > 0 for j, n in finite
            ),
        },
    }
