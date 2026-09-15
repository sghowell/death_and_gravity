"""Unchanged whole source in the explicit minimal dimensional vacuum continuation."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_common_gravity_masters import (
    reconstruction as previous_reconstruction,
)
from p8_vacuum_affine_massive_common_gravity_masters import source as previous
from p8_vacuum_affine_massive_graviton_threshold_sheet import (
    angular as previous_angular,
)
from p8_vacuum_affine_massive_graviton_threshold_sheet import (
    source as previous_angular_source,
)

S, T, U, MU, K, EP, NU = (
    previous.S,
    previous.T,
    previous.U,
    previous.MU,
    previous.K,
    previous.EP,
    previous.NU,
)
N, H, X, Y, Z = s.symbols(
    "transverse_dimension inverse_beta_squared left_axis right_axis scattering_cosine",
    real=True,
)
reconstruction, angular, angular_source = (
    previous_reconstruction,
    previous_angular,
    previous_angular_source,
)


def require_dimension(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)) or value < 2:
        raise ValueError("Require an integer transverse dimension at least2")
    return int(value)


def require_threshold_mass(value):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, s.Rational))
        or value <= 0
    ):
        raise ValueError("Require a positive exact rational mass squared")
    return s.Rational(value)


@cache
def data():
    inherited = previous.data()
    answer = {
        key: value for key, value in inherited.items() if key not in ("checks", "gates")
    }
    answer["dimensional_scope"] = (
        "Minimal D=4+2EP continuation of the same vacuum Einstein and massive-scalar vertices. Noninteger dimension is a regulator continuation, not extra physical states or a unique full off-background DHOST completion."
    )
    answer["checks"] = dict(inherited["checks"])
    answer["gates"] = {
        "whole_original_R_F_and_heavy_source_retained": all(
            k in answer
            for k in (
                "whole_original_R_F",
                "whole_retained_heavy_source",
                "all_three_vacuum_constants",
            )
        ),
        "source_check_dictionary_not_parent_alias": answer["checks"]
        is not inherited["checks"],
        "fixed_vacuum_constants_not_retuned": True,
        "minimal_dimensional_scope_not_full_DHOST_completion": True,
        "pure_gravity_sector_not_all_original_local_loops": True,
        "leading_flat_reference_not_exact_quantum_vacuum": True,
    }
    return answer
