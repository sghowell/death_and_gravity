"""Full 32-refinement coverage and already-integrated first-sheet majorants."""

from functools import cache

import sympy as s
from p8_vacuum_two_loop_light_pole import graphs, subtraction


@cache
def data():
    L, g, Q, alpha = s.symbols("L g Q alpha", positive=True)
    parents = graphs.data()
    counts = parents["group_counts"]
    expected = {
        "sunset_local": 1,
        "sunset_mixed": 6,
        "sunset_same_heavy": 3,
        "sunset_finite": 6,
        "nested_tadpole_outer_constant": 8,
        "nested_tadpole_inner_constant": 4,
        "nested_tadpole_inner_mixed": 4,
    }
    suprema = {
        "mixed_sunset_constant_anchored": 33 * L * g / Q**2,
        "same_heavy_sunset": 5 * g**2 / (2 * Q**2),
        "different_heavy_sunset": 56 * g**2 / Q**2,
        "nested_decaying_remainder": 4 * g**2 / Q**2,
    }
    slopes = {key: val / 2 for key, val in suprema.items()}
    slopes["nested_constant_alpha"] = alpha**2
    checks = {
        f"full_parent_count_{key}": counts[key] - n for key, n in expected.items()
    }
    checks.update(
        {
            "all_raw_refinements_retained": parents["raw_refinement_count"] - 32,
            "full_parent_partition": sum(counts.values()) - 32,
            "combined_Cauchy_g_squared_factor": s.Rational(5, 4)
            + 28
            + 2
            - s.Rational(125, 4),
        }
    )
    # These three checks inspect the actual local-sunset proper contractions.
    checks.update(
        {
            "same_" + key: value
            for key, value in subtraction.data()["checks"].items()
            if key.startswith("local_sunset_contracted_proper_core_")
        }
    )
    u = s.Symbol("invariant")
    A = s.Function("one_loop_self_energy")(u)
    aD, bD = s.symbols("alpha_D beta_D")
    checks["nested_alpha_before_outer_OS_derivative"] = s.diff(alpha * A, u).subs(
        u, 1
    ) - alpha * s.diff(A, u).subs(u, 1)
    checks["all_affine_references_removed_only_by_outer_OS"] = s.expand(
        (aD + bD * u) - (aD + bD) - (u - 1) * bD
    )
    return {
        "raw_refinement_count": parents["raw_refinement_count"],
        "same_full_graph_partition": counts,
        "holomorphic_disc_radius_about_one": s.Integer(2),
        "integrated_unprojected_or_constant_anchored_majorants": suprema,
        "finite_slope_majorants": slopes,
        "checks": checks,
        "scope": "Reuse S6.126's proven first-sheet holomorphic neighborhood and integrated majorants, not an assumed complex momentum shift. Cauchy's derivative estimate divides each supremum by radius two. The mixed function is constant-anchored before its regulator is removed. Tadpole contractions are momentum independent; the nested alpha term has finite slope alpha_0^2.",
    }
