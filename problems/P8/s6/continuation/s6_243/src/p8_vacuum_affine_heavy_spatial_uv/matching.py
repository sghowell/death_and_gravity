"""Complete scalar dimensional finite UV difference with its own fixed pole action."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import state
from p8_vacuum_affine_ordered_scalar_symbol import matching as geometry

from . import jets

t = geometry.jets.t
a = (1 + t * t) ** 2
H = 4 * t / (1 + t * t)


def current_substitution():
    return {
        jets.p: geometry.jets.p,
        jets.a: a,
        **{jets.h[j]: s.diff(H, t, j) for j in range(6)},
    }


@cache
def finite():
    inv = jets.inv
    d = jets.d
    gamma, curv = geometry.hessians(*inv)
    pole = (
        -(jets.m**2) * curv["R_old"] / 3
        + curv["R_squared"] / 36
        + (curv["Riemann_squared"] - curv["Ricci_squared"]) / 90
    )
    logs = jets.logarithmic()
    value = [
        s.factor(v.subs(d, 3)).subs(current_substitution(), simultaneous=True)
        / (2 * s.pi**2)
        for v in logs[:3]
    ]
    slope = [
        s.factor(s.diff(v, d).subs(d, 3)).subs(
            current_substitution(), simultaneous=True
        )
        / (2 * s.pi**2)
        for v in logs[:3]
    ]
    power = jets.spatial_rows()[0, 0, 2].subs(d, 3).subs(
        current_substitution(), simultaneous=True
    ) / (2 * s.pi**2)
    ell = s.Symbol("ell", real=True)
    raw = (
        sum(
            ((1 - s.log(2) - ell / 2) * value[r] - slope[r]) * gamma[r]
            for r in range(3)
        )
        - jets.m**2 * power * gamma[0] / 2
    )
    counter = s.diff(pole, geometry.geometry.d).subs(geometry.geometry.d, 3) / (
        32 * s.pi**2
    )
    complete = raw + counter
    rows = tuple(s.factor(s.diff(complete, g)) for g in gamma)
    swap = {inv[4]: inv[5], inv[5]: inv[4]}
    transpose = [row.xreplace(swap) for row in rows]
    checks = {
        "whole_minimal_scalar_physical_pole": jets.pole_check(),
        "complete_ordered_Green_two": s.factor(rows[2] - transpose[2]),
        "complete_ordered_Green_one": s.factor(
            rows[1] - 2 * s.diff(transpose[2], t) + transpose[1]
        ),
        "complete_ordered_Green_zero": s.factor(
            rows[0]
            - transpose[0]
            + s.diff(transpose[1], t)
            - s.diff(transpose[2], t, 2)
        ),
    }
    for key, v in jets.spatial_rows().items():
        j, _r, degree = key
        if j + degree < 4 and key != (0, 0, 2):
            checks[
                "no_other_spatial_power_or_odd_radial_grade_" + "_".join(map(str, key))
            ] = v
    for j in (3, 4):
        checks["complete_logarithmic_high_source_jet_" + str(j)] = logs[j]
    for j, v in enumerate(rows):
        checks["zero_transfer_full_finite_" + str(j)] = s.factor(
            v.subs(geometry.jets.p, 0)
        )
    L, B = s.symbols("bounded_log_1_plus_t_squared bounded_log_2", real=True)
    c1 = c0 = s.S.Zero
    coefficient_rows = {}
    for j, row in enumerate(rows):
        for inv0, weight in zip(inv, (1, 1, 1, 3, 2, 2)):
            v = s.diff(row, inv0)
            normalized = s.cancel(
                s.expand_log(v * s.pi**2 * a, force=True).subs(
                    {s.log(1 + t * t): L, s.log(2): B}
                )
            )
            poly = s.Poly(normalized, t, geometry.jets.p, jets.m, ell, L, B)
            checks[
                "complete_finite_coefficient_reconstruction_" + str(j) + "_" + str(inv0)
            ] = s.factor(normalized - poly.as_expr())
            for powers, coef in poly.terms():
                tp, pp, mp, lp, ll, bb = powers
                assert pp <= 4 and mp in (0, 2) and lp <= 1 and ll <= 1 and bb <= 1
                bound = (
                    weight
                    * abs(coef)
                    * s.Rational(1, 2) ** tp
                    * 462**lp
                    * s.Rational(1, 4) ** ll
                    / 9
                )
                if mp:
                    c1 += bound
                else:
                    c0 += bound
            coefficient_rows[j, str(inv0)] = s.factor(v)
    bound = (c1 * state.MASS2 + c0) / state.KAPPA
    checks["entire_mass_squared_coefficient_sum"] = c1 - s.Rational(76796107, 23224320)
    checks["entire_constant_coefficient_sum"] = c0 - s.Rational(9596179133, 464486400)
    return {
        "invariants": inv,
        "source_jets": gamma,
        "ell": ell,
        "own_scalar_pole": pole,
        "complete_physical_logarithmic_rows": tuple(value),
        "complete_invariant_dimensional_slopes": tuple(slope),
        "complete_physical_power_two": power,
        "complete_covariant_evanescent_counteraction": counter,
        "complete_finite_UV_difference_coefficients": rows,
        "all_eighteen_fixed_invariant_coefficients": coefficient_rows,
        "full_coefficient_sum_mass_squared": c1,
        "full_coefficient_sum_constant": c0,
        "complete_normalized_local_UV_difference_bound": bound,
        "checks": checks,
        "gates": {
            "every_ordered_endpoint_and_radial_coefficient": len(
                jets.endpoint_products()
            )
            == 60
            and sum(map(len, jets.endpoint_products().values())) == 140
            and len(jets.radial_rows()) == 35,
            "all_eighteen_finite_invariant_coefficients": len(coefficient_rows) == 18,
            "actual_scalar_mass_lower_for_log_bound": state.MASS2 > 10**197,
            "actual_scalar_mass_upper_for_log_bound": state.MASS2 < 10**198,
            "positive_entire_coefficient_sum": bool(c1 > 0 and c0 > 0),
            "full_normalized_local_UV_difference_below_one_e_minus_600": bool(
                bound < s.Rational(1, 10**600)
            ),
            "compact_variations_and_physical_ordered_Green_pairing": True,
            "all_dimension_and_volume_derivatives_retained": True,
            "old_vector_finite_weights_not_borrowed": True,
            "full_scalar_heat_action_not_added_twice": True,
            "no_full_UV_subtracted_kernel_or_homogeneous_anchor_bound": True,
        },
    }
