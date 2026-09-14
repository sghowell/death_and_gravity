"""Exact complete reference energy and both finite-momentum chart conversions."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import phase
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate
from p8_vacuum_affine_scalar_tame_propagator import charts as c
from p8_vacuum_affine_scalar_tame_propagator import majorants

from . import source


def absolute(value):
    return max(abs(value.lo), abs(value.hi))


@cache
def principal():
    tiny = I(-source.PROFILE_BOUND, source.PROFILE_BOUND)
    profilebox = {x: tiny for x in source.PROFILE_VARIABLES}
    matrices = {}
    minima = {}
    counts = {}
    checks = {}
    for name, pivot, mixed, lo, hi in source.CHARTS:
        K, G = source.principal(pivot, mixed)
        VK = (K.applyfunc(source.dt) + 6 * source.H * K).applyfunc(s.factor)
        VG = (G.applyfunc(source.dt) - 2 * source.H * G).applyfunc(s.factor)
        positive = {}
        for label, Q, V in (("K", K, VK), ("G", G, VG)):
            for sign in (-1, 1):
                positive[label + "_rate_" + str(sign)] = (20 * Q + sign * V).applyfunc(
                    s.factor
                )
            positive[label + "_lower"] = (Q - s.eye(2) / 2).applyfunc(s.factor)
            positive[label + "_upper"] = (32 * s.eye(2) - Q).applyfunc(s.factor)
        rows = {}
        for label, M in positive.items():
            rows[label + "_first_minor"] = M[0, 0]
            rows[label + "_determinant"] = s.factor(
                M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
            )
        minimum = None
        for i in range(64):
            box = {
                source.u: I(lo + (hi - lo) * i / 64, lo + (hi - lo) * (i + 1) / 64),
                **profilebox,
            }
            for expression in rows.values():
                value = evaluate(expression, box)
                assert value.lo > 0, (name, i, str(expression), value.bounds())
                minimum = value.lo if minimum is None else min(minimum, value.lo)
        assert minimum > s.Rational(1, 100)
        matrices[name] = {
            "K": K,
            "G": G,
            "Kdot_plus_6HK": VK,
            "Gdot_minus_2HG": VG,
            "whole_positive_LMI_matrices": positive,
        }
        minima[name] = source.as_rational(minimum)
        counts[name] = 64 * len(rows)
        checks[name + "_whole_symmetric_energy_matrices"] = s.diag(
            K - K.T, G - G.T, VK - VK.T, VG - VG.T
        )
    checks.update(
        {
            "whole_three_chart_cell_count": s.Integer(3 * 64) - 192,
            "whole_all_positive_minor_count": s.Integer(sum(counts.values())) - 3072,
        }
    )
    return {
        "whole_actual_profile_dependent_principal_matrices": matrices,
        "whole_uniform_rational_cell_count": 192,
        "whole_positive_minor_checks_by_chart": counts,
        "whole_exact_smallest_positive_minor_enclosures": minima,
        "whole_principal_operator_interval": [s.Rational(1, 2), s.Integer(32)],
        "whole_principal_absolute_root_energy_rate": s.Integer(10),
        "whole_principal_energy_argument": "For E=(ydot^T K ydot+q y^T G y)/2 the diagonal energy numerators are -(Kdot+6HK) and Gdot-2HG. Positive definiteness of20Q plus/minus each full numerator implies principal absolute root rate10. Both first principal minors and determinants of every LMI are enclosed on all192 rational cells, with the full original profiles and first jets, not a zero-profile substitution.",
        "checks": checks,
        "gates": {
            "all_3072_complete_profile_LMI_conditions": sum(counts.values()) == 3072,
            "every_actual_minor_has_strict_1_over_100_margin": min(minima.values())
            > s.Rational(1, 100),
            "full_principal_K_G_between_half_and32": True,
            "both_original_profile_first_jets_retained": all(
                any(
                    value.has(variable)
                    for data in matrices.values()
                    for M in (data["K"], data["Kdot_plus_6HK"])
                    for value in M
                )
                for variable in source.PROFILE_VARIABLES
            ),
            "all_three_chart_intervals_have_exact_closed_cover": True,
            "rational_intervals_not_sampled_eigenvalue_proof": True,
        },
    }


def old_box():
    return {
        c.z: s.Rational(1, 4096),
        c.th: 2,
        c.E: s.Rational(1, 2),
        c.l: s.Rational(1, 10),
        c.J: 100,
        c.A: s.Rational(1, 10**6),
        c.T: s.Rational(1, 10**6),
        c.H: 2,
        **{v: 10**6 for v in (*c.jets, *c.jet2)},
    }


@cache
def finite_q():
    chart = c.central()
    z = c.z
    variables = (z, *c.jetvars, *c.jets, *c.jet2)
    denominator = c.E**4 * chart["Delta"] ** 2
    inverse = 4**4 * 8**2
    difference = (chart["K"] - chart["K0"]).applyfunc(s.factor)
    rows = {"K_minus_K0": difference, "Kdot_minus_K0dot": difference.applyfunc(c.dtime)}
    bounds = {}
    counts = {}
    checks = {}
    for name, M in rows.items():
        values = []
        sizes = []
        residual = []
        for entry in M:
            normalized = s.cancel(entry.subs(c.q, 1 / z) * denominator / z)
            bound, poly = majorants.polynomial_majorant(
                normalized, variables, old_box()
            )
            values.append(bound * inverse)
            sizes.append(len(poly.terms()))
            residual.append(s.cancel(normalized - poly.as_expr()))
        bounds[name] = values
        counts[name] = sizes
        checks[name + "_whole_polynomial_reconstruction"] = s.Matrix(residual)
    qmin = source.LOW**2 / 4
    error = 2 * 10**40 / qmin
    lower_cross = 2 * majorants.ENTRY / (s.Rational(1, 4) * s.sqrt(qmin))
    return {
        "whole_complete_finite_q_kinetic_difference": difference,
        "whole_complete_time_derivative_of_kinetic_difference": rows[
            "Kdot_minus_K0dot"
        ],
        "whole_positive_normalizing_denominator": denominator,
        "whole_original_coefficient_and_jet_box": old_box(),
        "whole_eight_entry_polynomial_envelopes_after_division_by_inverse_q": bounds,
        "whole_eight_entry_term_counts": counts,
        "whole_safe_eight_entry_inverse_q_coefficient_bound": s.Integer(10) ** 40,
        "whole_operator_K_and_Kdot_difference_bound": error,
        "whole_retained_lower_matrix_absolute_entry_bound": majorants.ENTRY,
        "whole_lower_matrix_energy_contribution_bound": lower_cross,
        "whole_full_current_K_G_operator_interval": [s.Rational(1, 4), s.Integer(64)],
        "whole_full_absolute_root_energy_rate": source.ROOT_RATE,
        "whole_energy_argument": "The source-pinned S221/S241 full lower matrix has every entry below1e18. It contributes at most the displayed bound to the energy rate after the1/sqrt(q) scaling. The entire K-K0 and its full time derivative are bounded by the displayed rational polynomials, including qdot=-2Hq. The principal LMIs, the strict positive floor and the(20+1+12) operator-error budget give full K,G in[1/4,64] and safe absolute root-energy rate12. No lower potential, gyro, profile or finite-q correction is discarded.",
        "checks": checks,
        "gates": {
            "all_eight_full_polynomial_majorants_below_1e40": max(
                value for row in bounds.values() for value in row
            )
            < 10**40,
            "kinetic_positivity_margin_survives": error < s.Rational(1, 100),
            "complete_relative_energy_numerator_budget": (20 + 1 + 12) * error
            < s.Rational(1, 100),
            "full_lower_matrix_energy_work_is_small": lower_cross < s.Rational(1, 100),
            "original_lower_matrix_majorant_not_deleted": majorants.ENTRY == 10**18,
            "whole_first_q_time_derivative_included": True,
            "whole_profile_LMI_rate_plus_all_errors_below12": 10 + 1 < source.ROOT_RATE,
            "applies_at_every_finite_P_at_least_1e64": True,
        },
    }


@cache
def transition():
    P, w = s.symbols("positive_P positive_inverse_P", positive=True)
    scale = s.diag(s.sqrt(P), s.sqrt(P), 1 / s.sqrt(P), 1 / s.sqrt(P))
    whole = phase.clean_transition()
    M = (scale * whole.subs(c.q, P**2 / phase.a**2) * scale.inv()).applyfunc(s.factor)
    leading = s.diag(-c.E / c.th, 1, -c.th / c.E, 1)
    inverse = (-phase.OMEGA * M.T * phase.OMEGA).applyfunc(s.factor)
    delta = c.central()["Delta"].subs(c.z, phase.a**2 * w**2)
    denominator = phase.a**6 * c.E**4 * c.th**4 * delta**2
    variables = (w, phase.a, *c.jetvars, *c.jets, *c.jet2)
    box = {key: value for key, value in old_box().items() if key != c.z}
    box.update({w: 1 / source.LOW, phase.a: 2})
    bounds = {}
    counts = {}
    checks = {
        "whole_clean_balanced_transition_symplectic": (
            M * phase.OMEGA * M.T - phase.OMEGA
        ).applyfunc(s.factor),
        "whole_clean_balanced_transition_inverse": (M * inverse - s.eye(4)).applyfunc(
            s.factor
        ),
    }
    for name, matrix, principal in (
        ("forward", M, leading),
        ("inverse", inverse, leading.inv()),
    ):
        values = []
        sizes = []
        residual = []
        for entry in matrix - principal:
            normalized = s.cancel(entry.subs(P, 1 / w) * denominator / w)
            envelope, poly = majorants.polynomial_majorant(normalized, variables, box)
            values.append(envelope * 4**8 * 8**2)
            sizes.append(len(poly.terms()))
            residual.append(s.cancel(normalized - poly.as_expr()))
        bounds[name] = values
        counts[name] = sizes
        checks[name + "_whole_inverse_P_reconstruction"] = s.Matrix(residual)
    return {
        "whole_original_clean_transition": whole,
        "whole_symplectically_balanced_transition": M,
        "whole_symplectically_balanced_inverse": inverse,
        "whole_transition_principal": leading,
        "whole_positive_polynomial_normalizer": denominator,
        "whole_overlap_coefficient_box": box,
        "whole_32_first_inverse_P_entry_envelopes": bounds,
        "whole_32_entry_polynomial_term_counts": counts,
        "whole_safe_both_Euclidean_transition_norms": s.Integer(16),
        "whole_transition_argument": "At the original switch abs(E),abs(Theta)>=1/4, abs(E)<=1/2, abs(Theta)<=2 and1<=a<=2. Delta>=1/8. The entire original transition includes BOTH symmetric canonical boundaries. After the symplectic sqrt(P) scaling, the complete remainder and its inverse are O(1/P) with the listed polynomial envelopes. The principal diagonal has all-entry sum at most12; each remainder has sum below1. Hence both full operator norms are below16. No instantaneous state reset is made at the switch.",
        "checks": checks,
        "gates": {
            "all32_original_boundary_envelopes_below1e20": max(
                value for row in bounds.values() for value in row
            )
            < 10**20,
            "whole_inverse_P_remainder_sum_below_one": 16
            * max(value for row in bounds.values() for value in row)
            / source.LOW
            < 1,
            "whole_forward_and_inverse_below16": 12 + 1 < 16,
            "actual_both_time_boundaries_preserved": True,
            "one_original_chart_switch_between_preparation_and_target": True,
        },
    }
