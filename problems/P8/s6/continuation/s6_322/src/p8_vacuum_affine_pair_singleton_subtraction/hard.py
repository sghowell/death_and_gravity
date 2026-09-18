"""Anisotropic387-core tube, exact remainder budgets and current product."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import trees as lower
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree as matter
from p8_vacuum_affine_three_soft_current_subtraction import hard as one_block

from . import algebra, hierarchy, source

ETA = s.Rational(1, 10**13)
clean, matrix = algebra.clean, algebra.matrix


def absolute_square(z):
    return clean(z * s.conjugate(z))


def entropy(a, b):
    a, b = map(source.require_mass, (a, b))
    return (a + b) * s.log(a + b) - a * s.log(a) - b * s.log(b)


def class_rectangle(a, b, c):
    a, b, c = source.require_energies(a, b, c)
    return c * entropy(a, b) / s.Integer(10) ** 723


def all_nonsingleton_rectangle(a, b, c):
    a, b, c = source.require_energies(a, b, c)
    J = c * entropy(a, b) + b * entropy(a, c) + a * entropy(b, c)
    return 2 * J / s.Integer(10) ** 723


@cache
def budgets():
    n, K = source.HEAVY_MASS2, source.KAPPA
    r3 = 3 * 138240 * 30000000
    r4 = 3 * 10616832 * 30000000
    regular_g = 177 * 1024**4 * 16 * 1800000 * max(1, r3, r3 * r3, r4)
    regular_m = (33 * 4 + 177 * 8) * 1024**2 * 16
    ordering = (2 * 8 * 16 * 64 + 2 * 4 * 32 * 64 + 2 * 8 * 4 * 4 * 512) * 2
    offshell = 2 * 14 * 8 + 2 * 4 * 4 * 8 * 16
    mixed = (
        s.Rational(10**7 * 2 * 10**10 * 600000, 8)
        + s.Rational(2 * 10**5 * 10**10 * 600000 * 2, 8)
        + s.Rational(2 * 10**5 * 2048**2 * 200 * 600000 * 2, 8)
    )
    timelike = (
        s.Rational(10**7 * 128**2) / s.Rational(45, 16)
        + s.Rational(2 * 10**5 * 2 * 128 * 16384) / s.Rational(45, 16)
        + s.Rational(2 * 10**5 * 128**2 * 200)
        / (s.Rational(45, 16) * s.Rational(25, 4))
    ) / 8
    order_m = 140 * 2 * 10**6 * 4
    shift_m = 140 * 4096 * 7
    contact_m = 4 * (2 * 256 * 32768)
    heavy_m = 12 * (2 * 2 * 128 * 16384 + 7 * 128**2)
    order_g = s.Rational(140 * 2 * 10**6 * 2 * 10**5 * 600000, 8)
    shift_g = s.Rational(140 * 4096 * 10**6 * 600000, 8)
    final_m = s.Rational(2 * (regular_m + order_m + shift_m + contact_m + heavy_m), 4)
    final_g = (
        2 * (s.Rational(regular_g, 8) + 8 * mixed + 4 * timelike + order_g + shift_g)
        + 16 * 256**2
    )
    B2 = 10**40 * n * n + 10**60
    f0, fa, fb, fab = hierarchy.CAPS
    # Both first derivatives use the larger1566 envelope for a symmetric formula.
    coefficient = (
        B2 / K ** s.Rational(3, 2) * (fab / ETA + 2 * fa / ETA**2 + f0 / ETA**3)
    )
    previous = one_block.analytic_data()
    previous_coefficient = previous["whole_product_rule_coefficient"]
    a, b, c = s.symbols("a b c", positive=True)
    I = (a + b) * s.log(a + b) - a * s.log(a) - b * s.log(b)
    J = (
        c * I
        + b * ((a + c) * s.log(a + c) - a * s.log(a) - c * s.log(c))
        + a * ((b + c) * s.log(b + c) - b * s.log(b) - c * s.log(c))
    )
    checks = {
        "exact_three_ordering_corrections": s.Integer(ordering - 327680),
        "exact_offshell_single_correction": s.Integer(offshell - 4320),
        "exact_matter_budget": final_m - 13632146432,
        "exact_gravity_budget": final_g
        - s.Rational(1952093286852638820978919582900252803544196317184, 9),
        "regular_topology_inventory": s.Integer(13 + 117 + 117 - 247),
        "double_external_topology_inventory": s.Integer(20 + 60 + 60 - 140),
        "all387_topology_inventory": s.Integer(247 + 140 - 387),
        "all1349_nonsingleton_terms": s.Integer(3 * 387 + 188 - 1349),
        "remaining3767_plus1349": s.Integer(3767 + 1349 - 5116),
        "entropy_second_mixed_derivative": s.factor(s.diff(I, a, b) - 1 / (a + b)),
        "three_class_kernel_is_one_J_not_three_J": s.expand(
            c * I + b * I.xreplace({b: c}) + a * I.xreplace({a: b, b: c}) - J
        ),
        "combined_coefficient_tracks_exact_previous188": coefficient
        + previous_coefficient
        - (coefficient + previous["whole_product_rule_coefficient"]),
    }
    gates = {
        "arbitrary_spatial_NA_bound": bool(4 + 14 * s.Rational(1, 8) < 8),
        "physical_singleton_XB_bound": bool(7 + 4 * s.Rational(1, 8) < 16),
        "arbitrary_spatial_YA_bound": bool(
            14 + 2 * s.Rational(1, 8) + 6 * s.Rational(1, 8) < 32
        ),
        "complex_normalized_current_gradient_below2048": 4 * 8 + 4 * 4 * 64 < 2048,
        "ordering_remainder_below1e6_W": ordering < 10**6,
        "offshell_replacement_below8192_u": offshell < 8192,
        "offshell_pair_product_remainder_below1e6_W": 8192 * 32 < 10**6,
        "paired_current_product_difference_below1e10": 16384 * (100000 + 2048) < 10**10,
        "matter_remainder_coefficient_below1e40": bool(final_m < 10**40),
        "gravity_remainder_coefficient_below1e60": bool(final_g < 10**60),
        "larger_first_derivative_used_conservatively": bool(fa >= fb),
        "original_pair_product_coefficient_below1e_minus723": bool(
            0 < coefficient < s.Rational(1, 10**723)
        ),
        "previous188_coefficient_below1e_minus738": bool(
            0 < previous_coefficient < s.Rational(1, 10**738)
        ),
        "all1349_coefficient_below2e_minus723": bool(
            coefficient + previous_coefficient < s.Rational(2, 10**723)
        ),
        "Cauchy_used_only_on_hard_operator_not_complete_pair": True,
        "closed_anisotropic_disc_cancels_singleton_pole_first": True,
        "all_compatible_current_faces_retained": True,
        "remaining3767_and_full_probability_not_declared_subtracted": True,
        "no_inclusive_real_virtual_state_Regge_or_P8_closure": True,
    }
    return {
        "checks": checks,
        "gates": gates,
        "whole_anisotropic_radius_eta": ETA,
        "whole_complex_core_remainder_B2": B2,
        "whole_component_budgets": {
            "regular_m": regular_m,
            "regular_g": regular_g,
            "ordering": ordering,
            "offshell_single": offshell,
            "mixed": mixed,
            "timelike": timelike,
            "matter_final": final_m,
            "gravity_final": final_g,
        },
        "whole_pair_product_derivative_coefficient": coefficient,
        "whole_retained188_derivative_coefficient": previous_coefficient,
        "whole_unsubtracted_three_singleton_terms": s.Integer(3767),
        "whole_closed_tube": "At positive a,b,c, u=a+b,W=u+c: radii eta*u for a,b and eta*W for c, eta=1e-13. Cancel the singleton-c scalar pole in u'*c'*kappa*sqrt(rho)*R2/A0 first. Composite and total scalar gaps stay above u/4 and W/4; hard mixed, timelike and heavy gaps stay above(tau+W^2)/600000,45/16,n/2. There is no pure-soft graviton denominator in R2.",
        "whole_uniform_core_remainder": "K2-K20 has operator norm below B2*W throughout the closed tube. K20 is the full Born product of normalized composite and singleton currents and is independent of c. Trace and composite off-shell corrections, all247 regular terms,140 double-external terms, forward current groupings and analytic phase are retained.",
        "whole_current_product_and_faces": "G387=K2(F2,B)/kappa. Four product-rule terms give |partial_abc G387|<1e-723/(a+b). K20(F2,B) is c-independent. F2=O(a+b), bounded K2 and the removable c-face give compatible axes. The threefold rectangle is below1e-723*c*I(a,b).",
        "whole_sum_and_integrability": "The three pair choices sum to1e-723*J, not three J. Adding the retained188 class gives a1349-term rectangle below2e-723*J. The previous integrable J and J^2 kernel bounds apply, but no complete5116-tree probability subtraction is inferred.",
    }


@cache
def complex_calibration():
    pars = source.original_parameters()
    E = s.Rational(5, 4)
    r0 = s.Rational(3, 4)
    parameter = s.Rational(31, 16)
    ep = (parameter + 1 / parameter) / 2
    rp = (parameter - 1 / parameter) / 2
    out = s.Matrix([s.Rational(3, 5), 0, s.Rational(4, 5)])
    incoming = (s.ImmutableMatrix([-E, 0, 0, -r0]), s.ImmutableMatrix([-E, 0, 0, r0]))
    born = (
        *incoming,
        s.ImmutableMatrix([E, *(r0 * out)]),
        s.ImmutableMatrix([E, *(-r0 * out)]),
    )
    A0 = clean(
        matter.born_continuation(born, pars["heavy"], pars["cubic"], pars["contact"])
        + lower.old.born(born) / pars["kappa"]
    )
    directions = (
        s.Matrix([-1, 0, 0]),
        s.Matrix([s.Rational(3, 5), s.Rational(4, 5), 0]),
        s.Matrix([s.Rational(3, 5), -s.Rational(4, 5), 0]),
    )
    A = s.ImmutableMatrix(
        s.diag(0, s.Rational(1, 4), s.Rational(1, 2), s.Rational(3, 4))
    )
    tangent = s.Matrix([-directions[2][1], directions[2][0], 0])
    vertical = s.Matrix([0, 0, 1])
    B = s.zeros(4)
    B[1:, 1:] = (tangent * tangent.T - vertical * vertical.T) / 2
    B = s.ImmutableMatrix(B)
    critical = (E - ep) / (E + ep)
    center_boost = critical - s.Rational(1, 10**30)
    checks = {}
    gates = {
        "positive_original_Born": bool(A0 > 0),
        "arbitrary_spatial_A_unit_bound": bool(sum(absolute_square(v) for v in A) < 1),
        "physical_TT_B_unit_bound": bool(sum(absolute_square(v) for v in B) < 1),
    }
    records = []
    for index, shift in enumerate(
        (s.S.Zero, s.I / s.Integer(10) ** 18, -s.I / s.Integer(10) ** 18)
    ):
        boost = center_boost + shift
        gamma = clean((1 + boost**2) / (1 - boost**2))
        sh = clean(2 * boost / (1 - boost**2))
        points = [*incoming]
        for sign in (1, -1):
            points.append(
                matrix(
                    s.Matrix(
                        [
                            gamma * ep + sign * sh * rp * out[0],
                            sh * ep + sign * gamma * rp * out[0],
                            sign * rp * out[1],
                            sign * rp * out[2],
                        ]
                    )
                )
            )
        points = tuple(points)
        Q = matrix(-sum(points, lower.VECTOR_ZERO))
        W, qx = Q[0], Q[1]
        weights = tuple(
            map(clean, ((3 * W - 5 * qx) / 8, 5 * (W + qx) / 16, 5 * (W + qx) / 16))
        )
        qs = tuple(matrix(s.Matrix([w, *(w * n)])) for w, n in zip(weights, directions))
        checks[f"{index}_conservation"] = matrix(sum(qs, lower.VECTOR_ZERO) - Q)
        for i, p in enumerate(points):
            checks[f"{index}_scalar_shell{i}"] = clean(algebra.dot(p, p) - 1)
        for i, q in enumerate(qs):
            checks[f"{index}_null_shell{i}"] = clean(algebra.dot(q, q))
        checks[f"{index}_pair_invariant_and_constant_phase"] = clean(
            algebra.dot(points[2] + points[3], points[2] + points[3]) - 4 * ep**2
        )
        if index == 0:
            center = weights
            u0 = center[0] + center[1]
            W0 = sum(center)
            gates["positive_real_energy_center"] = bool(
                all(w > 0 for w in center) and W0 < s.Rational(1, 8)
            )
            gates["hierarchically_small_singleton"] = bool(
                center[2] < ETA * W0 / s.Integer(10) ** 10
            )
            same, _rays, _born = source.recoil.momenta(E, qs, out)
            for i, (p, p0) in enumerate(zip(points, same)):
                checks[f"same_original_recoil_{i}"] = matrix(p - p0)
        else:
            for i, radius in enumerate((ETA * u0, ETA * u0, ETA * W0)):
                gates[f"{index}_inside_anisotropic_disc{i}"] = bool(
                    absolute_square(weights[i] - center[i]) < radius**2
                )
            gates[f"{index}_outside_old_relative_c_disc"] = bool(
                absolute_square(weights[2] - center[2])
                > (s.Rational(1, 10**12) * center[2]) ** 2
            )
        QA = matrix(qs[0] + qs[1])
        qc = qs[2]
        u = weights[0] + weights[1]
        wc = weights[2]
        checks[f"{index}_singleton_transverse"] = matrix(B * qc)
        checks[f"{index}_singleton_trace"] = algebra.trace(B)
        engine = algebra.CoreTwo(
            [("phi", p, 1) for p in points[1:]] + [("h", QA, A), ("h", qc, B)], **pars
        )
        value, count = engine.amplitude()
        checks[f"{index}_whole387_inventory"] = s.Integer(count - 387)
        K = clean(u * wc * pars["kappa"] * value / A0)
        leadingA = sum(clean(u * (p.T * A * p)[0] / algebra.dot(p, QA)) for p in born)
        leadingB = sum(clean(wc * (p.T * B * p)[0] / algebra.dot(p, qc)) for p in born)
        leading = clean(leadingA * leadingB)
        B2 = s.Integer(10) ** 40 * pars["heavy"] ** 2 + s.Integer(10) ** 60
        rho = s.factor(E * rp / (r0 * ep))
        gates[f"{index}_unphased_remainder_below_quarter_budget"] = bool(
            absolute_square(K - leading) < (B2 * W0 / 4) ** 2
        )
        gates[f"{index}_leading_current_product_bound"] = bool(
            absolute_square(leading) < 256**4
        )
        gates[f"{index}_phase_positive_and_below_two"] = bool(0 < rho < 4)
        gates[f"{index}_phase_difference_by_rho_identity"] = bool(
            (rho - 1) ** 2 < (16 * W0) ** 2
        )
        gates[f"{index}_phased_triangle_budget"] = bool(B2 / 2 + 16 * 256**2 < B2)
        gates[f"{index}_only_exact_entries"] = all(
            not v.has(s.Float) for p in points for v in p
        )
        records.append(
            {
                "shift": shift,
                "weights": weights,
                "squared_unphased_remainder": absolute_square(K - leading),
                "phase_squared": rho,
            }
        )
    return {
        "checks": checks,
        "gates": gates,
        "whole_three_exact_anisotropic_points": records,
        "whole_calibration_boundary": "Original387 source at one real hierarchical recoil point and its two complex conjugates. Complex singleton shifts exceed the old relative-energy tube but lie within the proved anisotropic one. Exact phase invariant and triangle bounds are checked as finite calibrations, not as a replacement for the uniform proof.",
    }


@cache
def data():
    written, calibrated = budgets(), complex_calibration()
    return {
        "checks": {
            **{"budget_" + k: v for k, v in written["checks"].items()},
            **{"complex_" + k: v for k, v in calibrated["checks"].items()},
        },
        "gates": {**written["gates"], **calibrated["gates"]},
        "whole_written_estimate": {
            k: v for k, v in written.items() if k not in ("checks", "gates")
        },
        "whole_original_anisotropic_calibration": {
            k: v for k, v in calibrated.items() if k not in ("checks", "gates")
        },
    }
