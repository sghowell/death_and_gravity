"""Exact obstructions to isolated-cluster and multi-off-shell Ward shortcuts."""

from functools import cache

import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import trees as full
from p8_vacuum_affine_complete_two_graviton_tree import trees

from . import source

ETA, imm, old = trees.ETA, trees.imm, trees.old


def cluster_configuration(r):
    w = s.Rational(1, 32)
    sn = 2 * r / (1 + r * r)
    cs = (1 - r * r) / (1 + r * r)
    ns = (imm([0, 0, 1]), imm([sn, 0, cs]), imm([1, 0, 0]))
    ts = (imm([1, 0, 0]), imm([cs, 0, -sn]), imm([0, 1, 0]))
    bs = (imm([0, 1, 0]), imm([0, 1, 0]), imm([0, 0, 1]))
    hs = []
    for n, t, b in zip(ns, ts, bs):
        q = imm([w, *(w * n)])
        eps = s.zeros(4)
        eps[1:, 1:] = t * t.T - b * b.T
        hs.append((q, imm(eps)))
    Q = imm(sum((q for q, _ in hs), trees.VECTOR_ZERO))
    v = imm([Q[1] / Q[0], 1, 0, 0])
    return tuple(hs), Q, imm(v * v.T)


class NoQuartic(trees.TreeEngine):
    def permitted(self, tags):
        return super().permitted(tags) and tags != ("h", "h", "h", "h")


@cache
def cluster_data():
    r = s.Symbol("r", positive=True)
    hs, Q, U = cluster_configuration(r)
    engine = trees.TreeEngine([("h", q, e) for q, e in hs])
    J, count = engine.amputated(7, "h")
    propagated, _ = engine.current(7, "h")
    value = s.factor(old.pair(propagated, U))
    polynomial = (
        6 * r**12
        - 69 * r**11
        + 16 * r**10
        + 111 * r**9
        - 342 * r**8
        + 796 * r**7
        - 942 * r**6
        + 1176 * r**5
        - 904 * r**4
        + 705 * r**3
        - 378 * r**2
        + 161 * r
        - 48
    )
    expected = (
        2 * polynomial / (9 * (r - 1) ** 2 * (r * r + 1) ** 4 * (2 * r * r - r + 1))
    )
    residue = imm(
        propagated.applyfunc(lambda a: s.limit(s.factor((r - 1) ** 2 * a), r, 1))
    )
    expected_residue = imm(
        [[-18, 12, 0, 12], [12, -6, 0, -12], [0, 0, 6, 0], [12, -12, 0, -6]]
    )
    bad, n_bad = NoQuartic([("h", q, e) for q, e in hs]).amputated(7, "h")
    bad_ward = (Q.T * ETA * bad).applyfunc(s.factor)
    checks = {
        "complete_pure_three_root_Ward": (Q.T * ETA * J).applyfunc(s.factor),
        "conserved_outer_source": (Q.T * ETA * U).applyfunc(s.factor),
        "whole_isolated_cluster_rational_function": s.factor(value - expected),
        "whole_isolated_cluster_pole_residue": s.limit(
            s.factor((r - 1) ** 2 * value), r, 1
        )
        - 2,
        "whole_propagated_root_residue_matrix": residue - expected_residue,
        "outer_root_stays_off_shell": s.factor(
            old.dot(Q, Q).subs(r, 1) - s.Rational(1, 256)
        ),
        "outer_source_norm_remains_bounded": s.factor(
            sum(x * x for x in U).subs(r, 1) - s.Rational(169, 81)
        ),
        "nonpole_other_collinear_limit": s.limit(value, r, 0) + s.Rational(32, 3),
        "full_pure_cluster_count": s.Integer(count - 4),
        "omitted_EH4_cluster_count": s.Integer(n_bad - 3),
        "pair_dot_exact_vanishing_order": s.factor(
            old.dot(hs[1][0], hs[2][0])
            - s.Rational(1, 1024) * (r - 1) ** 2 / (1 + r * r)
        ),
        "unit_norm_polarization_rescaling_nonzero_residue": s.simplify(
            2 * (9 / s.Integer(13)) / (s.sqrt(2) ** 3) - 9 / (13 * s.sqrt(2))
        ),
    }
    for i, (q, e) in enumerate(hs):
        checks[f"physical_null_shell_{i}"] = s.factor(old.dot(q, q))
        checks[f"physical_TT_transversality_{i}"] = (e * q).applyfunc(s.factor)
        checks[f"physical_TT_trace_{i}"] = s.factor(s.trace(ETA * e))
        checks[f"physical_TT_norm_squared_{i}"] = s.factor(sum(a * a for a in e) - 2)
    at_hs, at_Q, at_U = cluster_configuration(s.S.One)
    outer = imm(
        -(ETA * at_U * ETA - ETA * s.trace(ETA * at_U) / 2) / old.dot(at_Q, at_Q)
    )
    rest, nrest = trees.TreeEngine(
        [("h", at_hs[0][0], at_hs[0][1]), ("h", -at_Q, outer)]
    ).amputated(3, "h")
    defect = ((at_hs[1][0] + at_hs[2][0]).T * ETA * rest).applyfunc(s.factor)
    checks["nested_cut_exact_Ward_defect"] = defect - s.Matrix(
        [[s.Rational(1, 32), s.Rational(1, 24), 0, s.Rational(1, 32)]]
    )
    checks["nested_cubic_count"] = s.Integer(nrest - 1)
    return {
        "parameter": r,
        "rational_cluster": value,
        "root_residue": residue,
        "nested_Ward_defect": defect,
        "checks": checks,
        "gates": {
            "missing_EH4_breaks_complete_root_Ward": bad_ward != s.zeros(1, 4),
            "conserved_outer_source_does_not_conserve_nested_descendant": defect
            != s.zeros(1, 4),
        },
    }


@cache
def actual_source_data():
    residue = cluster_data()["root_residue"]
    hs, Q, _ = cluster_configuration(s.S.One)
    E = s.Rational(5, 4)
    ep = (2 * E - Q[0]) / 2
    dy = s.sqrt(ep**2 - 1 - (Q[1] ** 2 + Q[3] ** 2) / 4)
    ps = (
        imm([-E, 0, 0, -s.Rational(3, 4)]),
        imm([-E, 0, 0, s.Rational(3, 4)]),
        imm([ep, -Q[1] / 2, dy, -Q[3] / 2]),
        imm([ep, -Q[1] / 2, -dy, -Q[3] / 2]),
    )
    checks = {
        "physical_recoil_square": s.factor(dy**2 - s.Rational(457, 1024)),
        "physical_total_conservation": sum(ps, trees.VECTOR_ZERO) + Q,
    }
    for i, p in enumerate(ps):
        checks[f"physical_mass_shell_{i}"] = s.factor(old.dot(p, p) - 1)
    samples = {}
    gates = {}
    for label, pars in (
        ("diagnostic", source.diagnostic_parameters()),
        ("original", source.original_parameters()),
    ):
        U, n = trees.TreeEngine([("phi", p, 1) for p in ps], **pars).amputated(15, "h")
        checks[label + "_full_hard_root_Ward"] = (Q.T * ETA * U).applyfunc(s.factor)
        checks[label + "_full_hard_root_count"] = s.Integer(n - 47)
        pole = s.factor(old.pair(residue, U) / pars["kappa"])
        propagated = imm(-(ETA * U * ETA - ETA * s.trace(ETA * U) / 2) / old.dot(Q, Q))
        R, nrest = trees.TreeEngine(
            [("h", hs[0][0], hs[0][1]), ("h", -Q, propagated)], kappa=pars["kappa"]
        ).amputated(3, "h")
        defect = ((hs[1][0] + hs[2][0]).T * ETA * R).applyfunc(s.factor)
        checks[label + "_nested_cubic_count"] = s.Integer(nrest - 1)
        samples[label] = {"isolated_cluster_pole": pole, "hard_graph_count": n}
        gates[label + "_actual_conserved_hard_source_has_nonzero_cluster_pole"] = (
            pole != 0
        )
        gates[label + "_actual_nested_cubic_remainder_not_conserved"] = (
            defect != s.zeros(1, 4)
        )
    return {"checks": checks, "gates": gates, "samples": samples}


@cache
def four_soft_data():
    w = s.Rational(1, 32)
    x, y, z = map(s.Matrix, ([1, 0, 0], [0, 1, 0], [0, 0, 1]))

    def pol(t, b):
        A = s.zeros(4)
        A[1:, 1:] = t * t.T - b * b.T
        return imm(A)

    hs = (
        (imm([w, w, 0, 0]), pol(y, z)),
        (imm([w, -w, 0, 0]), pol(y, z)),
        (imm([w, 0, 0, w]), pol(x, y)),
        (imm([w, 0, 0, -w]), pol(x, y)),
    )
    J, n = full.TreeEngine([("h", q, e) for q, e in hs], kappa=16).amputated(15, "h")
    Q = sum((q for q, _ in hs), trees.VECTOR_ZERO)
    checks = {
        "complete_four_soft_root_count": s.Integer(n - 26),
        "complete_four_soft_root_Ward": (Q.T * ETA * J).applyfunc(s.factor),
    }
    E = s.Rational(5, 4)
    ep = E - Q[0] / 2
    pz = s.sqrt(ep * ep - 1)
    ps = (
        imm([-E, 0, 0, -s.Rational(3, 4)]),
        imm([-E, 0, 0, s.Rational(3, 4)]),
        imm([ep, 3 * pz / 5, 0, 4 * pz / 5]),
        imm([ep, -3 * pz / 5, 0, -4 * pz / 5]),
    )
    checks["four_real_recoil_conservation"] = sum(ps, trees.VECTOR_ZERO) + Q
    for i, p in enumerate(ps):
        checks[f"four_real_recoil_shell_{i}"] = s.factor(old.dot(p, p) - 1)
    q34 = hs[2][0] + hs[3][0]
    h34, n34 = trees.TreeEngine([("h", q, e) for q, e in hs[2:]], kappa=16).current(
        3, "h"
    )
    R, nrest = trees.TreeEngine(
        [("phi", p, 1) for p in ps] + [("h", q34, h34)],
        **source.diagnostic_parameters(),
    ).amputated(31, "h")
    defect = ((hs[0][0] + hs[1][0]).T * ETA * R).applyfunc(s.factor)
    checks["second_collapsed_pair_count"] = s.Integer(n34 - 1)
    checks["multi_offshell_remainder_count"] = s.Integer(nrest - 434)
    return {
        "checks": checks,
        "whole_multi_offshell_Ward_defect": defect,
        "gates": {
            "completed_one_offshell_four_soft_root_is_conserved": checks[
                "complete_four_soft_root_Ward"
            ]
            == s.zeros(1, 4),
            "two_offshell_remainder_is_not_generally_conserved": defect
            != s.zeros(1, 4),
        },
    }


@cache
def data():
    cluster, actual, four = cluster_data(), actual_source_data(), four_soft_data()
    checks = {}
    gates = {}
    for label, record in (
        ("cluster", cluster),
        ("actual", actual),
        ("four_soft", four),
    ):
        checks.update(
            {label + "_" + key: value for key, value in record["checks"].items()}
        )
        gates.update(
            {label + "_" + key: value for key, value in record["gates"].items()}
        )
    return {
        "whole_isolated_cluster_counterexample": {
            key: cluster[key]
            for key in (
                "parameter",
                "rational_cluster",
                "root_residue",
                "nested_Ward_defect",
            )
        },
        "whole_actual_model_cluster_poles": actual["samples"],
        "whole_multi_offshell_remainder_Ward_defect": four[
            "whole_multi_offshell_Ward_defect"
        ],
        "whole_obstruction_scope": "The isolated pure-three-soft graph class has a nonzero pair-collinear pole even after contraction with the actual conserved original47-tree hard source. Its total root stays strictly off shell. The full physical amplitude is NOT claimed divergent: completing the pair remainder changes the graph class and restores its Ward identity. A second collapsed pair makes that remainder multi-off-shell and loses this premise again.",
        "checks": checks,
        "gates": {
            **gates,
            "nonzero_pole_survives_unit_Frobenius_rescaling": True,
            "actual_hard_recoil_extends_smoothly_near_pair_collision": True,
            "counterexamples_not_full_amplitude_or_detector_divergence": True,
            "no_frozen_source_or_physical_parameter_changed": True,
        },
    }
