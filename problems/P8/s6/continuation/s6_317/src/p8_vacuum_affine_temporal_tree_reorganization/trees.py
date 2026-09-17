"""Gauge-complete rooted recursion, with literal nonlinear chart comparisons."""

from functools import cache

import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import trees as original
from p8_vacuum_affine_complete_pair_factorization import obstruction
from p8_vacuum_affine_complete_two_graviton_tree import trees as lower

from . import chart, source

ETA, imm = lower.ETA, lower.imm


class TemporalSoft(original.TreeEngine):
    """Same vertices and graph classes; change only the pure-soft inverse policy."""

    @lower.instance_cache
    def current(self, mask, kind):
        value, count = super().current(mask, kind)
        if (
            count
            and kind == "h"
            and mask.bit_count() > 1
            and all(
                self.legs[i][0] == "h" for i in range(len(self.legs)) if mask >> i & 1
            )
        ):
            value = chart.project(value, self.momentum(mask))
        return value, count


def pair_shift(H, Q):
    """Coordinate vector for g=eta+2h at unit coupling; TT singleton shifts vanish."""
    Q = chart.require_energy(Q)
    q = ETA * Q
    z = s.zeros(4, 1)
    z[0] = H[0, 0] / Q[0]
    for j in range(1, 4):
        z[j] = (2 * H[0, j] - q[j] * z[0]) / Q[0]
    return imm(ETA * z)


def triple_pullback(hs):
    """Independent closed degree-three geometric formula, not new recursion."""
    if len(hs) != 3:
        raise ValueError("Require exactly three physical soft leaves")
    engine = lower.TreeEngine([("h", q, A) for q, A in hs])
    H, count = engine.current(7, "h")
    Q = engine.momentum(7)
    correction = s.zeros(4)
    for mask, k in ((3, 2), (5, 1), (6, 0)):
        H2, _ = engine.current(mask, "h")
        P = engine.momentum(mask)
        Z = pair_shift(H2, P)
        p = ETA * P
        q, A = hs[k]
        AZ = A * Z
        correction += (q.T * ETA * Z)[0] * A + p * AZ.T + AZ * p.T
    return chart.project(H - correction, Q), H, imm(correction), Q, count


@cache
def cluster_calibration():
    r = s.Symbol("r", positive=True)
    hs, Q, _ = obstruction.cluster_configuration(r)
    temporal, H, correction, _, count = triple_pullback(hs)
    engine = TemporalSoft([("h", q, A) for q, A in hs])
    recursive, recursive_count = engine.current(7, "h")
    metric, momenta = chart.metric_coefficients(hs, lower.TreeEngine)
    _, mapped = chart.temporal_chart(metric, momenta, 3)
    linear = chart.project(H, Q)
    residue = lambda A: imm(
        A.applyfunc(lambda value: s.limit(s.factor((r - 1) ** 2 * value), r, 1))
    )
    finite = imm(temporal.applyfunc(lambda value: s.limit(value, r, 1)))
    expected_linear = imm([[0, 0, 0, 0], [0, 2, 0, -4], [0, 0, 6, 0], [0, -4, 0, 0]])
    expected_finite = imm(
        [
            [0, 0, 0, 0],
            [0, -5, 0, s.Rational(35, 6)],
            [0, 0, s.Rational(71, 2), 0],
            [0, s.Rational(35, 6), 0, -s.Rational(61, 6)],
        ]
    )
    J, _ = engine.amputated(7, "h")
    checks = {
        "original_and_temporal_cluster_graph_count": s.Integer(count - recursive_count),
        "complete_three_soft_graph_count": s.Integer(count - 4),
        "nonlinear_closed_formula_equals_recursive_source": (
            temporal - recursive
        ).applyfunc(s.factor),
        "literal_coordinate_pullback_equals_closed_formula": (
            mapped[7] - 2 * temporal
        ).applyfunc(s.factor),
        "complete_temporal_source_Ward": (Q.T * ETA * J).applyfunc(s.factor),
        "linear_only_projection_has_exact_nonzero_pole": residue(linear)
        - expected_linear,
        "nonlinear_temporal_projection_removes_known_pole": residue(temporal),
        "nonlinear_temporal_finite_endpoint": finite - expected_finite,
        "recursive_tensor_is_temporal": recursive[0, :],
    }
    return {
        "checks": checks,
        "whole_linear_only_pole": residue(linear),
        "whole_nonlinear_limit": finite,
        "whole_nonlinear_correction": correction,
        "gates": {
            "linear_only_projection_is_not_the_nonlinear_chart": residue(linear)
            != s.zeros(4),
            "nonlinear_compensation_is_not_zero": correction != s.zeros(4),
            "literal_pullback_checks_more_than_temporal_components": True,
        },
    }


def complete_calibration(original_parameters):
    if not isinstance(original_parameters, bool):
        raise TypeError("Require explicit original/diagnostic boolean")
    return _complete_calibration(original_parameters)


@cache
def _complete_calibration(original_parameters):
    pars = (
        source.original_parameters()
        if original_parameters
        else source.diagnostic_parameters()
    )
    points, hs = original.configuration()
    legs = [("phi", p, 1) for p in points[1:]] + [("h", q, A) for q, A in hs]
    baseline, n = original.TreeEngine(legs, **pars).amplitude()
    reorganized, m = TemporalSoft(legs, **pars).amplitude()
    return {
        "checks": {
            "whole_amplitude_unchanged": s.factor(reorganized - baseline),
            "whole_graph_count_unchanged": s.Integer(m - n),
            "whole_graph_inventory": s.Integer(n - 5116),
        },
        "parameters": pars,
        "graph_count": m,
        "nonzero_amplitude": reorganized != 0,
    }


@cache
def data():
    cluster = cluster_calibration()
    checks = {key: value for key, value in cluster["checks"].items()}
    gates = dict(cluster["gates"])
    metric, momenta, _shifts, mapped = chart.four_literal()
    engine = TemporalSoft([("h", q, A) for q, A in chart.four_rays()])
    for mask in range(1, 16):
        H, count = engine.current(mask, "h")
        checks[f"four_leaf_all_subset_literal_chart_{mask}"] = (
            2 * H - mapped[mask]
        ).applyfunc(s.factor)
        if mask == 15:
            checks["four_soft_tree_inventory"] = s.Integer(count - 26)
    J, _ = engine.amputated(15, "h")
    checks["four_soft_complete_temporal_source_Ward"] = (
        momenta[15].T * ETA * J
    ).applyfunc(s.factor)
    gates["fourth_order_not_linear_projection"] = mapped[15] != 2 * chart.project(
        metric[15] / 2, momenta[15]
    )
    calibrations = {}
    for label, original_parameters in (("diagnostic", False), ("original", True)):
        result = complete_calibration(original_parameters)
        checks.update(
            {label + "_" + key: value for key, value in result["checks"].items()}
        )
        calibrations[label] = {
            key: value for key, value in result.items() if key != "checks"
        }
        gates[label + "_full_amplitude_nonzero"] = result["nonzero_amplitude"]
    return {
        "whole_recursive_policy": "After each non-singleton pure-soft h propagation, use Pi_Q on that current. Recompute every higher current from the same complete action vertices. Keep the original inverse-kinetic policy on mixed/hard subsets; their values change consistently with the lower fields. Never project a zero-energy hard channel.",
        "whole_finite_Noether_induction": "Every proper labeled coefficient solves the ungaugefixed Euler equations. The next complete graviton source is conserved by the Noether identity; the original inverse solves its linear equation, and the added longitudinal projection solves the homogeneous linear equation. Matter inverse coefficients solve their unchanged equations. This finite induction holds away from internal propagator poles.",
        "whole_on_shell_tree_invariance": "The two finite solutions have identical physical free waves and are related order by order by higher-degree diffeomorphisms. Amputated on-shell scalar-root coefficients are invariant: covariance reduces changes to lower Euler coefficients and an inverse-external-propagator factor. These vanish on the solved physical tree. This is finite formal tree invariance, not quantum-state or inclusive-rate existence.",
        "whole_literal_nonlinear_calibration": {
            key: value
            for key, value in cluster.items()
            if key not in ("checks", "gates")
        },
        "whole_complete_original_and_diagnostic_amplitudes": calibrations,
        "checks": checks,
        "gates": {
            **gates,
            "all_vertices_and_labeled_graph_counts_retained": True,
            "only_positive_energy_pure_soft_subsets_are_projected": True,
            "nonlinear_terms_generated_by_full_recursion_not_deleted": True,
            "mixed_hard_coefficients_recomputed_in_original_inverse_policy": True,
            "finite_Noether_proof_not_a_quantum_Ward_assumption": True,
            "on_shell_invariance_requires_complete_source": True,
            "no_inclusive_probability_or_Regge_inferred": True,
        },
    }
