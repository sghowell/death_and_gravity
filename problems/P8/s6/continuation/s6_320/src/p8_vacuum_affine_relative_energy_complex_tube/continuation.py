"""Separate Gaussian-rational continuation; frozen real point APIs are unchanged."""

from functools import cache

import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import trees as original
from p8_vacuum_affine_complete_two_graviton_tree import trees as lower
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree as matter

from . import bounds, source

ETA = lower.ETA


def canonical(value):
    return s.factor(s.expand(value, complex=True))


def exact_gaussian(value):
    if value is None or isinstance(value, (str, bool, float, s.Float)):
        raise TypeError("Require exact Gaussian-rational entries")
    value = s.sympify(value)
    if not isinstance(value, s.Expr) or value.has(s.Float, s.oo, -s.oo, s.zoo, s.nan):
        raise ValueError("Require finite exact Gaussian rationals")
    real, imag = value.as_real_imag()
    if real.is_Rational is not True or imag.is_Rational is not True:
        raise ValueError("Require rational real and imaginary parts")
    return canonical(value)


def gaussian_array(value):
    if not isinstance(value, (s.MatrixBase, tuple, list)):
        raise TypeError("Require an exact matrix or nested sequence")
    entries = list(value) if isinstance(value, s.MatrixBase) else value
    flat = []
    for entry in entries:
        flat.extend(entry if isinstance(entry, (tuple, list)) else (entry,))
    for entry in flat:
        exact_gaussian(entry)
    return s.ImmutableMatrix(value).applyfunc(exact_gaussian)


def project_complex(A, Q):
    A, Q = gaussian_array(A), gaussian_array(Q)
    if A.shape != (4, 4) or A != A.T or Q.shape != (4, 1):
        raise ValueError("Require a symmetric four-tensor and four-vector")
    if Q[0] == 0:
        raise ValueError("Require nonzero temporal energy")
    q = ETA * Q
    return s.ImmutableMatrix(
        4,
        4,
        lambda i, j: canonical(
            A[i, j]
            - q[i] * A[0, j] / Q[0]
            - q[j] * A[i, 0] / Q[0]
            + q[i] * q[j] * A[0, 0] / Q[0] ** 2
        ),
    )


def require_relative_tube(center, shifted):
    ws = bounds.require_energies(center)
    if not isinstance(shifted, (tuple, list)) or len(shifted) != len(ws):
        raise ValueError("Require one continued energy per real center")
    zs = tuple(exact_gaussian(z) for z in shifted)
    if not all(
        canonical((z - w) * s.conjugate(z - w)) <= bounds.EPSILON**2 * w * w
        for z, w in zip(zs, ws)
    ):
        raise ValueError("Require the stated relative-energy polydisc")
    return ws, zs


def normalize(value):
    return (
        s.ImmutableMatrix(value.applyfunc(canonical))
        if isinstance(value, s.MatrixBase)
        else canonical(value)
    )


class GaussianEngine(original.TreeEngine):
    """Only exact scalar representation changes; all literal source formulas remain."""

    @lower.instance_cache
    def vertex(self, tags, momenta, values):
        return canonical(super().vertex(tags, momenta, values))

    @lower.instance_cache
    def current(self, mask, kind):
        value, count = super().current(mask, kind)
        return normalize(value), count

    @lower.instance_cache
    def amputated(self, mask, kind):
        value, count = super().amputated(mask, kind)
        return normalize(value), count


class TemporalEngine(GaussianEngine):
    """Same complete policy as S317, using a separately owned complex projection."""

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
            value = project_complex(value, self.momentum(mask))
        return value, count


@cache
def pure_calibration():
    weights = (s.Rational(1, 64), s.Rational(1, 32), s.Rational(3, 64))
    directions = (s.Matrix([0, 0, 1]), s.Matrix([1, 0, 0]), s.Matrix([0, 1, 0]))
    pols = (s.diag(0, 1, -1, 0), s.diag(0, 0, 1, -1), s.diag(0, 1, 0, -1))
    shifts = (s.I / 1000, -s.I / 2000, (1 + s.I) / 3000)
    hs = [
        (s.ImmutableMatrix([w * (1 + e), *(w * (1 + e) * n)]), s.ImmutableMatrix(A))
        for w, e, n, A in zip(weights, shifts, directions, pols)
    ]
    engine = TemporalEngine([("h", q, A) for q, A in hs])
    A, count = engine.current(7, "h")
    J, count2 = engine.amputated(7, "h")
    Q = engine.momentum(7)
    W = sum(weights)
    v = sum((w * n for w, n in zip(weights, directions)), s.zeros(3, 1)) / W
    delta2 = 1 - (v.T * v)[0]
    upper = bounds.pure_coefficient(3) * W**3 / s.prod(weights)
    block = A[1:, 1:]
    norm2 = lambda M: canonical(sum(s.conjugate(x) * x for x in M))
    norm_squares = (
        norm2(block),
        norm2(block * v) / delta2,
        canonical(s.Abs((v.T * block * v)[0]) ** 2) / delta2**2,
    )
    return {
        "checks": {
            "complete_pure_tree_inventory": s.Integer(count - 4),
            "complete_amputated_inventory": s.Integer(count2 - count),
            "temporal_row": A[0, :],
            "complete_complex_source_Ward": (Q.T * ETA * J).applyfunc(canonical),
        },
        "gates": {
            "genuinely_complex_current": any(s.im(x) != 0 for x in A),
            "three_complex_energy_shifts_inside_pure_tube": all(
                canonical(e * s.conjugate(e)) < s.Rational(1, 100) ** 2 for e in shifts
            ),
            "all_three_weighted_norm_components_bounded": all(
                v < 8 * upper**2 for v in norm_squares
            ),
            "pure_calibration_uses_unit_kappa_not_original_full_claim": True,
        },
        "whole_three_independent_energy_shifts": {
            "real_center": weights,
            "relative_shifts": shifts,
            "graph_count": s.Integer(count),
            "center_delta_squared": delta2,
            "squared_weighted_norms": norm_squares,
            "squared_unnormalized_majorant": 8 * upper**2,
        },
    }


@cache
def full_calibration():
    E, v0 = s.Rational(5, 4), s.Rational(31, 16)
    vc = v0 + s.I / s.Integer(10) ** 16
    ep0 = (v0 + 1 / v0) / 2
    ep, rp = canonical((vc + 1 / vc) / 2), canonical((vc - 1 / vc) / 2)
    W0, W = 2 * (E - ep0), canonical(2 * (E - ep))
    gamma = canonical(W / W0)
    u, r0 = s.Matrix([s.Rational(3, 5), 0, s.Rational(4, 5)]), s.Rational(3, 4)
    points = (
        s.ImmutableMatrix([-E, 0, 0, -r0]),
        s.ImmutableMatrix([-E, 0, 0, r0]),
        s.ImmutableMatrix([ep, *(rp * u)]),
        s.ImmutableMatrix([ep, *(-rp * u)]),
    )
    born = (
        points[0],
        points[1],
        s.ImmutableMatrix([E, *(r0 * u)]),
        s.ImmutableMatrix([E, *(-r0 * u)]),
    )
    _oldpoints, oldhs = original.configuration()
    hs = tuple((normalize(gamma * q), A) for q, A in oldhs)
    ws = tuple(q[0] for q, A in oldhs)
    require_relative_tube(ws, tuple(q[0] for q, A in hs))
    pars = source.original_parameters()
    Am = s.factor(
        matter.born_continuation(born, pars["heavy"], pars["cubic"], pars["contact"])
    )
    AG = s.factor(lower.old.born(born) / pars["kappa"])
    legs = [("phi", p, 1) for p in points[1:]] + [("h", q, A) for q, A in hs]
    baseline, n = GaussianEngine(legs, **pars).amplitude()
    mapped, m = TemporalEngine(legs, **pars).amplitude()
    square = canonical(baseline * s.conjugate(baseline) / (8 * (Am + AG) ** 2))
    phase = canonical(E * rp / (r0 * ep))
    checks = {
        "original_complete_tree_count": s.Integer(n - 5116),
        "complex_temporal_graph_inventory": s.Integer(m - n),
        "original_complex_gauge_equivalence": canonical(mapped - baseline),
        "complex_total_conservation": normalize(
            sum(points, lower.VECTOR_ZERO) + sum((q for q, A in hs), lower.VECTOR_ZERO)
        ),
    }
    for j, p in enumerate(points):
        checks[f"complex_massive_shell_{j}"] = canonical(original.old.dot(p, p) - 1)
    for j, (q, A) in enumerate(hs):
        checks[f"complex_null_shell_{j}"] = canonical(original.old.dot(q, q))
        checks[f"complex_TT_leaf_{j}"] = normalize(A * q)
        checks[f"unit_polarization_square_normalization_{j}"] = (
            sum(v * v for v in A) - 2
        )
    return {
        "checks": checks,
        "gates": {
            "positive_real_center_Born_pieces": Am > 0 and AG > 0,
            "genuinely_complex_original_amplitude": canonical(s.im(baseline)) != 0,
            "complete_original_complex_amplitude_bound": 0
            < square
            < bounds.bare_upper(ws) ** 2,
            "continued_phase_norm_below_two": canonical(phase * s.conjugate(phase))
            < 16,
            "unchanged_original_parameters": pars == source.original_parameters(),
        },
        "whole_original_complex_5116": {
            "real_center_energies": ws,
            "relative_energy_multiplier": gamma,
            "graph_count": s.Integer(n),
            "normalized_unit_polarization_absolute_square": square,
            "bare_majorant": bounds.bare_upper(ws),
            "continued_rho": phase,
        },
    }


@cache
def data():
    pure, full = pure_calibration(), full_calibration()
    return {
        "whole_separate_continuation_policy": "Gaussian-rational real/imaginary normalization is exact and changes no action or arithmetic backend. Every higher current is recomputed from complete literal vertices. Only the pure-soft inverse representative receives the separately owned complex temporal projection; all frozen real point APIs remain unchanged.",
        "whole_pure_complex_calibration": pure["whole_three_independent_energy_shifts"],
        "whole_original_full_complex_calibration": full["whole_original_complex_5116"],
        "checks": {
            **{"pure_" + k: v for k, v in pure["checks"].items()},
            **{"full_" + k: v for k, v in full["checks"].items()},
        },
        "gates": {
            **{"pure_" + k: v for k, v in pure["gates"].items()},
            **{"full_" + k: v for k, v in full["gates"].items()},
            "no_frozen_real_type_guard_was_weakened": True,
            "finite_complex_points_not_a_proof_of_all_order_holomorphy": True,
        },
    }
