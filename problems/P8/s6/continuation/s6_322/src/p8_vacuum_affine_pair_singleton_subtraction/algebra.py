"""Literal off-shell scalar identities and independent original387-class calibration."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import trees as lower
from p8_vacuum_affine_relative_energy_complex_tube import continuation as analytic
from p8_vacuum_affine_two_real_soft_overlap import external

from . import source

clean = analytic.canonical


def matrix(A):
    return s.ImmutableMatrix(A.applyfunc(clean))


def dot(p, q):
    return (p.T * lower.ETA * q)[0]


def trace(A):
    return s.trace(lower.ETA * A)


def numerator(p, A, Q, mass=1):
    return (
        (p.T * A * p)[0]
        + (p.T * A * Q)[0]
        - trace(A) * (dot(p, p) - mass + dot(p, Q)) / 2
    )


class CoreTwo(analytic.GaussianEngine):
    """All387 hard trees, with nonsingleton pure-soft branches disabled."""

    marked = False
    marker = s.Symbol("light_scalar_propagator_marker", real=True)

    @lower.instance_cache
    def current(self, mask, kind):
        if mask.bit_count() > 1 and all(
            self.legs[i][0] == "h" for i in range(len(self.legs)) if mask >> i & 1
        ):
            return (lower.ZERO if kind == "h" else s.S.Zero), 0
        value, count = super().current(mask, kind)
        if self.marked and kind == "phi" and count and mask.bit_count() > 1:
            value = self.marker * value
        return value, count


@cache
def identities():
    p = s.Matrix(s.symbols("p0:4", real=True))
    Q = s.Matrix(s.symbols("q0:4", real=True))
    k = s.Matrix(s.symbols("k0:4", real=True))
    mass = s.Symbol("mass", real=True)
    entries = s.symbols("A0:10", real=True)
    A = s.zeros(4)
    at = 0
    for i in range(4):
        for j in range(i, 4):
            A[i, j] = A[j, i] = entries[at]
            at += 1
    shift = (
        2 * (p.T * A * k)[0]
        + (k.T * A * k)[0]
        + (k.T * A * Q)[0]
        - trace(A) * (2 * dot(p, k) + dot(k, k) + dot(k, Q)) / 2
    )
    den = lambda q: dot(p, q) + dot(q, q) / 2
    checks = {
        "generic_arbitrary_tensor_literal_vertex": s.expand(
            lower.scalar_vertex(
                s.ImmutableMatrix(p),
                s.ImmutableMatrix(-p - Q),
                (s.ImmutableMatrix(A),),
                mass,
            )
            + 2 * numerator(p, A, Q, mass)
        ),
        "generic_offshell_shift_including_trace": s.expand(
            numerator(p + k, A, Q, mass) - numerator(p, A, Q, mass) - shift
        ),
        "generic_total_denominator": s.expand(den(Q + k) - den(Q) - den(k) - dot(Q, k)),
    }
    old, controls = external.general_identities()
    checks.update({"ordering_" + name: value for name, value in old.items()})
    return {
        "checks": checks,
        "gates": {
            "all_ordering_negative_controls_retained": all(controls.values()),
            "no_null_TT_or_on_shell_premise_for_generic_vertex": True,
        },
        "whole_generic_numerator": "N_A(p)=p.A.p+p.A.Q-tr(eta*A)*(p^2-m+p.Q)/2. Literal Phi2-h vertex(p,-p-Q)=-2N_A. The full shifted numerator retains trace, virtual-mass and nontransverse terms.",
        "whole_ordering_controls": controls,
    }


@cache
def calibration():
    pars = source.original_parameters()
    E = s.Rational(5, 4)
    r0 = s.Rational(3, 4)
    parameter = s.Rational(31, 16)
    ep = (parameter + 1 / parameter) / 2
    rp = (parameter - 1 / parameter) / 2
    out = s.Matrix([s.Rational(3, 5), 0, s.Rational(4, 5)])
    points = (
        s.ImmutableMatrix([-E, 0, 0, -r0]),
        s.ImmutableMatrix([-E, 0, 0, r0]),
        s.ImmutableMatrix([ep, *(rp * out)]),
        s.ImmutableMatrix([ep, *(-rp * out)]),
    )
    W = 2 * (E - ep)
    weights = (3 * W / 8, 5 * W / 16, 5 * W / 16)
    directions = (
        s.Matrix([-1, 0, 0]),
        s.Matrix([s.Rational(3, 5), s.Rational(4, 5), 0]),
        s.Matrix([s.Rational(3, 5), -s.Rational(4, 5), 0]),
    )
    hs = []
    for w, n in zip(weights, directions):
        U = s.Matrix([-n[1], n[0], 0])
        V = s.Matrix([0, 0, 1])
        A = s.zeros(4)
        A[1:, 1:] = U * U.T - V * V.T
        hs.append((s.ImmutableMatrix([w, *(w * n)]), s.ImmutableMatrix(A)))
    checks = {
        "total_conservation": sum(points, lower.VECTOR_ZERO)
        + sum((q for q, A in hs), lower.VECTOR_ZERO)
    }
    recoil_points, _rays, _born = source.recoil.momenta(E, tuple(q for q, A in hs), out)
    for i, (p, q) in enumerate(zip(points, recoil_points)):
        checks[f"same_original_recoil_{i}"] = matrix(p - q)
    pair_engine = analytic.TemporalEngine([("h", q, A) for q, A in hs[:2]], **pars)
    pair, pair_count = pair_engine.current(3, "h")
    QA = hs[0][0] + hs[1][0]
    qc, B = hs[2]
    checks.update(
        {
            "pair_inventory": s.Integer(pair_count - 1),
            "temporal_pair_time_row": pair[0, :],
            "singleton_null": dot(qc, qc),
        }
    )

    def core(ps):
        return clean(
            external.core(ps, "C", pars["heavy"], pars["cubic"], pars["contact"])
            + external.core(ps, "H", pars["heavy"], pars["cubic"], pars["contact"])
            + external.core(ps, "GR", pars["heavy"], pars["cubic"], pars["contact"])
            / pars["kappa"]
        )

    def explicit_double(A):
        fields = (A, B)
        momenta = (QA, qc)
        nums = [[numerator(p, H, q) for p in points] for H, q in zip(fields, momenta)]
        dens = [[dot(p, q) + dot(q, q) / 2 for p in points] for q in momenta]
        answer = s.S.Zero
        for i in range(4):
            for j in range(4):
                if i == j:
                    continue
                shifted = list(points)
                shifted[i] = matrix(shifted[i] + QA)
                shifted[j] = matrix(shifted[j] + qc)
                answer += (
                    nums[0][i] * nums[1][j] * core(shifted) / (dens[0][i] * dens[1][j])
                )
        for i, p in enumerate(points):
            Ad, Bd = dens[0][i], dens[1][i]
            Z = dot(QA, qc)
            after_B = numerator(p + QA, B, qc)
            after_A = numerator(p + qc, A, QA)
            factor = nums[0][i] * after_B / (Ad * (Ad + Bd + Z)) + nums[1][
                i
            ] * after_A / (Bd * (Ad + Bd + Z))
            shifted = list(points)
            shifted[i] = matrix(p + QA + qc)
            answer += factor * core(shifted)
        return clean(answer / pars["kappa"])

    for label, A in (
        ("actual_complete_pair", pair),
        ("arbitrary_spatial_trace", s.ImmutableMatrix(s.diag(0, 1, 2, 3))),
    ):
        engine = CoreTwo(
            [("phi", p, 1) for p in points[1:]] + [("h", QA, A), ("h", qc, B)], **pars
        )
        engine.marked = True
        value, count = engine.amplitude()
        polynomial = s.Poly(value, engine.marker)
        checks[label + "_whole387_inventory"] = s.Integer(count - 387)
        checks[label + "_maximum_light_degree"] = s.Integer(polynomial.degree() - 2)
        checks[label + "_independent140_double_external"] = clean(
            polynomial.coeff_monomial(engine.marker**2) - explicit_double(A)
        )
        if label == "actual_complete_pair":
            actual_complete = clean(value.subs(engine.marker, 1))

    class CompletePairPlusSingleton(analytic.GaussianEngine):
        @lower.instance_cache
        def current(self, mask, kind):
            if all(
                self.legs[i][0] == "h" for i in range(len(self.legs)) if mask >> i & 1
            ):
                if mask == 24 and kind == "h":
                    return pair, 1
                if mask == 32:
                    return super().current(mask, kind)
                return (lower.ZERO if kind == "h" else s.S.Zero), 0
            return super().current(mask, kind)

    restricted = CompletePairPlusSingleton(
        [("phi", p, 1) for p in points[1:]] + [("h", q, A) for q, A in hs], **pars
    )
    value, count = restricted.amplitude()
    checks["independent_temporal387_inventory"] = s.Integer(count - 387)
    checks["collapsed_core_equals_independent_temporal_class"] = clean(
        value - actual_complete
    )
    return {
        "checks": checks,
        "gates": {
            "physical_recoil_and_energy_domain": bool(0 < W < s.Rational(1, 8)),
            "actual_composite_momentum_timelike": bool(dot(QA, QA) > 0),
            "actual_complete_pair_has_nonzero_trace": bool(trace(pair) != 0),
            "actual_complete_pair_not_transverse": matrix(pair * QA) != s.zeros(4, 1),
            "original_couplings_used_without_diagnostic_substitution": pars
            == source.original_parameters(),
            "no_multioffshell_Ward_premise": True,
        },
        "whole_original387_calibration": "Two independent marked387 sums, for the actual complete temporal pair and an arbitrary traceful spatial tensor; full degree-two coefficient equals the explicit140 double-external formula. An independently filtered temporal387 class equals the collapsed composite-pair core.",
        "whole_marked_calibration_count": 2,
        "whole_temporal_class_inventory": s.Integer(387),
    }


@cache
def data():
    generic, actual = identities(), calibration()
    return {
        "checks": {
            **{"generic_" + k: v for k, v in generic["checks"].items()},
            **{"original_" + k: v for k, v in actual["checks"].items()},
        },
        "gates": {**generic["gates"], **actual["gates"]},
        "whole_generic_identities": {
            k: v for k, v in generic.items() if k not in ("checks", "gates")
        },
        "whole_original_calibration": {
            k: v for k, v in actual.items() if k not in ("checks", "gates")
        },
    }
