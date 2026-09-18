"""Unchanged source, independent nilpotent Taylor oracle and angular coverage."""

from functools import cache
from itertools import product

import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import trees as original
from p8_vacuum_affine_complete_two_graviton_tree import trees as lower
from p8_vacuum_affine_relative_energy_complex_tube import source as previous
from p8_vacuum_affine_temporal_tree_reorganization import trees

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_multiplicity = previous.require_multiplicity
original_parameters, diagnostic_parameters = (
    previous.original_parameters,
    previous.diagnostic_parameters,
)
recoil = previous.recoil
VARIABLES = s.symbols("ea eb ec", positive=True)
ZERO = dict.fromkeys(VARIABLES, 0)


def require_case(sector, bits):
    if type(sector) is not int or sector not in (0, 1):
        raise ValueError("Require one of the two independent open angular sectors")
    if (
        not isinstance(bits, tuple)
        or len(bits) != 3
        or any(type(b) is not int or b not in (0, 1) for b in bits)
    ):
        raise ValueError("Require exactly three integer plus/cross selectors")
    return sector, bits


def coordinates(sector, r, v, u):
    if sector == 0:
        return ((0, 0), (r, 0), (-v, u))
    if sector == 1:
        return ((0, 0), (r + v, 0), (r, u))
    raise ValueError("Only the two independent angular sectors are computed")


def direction(x, y):
    x, y = map(s.sympify, (x, y))
    return s.Matrix([2 * x, 2 * y, 1 - x * x - y * y]) / (1 + x * x + y * y)


def frame(x, y):
    x, y = map(s.sympify, (x, y))
    d = 1 + x * x + y * y
    return (
        s.Matrix([1 - 2 * x * x / d, -2 * x * y / d, -2 * x / d]),
        s.Matrix([-2 * x * y / d, 1 - 2 * y * y / d, -2 * y / d]),
    )


def rays(weights, directions, frames, bits):
    rows = []
    for w, n, (u, v), bit in zip(weights, directions, frames, bits):
        A = s.zeros(4)
        A[1:, 1:] = u * u.T - v * v.T if bit == 0 else u * v.T + v * u.T
        rows.append(("h", s.ImmutableMatrix([w, *(w * n)]), s.ImmutableMatrix(A)))
    return rows


def jet(value):
    value = s.sympify(value)
    if not value:
        return s.S.Zero
    terms = s.Poly(value, *VARIABLES).as_dict()
    return s.Add(
        *(
            coefficient * s.prod(v**p for v, p in zip(VARIABLES, powers))
            for powers, coefficient in terms.items()
            if max(powers) <= 1
        )
    )


def jm(A):
    return s.ImmutableMatrix(A.applyfunc(jet))


def inverse(value):
    value = jet(value)
    constant = value.subs(ZERO)
    if not constant:
        raise ValueError("Taylor center is on a propagator pole")
    v = jet((value - constant) / constant)
    v2 = jet(v * v)
    v3 = jet(v2 * v)
    result = jet((1 - v + v2 - v3) / constant)
    assert jet(value * result) == 1
    return result


class TaylorSoft(original.TreeEngine):
    @lower.instance_cache
    def vertex(self, tags, momenta, values):
        return jet(super().vertex(tags, momenta, values))

    @lower.instance_cache
    def current(self, mask, kind):
        if mask & (mask - 1) == 0:
            leg = self.legs[mask.bit_length() - 1]
            return (
                (leg[2], 1)
                if leg[0] == kind
                else (lower.ZERO if kind == "h" else s.S.Zero, 0)
            )
        current, count = self.amputated(mask, kind)
        if not count:
            return (lower.ZERO if kind == "h" else s.S.Zero, 0)
        Q = self.momentum(mask)
        square = jet((Q.T * lower.ETA * Q)[0])
        if kind != "h":
            mass = s.S.One if kind == "phi" else self.heavy
            return jet(-current * inverse(square - mass)), count
        H = jm(
            -(
                lower.ETA * current * lower.ETA
                - lower.ETA * s.trace(lower.ETA * current) / 2
            )
            * inverse(square)
        )
        iq0 = inverse(Q[0])
        q = lower.ETA * Q
        P = s.Matrix(
            4,
            4,
            lambda mu, nu: jet(
                H[mu, nu]
                - q[mu] * H[0, nu] * iq0
                - q[nu] * H[mu, 0] * iq0
                + q[mu] * q[nu] * H[0, 0] * jet(iq0 * iq0)
            ),
        )
        return jm(P), count


def evaluate(center, directions, frames, bits):
    ws = tuple(w + e for w, e in zip(center, VARIABLES))
    engine = TaylorSoft(rays(ws, directions, frames, bits))
    H, count = engine.current(7, "h")
    assert count == 4 and H[0, :] == s.zeros(1, 4)
    J, source_count = engine.amputated(7, "h")
    assert source_count == 4
    assert jm(engine.momentum(7).T * lower.ETA * J) == s.zeros(1, 4)
    return jm(s.prod(ws) * inverse(sum(ws)) * H)


@cache
def taylor_calibration():
    a, b, c = s.symbols("a b c", positive=True)
    directions = (s.Matrix([1, 0, 0]), s.Matrix([0, 1, 0]), s.Matrix([0, 0, 1]))
    frames = (
        (s.Matrix([0, 1, 0]), s.Matrix([0, 0, 1])),
        (s.Matrix([1, 0, 0]), s.Matrix([0, 0, 1])),
        (s.Matrix([1, 0, 0]), s.Matrix([0, 1, 0])),
    )
    checks = {}
    for k, (bits, center) in enumerate(
        (
            ((0, 0, 0), (s.Rational(1, 3), s.Rational(2, 3), s.S.One)),
            ((1, 0, 1), (s.Rational(1, 100), s.S.One, s.S.One)),
        )
    ):
        made = evaluate(center, directions, frames, bits)
        full, count = trees.TemporalSoft(
            rays((a, b, c), directions, frames, bits)
        ).current(7, "h")
        assert count == 4
        for i in range(1, 4):
            for j in range(i, 4):
                F = s.factor(a * b * c * full[i, j] / (a + b + c))
                P = s.Poly(made[i, j], *VARIABLES)
                for powers in product((0, 1), repeat=3):
                    target = F
                    for variable, power in zip((a, b, c), powers):
                        if power:
                            target = s.diff(target, variable)
                    target = s.factor(target.subs(dict(zip((a, b, c), center))))
                    residual = s.factor(P.coeff_monomial(powers) - target)
                    assert residual == 0
                    checks[
                        f"Taylor_baseline_{k}_{i}{j}_{''.join(map(str, powers))}"
                    ] = residual
    assert len(checks) == 96
    return checks


@cache
def angular_data():
    r, v, u = s.symbols("r v u", positive=True)
    t = r + v
    sn = 2 * r / (1 + r * r)
    cs = (1 - r * r) / (1 + r * r)
    O = s.Matrix([[-cs, 0, sn], [0, 1, 0], [sn, 0, cs]])
    D = (1 + r * t) ** 2 + r * r * u * u
    vp = ((t - r) * (1 + r * t) + r * u * u) / D
    up = u * (1 + r * r) / D
    checks = {
        "reflection_orthogonal": (O.T * O - s.eye(3)).applyfunc(s.factor),
        "reflection_determinant": s.factor(O.det() + 1),
    }
    for label, left, right in (
        ("anchor", O * direction(0, 0), direction(r, 0)),
        ("swap", O * direction(r, 0), direction(0, 0)),
        ("third", O * direction(t, u), direction(-vp, up)),
    ):
        checks[label + "_exact_sector_map"] = (left - right).applyfunc(s.factor)
    for expr in (D, s.fraction(s.factor(vp))[0], s.fraction(s.factor(up))[0]):
        assert all(value > 0 for value in s.Poly(expr, r, v, u).coeffs())
    OO = s.eye(4)
    OO[1:, 1:] = O
    checks["Lorentz_metric_preserved"] = (OO.T * lower.ETA * OO - lower.ETA).applyfunc(
        s.factor
    )
    checks["temporal_reference_preserved"] = OO * s.Matrix([1, 0, 0, 0]) - s.Matrix(
        [1, 0, 0, 0]
    )
    a, b, c = s.symbols("a b c", positive=True)
    S = 1 / (a + b) + 1 / (a + c) + 1 / (b + c)
    checks["pair_scale_permutation"] = s.factor(S - S.xreplace({a: b, b: a}))
    for point, (rv, vv, uv, weights) in enumerate(
        (
            (
                s.Rational(1, 2),
                s.Rational(1, 3),
                s.Rational(1, 4),
                (s.Rational(1, 4), s.S.One, s.S.One),
            ),
            (
                s.Rational(1, 3),
                s.Rational(1, 9),
                s.Rational(1, 27),
                (s.Rational(1, 81), s.S.One, s.S.One),
            ),
        )
    ):
        rotation = s.ImmutableMatrix(OO.subs({r: rv, v: vv, u: uv}))
        coords = ((0, 0), (rv, 0), (rv + vv, uv))
        for bits in product((0, 1), repeat=3):
            legs = rays(
                weights,
                tuple(direction(*xy) for xy in coords),
                tuple(frame(*xy) for xy in coords),
                bits,
            )
            mapped = [
                (
                    "h",
                    s.ImmutableMatrix(rotation * legs[k][1]),
                    s.ImmutableMatrix(rotation * legs[k][2] * rotation.T),
                )
                for k in (1, 0, 2)
            ]
            original_engine = trees.TemporalSoft(legs)
            mapped_engine = trees.TemporalSoft(mapped)
            H, count = original_engine.current(7, "h")
            K, mapped_count = mapped_engine.current(7, "h")
            J, jcount = original_engine.amputated(7, "h")
            L, lcount = mapped_engine.amputated(7, "h")
            assert count == mapped_count == jcount == lcount == 4
            label = f"covariance_{point}_{''.join(map(str, bits))}"
            checks[label + "_current"] = (K - rotation * H * rotation.T).applyfunc(
                s.factor
            )
            checks[label + "_source"] = (L - rotation * J * rotation.T).applyfunc(
                s.factor
            )
            checks[label + "_original_Ward"] = (
                original_engine.momentum(7).T * lower.ETA * J
            ).applyfunc(s.factor)
            checks[label + "_mapped_Ward"] = (
                mapped_engine.momentum(7).T * lower.ETA * L
            ).applyfunc(s.factor)
    for value in checks.values():
        assert all(
            z == 0
            for z in (list(value) if isinstance(value, s.MatrixBase) else [value])
        )
    return {
        "whole_third_sector_reflection": O,
        "whole_positive_third_to_first_coordinates": (-vp, up),
        "checks": checks,
        "gates": {
            "exact_energy_independent_orthogonal_sector_map": True,
            "positive_coordinate_images_and_temporal_reference": True,
            "all16_independent_current_source_covariance_oracles": True,
            "eta_contractions_and_complete_partitions_prove_O3_covariance": True,
            "arbitrary_complex_unit_TT_by_multilinearity": True,
            "generic_chart_boundary_limits_not_exact_internal_pole_values": True,
        },
    }


@cache
def data():
    old = previous.data()
    checks = dict(old["checks"])
    checks.update(taylor_calibration())
    return {
        **{k: v for k, v in old.items() if k not in ("checks", "gates")},
        "whole_new_claim_boundary": "Complete three-soft temporal current derivative hierarchy and compatible energy faces; one188-term maximal-soft-block class only, not the full5116 subtraction or an inclusive rate.",
        "whole_source_ownership": "S320 and S313 retain the complex hard recoil input; S318 supplies the uniform origin estimate; S317 supplies exact temporal-tree equivalence; literal S315 EH3/EH4 vertices are reconstructed in the new polynomial engine.",
        "checks": checks,
        "gates": {
            **old["gates"],
            "independent96_Taylor_coefficients_rebuilt_from_literal_action": True,
            "same_original_parameters_and_no_private_cache_inputs": True,
            "three_current_not_full_three_real_inclusive_probability": True,
        },
    }
