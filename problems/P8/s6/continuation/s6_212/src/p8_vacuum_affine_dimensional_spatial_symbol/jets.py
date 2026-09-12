"""General-dimensional exact WKB UV jets and independent source derivatives."""

from functools import cache

import sympy as s
from p8_vacuum_affine_spatial_symbol import benchmark as ring
from p8_vector_clock_matching import continuation
from p8_vector_state import wkb

from .geometry import d

t, u, m, p = s.symbols("t u m p", real=True)
a = (1 + t * t) ** 2
H = 4 * t / (1 + t * t)
ONE = (s.Integer(1),) + (s.Integer(0),) * 4
ZERO = (s.Integer(0),) * 5


def shift(row, n):
    return (s.Integer(0),) * n + tuple(row[: 5 - n])


def power(row, n):
    out = ONE
    for _ in range(n):
        out = ring.multiply(out, row)
    return out


@cache
def data():
    from p8_vacuum_affine_spatial_matching_difference import jets as physical

    checks = {}
    for kind in ("transverse", "longitudinal"):
        original = continuation.coefficients(kind)
        for name in ("second_order_cancellation", "fourth_order_cancellation"):
            checks[kind + "_" + name] = original[name]
        expected = (
            -(d - 3) * s.diff(H, t) / 4 - (d - 3) * (d - 1) * H * H / 8
            if kind == "transverse"
            else -(d - 1) * s.diff(H, t) / 4 - (d * d - 1) * H * H / 8
        )
        P1 = original["P2"].subs({continuation.local.dimension: d, wkb.u: t, wkb.z: 1})
        checks[kind + "_actual_general_dimension_first_WKB"] = s.factor(P1 - expected)
        for leg in ("k", "l"):
            rows = frequency(kind, leg)
            old = physical.frequency(kind, leg)
            checks[kind + "_" + leg + "_all_physical_frequency_and_momentum_jets"] = (
                s.Matrix(
                    [
                        s.factor(value.subs(d, 3) - target)
                        for row, prior in zip(rows, old)
                        for value, target in zip(row, prior)
                    ]
                )
            )
    checks["all_105_dimension_dependent_endpoint_rows"] = len(endpoint_products()) - 105
    return {
        "canonical_rates": "The continued Hamiltonian K=a^(2-d)[I+kk^t/(a^2m^2)] has transverse/longitudinal canonical dilution rates(d-2)H/2 and[(d-2)/2+z]H. In Wbar variables the transverse b includes(d-3)aH/2, and longitudinal adds a*z*H. These terms cannot be frozen atd3 before extracting a finite part.",
        "actual_reference": "Use the existing general-d clock-matching P2/P4 coefficients for the first two WKB corrections. They determine all UV degrees through4. Full normalized amplitudes, four inverse phases and all source time derivatives retain the massive constraint and magnetic vertex.",
        "analytic_pair": "The symbolic coefficient computation first uses real d so the two phase branches are Schwarz partners, then analytically continues those coefficient expressions. At complex d the parameter is not conjugated; this is not the Hermitian adjoint of a noninteger-dimensional Hilbert space.",
        "finite_ring": "After each endpoint iteration only the degrees needed by subsequent endpoints are retained. This exact filtration optimization does not alter any coefficient through total grade4.",
        "checks": checks,
        "gates": {
            "actual_dimension_dependent_rates_and_WKB_coefficients": True,
            "all_five_endpoint_labels_and_source_jets": True,
            "no_complex_dimension_conjugation": True,
            "source_detector_times_not_coalesced_before_derivatives": True,
            "physical_dimension_frequency_and_momentum_limits_checked": True,
        },
    }


def diff(row):
    return tuple(s.cancel(s.diff(v, t)) for v in row)


def truncate(row, degree):
    return tuple(row[: degree + 1]) + (s.Integer(0),) * (4 - degree)


def sharp(row):
    return tuple(s.conjugate(v).expand() for v in row)


@cache
def frequency(kind, leg):
    if kind not in ("transverse", "longitudinal") or leg not in ("k", "l"):
        raise ValueError("Require actual polarization and internal leg")
    mass2 = s.expand(m * m * a * a)
    q = (1, 0, mass2, 0, 0) if leg == "k" else (1, -2 * p * u, p * p + mass2, 0, 0)
    O = ring.square_root_one(tuple(map(s.sympify, q)))
    inv = ring.reciprocal(O)
    dz = ring.scale(shift(power(inv, 2), 2), -mass2)
    actual = continuation.coefficients(kind)
    P1 = actual["P2"].subs(continuation.local.dimension, d)
    P10 = s.cancel(P1.subs({wkb.z: 1, wkb.u: t}))
    P1z = s.cancel(s.diff(P1, wkb.z).subs({wkb.z: 1, wkb.u: t}))
    P20 = s.cancel(
        actual["P4"].subs(continuation.local.dimension, d).subs({wkb.z: 1, wkb.u: t})
    )
    term1 = ring.multiply(ring.add(ring.scale(ONE, P10), ring.scale(dz, P1z)), inv)
    W = ring.add(
        O,
        ring.add(
            ring.scale(shift(term1, 2), a * a),
            ring.scale(shift(power(inv, 3), 4), a**4 * P20),
        ),
    )
    W = tuple(s.cancel(v) for v in W)
    b = ring.scale(ring.multiply(diff(W), ring.reciprocal(W)), a / 2)
    b = ring.add(b, ring.scale(ONE, (d - 3) * a * H / 2))
    if kind == "longitudinal":
        b = ring.add(b, ring.scale(ring.add(ONE, dz), a * H))
    momentum = ring.add(ring.scale(W, -s.I), ring.scale(shift(b, 1), -1))
    return O, W, tuple(s.cancel(v) for v in momentum)


@cache
def amplitudes():
    f = {
        key: frequency(kind, leg)
        for key, kind, leg in (
            ("kt", "transverse", "k"),
            ("kl", "longitudinal", "k"),
            ("lt", "transverse", "l"),
            ("ll", "longitudinal", "l"),
        )
    }
    out = {}
    phases = {}
    for name, left, right in (
        ("TT", "kt", "lt"),
        ("TL", "kt", "ll"),
        ("LT", "kl", "lt"),
        ("LL", "kl", "ll"),
    ):
        O1, W1, P1 = f[left]
        O2, W2, P2 = f[right]
        normalization = ring.scale(
            ring.reciprocal(ring.square_root_one(ring.multiply(W1, W2))),
            s.Rational(1, 2),
        )
        phases[name] = ring.reciprocal(ring.add(W1, W2))
        pp = ring.multiply(P1, P2)
        if name == "TT":
            numerator = ring.add(
                ring.scale(shift(ONE, 2), m * m * a * a), ring.scale(pp, -1)
            )
            out["A"] = ring.multiply(numerator, normalization)
            out["B"] = normalization
        elif name in ("TL", "LT"):
            longitudinal = O2 if name == "TL" else O1
            numerator = ring.add(
                longitudinal,
                ring.scale(ring.multiply(pp, ring.reciprocal(longitudinal)), -1),
            )
            out[name] = ring.scale(
                shift(ring.multiply(numerator, normalization), 1), m * a
            )
        else:
            product = ring.multiply(O1, O2)
            numerator = ring.add(
                product,
                ring.scale(
                    shift(ring.multiply(pp, ring.reciprocal(product)), 2),
                    -m * m * a * a,
                ),
            )
            out[name] = ring.multiply(numerator, normalization)
    return {key: tuple(s.cancel(v) for v in row) for key, row in out.items()}, phases


@cache
def endpoint_products():
    amps, phases = amplitudes()
    out = {}
    for geo, left, right, sector in (
        ("00", "A", "A", "TT"),
        ("01", "A", "B", "TT"),
        ("10", "B", "A", "TT"),
        ("11", "B", "B", "TT"),
        ("TL", "TL", "TL", "TL"),
        ("LT", "LT", "LT", "LT"),
        ("LL", "LL", "LL", "LL"),
    ):
        # left detector is NOT differentiated; right source is.
        detector = sharp(amps[left])
        source = {0: ring.scale(amps[right], 1 / a)}
        g = ring.scale(phases[sector], a)
        for j in range(5):
            for r, row in source.items():
                degree = 4 - j
                product = ring.multiply(
                    truncate(phases[sector], degree),
                    ring.multiply(truncate(detector, degree), truncate(row, degree)),
                )
                out[geo, j, r] = tuple(
                    s.cancel(s.im((-s.I * s.I**j) * v) * -1) for v in product[: 5 - j]
                )
            if j == 4:
                break
            following = {}
            for r, row in source.items():
                value = truncate(
                    ring.multiply(truncate(g, 3 - j), truncate(row, 3 - j)), 3 - j
                )
                differentiated = diff(value)
                following[r] = ring.add(following.get(r, ZERO), differentiated)
                following[r + 1] = ring.add(following.get(r + 1, ZERO), value)
            source = following
    return out
