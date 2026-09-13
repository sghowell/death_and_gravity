"""Full scalar four-jet ring, ordered endpoints and six-invariant geometry."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import quantum

d, p, u, m = s.symbols("d p u m", real=True)
a = s.Symbol("a", positive=True)
h = s.symbols("H0:6", real=True)
inv = s.symbols("T V W S UD UG", real=True)
ONE = (s.S.One, s.S.Zero, s.S.Zero, s.S.Zero, s.S.Zero)
ZERO = (s.S.Zero,) * 5


def add(x, y):
    return tuple(s.expand(x[j] + y[j]) for j in range(5))


def scale(x, c):
    return tuple(s.expand(c * v) for v in x)


def mul(x, y):
    return tuple(s.expand(sum(x[k] * y[j - k] for k in range(j + 1))) for j in range(5))


def reciprocal(x):
    out = [1 / x[0]]
    for j in range(1, 5):
        out.append(s.expand(-sum(x[k] * out[j - k] for k in range(1, j + 1)) / x[0]))
    return tuple(out)


def sqrtone(x):
    assert x[0] == 1
    out = [s.S.One]
    for j in range(1, 5):
        out.append(s.expand((x[j] - sum(out[k] * out[j - k] for k in range(1, j))) / 2))
    return tuple(out)


def shift(x, n):
    return (s.S.Zero,) * n + tuple(x[: 5 - n])


def power(x, n):
    out = ONE
    for _ in range(n):
        out = mul(out, x)
    return out


def derivative(value):
    if s.sympify(value).has(h[-1]):
        raise ValueError("The represented time jet does not supply H6")
    return s.expand(
        a * h[0] * s.diff(value, a)
        + sum(h[j + 1] * s.diff(value, h[j]) for j in range(5))
    )


def diff(x):
    return tuple(derivative(v) for v in x)


def truncate(x, k):
    return tuple(x[: k + 1]) + (s.S.Zero,) * (4 - k)


def sharp(x):
    return tuple(s.expand(s.conjugate(v)) for v in x)


@cache
def frequency(leg):
    if not isinstance(leg, str) or leg not in ("k", "l"):
        raise ValueError("Retain one of the two stated internal scalar legs")
    mass2 = m * m * a * a
    square = ONE if leg == "k" else (s.S.One, -2 * p * u, p * p, s.S.Zero, s.S.Zero)
    square = add(square, (0, 0, mass2, 0, 0))
    O = sqrtone(square)
    Oinv = reciprocal(O)
    dz = scale(shift(power(Oinv, 2), 2), -mass2)
    ad = quantum.adiabatic()
    sub = {quantum.D: d, **dict(zip(quantum.h, h))}
    P20 = s.expand(ad["P2"].subs(quantum.z, 1).subs(sub, simultaneous=True))
    P2z = s.expand(
        s.diff(ad["P2"], quantum.z).subs(quantum.z, 1).subs(sub, simultaneous=True)
    )
    P40 = s.expand(ad["P4"].subs(quantum.z, 1).subs(sub, simultaneous=True))
    W = add(
        O,
        add(
            scale(shift(mul(add(scale(ONE, P20), scale(dz, P2z)), Oinv), 2), a * a),
            scale(shift(power(Oinv, 3), 4), a**4 * P40),
        ),
    )
    rate = scale(add(mul(diff(W), reciprocal(W)), scale(ONE, (d - 1) * h[0])), a / 2)
    momentum = add(scale(W, -s.I), scale(shift(rate, 1), -1))
    return O, W, momentum


@cache
def amplitudes():
    _, W1, P1 = frequency("k")
    _, W2, P2 = frequency("l")
    norm = scale(reciprocal(sqrtone(mul(W1, W2))), s.Rational(1, 2))
    dot = (-s.S.One, p * u, 0, 0, 0)
    numerator = add(
        scale(mul(P1, P2), -1), add(scale(dot, -1), (0, 0, m * m * a * a, 0, 0))
    )
    A = scale(mul(numerator, norm), s.Rational(1, 2))
    return {"trace": A, "gradient": norm}, reciprocal(add(W1, W2))


@cache
def endpoint_products():
    amps, phase = amplitudes()
    out = {}
    for geo, left, right in (
        ("00", "trace", "trace"),
        ("01", "trace", "gradient"),
        ("10", "gradient", "trace"),
        ("11", "gradient", "gradient"),
    ):
        detector = sharp(amps[left])
        source = {0: scale(amps[right], 1 / a)}
        g = scale(phase, a)
        for j in range(5):
            degree = 4 - j
            for r, row in source.items():
                product = mul(
                    truncate(phase, degree),
                    mul(truncate(detector, degree), truncate(row, degree)),
                )
                out[geo, j, r] = tuple(
                    s.expand(s.im(s.I * (-s.I) ** j * v)) for v in product[: 5 - j]
                )
            if j == 4:
                break
            following = {}
            for r, row in source.items():
                value = truncate(mul(truncate(g, 3 - j), truncate(row, 3 - j)), 3 - j)
                following[r] = add(following.get(r, ZERO), diff(value))
                following[r + 1] = add(following.get(r + 1, ZERO), value)
            source = following
    return out


@cache
def geometries():
    T, V, W, S, UD, UG = inv
    y = s.Symbol("y", real=True)
    rr = 1 - u * u
    b = d - 1
    c = u * (y - u)
    cA = -rr / b
    cB = rr / b + c
    return y, {
        "00": S,
        "01": s.expand(cA * S + cB * UD),
        "10": s.expand(cA * S + cB * UG),
        "11": s.expand(
            rr**2 * (2 * T - 4 * V + S - UD - UG + 3 * W) / (b * (b + 2))
            + rr * (y - 2 * u) ** 2 * (V - W) / b
            - rr * c * (UD + UG - 2 * W) / b
            + c * c * W
        ),
    }


def sphere(value):
    polynomial = s.Poly(s.expand(value), u)
    result = s.S.Zero
    for (j,), c in polynomial.terms():
        if j % 2 == 0:
            result += c * s.rf(s.Rational(1, 2), j // 2) / s.rf(d / 2, j // 2)
    return s.factor(result)


@cache
def radial_rows():
    y, geo = geometries()
    end = endpoint_products()
    out = {}
    for j in range(5):
        for r in range(j + 1):
            for degree in range(5 - j):
                value = 0
                for name in geo:
                    gy = s.Poly(geo[name], y)
                    for (gg,), coeff in gy.terms():
                        if gg <= degree:
                            value += coeff * p**gg * end[name, j, r][degree - gg]
                out[j, r, degree] = sphere(value)
    return out


@cache
def spatial_rows():
    return {
        key: s.factor(value - value.subs(p, 0)) for key, value in radial_rows().items()
    }


@cache
def logarithmic():
    rows = spatial_rows()
    return tuple(
        s.factor(sum(rows[j, r, 4 - j] for j in range(r, 5))) for r in range(5)
    )


def pole_check():
    from p8_vacuum_affine_ordered_scalar_symbol import matching as geometry

    _gg, curv = geometry.hessians(*inv)
    oldd = geometry.geometry.d
    oldt = geometry.jets.t
    # Compare at the actual reference only AFTER independent generic jet extraction.
    aa = (1 + oldt * oldt) ** 2
    HH = 4 * oldt / (1 + oldt * oldt)
    sub = {
        d: oldd,
        p: geometry.jets.p,
        a: aa,
        **{h[j]: s.diff(HH, oldt, j) for j in range(6)},
    }
    source = s.symbols("Gamma0 Gamma1 Gamma2 Gamma3 Gamma4", real=True)
    mode = (
        sum(logarithmic()[j] * source[j] for j in range(5))
        .subs(sub, simultaneous=True)
        .subs(oldd, 3)
    )
    pole = (
        -m * m * curv["R_old"] / 3
        + curv["R_squared"] / 36
        + (curv["Riemann_squared"] - curv["Ricci_squared"]) / 90
    ).subs(oldd, 3)
    return s.factor(pole - 16 * mode)
