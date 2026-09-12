"""Independent literal curvature, integer geometry and continued finite-part checks."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_dimensional_spatial_symbol import (
    audit,
    density,
    geometry,
    jets,
    matching,
)
from p8_vacuum_affine_spatial_symbol import sectors
from p8_vector_clock_matching import continuation
from p8_vector_hadamard import series as physical_series
from p8_vector_state import wkb

A = s.symbols("a0:5", positive=True)
Dtime = s.symbols("D0:5", real=True)
Gtime = s.symbols("G0:5", real=True)
p = s.Symbol("p", real=True)
ZERO = (s.Integer(0),) * 4


def add(*rows):
    return tuple(s.expand(sum(row[j] for row in rows)) for j in range(4))


def scale(row, value):
    return tuple(s.expand(v * value) for v in row)


def mul(a, b):
    if a == ZERO or b == ZERO:
        return ZERO
    return (
        s.expand(a[0] * b[0]),
        s.expand(a[0] * b[1] + a[1] * b[0]),
        s.expand(a[0] * b[2] + a[2] * b[0]),
        s.expand(a[0] * b[3] + a[1] * b[2] + a[2] * b[1] + a[3] * b[0]),
    )


def scalar(value):
    return (s.sympify(value), 0, 0, 0)


def time_derivative(value):
    return s.expand(
        sum(
            s.diff(value, jets[j]) * jets[j + 1]
            for jets in (A, Dtime, Gtime)
            for j in range(4)
        )
    )


def matrix_product(a, b):
    return [
        [add(*(mul(a[i][k], b[k][j]) for k in range(len(b)))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def trace_square(matrix):
    return add(
        *(
            mul(matrix[i][j], matrix[j][i])
            for i in range(len(matrix))
            for j in range(len(matrix))
        )
    )


def literal_curvatures(D, G):
    dim = D.rows
    N = dim + 1
    a, a1, a2 = A[:3]
    S = (D * G + G * D) / 2

    def e(i, j, r=0, z=0):
        if i == 0 or j == 0:
            return scalar(-1) if i == j == 0 and r == z == 0 else ZERO
        i, j = i - 1, j - 1
        if z == 0:
            mixed = sum(
                s.binomial(r, k) * Dtime[k] * Gtime[r - k] for k in range(r + 1)
            )
            return (
                s.Integer(i == j) if r == 0 else 0,
                Dtime[r] * D[i, j],
                Gtime[r] * G[i, j],
                mixed * S[i, j],
            )
        return (
            0,
            (-s.I * p) ** z * Dtime[r] * D[i, j],
            (s.I * p) ** z * Gtime[r] * G[i, j],
            0,
        )

    def metric(i, j, r=0, z=0):
        if r == 0:
            return scale(e(i, j, 0, z), a * a)
        if r == 1:
            return add(scale(e(i, j, 0, z), 2 * a * a1), scale(e(i, j, 1, z), a * a))
        if r == 2:
            return add(
                scale(e(i, j, 0, z), 2 * (a1 * a1 + a * a2)),
                scale(e(i, j, 1, z), 4 * a * a1),
                scale(e(i, j, 2, z), a * a),
            )
        raise ValueError("Only metric jets through two")

    @cache
    def derivative(i, j, indices):
        if any(k not in (0, dim) for k in indices):
            return ZERO
        return metric(i, j, indices.count(0), indices.count(dim))

    inverse = []
    for i in range(N):
        row = []
        for j in range(N):
            if i == 0 or j == 0:
                value = scalar(-1 / a**2) if i == j == 0 else ZERO
            else:
                value = scale(
                    (
                        s.Integer(i == j),
                        -Dtime[0] * D[i - 1, j - 1],
                        -Gtime[0] * G[i - 1, j - 1],
                        Dtime[0] * Gtime[0] * S[i - 1, j - 1],
                    ),
                    a**-2,
                )
            row.append(value)
        inverse.append(row)

    @cache
    def christoffel(i, j, k):
        return scale(
            add(
                derivative(i, k, (j,)),
                derivative(i, j, (k,)),
                scale(derivative(j, k, (i,)), -1),
            ),
            s.Rational(1, 2),
        )

    @cache
    def riemann(i, j, k, l):
        first = scale(
            add(
                derivative(i, l, (j, k)),
                derivative(j, k, (i, l)),
                scale(derivative(i, k, (j, l)), -1),
                scale(derivative(j, l, (i, k)), -1),
            ),
            s.Rational(1, 2),
        )
        terms = []
        for r in range(N):
            for v in range(N):
                if inverse[r][v] != ZERO:
                    terms.append(
                        mul(
                            inverse[r][v],
                            add(
                                mul(christoffel(r, j, k), christoffel(v, i, l)),
                                scale(
                                    mul(christoffel(r, j, l), christoffel(v, i, k)), -1
                                ),
                            ),
                        )
                    )
        return add(first, *terms)

    Ric = [
        [
            add(
                *(
                    mul(inverse[i][k], riemann(i, j, k, l))
                    for i in range(N)
                    for k in range(N)
                )
            )
            for l in range(N)
        ]
        for j in range(N)
    ]
    R = add(*(mul(inverse[i][j], Ric[i][j]) for i in range(N) for j in range(N)))
    Ric2 = trace_square(matrix_product(inverse, Ric))
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    RB = [[riemann(i, j, k, l) for k, l in pairs] for i, j in pairs]
    GB = [
        [
            add(
                mul(inverse[i][k], inverse[j][l]),
                scale(mul(inverse[i][l], inverse[j][k]), -1),
            )
            for k, l in pairs
        ]
        for i, j in pairs
    ]
    Riem2 = scale(trace_square(matrix_product(GB, RB)), 4)
    return {
        name: s.factor(a ** (dim + 1) * row[3])
        for name, row in (
            ("R_old", R),
            ("R_squared", mul(R, R)),
            ("Ricci_squared", Ric2),
            ("Riemann_squared", Riem2),
        )
    }


def spatial_operator(expression):
    expression = s.expand(expression - expression.subs(p, 0))
    out = 0
    for r in range(3):
        term = s.diff(expression, Dtime[r])
        for _ in range(r):
            term = -time_derivative(term)
        out += term
    return s.factor(out)


def channel(dim, name):
    if name == "tensor":
        D = s.zeros(dim)
        D[0, 0] = 1
        D[1, 1] = -1
    elif name == "vector":
        D = s.zeros(dim)
        D[0, dim - 1] = D[dim - 1, 0] = 1
    elif name == "scalar":
        D = s.diag(*([-1] * (dim - 1) + [dim - 1]))
    else:
        raise ValueError(name)
    return D


@pytest.mark.parametrize("dimension", (3, 4, 5, 6))
@pytest.mark.parametrize("name", ("tensor", "vector", "scalar", "noncommuting"))
def test_literal_full_metric_second_variation_and_continued_curvatures(dimension, name):
    if name == "noncommuting":
        D = channel(dimension, "tensor") + channel(dimension, "vector") / 3
        G = channel(dimension, "scalar") / 7 - channel(dimension, "vector") / 5
        assert D * G != G * D
    else:
        D = channel(dimension, name)
        G = D
    T = s.trace(D * G)
    V = (D[:, -1].T * G[:, -1])[0]
    W = D[-1, -1] * G[-1, -1]
    source, expected = matching.hessians(T, V, W)
    actual = literal_curvatures(D, G)
    mapping = {
        A[0]: jets.a,
        A[1]: jets.a * s.diff(jets.a, jets.t),
        A[2]: jets.a * s.diff(jets.a * s.diff(jets.a, jets.t), jets.t),
        p: jets.p,
        Gtime[0]: source[0],
        Gtime[1]: jets.a * source[1],
        Gtime[2]: jets.a**2 * source[2] + jets.a * s.diff(jets.a, jets.t) * source[1],
    }
    for key, value in actual.items():
        operator = spatial_operator(value)
        difference = s.factor(
            operator.subs(mapping) / jets.a - expected[key].subs(geometry.d, dimension)
        )
        assert difference == 0, (dimension, name, key, difference)


@cache
def continued_functions():
    rows = density.invariant_spatial_coefficients()
    inv, g, pole = matching.fixed_pole()
    power = sum(v * c for v, c in zip(rows[0, 0, 2], inv))
    logs = [
        sum(sum(rows[j, r, 4 - j][k] for j in range(r, 5)) * inv[k] for k in range(3))
        for r in range(3)
    ]
    functions = s.lambdify(
        (geometry.d, jets.t, jets.p, jets.m, *inv), [power, *logs], "mpmath", cse=True
    )
    counterterm = s.lambdify(
        (geometry.d, jets.t, jets.p, jets.m, *inv),
        [s.diff(pole, x) for x in g],
        "mpmath",
        cse=True,
    )
    _inv, _g, ell, finite = matching.finite_input()
    expected = s.lambdify(
        (jets.t, jets.p, jets.m, *inv),
        [v.subs(ell, 0) for v in finite],
        "mpmath",
        cse=True,
    )
    return functions, counterterm, expected


def continued_radial_finite_values(epsilon, time, transfer, invariants):
    fn, ct, _expected = continued_functions()
    dimension = 3 - 2 * epsilon
    mass = mp.mpf(1000)
    power, *logs = fn(dimension, time, transfer, mass, *invariants)
    counter = ct(dimension, time, transfer, mass, *invariants)
    measure = (
        2
        * mp.pi ** (dimension / 2)
        / mp.gamma(dimension / 2)
        / (2 * mp.pi) ** dimension
    )
    # Exactly the original MSbar convention, checked separately by its Beta integrals.
    scale_ms = mp.exp(mp.euler * epsilon) * (mass * mass / (4 * mp.pi)) ** epsilon
    out = []
    for r in range(3):
        radial = logs[r] * mass ** (-2 * epsilon) / (2 * epsilon)
        if r == 0:
            radial += power * mass ** (2 - 2 * epsilon) / (2 * epsilon - 2)
        out.append(measure * scale_ms * radial - counter[r] / (64 * mp.pi**2 * epsilon))
    return out


def direct_dimension_constant(time, transfer, invariants, radius="0.005", points=32):
    epsradius = mp.mpf(radius)
    values = [
        continued_radial_finite_values(
            epsradius * mp.exp(2j * mp.pi * i / points), time, transfer, invariants
        )
        for i in range(points)
    ]
    return [mp.fsum(row[r] for row in values) / points for r in range(3)]


@pytest.mark.parametrize("time", ("-0.5", "0", "0.25", "0.5"))
@pytest.mark.parametrize("transfer", ("13", "1000", "1000000"))
@pytest.mark.parametrize(
    "invariants",
    ((1, 0, 0), (1, s.Rational(1, 2), 0), (1, s.Rational(2, 3), s.Rational(2, 3))),
)
def test_independent_continued_radial_pole_subtraction_finite_part(
    time, transfer, invariants
):
    with mp.workdps(135):
        time, transfer = map(mp.mpf, (time, transfer))
        inv = tuple(mp.mpf(str(s.N(v, 140))) for v in invariants)
        actual = direct_dimension_constant(time, transfer, inv)
        expected = continued_functions()[2](time, transfer, 1000, *inv)
        for a, b in zip(actual, expected):
            assert abs(a - b) < mp.mpf("1e-55") * (1 + abs(b))


def test_finite_radial_constant_two_dimension_Cauchy_resolutions():
    with mp.workdps(135):
        inv = (mp.mpf(1), mp.mpf("0.4"), mp.mpf("0.2"))
        a = direct_dimension_constant(mp.mpf("0.25"), mp.mpf("1000"), inv, "0.005", 32)
        b = direct_dimension_constant(mp.mpf("0.25"), mp.mpf("1000"), inv, "0.01", 40)
        assert max(abs(x - y) / (1 + abs(x)) for x, y in zip(a, b)) < mp.mpf("1e-55")


@pytest.mark.parametrize("module", (geometry, jets, density, matching))
def test_complete_actual_dimensional_core(module):
    data = module.data()
    for name, value in data["checks"].items():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.cancel(v) == 0 for v in entries), name
    assert all(bool(v) for v in data["gates"].values())


def test_evanescent_odd_endpoint_and_scalar_basis_cannot_be_frozen():
    rows = density.invariant_dimension_jet()
    assert all(v[0] == 0 for v in rows[1, 0, 3])
    assert any(v[1] != 0 for v in rows[1, 0, 3])
    raw = density.first_dimension_jet()["scalar", 0, 0, 4][1]
    fixed = sum(
        v[1] * w for v, w in zip(rows[0, 0, 4], (1, s.Rational(2, 3), s.Rational(2, 3)))
    )
    assert s.factor(raw - fixed) != 0


def test_Euler_dimension_derivative_and_volume_log_not_deleted():
    inv, g, pole = matching.fixed_pole()
    _g, rows = matching.hessians(*inv)
    assert s.factor(rows["Euler"].subs(geometry.d, 3)) == 0
    assert s.factor(s.diff(rows["Euler"], geometry.d).subs(geometry.d, 3)) != 0
    _inv, _g, _ell, finite = matching.finite_input()
    assert s.factor(finite[1] - s.diff(finite[2], jets.t)) == 0
    physical_second = s.diff(pole, g[2]).subs(geometry.d, 3) / (32 * s.pi**2)
    wrong = [
        finite[r]
        - s.log(jets.a) * s.diff(pole, g[r]).subs(geometry.d, 3) / (32 * s.pi**2)
        for r in range(3)
    ]
    assert s.factor(wrong[1] - s.diff(wrong[2], jets.t) - jets.H * physical_second) == 0
    assert physical_second != 0


@pytest.mark.parametrize(
    "time", (s.Rational(-1, 2), 0, s.Rational(1, 4), s.Rational(1, 2))
)
@pytest.mark.parametrize("transfer", (0, 13, 1000, 10**9))
def test_fixed_finite_coefficient_norm_and_origin(time, transfer):
    with mp.workdps(90):
        values = continued_functions()[2](
            mp.mpf(str(time.evalf(95))) if isinstance(time, s.Expr) else time,
            transfer,
            1000,
            1,
            s.Rational(2, 3),
            s.Rational(2, 3),
        )
        bound = mp.mpf("1e5") * (1 + mp.mpf(transfer) ** 2) ** 2
        if transfer == 0:
            assert sum(abs(v) for v in values) == 0
        else:
            assert sum(abs(v) for v in values) < bound


def test_finite_UV_piece_not_complete_original_response_or_P8():
    assert matching.data()["gates"][
        "full_regulator_limit_and_current_assembly_not_claimed"
    ]
    assert matching.finite_coefficient_norm()["sum"] < 100000


@cache
def integer_dimension_geometry(dimension, name):
    coordinates = s.symbols("n0:" + str(dimension - 1), real=True)
    u, y = geometry.u, geometry.y
    k = s.Matrix([*coordinates, u])
    ell = -k + s.Matrix([0] * (dimension - 1) + [y])
    ell2 = 1 - 2 * u * y + y * y
    D = channel(dimension, name)
    norm = s.trace(D * D)
    identity = s.eye(dimension)
    P = identity - k * k.T
    Q = identity - ell * ell.T / ell2
    # Direct field-strength contraction, no candidate magnetic matrix.
    M = s.zeros(dimension)
    for c in range(dimension):
        Fk = s.Matrix(
            dimension,
            dimension,
            lambda i, j, c=c: k[i] * int(c == j) - k[c] * int(i == j),
        )
        Fl = s.Matrix(
            dimension,
            dimension,
            lambda i, j, c=c: ell[i] * int(c == j) - ell[c] * int(i == j),
        )
        M -= Fk.T * D * Fl
    rows = {
        "00": s.trace(P * D * Q * D),
        "01": s.trace(P * M * Q * D),
        "10": s.trace(P * D * Q * M.T),
        "11": s.trace(P * M * Q * M.T),
        "TL": (ell.T * D * P * D * ell)[0] / ell2,
        "LT": (k.T * D * Q * D * k)[0],
        "LL": (k.T * D * ell)[0] ** 2 / ell2,
    }
    actual = {}
    for key, value in rows.items():
        numerator, denominator = s.fraction(s.cancel(value / norm))
        assert not any(denominator.has(v) for v in coordinates)
        polynomial = s.Poly(numerator, *coordinates)
        total = 0
        for powers, c in polynomial.terms():
            if any(v % 2 for v in powers):
                continue
            degree = sum(powers) // 2
            moment = (
                (1 - u * u) ** degree
                * s.prod(s.rf(s.Rational(1, 2), v // 2) for v in powers)
                / s.rf(s.Rational(dimension - 1, 2), degree)
            )
            total += c * moment
        actual[key] = s.factor(total / denominator)
    return actual


@pytest.mark.parametrize("dimension", (3, 4, 5, 6))
@pytest.mark.parametrize("name", ("tensor", "vector", "scalar"))
def test_literal_integer_dimensional_field_strength_and_all_angular_geometry(
    dimension, name
):
    literal = integer_dimension_geometry(dimension, name)
    proposed = geometry.contractions(name)
    for key, value in literal.items():
        assert s.factor(value - proposed[key].subs(geometry.d, dimension)) == 0, (
            dimension,
            name,
            key,
        )


@pytest.mark.parametrize("dimension", (3, 4, 5, 6))
@pytest.mark.parametrize("power", range(5))
def test_exact_grouped_transverse_radius_moments(dimension, power):
    u = geometry.u
    actual = geometry.azimuth(geometry.Z ** (2 * power)).subs(geometry.d, dimension)
    expected = (
        (1 - u * u) ** power
        * s.rf(s.Rational(dimension - 3, 2), power)
        / s.rf(s.Rational(dimension - 1, 2), power)
    )
    assert s.factor(actual - expected) == 0


def test_grouped_radius_is_not_a_signed_cartesian_coordinate():
    with pytest.raises(ValueError, match="even power"):
        geometry.azimuth(geometry.Z)


@cache
def full_dimensional_WKB_coefficient(kind, order):
    if order == 0:
        return s.Integer(1)
    u, z = wkb.u, wkb.z
    U = continuation.coefficients(kind)["U"]
    lam = wkb.background()["lambda"]

    def time_slope(value):
        return s.factor(
            s.diff(value, u) + wkb.background()["z_prime"] * s.diff(value, z)
        )

    previous = [full_dimensional_WKB_coefficient(kind, n) for n in range(order)]
    inverse = [s.Integer(1)]
    for n in range(1, order):
        inverse.append(
            s.factor(-sum(previous[j] * inverse[n - j] for j in range(1, n + 1)))
        )
    DS = [s.Integer(0)] + [
        time_slope(previous[n]) - 2 * n * lam * previous[n] for n in range(1, order)
    ]
    rate = [lam] + [
        s.factor(sum(DS[j] * inverse[n - j] for j in range(1, n + 1)))
        for n in range(1, order)
    ]
    value = (
        -sum(previous[j] * previous[order - j] for j in range(1, order))
        - time_slope(rate[-1]) / 2
        + (order - 1) * lam * rate[-1]
        + sum(rate[j] * rate[order - 1 - j] for j in range(order)) / 4
        - (U if order == 1 else 0)
    )
    return s.factor(value / 2)


@cache
def dimensional_coefficient_function(kind):
    expressions = []
    for order in range(1, 5):
        value = full_dimensional_WKB_coefficient(kind, order)
        expressions.extend((value, s.diff(value, wkb.u), s.diff(value, wkb.z)))
    return s.lambdify(
        (wkb.u, wkb.z, continuation.local.dimension), expressions, "mpmath", cse=True
    )


@cache
def amplitude_function():
    a, mass = s.symbols("a mass")
    point = {"a": a}
    args = [a, mass]
    for key in ("kt", "kl", "lt", "ll"):
        point[key] = s.symbols("f_" + key + " p_" + key + " o_" + key)
        args.extend(point[key])
    expressions = sectors.scalar_amplitudes(point, mass)
    return s.lambdify(args, list(expressions.values()), "mpmath", cse=True)


def scalar_values(point, mass):
    args = [point["a"], mass]
    for key in ("kt", "kl", "lt", "ll"):
        args.extend(point[key])
    return amplitude_function()(*args)


def scaled_point(u, x, k, l, phase_sign=-1, orders=4, dimension=3):
    a = (1 + u * u) ** 2
    H = 4 * u / (1 + u * u)
    point = {"a": a}
    frequencies = {}
    for prefix, v in (("k", k), ("l", l)):
        omega = mp.sqrt((v.T * v)[0] / a**2 + 1000**2 * x * x)
        z = 1 - 1000**2 * x * x / omega**2
        domega = -H * z * omega
        dz = -2 * H * z * (1 - z)
        for suffix, kind in (("t", "transverse"), ("l", "longitudinal")):
            coeff = dimensional_coefficient_function(kind)(u, z, dimension)
            W, dW = omega, domega
            for q in range(1, orders + 1):
                P, du, dzP = coeff[3 * (q - 1) : 3 * q]
                W += P * x ** (2 * q) * omega ** (1 - 2 * q)
                dW += x ** (2 * q) * (
                    (du + dzP * dz) * omega ** (1 - 2 * q)
                    + (1 - 2 * q) * P * omega ** (-2 * q) * domega
                )
            f = 1 / mp.sqrt(2 * W)
            rate = (
                (dimension - 2) * H / 2
                if suffix == "t"
                else ((dimension - 2) / 2 + z) * H
            )
            p = (phase_sign * 1j * W - x * (dW / (2 * W) + rate)) * f
            point[prefix + suffix] = (f, p, omega)
            frequencies[prefix + suffix] = W
    return point, frequencies


def four_inverse_phases(freq):
    return [
        1 / (freq[a] + freq[b])
        for a, b in (("kt", "lt"), ("kt", "ll"), ("kl", "lt"), ("kl", "ll"))
    ]


GEOS = ("00", "01", "10", "11", "TL", "LT", "LL")
PAIRING = ((0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0), (2, 2, 1), (3, 3, 2), (4, 4, 3))


def mode_endpoint_rows(x, time, angle, transfer, dimension, samples=20):
    n = mp.matrix([mp.sqrt(1 - angle * angle), 0, angle])
    k = n
    ell = -n + mp.matrix([0, 0, x * transfer])
    detector, df = scaled_point(time, x, k, ell, 1, 4, dimension)
    detector_amps = scalar_values(detector, 1000 * x)
    gd = four_inverse_phases(df)
    clock_radius = mp.mpf("1e-5")
    phases = [mp.exp(2j * mp.pi * i / samples) for i in range(samples)]
    source_amps = []
    source_g = []
    for phase in phases:
        source, freq = scaled_point(
            time + clock_radius * phase, x, k, ell, -1, 4, dimension
        )
        source_amps.append(scalar_values(source, 1000 * x))
        source_g.append(four_inverse_phases(freq))
    ampjet = [
        [
            mp.fsum(source_amps[i][a] * phases[i] ** (-r) for i in range(samples))
            / samples
            / clock_radius**r
            for r in range(5)
        ]
        for a in range(5)
    ]
    gjet = [
        [
            mp.fsum(source_g[i][a] * phases[i] ** (-r) for i in range(samples))
            / samples
            / clock_radius**r
            for r in range(5)
        ]
        for a in range(4)
    ]
    out = {}
    for geo, (left, right, sector) in zip(GEOS, PAIRING):
        g = gjet[sector]
        for r in range(5):
            row = [
                ampjet[right][d - r] / mp.factorial(r) if d >= r else mp.mpc(0)
                for d in range(5)
            ]
            for j in range(5):
                if r <= j:
                    out[geo, j, r] = (
                        (-1j) * 1j**j * gd[sector] * detector_amps[left] * row[0]
                    )
                row = [
                    (d + 1) * mp.fsum(g[h] * row[d + 1 - h] for h in range(d + 2))
                    for d in range(len(row) - 1)
                ]
    return out


@cache
def complete_dimensional_endpoint_function():
    rows = jets.endpoint_products()
    keys = [(key, degree) for key, row in rows.items() for degree in range(len(row))]
    fn = s.lambdify(
        (geometry.d, jets.t, jets.u, jets.p, jets.m),
        [rows[key][degree] for key, degree in keys],
        "mpmath",
        cse=True,
    )
    return keys, fn


@pytest.mark.parametrize("dimension", ("2.9", "3.1", "4"))
@pytest.mark.parametrize("time", ("-0.5", "0.25"))
@pytest.mark.parametrize("transfer", ("13", "1000000"))
def test_independent_full_four_order_dimensional_modes_all_source_jets(
    dimension, time, transfer
):
    with mp.workdps(135):
        dimension, time, transfer = map(mp.mpf, (dimension, time, transfer))
        angle = mp.mpf("0.37")
        radius = mp.mpf("1e-7") / (1 + transfer / 1000)
        points = 24
        phases = [mp.exp(2j * mp.pi * i / points) for i in range(points)]
        values = [
            mode_endpoint_rows(radius * z, time, angle, transfer, dimension)
            for z in phases
        ]
        keys, fn = complete_dimensional_endpoint_function()
        expected = fn(dimension, time, angle, transfer, 1000)
        for (key, degree), target in zip(keys, expected):
            actual = -mp.im(
                mp.fsum(values[i][key] * phases[i] ** (-degree) for i in range(points))
                / points
                / radius**degree
            )
            assert abs(actual - target) < mp.mpf("1e-32") * (1 + abs(target)), (
                key,
                degree,
                actual,
                target,
            )


@pytest.mark.parametrize("kind", ("transverse", "longitudinal"))
@pytest.mark.parametrize("order", range(1, 5))
def test_full_dimensional_reference_reduces_to_actual_four_order_physical_WKB(
    kind, order
):
    assert (
        s.factor(
            full_dimensional_WKB_coefficient(kind, order).subs(
                continuation.local.dimension, 3
            )
            - physical_series.coefficient(kind, order)
        )
        == 0
    )


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_all_scoped_dimensional_spatial_finite_residuals(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(v == 0 for v in entries)


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_all_unsupported_dimensional_spatial_finite_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_exact_dimensional_spatial_finite_audit_counts():
    assert len(audit.residuals()) == 73
    assert audit.scalar_entry_count() == 379
    assert len(audit.gates()) == 41
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 134


def test_original_frontier_not_closed_by_dimensional_UV_finite_coefficient():
    assert audit.frontier() == audit.previous.frontier()
    assert len(audit.matching()) == len(audit.previous.matching()) + 1
    assert audit.gates()["full_dimensional_limit_and_response_assembly_still_required"]
