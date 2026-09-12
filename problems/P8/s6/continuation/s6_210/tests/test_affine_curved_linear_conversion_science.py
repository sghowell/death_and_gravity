"""Independent actual full-W8 second-grade and proper-time coefficient checks."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_cd_metric_noise import stress
from p8_vacuum_affine_curved_linear_conversion import (
    angular,
    audit,
    bounds,
    density,
    jets,
)
from p8_vacuum_affine_spatial_symbol import sectors
from p8_vacuum_affine_subleading_band_conversion import flat as old_flat
from p8_vector_hadamard import series
from p8_vector_state import wkb


def cross(a, b):
    return mp.matrix(
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    )


def real_frame(n):
    index = min(range(3), key=lambda j: abs(n[j]))
    axis = mp.matrix([1 if j == index else 0 for j in range(3)])
    first = axis - n * (n.T * axis)[0]
    first /= mp.norm(first)
    return first, cross(n, first)


def joint_frame(k, p):
    magnitude = mp.norm(k)
    n0 = -k / magnitude
    ell = -k + p
    rho = mp.sqrt((ell.T * ell)[0])
    n = ell / rho
    skew = n * n0.T - n0 * n.T
    R = mp.eye(3) + skew + skew * skew / (1 + (n.T * n0)[0])
    first, second = real_frame(n0)
    return ell, rho, n, R, R * first, R * second


@cache
def coefficient_function(kind):
    expressions = []
    for order in range(1, 5):
        P = series.coefficient(kind, order)
        expressions.extend((P, s.diff(P, wkb.u), s.diff(P, wkb.z)))
    return s.lambdify((wkb.u, wkb.z), expressions, "mpmath", cse=True)


def W_data(kind, u, ell):
    a = (1 + u * u) ** 2
    H = 4 * u / (1 + u * u)
    omega = mp.sqrt(1000**2 + (ell.T * ell)[0] / a**2)
    z = 1 - 1000**2 / omega**2
    omega_prime = -H * z * omega
    z_prime = -2 * H * z * (1 - z)
    values = coefficient_function(kind)(u, z)
    W = omega
    derivative = omega_prime
    for order in range(1, 5):
        P, du, dz = values[3 * (order - 1) : 3 * order]
        W += P * omega ** (1 - 2 * order)
        derivative += (du + dz * z_prime) * omega ** (1 - 2 * order) + (
            1 - 2 * order
        ) * P * omega ** (-2 * order) * omega_prime
    return a, H, omega, z, W, derivative


def physical_modes(u, ell, n, e1, e2, phase_sign=-1):
    result = []
    for kind, e in (("transverse", e1), ("transverse", e2), ("longitudinal", n)):
        a, H, omega, z, W, derivative = W_data(kind, u, ell)
        f = 1 / mp.sqrt(2 * W)
        rate = H / 2 if kind == "transverse" else (mp.mpf("0.5") + z) * H
        momentum = (phase_sign * 1j * W - derivative / (2 * W) - rate) * f
        if kind == "transverse":
            E = momentum * e
            B = -phase_sign * 1j * cross(ell, e) * f / a
            V0 = mp.mpf(0)
            V = 1000 * f * e
        else:
            E = 1000 / omega * momentum * n
            B = mp.zeros(3, 1)
            V0 = phase_sign * 1j * mp.sqrt((ell.T * ell)[0]) / (a * omega) * momentum
            V = omega * f * n
        result.append((mp.matrix([*E, *B, V0, *V]), W))
    return result


@cache
def complete_stress_matrices():
    def real_number(x):
        n, d = s.fraction(x)
        return mp.mpf(int(n)) / int(d)

    return [
        mp.matrix([[real_number(M[i, j]) for j in range(10)] for i in range(10)])
        for M in stress.data()["quadratic_matrices"].values()
    ]


def contracted_spatial_matrix(tensor):
    matrices = complete_stress_matrices()
    return sum(
        (
            tensor[i, j] * matrices[4 * (i + 1) + j + 1]
            for i in range(3)
            for j in range(3)
        ),
        mp.zeros(10),
    )


@cache
def geometric_function():
    kv = s.Matrix(s.symbols("k0:3"))
    lv = s.Matrix(s.symbols("l0:3"))
    D = s.Matrix(3, 3, s.symbols("D0:9"))
    G = s.Matrix(3, 3, s.symbols("G0:9"))
    expressions = sectors.geometric_factors(kv, lv, D, G)
    return s.lambdify(
        [*kv, *lv, *D, *G], list(expressions.values()), "mpmath", cse=True
    )


@cache
def amplitude_function():
    a, m = s.symbols("a m")
    point = {"a": a}
    args = [a, m]
    for key in ("kt", "kl", "lt", "ll"):
        point[key] = s.symbols("f_" + key + " p_" + key + " o_" + key)
        args.extend(point[key])
    expr = sectors.scalar_amplitudes(point, m)
    return s.lambdify(args, list(expr.values()), "mpmath", cse=True)


def scalar_values(point, mass):
    args = [point["a"], mass]
    for key in ("kt", "kl", "lt", "ll"):
        args.extend(point[key])
    return amplitude_function()(*args)


def projected_products(k, l, D, G, source, detector, x):
    c = geometric_function()(*k, *l, *D, *G)
    A = scalar_values(detector, 1000 * x)
    B = scalar_values(source, 1000 * x)
    return [
        A[0] * B[0] * c[0]
        + A[0] * B[1] * c[1]
        + A[1] * B[0] * c[2]
        + A[1] * B[1] * c[3],
        A[2] * B[2] * c[4],
        A[3] * B[3] * c[5],
        A[4] * B[4] * c[6],
    ]


def scaled_point(u, x, k, l, phase_sign=-1, orders=4):
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
            coeff = coefficient_function(kind)(u, z)
            W, dW = omega, domega
            for q in range(1, orders + 1):
                P, du, dzP = coeff[3 * (q - 1) : 3 * q]
                W += P * x ** (2 * q) * omega ** (1 - 2 * q)
                dW += x ** (2 * q) * (
                    (du + dzP * dz) * omega ** (1 - 2 * q)
                    + (1 - 2 * q) * P * omega ** (-2 * q) * domega
                )
            f = 1 / mp.sqrt(2 * W)
            rate = H / 2 if suffix == "t" else (mp.mpf("0.5") + z) * H
            p = (phase_sign * 1j * W - x * (dW / (2 * W) + rate)) * f
            point[prefix + suffix] = (f, p, omega)
            frequencies[prefix + suffix] = W
    return point, frequencies


def scaled_vectors(momentum, e1, e2, point, prefix, x, phase_sign=-1):
    rho = mp.sqrt((momentum.T * momentum)[0])
    n = momentum / rho
    f, p, _omega = point[prefix + "t"]
    fl, pl, ol = point[prefix + "l"]
    out = []
    for e in (e1, e2):
        E = p * e
        B = -phase_sign * 1j * cross(momentum, e) * f / point["a"]
        out.append(mp.matrix([*E, *B, 0, *(1000 * x * f * e)]))
    E = 1000 * x * pl * n / ol
    out.append(
        mp.matrix(
            [
                *E,
                0,
                0,
                0,
                phase_sign * 1j * rho * pl / (point["a"] * ol),
                *(ol * fl * n),
            ]
        )
    )
    return out


def full_products(k, l, e1, e2, f1, f2, D, G, source, detector, x):
    sk = scaled_vectors(k, e1, e2, source, "k", x)
    sl = scaled_vectors(l, f1, f2, source, "l", x)
    dk = scaled_vectors(k, e1, e2, detector, "k", x, 1)
    dl = scaled_vectors(l, f1, f2, detector, "l", x, 1)
    MD, MG = contracted_spatial_matrix(D), contracted_spatial_matrix(G)
    return [
        (dk[i].T * MD * dl[j])[0] * (sk[i].T * MG * sl[j])[0]
        for i in range(3)
        for j in range(3)
    ]


def grouping(values):
    return [
        sum(values[3 * i + j] for i in range(2) for j in range(2)),
        values[2] + values[5],
        values[6] + values[7],
        values[8],
    ]


def four_inverse_phases(freq):
    return [
        1 / (freq[a] + freq[b])
        for a, b in (("kt", "lt"), ("kt", "ll"), ("kl", "lt"), ("kl", "ll"))
    ]


def tensors():
    D = mp.matrix([[2, 1, -1], [1, -3, 2], [-1, 2, 1]]) / 7
    G = mp.matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 17
    return D, G


def mp_number(value):
    n, d = s.fraction(value)
    return mp.mpf(int(n)) / int(d)


def analytic_frame(v, n0):
    n = v / mp.sqrt((v.T * v)[0])
    skew = n * n0.T - n0 * n.T
    R = mp.eye(3) + skew + skew * skew / (1 + (n.T * n0)[0])
    e1, e2 = real_frame(n0)
    return R * e1, R * e2


def flat_point(x, k, l, sign):
    point = {"a": mp.mpf(1)}
    freq = {}
    for prefix, v in (("k", k), ("l", l)):
        omega = mp.sqrt((v.T * v)[0] + 1000**2 * x * x)
        f = 1 / mp.sqrt(2 * omega)
        for suffix in ("t", "l"):
            point[prefix + suffix] = (f, sign * 1j * omega * f, omega)
            freq[prefix + suffix] = omega
    return point, freq


def fixtures(index):
    if index == 0:
        return tensors()
    D = mp.matrix([[1, 2, 0], [2, -3, 1], [0, 1, 2]]) / 11
    G = mp.matrix([[2, 1, -2], [1, 1, 2], [-2, 2, -3]]) / 13
    return D, G


def leading_point(n, D, G, a):
    return (
        4 * mp.fsum((D * G)[i, i] for i in range(3))
        - 4 * (n.T * (D * G + G * D) * n)[0]
        + 3 * (n.T * D * n)[0] * (n.T * G * n)[0]
    ) / (8 * a)


def current_point(x, time, n, P, D, G, method="projected", flat_comparison=False):
    k = n + x * P / 2
    ell = -n + x * P / 2
    a = (1 + time * time) ** 2
    argument = a * x if flat_comparison else x
    if flat_comparison:
        source, freq = flat_point(argument, k, ell, -1)
        detector, _ = flat_point(argument, k, ell, 1)
    else:
        source, freq = scaled_point(time, x, k, ell, -1, 4)
        detector, _ = scaled_point(time, x, k, ell, 1, 4)
    if method == "projected":
        values = projected_products(k, ell, D, G, source, detector, argument)
    elif method == "full":
        e1, e2 = analytic_frame(k, n)
        f1, f2 = analytic_frame(ell, -n)
        values = grouping(
            full_products(k, ell, e1, e2, f1, f2, D, G, source, detector, argument)
        )
    else:
        raise ValueError("Unknown independent physical route")
    g = four_inverse_phases(freq)
    value = mp.fsum(values[i] * g[i] for i in range(4))
    return value / a if flat_comparison else value


@cache
def actual_point_coefficients(
    time, direction, transfer_scale, tensor_index, method="projected"
):
    with mp.workdps(125):
        t = mp.mpf(time)
        n = mp.matrix(direction)
        n /= mp.norm(n)
        P = mp.mpf(transfer_scale) * mp.matrix([11, -7, 5])
        D, G = fixtures(tensor_index)
        radius = mp.mpf("1e-6") / (1 + mp.norm(P) / 1000)
        points = 24
        phases = [mp.exp(2j * mp.pi * i / points) for i in range(points)]
        actual = [current_point(radius * z, t, n, P, D, G, method) for z in phases]
        flat_values = [
            current_point(radius * z, t, n, P, D, G, method, True) for z in phases
        ]
        return {
            (kind, d): mp.re(
                mp.fsum(values[i] * phases[i] ** (-d) for i in range(points))
                / points
                / radius**d
            )
            for kind, values in (("actual", actual), ("flat", flat_values))
            for d in (0, 2)
        }


@pytest.mark.parametrize("module", (jets, density, angular, bounds))
def test_complete_actual_curved_linear_exact_core(module):
    data = module.data()
    for name, value in data["checks"].items():
        if isinstance(value, s.MatrixBase):
            assert all(s.cancel(v) == 0 for v in value), name
        else:
            assert s.cancel(value) == 0, name
    assert all(data["gates"].values())


@pytest.mark.parametrize("time", ("-0.5", "0", "0.25", "0.49", "0.5"))
@pytest.mark.parametrize("direction", ((3, 0, 4), (1, 2, 2), (0, 0, 1)))
@pytest.mark.parametrize("transfer_scale", ("0", "1", "1e6"))
@pytest.mark.parametrize("tensor_index", (0, 1))
def test_actual_full_W8_curved_source_value_second_grade(
    time, direction, transfer_scale, tensor_index
):
    with mp.workdps(110):
        t = mp.mpf(time)
        a = (1 + t * t) ** 2
        cL = -2 * (7 * t * t + 1) / (1 + t * t) ** 2
        n = mp.matrix(direction)
        n /= mp.norm(n)
        D, G = fixtures(tensor_index)
        values = actual_point_coefficients(
            time, direction, transfer_scale, tensor_index
        )
        curvature = -3 * a * cL * (n.T * D * n)[0] * (n.T * G * n)[0] / 8
        assert abs(values["actual", 2] - values["flat", 2] - curvature) < mp.mpf(
            "1e-35"
        )
        expected0 = leading_point(n, D, G, a)
        assert abs(values["actual", 0] - expected0) < mp.mpf("1e-43")
        assert abs(values["flat", 0] - expected0) < mp.mpf("1e-43")


@pytest.mark.parametrize("time", ("-0.5", "0.25", "0.5"))
@pytest.mark.parametrize("transfer_scale", ("1", "1e6"))
def test_independent_full_ten_readouts_all_nine_pairs_curved_and_mass_scaled_flat(
    time, transfer_scale
):
    with mp.workdps(110):
        projected = actual_point_coefficients(time, (3, 0, 4), transfer_scale, 0)
        full = actual_point_coefficients(time, (3, 0, 4), transfer_scale, 0, "full")
        for key in projected:
            assert abs(projected[key] - full[key]) < mp.mpf("1e-70")


def source_endpoint_rows(x, time, n, P, D, G, samples=16, averaged_channel=None):
    k = n + x * P / 2
    ell = -n + x * P / 2
    radius = mp.mpf("1e-5")
    detector, df = scaled_point(time, x, k, ell, 1, 4)
    gd = four_inverse_phases(df)
    phases = [mp.exp(2j * mp.pi * i / samples) for i in range(samples)]
    amplitudes = []
    inverse = []
    for phase in phases:
        source, freq = scaled_point(time + radius * phase, x, k, ell, -1, 4)
        if averaged_channel is None:
            amplitudes.append(projected_products(k, ell, D, G, source, detector, x))
        else:
            amplitudes.append(
                actual_averaged_products(
                    x, P[2], n[2], averaged_channel, source, detector
                )
            )
        inverse.append(four_inverse_phases(freq))
    result = mp.zeros(3, 3)
    for sector in range(4):
        g = [
            mp.fsum(inverse[i][sector] * phases[i] ** (-j) for i in range(samples))
            / samples
            / radius**j
            for j in range(3)
        ]
        amp = [
            mp.fsum(amplitudes[i][sector] * phases[i] ** (-j) for i in range(samples))
            / samples
            / radius**j
            for j in range(3)
        ]
        for r in range(3):
            row = [
                amp[d - r] / mp.factorial(r) if d >= r else mp.mpc(0) for d in range(3)
            ]
            for j in range(3):
                result[j, r] += (-1j) * 1j**j * gd[sector] * row[0]
                row = [
                    (d + 1) * mp.fsum(g[h] * row[d + 1 - h] for h in range(d + 2))
                    for d in range(len(row) - 1)
                ]
    return result


@pytest.mark.parametrize("time", ("-0.5", "-0.2", "0", "0.25", "0.5"))
@pytest.mark.parametrize("direction", ((3, 0, 4), (1, 2, 2), (0, 0, 1)))
@pytest.mark.parametrize("tensor_index", (0, 1))
def test_actual_complete_source_first_and_second_time_jets(
    time, direction, tensor_index
):
    with mp.workdps(115):
        t = mp.mpf(time)
        a = (1 + t * t) ** 2
        H = 4 * t / (1 + t * t)
        n = mp.matrix(direction)
        n /= mp.norm(n)
        P = mp.matrix([11, -7, 5])
        D, G = fixtures(tensor_index)
        rows = source_endpoint_rows(mp.mpf(0), t, n, P, D, G)
        f0 = leading_point(n, D, G, a)
        expected = (0, -a * a * f0 * H / 4, -a * a * f0 / 4)
        for r in range(3):
            assert abs(-mp.im(rows[2, r]) - expected[r]) < mp.mpf("1e-43")
        assert abs(mp.im(rows[1, 0])) < mp.mpf("1e-60")


@pytest.mark.parametrize("time", ("-0.5", "0", "0.25"))
@pytest.mark.parametrize("direction", ((3, 0, 4), (1, 2, 2)))
def test_actual_first_endpoint_degree_one_cancels_without_deleting_higher_slot(
    time, direction
):
    with mp.workdps(125):
        t = mp.mpf(time)
        n = mp.matrix(direction)
        n /= mp.norm(n)
        P = mp.matrix([11, -7, 5])
        D, G = fixtures(0)
        radius = mp.mpf("1e-6")
        points = 24
        phases = [mp.exp(2j * mp.pi * i / points) for i in range(points)]
        values = [source_endpoint_rows(radius * z, t, n, P, D, G) for z in phases]
        degree1 = mp.matrix(
            [
                -mp.im(
                    mp.fsum(values[i][1, r] * phases[i] ** (-1) for i in range(points))
                    / points
                    / radius
                )
                for r in range(2)
            ]
        )
        assert mp.norm(degree1) < mp.mpf("1e-43")
        if time == "0.25" and direction == (3, 0, 4):
            degree3 = -mp.im(
                mp.fsum(values[i][1, 0] * phases[i] ** (-3) for i in range(points))
                / points
                / radius**3
            )
            assert abs(degree3) > 1


@pytest.mark.parametrize(
    "time",
    (
        s.Rational(-1, 2),
        s.Rational(-1, 4),
        s.Integer(0),
        s.Rational(1, 4),
        s.Rational(1, 2),
    ),
)
@pytest.mark.parametrize(
    "p", (s.Integer(0), s.Integer(1), s.Integer(1000), s.Integer(10) ** 12)
)
@pytest.mark.parametrize("channel", ("tensor", "vector", "scalar"))
def test_all_actual_coefficient_norm_displays_and_origin(time, p, channel):
    inv = {
        "tensor": (1, 0, 0),
        "vector": (1, s.Rational(1, 2), 0),
        "scalar": (1, s.Rational(2, 3), s.Rational(2, 3)),
    }
    T, V, W = inv[channel]
    c3 = angular.cubic(time, p, T, V, W)
    c1 = angular.linear(time, p, 1000, T, V, W)
    with mp.workdps(85):
        value3 = abs(mp.mpf(str(s.N(c3, 85))))
        values1 = sum(abs(mp.mpf(str(s.N(value, 85)))) for value in c1)
        if p == 0:
            assert value3 == values1 == 0
        else:
            assert value3 < mp.mpf(p) / 100
            assert values1 < mp.mpf("1e4") * (1 + mp.mpf(p) ** 2) ** mp.mpf("1.5")


def test_actual_proper_time_green_identity_and_missing_first_derivative_control():
    t = s.symbols("t", real=True)
    a = (1 + t * t) ** 2
    bump = (1 - 4 * t * t) ** 4
    D = bump
    G = (1 + t * t) * bump
    full = lambda F: s.diff(a * s.diff(F, t), t)
    wrong = lambda F: a * s.diff(F, t, 2)
    left = s.integrate(
        s.expand(D * full(G) - G * full(D)), (t, -s.Rational(1, 2), s.Rational(1, 2))
    )
    missing = s.integrate(
        s.expand(D * wrong(G) - G * wrong(D)), (t, -s.Rational(1, 2), s.Rational(1, 2))
    )
    assert left == 0
    assert missing != 0


def test_coefficients_not_full_original_P8_closure():
    assert bounds.data()["gates"]["canonical_coefficients_not_cutoff_uniformity"]
    assert angular.data()["gates"]["no_finite_contact_or_covariant_matching_claim"]


@cache
def averaged_geometric_function(channel):
    u, y, values = old_flat.averaged_geometry(channel)
    return s.lambdify((u, y), list(values.values()), "mpmath", cse=True)


def actual_averaged_products(x, p, u, channel, source, detector):
    sigma = mp.sqrt(1 + p * u * x + p * p * x * x / 4)
    old_u = (u + p * x / 2) / sigma
    old_y = p * x / sigma
    c = list(averaged_geometric_function(channel)(old_u, old_y))
    c[1] *= sigma**2
    c[2] *= sigma**2
    c[3] *= sigma**4
    A = scalar_values(detector, 1000 * x)
    B = scalar_values(source, 1000 * x)
    return [
        A[0] * B[0] * c[0]
        + A[0] * B[1] * c[1]
        + A[1] * B[0] * c[2]
        + A[1] * B[1] * c[3],
        A[2] * B[2] * c[4],
        A[3] * B[3] * c[5],
        A[4] * B[4] * c[6],
    ]


def averaged_current_point(x, time, p, u, channel):
    n = mp.matrix([mp.sqrt(1 - u * u), 0, u])
    P = mp.matrix([0, 0, p])
    k = n + x * P / 2
    ell = -n + x * P / 2
    source, freq = scaled_point(time, x, k, ell, -1, 4)
    detector, _ = scaled_point(time, x, k, ell, 1, 4)
    pair = actual_averaged_products(x, p, u, channel, source, detector)
    g = four_inverse_phases(freq)
    return mp.fsum(pair[i] * g[i] for i in range(4))


def mp_channel(channel):
    D = old_flat.channel_tensor(channel)
    return mp.matrix(
        [[mp.mpf(str(s.N(D[i, j], 125))) for j in range(3)] for i in range(3)]
    )


@pytest.mark.parametrize("time", ("-0.5", "0.25", "0.5"))
@pytest.mark.parametrize("channel", ("tensor", "vector", "scalar"))
@pytest.mark.parametrize("transfer", ("13", "1000000"))
def test_actual_curved_centered_azimuth_reconstruction_against_all_physical_fields(
    time, channel, transfer
):
    with mp.workdps(115):
        t = mp.mpf(time)
        p = mp.mpf(transfer)
        u = mp.mpf("0.37")
        x = mp.mpf("1e-6") * (1 + mp.mpc(0, "0.4")) / (1 + p / 1000)
        P = mp.matrix([0, 0, p])
        D = mp_channel(channel)
        literal = mp.mpc(0)
        for i in range(8):
            phi = 2 * mp.pi * i / 8
            n = mp.matrix(
                [mp.sqrt(1 - u * u) * mp.cos(phi), mp.sqrt(1 - u * u) * mp.sin(phi), u]
            )
            literal += current_point(x, t, n, P, D, D, "full") / 8
        expected = averaged_current_point(x, t, p, u, channel)
        assert abs(literal - expected) < mp.mpf("1e-90")


@cache
def actual_integrated_linear_coefficients(time, transfer, channel, points=24):
    with mp.workdps(125):
        t = mp.mpf(time)
        p = mp.mpf(transfer)
        radius = mp.mpf("1e-6") / (1 + p / 1000)
        phases = [mp.exp(2j * mp.pi * i / points) for i in range(points)]
        nodes, weights = mp.gauss_quadrature(4, "legendre")
        P = mp.matrix([0, 0, p])
        D = mp_channel(channel)
        result = mp.matrix([0, 0, 0, 0])
        for z, w in zip(nodes, weights):
            u = (z + 1) / 2
            w /= 2
            values = [
                averaged_current_point(radius * phase, t, p, u, channel)
                for phase in phases
            ]
            f0 = mp.re(mp.fsum(values) / points)
            f2 = mp.re(
                mp.fsum(values[i] * phases[i] ** (-2) for i in range(points))
                / points
                / radius**2
            )
            n = mp.matrix([mp.sqrt(1 - u * u), 0, u])
            time_rows = source_endpoint_rows(
                mp.mpf(0), t, n, P, D, D, averaged_channel=channel
            )
            result[0] += -w * p * u * f0 / (4 * mp.pi**2)
            result[1] += (
                w * u * (p**3 * (3 - 5 * u * u) * f0 / 16 - p * f2 / 2) / (2 * mp.pi**2)
            )
            result[2] += -w * p * u * (-mp.im(time_rows[2, 1])) / (4 * mp.pi**2)
            result[3] += -w * p * u * (-mp.im(time_rows[2, 2])) / (4 * mp.pi**2)
        return result


@pytest.mark.parametrize("time", ("-0.5", "0", "0.25", "0.5"))
@pytest.mark.parametrize("transfer", ("13", "1000", "1000000"))
@pytest.mark.parametrize("channel", ("tensor", "vector", "scalar"))
def test_actual_full_curved_angular_integral_all_linear_source_jets(
    time, transfer, channel
):
    inv = {
        "tensor": (1, 0, 0),
        "vector": (1, s.Rational(1, 2), 0),
        "scalar": (1, s.Rational(2, 3), s.Rational(2, 3)),
    }
    with mp.workdps(115):
        T, V, W = inv[channel]
        t = s.Rational(time)
        p = s.Rational(transfer)
        expected = (angular.cubic(t, p, T, V, W),) + angular.linear(t, p, 1000, T, V, W)
        actual = actual_integrated_linear_coefficients(time, transfer, channel)
        for value, formula in zip(actual, expected):
            target = mp.mpf(str(s.N(formula, 115)))
            assert abs(value - target) < mp.mpf("1e-32")


@pytest.mark.parametrize("channel", ("tensor", "vector", "scalar"))
def test_actual_integrated_curved_coefficient_two_inverse_radius_resolutions(channel):
    with mp.workdps(115):
        a = actual_integrated_linear_coefficients("0.25", "1000", channel, 24)
        b = actual_integrated_linear_coefficients("0.25", "1000", channel, 32)
        assert mp.norm(a - b) < mp.mpf("1e-42")


@pytest.mark.parametrize("power", range(8))
def test_exact_four_node_angular_moments_through_degree_seven(power):
    u = s.symbols("u")
    sqrt30 = s.sqrt(30)
    roots = []
    for sign in (-1, 1):
        root = s.sqrt((15 + 2 * sign * sqrt30) / 35)
        weight = (18 - sign * sqrt30) / 36
        roots.extend(((root, weight), (-root, weight)))
    moment = sum(w * ((z + 1) / 2) ** power / 2 for z, w in roots)
    assert s.simplify(moment - s.integrate(u**power, (u, 0, 1))) == 0


def test_curvature_and_source_first_derivative_cannot_be_replaced_by_flat_values():
    t = s.Rational(1, 4)
    p = s.Integer(13)
    actual = angular.linear(t, p, 1000, 1, 0, 0)
    a, H, Hp, _c = jets.clock(t)
    missing_curvature = s.factor(
        actual[0] + p * a * (Hp + 2 * H * H) * 2 / (1024 * s.pi**2)
    )
    assert missing_curvature != actual[0]
    assert actual[1] != 0


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_all_named_scoped_actual_curved_coefficient_residuals(name):
    value = audit.residuals()[name]
    if isinstance(value, s.MatrixBase):
        assert all(v == 0 for v in value)
    else:
        assert value == 0


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_all_unsupported_actual_curved_coefficient_scope_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_exact_actual_curved_coefficient_audit_counts():
    assert len(audit.residuals()) == 49 and audit.scalar_entry_count() == 92
    assert len(audit.gates()) == 40 and all(audit.gates().values())
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 132


def test_original_frontier_not_closed_by_evaluated_cutoff_coefficient():
    assert len(audit.frontier()) == 9
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert (
        "NOT_FULL_ONE_BALL_CONTACT_COVARIANT_MATCHING_OR_V_G_B" in audit.ITEM["status"]
    )
