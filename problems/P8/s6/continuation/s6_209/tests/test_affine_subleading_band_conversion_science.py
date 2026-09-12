"""Independent original-shell, full-field flat and actual curved symmetry checks."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_cd_metric_noise import stress
from p8_vacuum_affine_spatial_symbol import extraction, sectors
from p8_vacuum_affine_subleading_band_conversion import (
    angular,
    audit,
    centered,
    flat,
    geometry,
)
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


@pytest.mark.parametrize("module", (geometry, angular, flat, centered))
def test_complete_exact_conversion_core(module):
    packet = module.data()
    for name, value in packet["checks"].items():
        if isinstance(value, s.MatrixBase):
            assert all(s.cancel(v) == 0 for v in value), name
        else:
            assert s.cancel(value) == 0, name
    assert all(packet["gates"].values())


@pytest.mark.parametrize("q", range(5))
@pytest.mark.parametrize("degree", range(9))
@pytest.mark.parametrize("epsilon", ("0.2", "0.001", "1e-12"))
def test_exact_original_grazing_shell_all_angular_powers(q, degree, epsilon):
    with mp.workdps(100):
        e = mp.mpf(epsilon)

        def radial(u):
            R = e * u + mp.sqrt(1 - e * e * (1 - u * u))
            return -mp.log(R) if q == 4 else (1 - R ** (4 - q)) / (4 - q)

        actual = 2 * mp.pi * mp.quad(lambda u: u**degree * radial(u), [-1, 0, e / 2])
        u = s.symbols("u")
        approximation = mp.mpf(0)
        if q < 4:
            approximation = mp.fsum(
                mp_number(angular.shape_coefficient(u**degree, u, q, h) / s.pi)
                * mp.pi
                * e**h
                for h in range(1, 5 - q)
            )
        bound = 8192 * (1 + degree * degree) * e ** (5 - q)
        assert abs(actual - approximation) < bound
        if degree == 0:
            assert actual > 0


@pytest.mark.parametrize("q", range(5))
@pytest.mark.parametrize("epsilon", ("0.2", "0.01", "1e-8"))
def test_independent_radial_cap_and_original_angular_shell(q, epsilon):
    with mp.workdps(90):
        e = mp.mpf(epsilon)
        cap = mp.quad(
            lambda r: mp.pi / e * r ** (2 - q) * ((r + e) ** 2 - 1), [1 - e, 1]
        )
        angular_integral = (
            2
            * mp.pi
            * mp.quad(
                lambda u: (
                    -mp.log(e * u + mp.sqrt(1 - e * e * (1 - u * u)))
                    if q == 4
                    else (1 - (e * u + mp.sqrt(1 - e * e * (1 - u * u))) ** (4 - q))
                    / (4 - q)
                ),
                [-1, 0, e / 2],
            )
        )
        assert abs(cap - angular_integral) < mp.mpf("1e-70") * e


def test_omitting_positive_grazing_changes_cubic_and_finite_shapes():
    u = s.symbols("u")
    wrong3 = 2 * s.pi * s.integrate(geometry.hemisphere_polynomial(0, 3, u), (u, -1, 0))
    wrong4 = (
        2 * s.pi * s.integrate(u * geometry.hemisphere_polynomial(0, 4, u), (u, -1, 0))
    )
    assert s.simplify(angular.shape_coefficient(1, u, 0, 3) - wrong3) == s.pi / 4
    assert s.simplify(angular.shape_coefficient(u, u, 0, 4) - wrong4) == s.pi / 24


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


@cache
def actual_flat_coefficients(channel, cosine, transfer, points=24):
    with mp.workdps(125):
        u = mp.mpf(cosine)
        P = mp.matrix([0, 0, mp.mpf(transfer)])
        D = mp.matrix(
            [
                [
                    mp_number(value)
                    if not value.has(s.sqrt(2), s.sqrt(6))
                    else mp.mpf(str(s.N(value, 125)))
                    for value in flat.channel_tensor(channel).row(i)
                ]
                for i in range(3)
            ]
        )
        radius = mp.mpf("1e-6")
        phases = [mp.exp(2j * mp.pi * j / points) for j in range(points)]
        values = []
        for phase in phases:
            x = radius * phase
            totals = [mp.mpc(0) for _ in range(5)]
            for azimuth in range(8):
                phi = 2 * mp.pi * azimuth / 8
                n = mp.matrix(
                    [
                        mp.sqrt(1 - u * u) * mp.cos(phi),
                        mp.sqrt(1 - u * u) * mp.sin(phi),
                        u,
                    ]
                )
                k = n
                ell = -n + x * P
                e1, e2 = real_frame(n)
                f1, f2 = analytic_frame(ell, -n)
                source, freq = flat_point(x, k, ell, -1)
                detector, _ = flat_point(x, k, ell, 1)
                pair = mp.fsum(
                    full_products(k, ell, e1, e2, f1, f2, D, D, source, detector, x)
                )
                g = 1 / (freq["kt"] + freq["lt"])
                for j in range(5):
                    totals[j] += (1, 0, -1, 0, 1)[j] * pair * g ** (j + 1) / 8
            values.append(totals)
        return {
            (j, d): mp.fsum(values[i][j] * phases[i] ** (-d) for i in range(points))
            / points
            / radius**d
            for j, d in extraction.slots()
        }


@pytest.mark.parametrize("channel", ("tensor", "vector", "scalar"))
@pytest.mark.parametrize("cosine", ("-0.8", "0.2", "0.91"))
@pytest.mark.parametrize("transfer", ("13", "-7"))
@pytest.mark.parametrize("slot", extraction.slots())
def test_full_nine_pair_flat_fields_all_angular_UV_coefficients(
    channel, cosine, transfer, slot
):
    with mp.workdps(115):
        j, d = slot
        u, m, p, rows = flat.coefficients()
        expected = mp_number(
            rows[channel, j][d].subs(
                {u: s.Rational(cosine), m: 1000, p: s.Rational(transfer)}
            )
        )
        actual = actual_flat_coefficients(channel, cosine, transfer)[slot]
        assert abs(actual - expected) < mp.mpf("1e-48") * max(mp.mpf(1), abs(expected))


def centered_endpoint_rows(x, reverse=False, method="projected", samples=16):
    n = mp.matrix([mp.mpf(3) / 5, 0, mp.mpf(4) / 5])
    if reverse:
        n = -n
    P = mp.matrix([11, -7, 5])
    k = n + x * P / 2
    ell = -n + x * P / 2
    e1, e2 = analytic_frame(k, n)
    f1, f2 = analytic_frame(ell, -n)
    time = mp.mpf("0.25")
    radius = mp.mpf("1e-5")
    D, G = tensors()
    detector, df = scaled_point(time, x, k, ell, 1, 4)
    gd = four_inverse_phases(df)
    phases = [mp.exp(2j * mp.pi * j / samples) for j in range(samples)]
    amplitudes = []
    inverse = []
    for phase in phases:
        source, freq = scaled_point(time + radius * phase, x, k, ell, -1, 4)
        if method == "projected":
            amplitudes.append(projected_products(k, ell, D, G, source, detector, x))
        elif method == "full":
            amplitudes.append(
                grouping(
                    full_products(k, ell, e1, e2, f1, f2, D, G, source, detector, x)
                )
            )
        else:
            raise ValueError("Unknown independent actual route")
        inverse.append(four_inverse_phases(freq))
    result = mp.zeros(5, 5)
    for sector in range(4):
        g = [
            mp.fsum(inverse[i][sector] * phases[i] ** (-j) for i in range(samples))
            / samples
            / radius**j
            for j in range(5)
        ]
        amp = [
            mp.fsum(amplitudes[i][sector] * phases[i] ** (-j) for i in range(samples))
            / samples
            / radius**j
            for j in range(5)
        ]
        for r in range(5):
            row = [
                amp[d - r] / mp.factorial(r) if d >= r else mp.mpc(0) for d in range(5)
            ]
            for j in range(5):
                result[j, r] += (-1j) * 1j**j * gd[sector] * row[0]
                row = [
                    (d + 1) * mp.fsum(g[h] * row[d + 1 - h] for h in range(d + 2))
                    for d in range(len(row) - 1)
                ]
    return result


@pytest.mark.parametrize("x", ("0.0001", "0.00000001", "complex"))
def test_actual_full_curved_centered_all_pairs_and_all_source_jets(x):
    with mp.workdps(115):
        x = mp.mpf("0.00001") * (1 + mp.mpc(0, "0.4")) if x == "complex" else mp.mpf(x)
        projected = centered_endpoint_rows(x)
        physical = centered_endpoint_rows(x, method="full")
        exchanged = centered_endpoint_rows(x, reverse=True)
        assert mp.norm(projected - physical) < mp.mpf("1e-70")
        assert mp.norm(projected - exchanged) < mp.mpf("1e-70")


@pytest.mark.parametrize("x", ("0.0001", "0.00001", "0.00000001"))
def test_actual_curved_Schwarz_parity_keeps_odd_endpoint_labels(x):
    with mp.workdps(115):
        x = mp.mpf(x)
        plus = centered_endpoint_rows(x)
        minus = centered_endpoint_rows(-x)
        for j in range(5):
            for r in range(j + 1):
                assert abs(mp.im(minus[j, r]) - (-1) ** j * mp.im(plus[j, r])) < mp.mpf(
                    "1e-65"
                )
        assert abs(mp.im(plus[1, 0])) > mp.mpf("1e-30")
        assert abs(mp.im(plus[3, 0])) > mp.mpf("1e-30")


@cache
def actual_centered_coefficients(points=24):
    with mp.workdps(125):
        radius = mp.mpf("1e-6")
        phases = [mp.exp(2j * mp.pi * j / points) for j in range(points)]
        values = [centered_endpoint_rows(radius * phase) for phase in phases]
        return {
            (j, d): mp.matrix(
                [
                    -mp.im(
                        mp.fsum(
                            values[i][j, r] * phases[i] ** (-d) for i in range(points)
                        )
                        / points
                        / radius**d
                    )
                    for r in range(j + 1)
                ]
            )
            for j, d in extraction.slots()
        }


@pytest.mark.parametrize("slot", extraction.slots())
def test_actual_curved_centered_UV_grades_two_resolutions(slot):
    with mp.workdps(115):
        a = actual_centered_coefficients(24)[slot]
        b = actual_centered_coefficients(32)[slot]
        assert mp.norm(a - b) < mp.mpf("1e-43")
        j, d = slot
        if (j + d) % 2:
            assert mp.norm(a) < mp.mpf("1e-43")
        if slot == (1, 1):
            assert mp.norm(a) < mp.mpf("1e-43")
        if slot == (1, 3):
            assert mp.norm(a) > mp.mpf("1e-10")


@pytest.mark.parametrize("degree", range(1, 9))
@pytest.mark.parametrize("order", (1, 2))
def test_exact_nonnegative_Chebyshev_derivative_expansions(degree, order):
    u = s.symbols("u")
    remainder = s.Poly(s.diff(s.chebyshevt(degree, u), u, order), u)
    coefficients = []
    while not remainder.is_zero:
        n = remainder.degree()
        T = s.Poly(s.chebyshevt(n, u), u)
        coefficient = remainder.LC() / T.LC()
        assert coefficient >= 0
        coefficients.append(coefficient)
        remainder -= coefficient * T
    expected = degree**2 if order == 1 else s.Rational(degree**2 * (degree**2 - 1), 3)
    assert sum(coefficients) == expected


@pytest.mark.parametrize(
    "value", (True, False, 1.0, s.Float(1), "1", None, s.Rational(1, 2), -1, 5)
)
def test_exact_shape_order_guards(value):
    with pytest.raises((TypeError, ValueError)):
        geometry.require_integer(value, 0, 4)


def test_complete_original_frontier_not_closed_by_conversion():
    assert centered.data()["gates"][
        "original_regulator_and_fixed_covariant_matching_unchanged"
    ]
    assert flat.converted()[2]["tensor", 0][0] == 0
    u, a = s.symbols("u a", positive=True)
    leading = angular.leading_average(u, a, 1, 0, 0)
    assert angular.shape_coefficient(leading, u, 0, 4) != 0


@cache
def cap_polynomial(degree):
    u = s.symbols("u")
    primitive = s.integrate(s.chebyshevt(degree, u), u)
    return s.lambdify(u, primitive - primitive.subs(u, -1), "mpmath", cse=True)


def complete_original_UV_shell(j, d, K, P):
    if P == 0:
        return mp.mpf(0)
    q = j + d
    primitive = cap_polynomial(d + 4)
    low = mp.mpf(1000) / K
    cuts = sorted(
        set(
            [low, mp.mpf(1)] + [v for v in [abs(K - P) / K, (K + P) / K] if low < v < 1]
        )
    )

    def integrand(v):
        r = K * v
        cap = (r * r + P * P - K * K) / (2 * r * P)
        if cap <= -1:
            return mp.mpf(0)
        return v ** (3 - q) * primitive(min(mp.mpf(1), cap))

    return -mp.quad(integrand, cuts) * K ** (4 - q) / (4 * mp.pi**2)


@pytest.mark.parametrize("slot", extraction.slots())
@pytest.mark.parametrize("cutoff", ("2000", "1000000"))
@pytest.mark.parametrize("fraction", ("0", "0.01", "0.25", "0.5", "1", "3"))
def test_all_original_two_leg_UV_slots_small_and_large_transfer_tail(
    slot, cutoff, fraction
):
    from p8_vacuum_affine_uniform_uv_remainder import domain, estimates

    with mp.workdps(90):
        j, d = slot
        q = j + d
        K = mp.mpf(cutoff)
        P = K * mp.mpf(fraction)
        U = 1000 + P
        actual = complete_original_UV_shell(j, d, K, P)
        u = s.symbols("u")
        polynomial = s.chebyshevt(d + 4, u)
        approximation = mp.mpf(0)
        if q < 4:
            approximation = -mp.fsum(
                K ** (4 - q - h)
                * P**h
                * mp_number(angular.shape_coefficient(polynomial, u, q, h) / s.pi)
                / (8 * mp.pi**2)
                for h in range(1, 5 - q)
            )
        coefficient = mp_number(estimates.row_bound(j) * domain.EPS ** (-d)) * U**d
        error = coefficient * abs(actual - approximation)
        if P == 0:
            assert error == 0
            return
        if P <= K / 4:
            local_bound = 8192 * 6121 * coefficient * P ** (5 - q) / K
        else:
            radial = K ** (4 - q) / (4 - q) if q < 4 else K / 1000
            shapes = (
                mp.fsum(K ** (4 - q - h) * P**h for h in range(1, 5 - q))
                if q < 4
                else 0
            )
            local_bound = coefficient * (radial + 64 * 6121 * shapes)
        assert error < local_bound
        assert error < mp.mpf("1e40") * (1 + P * P) ** 3 / K


def test_full_flat_grazing_changes_the_linear_coefficient():
    u, _m, p, rows = flat.coefficients()
    for channel in ("tensor", "vector", "scalar"):
        f0 = rows[channel, 0][0]
        correction = -(p**3) * f0.subs(u, 0) / (32 * s.pi**2)
        assert correction != 0
        expected = flat.converted()[2][channel, 0][1]
        assert s.factor(expected - correction) != expected


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_all_named_scoped_original_conversion_residuals(name):
    value = audit.residuals()[name]
    if isinstance(value, s.MatrixBase):
        assert all(v == 0 for v in value)
    else:
        assert value == 0


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_all_unsupported_original_conversion_scope_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_exact_complete_conversion_audit_counts():
    assert len(audit.residuals()) == 76
    assert audit.scalar_entry_count() == 181
    assert len(audit.gates()) == 44 and all(audit.gates().values())
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 131


def test_original_scope_matching_not_closed_by_centered_cancellation():
    assert len(audit.frontier()) == 9
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert (
        "NOT_CURVED_A1_EVALUATION_FULL_CONTACT_MATCHING_OR_V_G_B"
        in audit.ITEM["status"]
    )
