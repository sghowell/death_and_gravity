"""Independent actual inverse-radius domains, UV remainder and original-mask tails."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_cd_metric_noise import stress
from p8_vacuum_affine_spatial_symbol import extraction, sectors
from p8_vacuum_affine_uniform_uv_remainder import audit, domain, estimates, limit, tail
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


def uniform_configuration(transfer, direction=(3, 0, 4)):
    n = mp.matrix(direction)
    n /= mp.norm(n)
    P = mp.mpf(transfer) * mp.matrix([mp.mpf(2) / 3, mp.mpf(1) / 3, mp.mpf(2) / 3])
    rho = mp.mpf("0.01") / (1000 + mp.mpf(transfer))
    return n, P, rho


@pytest.mark.parametrize("transfer", ("0", "1", "1000", "1e8", "1e20"))
@pytest.mark.parametrize("time", ("-0.5", "0", "0.5"))
@pytest.mark.parametrize("angle", ("0", "0.7", "2.4"))
def test_actual_full_normalized_joint_frequency_domain(transfer, time, angle):
    with mp.workdps(100):
        n, P, rho = uniform_configuration(transfer)
        x = mp.mpf("0.99") * rho * mp.exp(1j * mp.mpf(angle))
        l, _, _, _, _, _ = joint_frame(n, x * P)
        u = mp.mpf(time) + mp.mpf("0.00009") * mp.exp(mp.mpc(0, "0.8"))
        point, freq = scaled_point(u, x, n, l)
        for key, W in freq.items():
            omega = point[key][2]
            assert mp.re(omega) > mp.mpf("0.98") / (mp.mpf(25) / 16)
            assert abs(omega) < 2
            assert abs(1 - 1000**2 * x * x / omega**2) < mp.mpf("1.01")
            assert abs(W / omega - 1) < mp.mpf("0.001")
            assert mp.re(W) > mp.mpf("0.5") and abs(W) < 3
        assert all(abs(g) < 1 for g in four_inverse_phases(freq))


@pytest.mark.parametrize("transfer", ("0", "1000", "1e8", "1e20"))
@pytest.mark.parametrize("direction", ((3, 0, 4), (1, 2, 2), (0, 0, 1)))
@pytest.mark.parametrize("phase_sign", (-1, 1))
def test_actual_all_ten_normalized_readouts_and_complex_frames(
    transfer, direction, phase_sign
):
    with mp.workdps(100):
        n, P, rho = uniform_configuration(transfer, direction)
        x = mp.mpf("0.99") * rho * mp.exp(mp.mpc(0, "1.1"))
        l, _, nl, R, f1, f2 = joint_frame(n, x * P)
        e1, e2 = real_frame(n)
        assert mp.norm(R.T * R - mp.eye(3)) < mp.mpf("1e-85")
        assert mp.norm(nl + n) < mp.mpf("0.025")
        assert all(mp.norm(v) < 2 for v in (nl, f1, f2))
        u = mp.mpf("0.25") + mp.mpf("0.00004") * mp.exp(mp.mpc(0, "0.8"))
        point, _ = scaled_point(u, x, n, l, phase_sign)
        for prefix, v, p1, p2 in (("k", n, e1, e2), ("l", l, f1, f2)):
            fields = scaled_vectors(v, p1, p2, point, prefix, x, phase_sign)
            assert all(
                abs(point[prefix + kind][0]) < 1 and abs(point[prefix + kind][1]) < 4
                for kind in ("t", "l")
            )
            assert all(mp.norm(vector) < 20 for vector in fields)


@pytest.mark.parametrize("transfer", ("0", "1000", "1e20"))
def test_complete_sixteen_component_stress_pairs_below_new_bound(transfer):
    with mp.workdps(95):
        n, P, rho = uniform_configuration(transfer)
        x = mp.mpf("0.98") * rho * mp.exp(mp.mpc(0, "0.9"))
        l, _, _, _, f1, f2 = joint_frame(n, x * P)
        e1, e2 = real_frame(n)
        point, _ = scaled_point(mp.mpf("0.3") + mp.mpc(0, "0.00004"), x, n, l)
        first = scaled_vectors(n, e1, e2, point, "k", x)
        second = scaled_vectors(l, f1, f2, point, "l", x)
        matrices = complete_stress_matrices()
        for a in first:
            for b in second:
                tensor = mp.matrix([(a.T * M * b)[0] for M in matrices])
                assert mp.norm(tensor) < 1000


@pytest.mark.parametrize("module", (domain, estimates, tail, limit))
def test_complete_exact_uniform_bound_core(module):
    d = module.data()
    for v in d["checks"].values():
        values = list(v) if isinstance(v, s.MatrixBase) else [v]
        assert all(s.simplify(x) == 0 for x in values)
    assert all(v is True for v in d["gates"].values())


def uniform_endpoint_rows(x, transfer, samples=16, orders=4, method="projected"):
    k, P, _ = uniform_configuration(transfer)
    l, _, _, _, f1, f2 = joint_frame(k, x * P)
    e1, e2 = real_frame(k)
    center = mp.mpf("0.25")
    radius = mp.mpf("1e-5")
    D, G = tensors()
    detector, df = scaled_point(center, x, k, l, 1, orders)
    gdet = four_inverse_phases(df)
    phases = [mp.exp(2j * mp.pi * i / samples) for i in range(samples)]
    amplitudes = []
    inverse_phases = []
    for phase in phases:
        source, freq = scaled_point(center + radius * phase, x, k, l, -1, orders)
        if method == "projected":
            values = projected_products(k, l, D, G, source, detector, x)
        elif method == "full":
            values = grouping(
                full_products(k, l, e1, e2, f1, f2, D, G, source, detector, x)
            )
        else:
            raise ValueError("Unknown independent comparison route")
        amplitudes.append(values)
        inverse_phases.append(four_inverse_phases(freq))
    result = mp.zeros(5, 5)
    for sector in range(4):
        g = [
            mp.fsum(
                inverse_phases[i][sector] * phases[i] ** (-j) for i in range(samples)
            )
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
        for source_jet in range(5):
            current = [
                amp[n - source_jet] / mp.factorial(source_jet)
                if n >= source_jet
                else mp.mpc(0)
                for n in range(5)
            ]
            for endpoint in range(5):
                result[endpoint, source_jet] += (
                    (-1j) * 1j**endpoint * gdet[sector] * current[0]
                )
                current = [
                    (n + 1) * mp.fsum(g[j] * current[n + 1 - j] for j in range(n + 2))
                    for n in range(len(current) - 1)
                ]
    return result


@cache
def uniform_coefficients(transfer, points=24):
    with mp.workdps(125):
        _, _, rho = uniform_configuration(transfer)
        radius = rho / 10
        phases = [mp.exp(2j * mp.pi * i / points) for i in range(points)]
        values = [uniform_endpoint_rows(radius * phase, transfer) for phase in phases]
        return {
            (j, d): mp.matrix(
                [
                    mp.fsum(values[i][j, r] * phases[i] ** (-d) for i in range(points))
                    / points
                    / radius**d
                    for r in range(j + 1)
                ]
            )
            for j, d in extraction.slots()
        }


@pytest.mark.parametrize("transfer", ("0", "1000", "1e20"))
@pytest.mark.parametrize("j,d", extraction.slots())
def test_actual_uniform_UV_coefficients_two_resolutions_and_all_jet_bounds(
    transfer, j, d
):
    with mp.workdps(110):
        _, _, rho = uniform_configuration(transfer)
        a = uniform_coefficients(transfer, 24)[j, d]
        b = uniform_coefficients(transfer, 32)[j, d]
        assert mp.norm(a - b) < mp.mpf("1e-32") * max(mp.mpf(1), mp.norm(b))
        D, G = tensors()
        norm = mp.norm(D) * mp.norm(G)
        for r in range(j + 1):
            bound = mp.mpf(str(estimates.jet_bound(j, r))) * rho ** (-d) * norm
            assert abs(a[r]) <= bound
        row = mp.mpf(str(estimates.row_bound(j))) * rho ** (-d) * norm
        assert mp.norm(a) <= row


@cache
def actual_uniform_remainders(transfer, where):
    with mp.workdps(125):
        _, _, rho = uniform_configuration(transfer)
        if where == "far4":
            x = rho / 4
        elif where == "far8":
            x = rho / 8
        elif where == "near":
            x = rho
        elif where == "mass":
            x = mp.mpf("0.001")
        else:
            raise ValueError("Unknown fixed domain point")
        actual = uniform_endpoint_rows(x, transfer)
        coeff = uniform_coefficients(transfer, 24)
        rows = []
        for j in range(5):
            row = mp.matrix(
                [
                    actual[j, r]
                    - sum(coeff[j, d][r] * x**d for d in extraction.retained_degrees(j))
                    for r in range(j + 1)
                ]
            )
            rows.append(x ** (j - 1) * row)
        return x, tuple(rows)


@pytest.mark.parametrize("transfer", ("0", "1000", "1e20"))
@pytest.mark.parametrize("where", ("far4", "far8"))
@pytest.mark.parametrize("j", range(5))
def test_actual_all_momentum_far_remainder_with_explicit_uniform_constant(
    transfer, where, j
):
    with mp.workdps(110):
        _, _, rho = uniform_configuration(transfer)
        x, rows = actual_uniform_remainders(transfer, where)
        D, G = tensors()
        norm = mp.norm(D) * mp.norm(G)
        bound = 2 * mp.mpf(str(estimates.row_bound(j))) * rho ** (j - 5) * x**4 * norm
        assert mp.norm(rows[j]) > x ** (j - 1) * mp.mpf("1e-60")
        assert mp.norm(rows[j]) < bound


@pytest.mark.parametrize("transfer", ("0", "1000", "1e20"))
@pytest.mark.parametrize("where", ("near", "mass"))
def test_actual_near_remainder_retains_full_row_and_every_subtracted_power(
    transfer, where
):
    with mp.workdps(110):
        _, _, rho = uniform_configuration(transfer)
        x, rows = actual_uniform_remainders(transfer, where)
        r = 1 / x
        D, G = tensors()
        norm = mp.norm(D) * mp.norm(G)
        nu = mp.sqrt(1000**2 + r * r / (mp.mpf(25) / 16) ** 2)
        bound = mp.mpf("4e27") * nu + sum(
            mp.mpf(str(estimates.row_bound(j))) * rho ** (-d) * x ** (j - 1 + d)
            for j, d in extraction.slots()
        )
        combined = mp.matrix([sum(rows[j][i] for j in range(i, 5)) for i in range(5)])
        assert mp.norm(combined) < bound * norm


def removed_fraction(r, P, K):
    if r > K:
        return mp.mpf(1)
    if r == 0:
        return mp.mpf(int(P > K))
    if P == 0:
        return mp.mpf(0)
    cosine = (r * r + P * P - K * K) / (2 * r * P)
    return max(mp.mpf(0), min(mp.mpf(1), (1 + cosine) / 2))


def regulator_configuration(transfer, fraction):
    P = mp.mpf(transfer)
    U = 1000 + P
    L = 200 * U
    K = max(mp.mpf(1000), mp.mpf(fraction) * L)
    return P, U, L, K


def finite_cuts(low, high, candidates):
    return sorted(set([low, high] + [x for x in candidates if low < x < high]))


@cache
def near_removed_moment(power, transfer, fraction):
    with mp.workdps(95):
        P, _U, L, K = regulator_configuration(transfer, fraction)
        lower = mp.mpf(1000) / L
        cuts = finite_cuts(lower, mp.mpf(1), [abs(K - P) / L, K / L, (K + P) / L])
        value = mp.quad(
            lambda v: v ** (3 - power) * removed_fraction(L * v, P, K), cuts
        )
        return L ** (4 - power) * value / (2 * mp.pi**2)


@pytest.mark.parametrize("power", range(5))
@pytest.mark.parametrize("transfer", ("0", "1000", "1e6"))
@pytest.mark.parametrize("fraction", ("0.001", "0.5", "1", "2"))
def test_complete_near_original_removed_union_power_and_log_bounds(
    power, transfer, fraction
):
    with mp.workdps(85):
        P, U, L, K = regulator_configuration(transfer, fraction)
        actual = near_removed_moment(power, transfer, fraction)
        full = (
            (L ** (4 - power) - 1000 ** (4 - power)) / (4 - power)
            if power < 4
            else mp.log(L / 1000)
        )
        full /= 2 * mp.pi**2
        assert 0 <= actual <= full * (1 + mp.mpf("1e-70"))
        assert actual <= 201 * U / K * full * (1 + mp.mpf("1e-70"))
        if P == 0 and K < L:
            exact = (
                (L ** (4 - power) - K ** (4 - power)) / (4 - power)
                if power < 4
                else mp.log(L / K)
            )
            exact /= 2 * mp.pi**2
            assert abs(actual - exact) < mp.mpf("1e-65") * max(mp.mpf(1), abs(exact))


@pytest.mark.parametrize("transfer", ("0", "1000", "1e6"))
@pytest.mark.parametrize("fraction", ("0.001", "0.5", "1", "2"))
def test_complete_far_original_removed_union_inverse_power_tail(transfer, fraction):
    with mp.workdps(90):
        P, _U, L, K = regulator_configuration(transfer, fraction)
        cuts = finite_cuts(
            mp.mpf(0), mp.mpf(1), [L / r for r in (abs(K - P), K, K + P) if r > L]
        )
        actual = mp.quad(
            lambda v: mp.mpf(1) if v == 0 else removed_fraction(L / v, P, K), cuts
        ) / (2 * mp.pi**2 * L)
        assert 0 < actual <= 1 / (2 * mp.pi**2 * L) * (1 + mp.mpf("1e-70"))
        assert actual <= 1 / (mp.pi**2 * K) * (1 + mp.mpf("1e-70"))


@pytest.mark.parametrize("transfer", ("0", "1000", "1e6"))
@pytest.mark.parametrize("fraction", ("0.001", "0.5", "1", "2"))
def test_complete_unexpanded_low_original_removed_union_tail(transfer, fraction):
    with mp.workdps(90):
        P, U, _L, K = regulator_configuration(transfer, fraction)
        A = mp.mpf(25) / 16
        m = mp.mpf(1000)
        cuts = finite_cuts(
            mp.mpf(0), mp.mpf(1), [r / m for r in (abs(K - P), K, K + P)]
        )
        weight = lambda v: v * v * mp.sqrt(1 + v * v / A**2)
        actual = mp.quad(lambda v: weight(v) * removed_fraction(m * v, P, K), cuts)
        full = mp.quad(weight, [0, 1])
        assert 0 <= actual <= full * (1 + mp.mpf("1e-70"))
        assert actual <= U / K * full * (1 + mp.mpf("1e-70"))


@pytest.mark.parametrize("power", range(7))
@pytest.mark.parametrize("P", (0, 1, 1000, 10**8, 10**20))
def test_same_six_spatial_derivatives_control_every_finite_and_tail_power(power, P):
    P = s.Integer(P)
    assert (domain.MASS + P) ** power <= estimates.spatial_weight(power) * (
        1 + P * P
    ) ** 3


BAD_ORDER = (
    True,
    False,
    1.0,
    s.Float(1),
    "1",
    None,
    s.oo,
    -s.oo,
    s.I,
    s.Symbol("j"),
    -1,
)


@pytest.mark.parametrize("value", BAD_ORDER + (7,))
def test_exact_spatial_weight_guard(value):
    estimates.spatial_weight(0)
    estimates.spatial_weight(1)
    with pytest.raises((TypeError, ValueError)):
        estimates.spatial_weight(value)


@pytest.mark.parametrize("value", BAD_ORDER + (5,))
def test_exact_complete_endpoint_row_guard(value):
    estimates.row_bound(0)
    estimates.row_bound(1)
    with pytest.raises((TypeError, ValueError)):
        estimates.row_bound(value)


@pytest.mark.parametrize("value", BAD_ORDER + (5,))
def test_exact_source_time_jet_guard(value):
    estimates.jet_bound(4, 0)
    estimates.jet_bound(4, 1)
    with pytest.raises((TypeError, ValueError)):
        estimates.jet_bound(4, value)


def test_actual_masks_not_replaced_and_repartition_anchor_retained():
    with mp.workdps(70):
        assert removed_fraction(mp.mpf(500), mp.mpf(2000), mp.mpf(1000)) == 1
        assert mp.mpf(500) < mp.mpf(1000)
    from p8_vacuum_affine_reference_state_prefactor import limit as original

    assert (
        original.constants()["complete_known_finite_curved_piece"]
        == 2 * s.Integer(10) ** 48
    )
    assert limit.KNOWN != s.Integer(5) * 10**54 + estimates.FINITE


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_all_named_uniform_remainder_residuals(name):
    v = audit.residuals()[name]
    values = list(v) if isinstance(v, s.MatrixBase) else [v]
    assert all(x == 0 for x in values)


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_all_unsupported_uniform_remainder_scopes(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_exact_uniform_remainder_audit_counts():
    assert len(audit.residuals()) == 34 and audit.scalar_entry_count() == 34
    assert len(audit.gates()) == 52 and len(audit.controls()) == 9
    assert audit.rejected_inputs() == 130


def test_original_frontier_not_closed_by_uniform_remainder():
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert all(v is True for v in audit.gates().values())
