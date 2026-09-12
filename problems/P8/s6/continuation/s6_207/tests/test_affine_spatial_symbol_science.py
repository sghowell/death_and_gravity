"""Independent physical pair, source-jet and full ultraviolet extraction tests."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_cd_metric_noise import stress
from p8_vacuum_affine_spatial_symbol import (
    audit,
    benchmark,
    extraction,
    sectors,
    symbol,
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


def conjugate(v):
    return mp.matrix([mp.conj(x) for x in v])


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


@pytest.mark.parametrize("time", ("-0.5", "0", "0.41"))
@pytest.mark.parametrize("direction", ((3, 0, 4), (1, 2, 2), (0, 0, 1)))
@pytest.mark.parametrize("angle", ("0", "0.7", "2"))
def test_actual_full_two_time_W8_projector_factorization(time, direction, angle):
    with mp.workdps(100):
        center = mp.mpf(time)
        u = center + mp.mpf("1e-5") * mp.exp(mp.mpc(0, "0.4"))
        x = mp.mpf("1e-5") * mp.exp(1j * mp.mpf(angle))
        k = mp.matrix(direction)
        k /= mp.norm(k)
        P = mp.matrix([11, -7, 5])
        l, _, _, _, f1, f2 = joint_frame(k, x * P)
        e1, e2 = real_frame(k)
        D, G = tensors()
        source, fs = scaled_point(u, x, k, l)
        detector, _ = scaled_point(center, x, k, l, 1)
        exact = grouping(full_products(k, l, e1, e2, f1, f2, D, G, source, detector, x))
        projected = projected_products(k, l, D, G, source, detector, x)
        for a, b in zip(exact, projected):
            assert abs(a - b) < mp.mpf("1e-85") * max(mp.mpf(1), abs(a))
        phases = four_inverse_phases(fs)
        assert abs(
            sum(a * g for a, g in zip(exact, phases))
            - sum(b * g for b, g in zip(projected, phases))
        ) < mp.mpf("1e-85")


@pytest.mark.parametrize("K", ("1e3", "1e6", "1e20"))
def test_scaled_point_matches_original_full_physical_W8_modes(K):
    with mp.workdps(100):
        K = mp.mpf(K)
        x = 1 / K
        u = mp.mpf("0.2")
        n = mp.matrix([mp.mpf(3) / 5, 0, mp.mpf(4) / 5])
        P = mp.matrix([11, -7, 5])
        ell, _rho, nl, _, f1, f2 = joint_frame(n, x * P)
        e1, e2 = real_frame(n)
        point, freq = scaled_point(u, x, n, ell)
        for pref, v, pol, nvec in (("k", n, (e1, e2), n), ("l", ell, (f1, f2), nl)):
            scaled = scaled_vectors(v, *pol, point, pref, x)
            actual = physical_modes(u, v / x, nvec, *pol)
            for i in range(3):
                assert mp.norm(scaled[i] - mp.sqrt(x) * actual[i][0]) < mp.mpf(
                    "1e-80"
                ) * max(1, mp.norm(scaled[i]))
                tag = pref + ("t" if i < 2 else "l")
                assert abs(freq[tag] - x * actual[i][1]) < mp.mpf("1e-80")


@pytest.mark.parametrize("module", (sectors, symbol, extraction, benchmark))
def test_complete_exact_symbol_core(module):
    packet = module.data()
    for value in packet["checks"].values():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.simplify(x) == 0 for x in entries)
    assert all(v is True for v in packet["gates"].values())


def complete_scaled_endpoint_rows(x, samples=16, orders=4, method="projected"):
    k = mp.matrix([mp.mpf(3) / 5, 0, mp.mpf(4) / 5])
    P = mp.matrix([11, -7, 5])
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
def actual_UV_coefficients(orders=4, points=24):
    with mp.workdps(125):
        radius = mp.mpf("1e-6")
        phases = [mp.exp(2j * mp.pi * i / points) for i in range(points)]
        values = [
            complete_scaled_endpoint_rows(radius * phase, orders=orders)
            for phase in phases
        ]
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


@pytest.mark.parametrize("x", ("0.00001", "0.00000001", "0.00001+0.000004j"))
def test_actual_curved_all_source_time_jets_two_independent_routes(x):
    with mp.workdps(110):
        x = (
            mp.mpc(x.replace("j", ""))
            if "+" not in x
            else mp.mpc("0.00001", "0.000004")
        )
        projected = complete_scaled_endpoint_rows(x, method="projected")
        full = complete_scaled_endpoint_rows(x, method="full")
        assert mp.norm(projected - full) < mp.mpf("1e-65") * max(
            mp.mpf(1), mp.norm(full)
        )
        for j in range(5):
            for r in range(j + 1, 5):
                assert abs(full[j, r]) < mp.mpf("1e-90")


@pytest.mark.parametrize("j,d", extraction.slots())
def test_actual_all_retained_curved_coefficients_two_resolutions(j, d):
    with mp.workdps(110):
        a = actual_UV_coefficients(4, 24)[j, d]
        b = actual_UV_coefficients(4, 32)[j, d]
        assert mp.norm(a - b) < mp.mpf("1e-35") * max(mp.mpf(1), mp.norm(b))


@pytest.mark.parametrize("j,d", extraction.slots())
def test_actual_all_retained_coefficients_independent_of_sixth_eighth_WKB_terms(j, d):
    with mp.workdps(110):
        full = actual_UV_coefficients(4, 24)[j, d]
        low = actual_UV_coefficients(2, 24)[j, d]
        assert mp.norm(full - low) < mp.mpf("1e-35") * max(mp.mpf(1), mp.norm(full))


def test_actual_curved_odd_endpoints_not_deleted_by_static_fixture():
    with mp.workdps(110):
        coefficients = actual_UV_coefficients(4, 24)
        first = [
            -mp.im(value)
            for d in extraction.retained_degrees(1)
            for value in coefficients[1, d]
        ]
        assert max(abs(v) for v in first) > 1
        # This fixture's third UV row cancels, but its finite endpoint does not.
        third = [
            -mp.im(value)
            for d in extraction.retained_degrees(3)
            for value in coefficients[3, d]
        ]
        assert max(abs(v) for v in third) < mp.mpf("1e-50")
        finite = complete_scaled_endpoint_rows(mp.mpf("0.0001"))
        for j in (1, 3):
            assert max(abs(mp.im(finite[j, r])) for r in range(j + 1)) > mp.mpf("1e-8")


@cache
def actual_remainder_rows(x):
    with mp.workdps(125):
        x = mp.mpf(x)
        actual = complete_scaled_endpoint_rows(x)
        coeff = actual_UV_coefficients(4, 24)
        result = []
        for j in range(5):
            row = mp.matrix(
                [
                    actual[j, r]
                    - sum(coeff[j, d][r] * x**d for d in extraction.retained_degrees(j))
                    for r in range(j + 1)
                ]
            )
            result.append(row * x ** (j - 1) / x**4)
        return tuple(result)


@pytest.mark.parametrize("j", range(5))
def test_actual_fixed_P_subtracted_remainder_has_integrable_fourth_power(j):
    with mp.workdps(110):
        first = actual_remainder_rows("1e-8")[j]
        second = actual_remainder_rows("5e-9")[j]
        assert mp.norm(second) > mp.mpf("1e-20")
        assert mp.norm(first - second) < mp.mpf("0.01") * mp.norm(second)


def test_source_detector_coincidence_order_control():
    with mp.workdps(80):
        t = mp.mpf("0.2")
        alpha = mp.mpf("1.5")
        beta = mp.mpf("0.7")
        correct = mp.diff(lambda u: mp.exp(alpha * u + beta * t), t)
        wrong = mp.diff(lambda u: mp.exp((alpha + beta) * u), t)
        assert abs(wrong - correct - beta * mp.exp((alpha + beta) * t)) < mp.mpf(
            "1e-70"
        )
        assert abs(wrong - correct) > mp.mpf("0.1")


def test_distinct_T_L_inverse_phases_cannot_be_replaced_by_one():
    with mp.workdps(90):
        x = mp.mpf("0.0001")
        t = mp.mpf("0.4")
        k = mp.matrix([mp.mpf(3) / 5, 0, mp.mpf(4) / 5])
        P = mp.matrix([11, -7, 5])
        l, _, _, _, _, _ = joint_frame(k, x * P)
        D, G = tensors()
        source, freq = scaled_point(t, x, k, l)
        detector, _ = scaled_point(t, x, k, l, 1)
        values = projected_products(k, l, D, G, source, detector, x)
        phases = four_inverse_phases(freq)
        actual = sum(a * g for a, g in zip(values, phases))
        wrong = sum(values) * phases[0]
        assert max(abs(g - phases[0]) for g in phases) > mp.mpf("1e-12")
        assert abs(actual - wrong) > mp.mpf("1e-16")


def flat_tensor(channel):
    if channel == "tensor":
        return mp.diag([1, -1, 0]) / mp.sqrt(2)
    if channel == "vector":
        return mp.matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / mp.sqrt(2)
    if channel == "scalar":
        return mp.diag([-1, -1, 2]) / mp.sqrt(6)
    raise ValueError("Unknown physical channel")


@cache
def direct_flat_coefficients(channel, transfer):
    with mp.workdps(110):
        transfer = mp.mpf(transfer)
        m = mp.mpf(1000)
        radius = mp.mpf("1e-6")
        points = 32
        phases = [mp.exp(2j * mp.pi * i / points) for i in range(points)]
        values = []
        for phase in phases:
            x = radius * phase
            k = mp.matrix([0, 0, 1])
            P = mp.matrix([0, 0, transfer])
            l, _, _, _, f1, f2 = joint_frame(k, x * P)
            e1, e2 = real_frame(k)
            source = {"a": mp.mpf(1)}
            detector = {"a": mp.mpf(1)}
            frequencies = {}
            for prefix, v in (("k", k), ("l", l)):
                omega = mp.sqrt((v.T * v)[0] + m * m * x * x)
                f = 1 / mp.sqrt(2 * omega)
                for suffix in ("t", "l"):
                    source[prefix + suffix] = (f, -1j * omega * f, omega)
                    detector[prefix + suffix] = (f, 1j * omega * f, omega)
                    frequencies[prefix + suffix] = omega
            D = flat_tensor(channel)
            F = sum(full_products(k, l, e1, e2, f1, f2, D, D, source, detector, x))
            inv = 1 / (frequencies["kt"] + frequencies["lt"])
            values.append([(1, 0, -1, 0, 1)[j] * F * inv ** (j + 1) for j in range(5)])
        return {
            (j, d): mp.fsum(values[i][j] * phases[i] ** (-d) for i in range(points))
            / points
            / radius**d
            for j, d in extraction.slots()
        }


@pytest.mark.parametrize("channel", ("tensor", "vector", "scalar"))
@pytest.mark.parametrize("transfer", ("13", "-7"))
@pytest.mark.parametrize("j,d", extraction.slots())
def test_all_exact_massive_flat_coefficients_against_full_field_pair_Cauchy(
    channel, transfer, j, d
):
    with mp.workdps(100):
        m, p, _, _, rows = benchmark.coefficients()
        expression = rows[channel, j][d].subs({m: 1000, p: int(transfer)})
        expected = mp.mpf(str(s.N(expression, 105)))
        actual = direct_flat_coefficients(channel, transfer)[j, d]
        assert abs(actual - expected) < mp.mpf("1e-45") * max(mp.mpf(1), abs(expected))


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
    5,
)


@pytest.mark.parametrize("value", BAD_ORDER)
def test_exact_endpoint_order_guards(value):
    extraction.retained_degrees(0)
    extraction.retained_degrees(1)
    with pytest.raises((TypeError, ValueError)):
        extraction.retained_degrees(value)
    with pytest.raises((TypeError, ValueError)):
        extraction.source_iterates(s.Symbol("F"), s.Symbol("g"), s.Symbol("t"), value)


@pytest.mark.parametrize(
    "v", (s.zeros(3, 1), s.Matrix([1, s.I, 0]), s.Matrix([1, 2]), s.Matrix([[1, 2, 3]]))
)
def test_bilinear_projector_domain_and_null_control(v):
    with pytest.raises(ValueError):
        sectors.projector(v)


def test_nonzero_constant_and_positive_root_branch_guards():
    with pytest.raises(ValueError):
        benchmark.reciprocal((0, 1, 0, 0, 0))
    with pytest.raises(ValueError):
        benchmark.square_root_one((-1, 0, 0, 0, 0))


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_all_named_scoped_symbol_residuals(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(x == 0 for x in entries)


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_all_unsupported_symbol_scope_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_exact_symbol_audit_counts():
    assert len(audit.residuals()) == 37 and audit.scalar_entry_count() == 119
    assert len(audit.gates()) == 39 and len(audit.controls()) == 9
    assert audit.rejected_inputs() == 129


def test_original_frontier_not_closed_by_finite_extraction():
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert all(v is True for v in audit.gates().values())
