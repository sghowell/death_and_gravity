"""Independent complete-W8 coefficients, source-time jets, contacts and pole tests."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_cd_metric_noise import stress
from p8_vacuum_affine_spatial_matching_difference import (
    audit,
    contact,
    density,
    jets,
    matching,
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


GEOS = ("00", "01", "10", "11", "TL", "LT", "LL")
PAIRING = ((0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0), (2, 2, 1), (3, 3, 2), (4, 4, 3))


@cache
def literal_expected_function():
    rows = jets.endpoint_products()
    keys = [(key, d) for key, row in rows.items() for d in range(len(row))]
    return keys, s.lambdify(
        (jets.t, jets.u, jets.p, jets.m),
        [rows[key][d] for key, d in keys],
        "mpmath",
        cse=True,
    )


def literal_endpoint_rows(x, time, angle, transfer, samples=20):
    n = mp.matrix([mp.sqrt(1 - angle * angle), 0, angle])
    k = n
    ell = -n + mp.matrix([0, 0, x * transfer])
    detector, df = scaled_point(time, x, k, ell, 1, 4)
    detector_amps = scalar_values(detector, 1000 * x)
    gd = four_inverse_phases(df)
    clock_radius = mp.mpf("1e-5")
    phases = [mp.exp(2j * mp.pi * i / samples) for i in range(samples)]
    source_amps = []
    source_g = []
    for phase in phases:
        source, freq = scaled_point(time + clock_radius * phase, x, k, ell, -1, 4)
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
def literal_coefficients(time, angle, transfer, points=24):
    with mp.workdps(135):
        time, angle, transfer = map(mp.mpf, (time, angle, transfer))
        radius = mp.mpf("1e-7") / (1 + transfer / 1000)
        phases = [mp.exp(2j * mp.pi * i / points) for i in range(points)]
        values = [
            literal_endpoint_rows(radius * z, time, angle, transfer) for z in phases
        ]
        keys, _fn = literal_expected_function()
        return {
            key: -mp.im(
                mp.fsum(
                    values[i][key[0]] * phases[i] ** (-key[1]) for i in range(points)
                )
                / points
                / radius ** key[1]
            )
            for key in keys
        }


@pytest.mark.parametrize("module", (jets, density, contact, matching))
def test_all_exact_core_identities_and_scope_gates(module):
    data = module.data()
    for name, v in data["checks"].items():
        entries = list(v) if isinstance(v, s.MatrixBase) else [v]
        assert all(s.cancel(x) == 0 for x in entries), name
    assert all(bool(v) for v in data["gates"].values())


@pytest.mark.parametrize("time", ("-0.5", "0", "0.25", "0.5"))
@pytest.mark.parametrize("angle", ("-0.6", "0.37"))
@pytest.mark.parametrize("transfer", ("13", "1000000"))
def test_independent_full_W8_all_endpoint_and_source_jets(time, angle, transfer):
    with mp.workdps(125):
        actual = literal_coefficients(time, angle, transfer)
        keys, fn = literal_expected_function()
        expected = fn(mp.mpf(time), mp.mpf(angle), mp.mpf(transfer), 1000)
        for key, value in zip(keys, expected):
            assert abs(actual[key] - value) < mp.mpf("1e-35") * (1 + abs(value)), (
                key,
                actual[key],
                value,
            )


@pytest.mark.parametrize("time", ("-0.5", "0.25", "0.5"))
@pytest.mark.parametrize("transfer", ("13", "1000000"))
def test_full_constrained_ten_field_nine_pair_azimuth_unshifted(time, transfer):
    with mp.workdps(115):
        t, p = mp.mpf(time), mp.mpf(transfer)
        angle = mp.mpf("0.37")
        x = mp.mpf("1e-6") * (1 + mp.mpc(0, "0.4")) / (1 + p / 1000)
        for channel in ("tensor", "vector", "scalar"):
            D = old_flat.channel_tensor(channel)
            D = mp.matrix(
                [[mp.mpf(str(s.N(v, 115))) for v in D.row(i)] for i in range(3)]
            )
            literal = mp.mpc(0)
            for index in range(8):
                phi = 2 * mp.pi * index / 8
                k = mp.matrix(
                    [
                        mp.sqrt(1 - angle * angle) * mp.cos(phi),
                        mp.sqrt(1 - angle * angle) * mp.sin(phi),
                        angle,
                    ]
                )
                ell = -k + mp.matrix([0, 0, x * p])
                e1, e2 = real_frame(k)
                f1, f2 = analytic_frame(ell, -k)
                source, freq = scaled_point(t, x, k, ell, -1, 4)
                detector, _ = scaled_point(t, x, k, ell, 1, 4)
                vals = grouping(
                    full_products(k, ell, e1, e2, f1, f2, D, D, source, detector, x)
                )
                g = four_inverse_phases(freq)
                literal += mp.fsum(v * h for v, h in zip(vals, g)) / 8
            ug, y, geo = old_flat.averaged_geometry(channel)
            c = s.lambdify((ug, y), list(geo.values()), "mpmath", cse=True)(
                angle, p * x
            )
            k = mp.matrix([mp.sqrt(1 - angle * angle), 0, angle])
            ell = -k + mp.matrix([0, 0, x * p])
            source, freq = scaled_point(t, x, k, ell, -1, 4)
            detector, _ = scaled_point(t, x, k, ell, 1, 4)
            A, B = scalar_values(detector, 1000 * x), scalar_values(source, 1000 * x)
            g = four_inverse_phases(freq)
            expected = mp.fsum(
                cg * A[left] * B[right] * g[sector]
                for cg, (left, right, sector) in zip(c, PAIRING)
            )
            assert abs(literal - expected) < mp.mpf("1e-80")


@pytest.mark.parametrize(
    "time", (s.Rational(-1, 2), 0, s.Rational(1, 4), s.Rational(1, 2))
)
@pytest.mark.parametrize("transfer", (0, 13, 1000, 10**9))
@pytest.mark.parametrize("channel", ("tensor", "vector", "scalar"))
def test_spatial_coefficient_bounds_and_exact_origin(time, transfer, channel):
    substitutions = {jets.t: time, jets.p: transfer, jets.m: 1000}
    f2 = density.quadratic(*density.CHANNELS[channel]).subs(substitutions)
    f4 = [
        value.subs(substitutions)
        for value in density.logarithmic(*density.CHANNELS[channel])
    ]
    with mp.workdps(90):
        lhs2 = abs(mp.mpf(str(s.N(f2, 90))))
        lhs4 = sum(abs(mp.mpf(str(s.N(v, 90)))) for v in f4)
        weight = (1 + mp.mpf(transfer) ** 2) ** 2
        if transfer == 0:
            assert lhs2 == lhs4 == 0
        else:
            assert lhs2 < weight / 200
            assert lhs4 < 10000 * weight


@pytest.mark.parametrize("transfer", ((0, 0, 0), (11, -7, 5), (10**6, 2, -4)))
def test_contact_cancels_only_spatial_difference_noncommuting_and_nonzero(transfer):
    D = s.Matrix([[2, 1, -1], [1, -3, 2], [-1, 2, 1]]) / 7
    G = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 17
    assert D * G != G * D
    P = s.Matrix(transfer)
    k = s.Matrix([3, -4, 5])
    full = contact.paired_contact(k, P, D, G, s.Rational(289, 256), 1000)
    zero = contact.paired_contact(k, s.zeros(3, 1), D, G, s.Rational(289, 256), 1000)
    assert full == zero
    assert s.trace(full) != 0


def test_fixed_pole_not_finite_matching_weights():
    assert s.Rational(13, 60) != -s.Rational(1, 30)
    assert s.Rational(1, 36) != -s.Rational(1, 18)
    assert matching.data()["gates"]["evanescent_finite_terms_still_required"]
    # One-ball quadratic coefficients are not proportional to the Einstein Hessian.
    f2 = [
        density.quadratic(*density.CHANNELS[ch])
        for ch in ("tensor", "vector", "scalar")
    ]
    assert f2[1] != 0
    assert matching.hessian_difference("vector")[1]["R_old"] == 0


def test_actual_odd_endpoint_is_nonzero_before_spatial_difference():
    rows = density.angular_coefficients()
    assert rows["tensor", 1, 0, 3] != 0
    assert density.spatial_difference()["tensor", 1, 0, 3] == 0


def test_exact_source_detector_green_form_and_missing_term_control():
    t = jets.t
    D = (1 - 4 * t * t) ** 4
    G = (1 + t * t) * D
    exact = s.expand(
        D * s.diff(jets.a * s.diff(G, t), t) - G * s.diff(jets.a * s.diff(D, t), t)
    )
    wrong = s.expand(jets.a * (D * s.diff(G, t, 2) - G * s.diff(D, t, 2)))
    assert s.integrate(exact, (t, -s.Rational(1, 2), s.Rational(1, 2))) == 0
    assert s.integrate(wrong, (t, -s.Rational(1, 2), s.Rational(1, 2))) != 0


@cache
def full_angular_expected_function():
    rows = density.angular_coefficients()
    keys = list(rows)
    return keys, s.lambdify(
        (jets.t, jets.p, jets.m), list(rows.values()), "mpmath", cse=True
    )


@cache
def exact_geometry_function(channel):
    u, y, geo = old_flat.averaged_geometry(channel)
    return s.lambdify((u, y), list(geo.values()), "mpmath", cse=True)


@cache
def independent_integrated_coefficients(time, transfer, points=24):
    with mp.workdps(135):
        time, p = mp.mpf(time), mp.mpf(transfer)
        radius = mp.mpf("1e-7") / (1 + p / 1000)
        phases = [mp.exp(2j * mp.pi * i / points) for i in range(points)]
        nodes, weights = mp.gauss_quadrature(5, "legendre")
        keys, _fn = full_angular_expected_function()
        integrated = {key: mp.mpc(0) for key in keys}
        for u, weight in zip(nodes, weights):
            values = [literal_endpoint_rows(radius * z, time, u, p) for z in phases]
            geometries = {
                ch: [exact_geometry_function(ch)(u, p * radius * z) for z in phases]
                for ch in ("tensor", "vector", "scalar")
            }
            for key in keys:
                ch, j, r, d = key
                summed = [
                    mp.fsum(
                        geometries[ch][i][q] * values[i][geo, j, r]
                        for q, geo in enumerate(GEOS)
                    )
                    for i in range(points)
                ]
                value = -mp.im(
                    mp.fsum(summed[i] * phases[i] ** (-d) for i in range(points))
                    / points
                    / radius**d
                )
                integrated[key] += weight * value / (4 * mp.pi**2)
        return integrated


@pytest.mark.parametrize("time", ("-0.5", "0.25", "0.5"))
@pytest.mark.parametrize("transfer", ("0", "1000", "1000000"))
def test_independent_all_angular_UV_coefficients_actual_four_W8(time, transfer):
    with mp.workdps(125):
        actual = independent_integrated_coefficients(time, transfer)
        keys, fn = full_angular_expected_function()
        expected = fn(mp.mpf(time), mp.mpf(transfer), 1000)
        for key, value in zip(keys, expected):
            assert abs(actual[key] - value) < mp.mpf("1e-32") * (1 + abs(value)), (
                key,
                actual[key],
                value,
            )


def test_independent_complete_angular_coefficients_two_Cauchy_resolutions():
    with mp.workdps(125):
        a = independent_integrated_coefficients("0.25", "1000", 24)
        b = independent_integrated_coefficients("0.25", "1000", 32)
        assert max(abs(a[k] - b[k]) / (1 + abs(a[k])) for k in a) < mp.mpf("1e-35")


@pytest.mark.parametrize("degree", range(10))
def test_independent_five_node_angular_moments(degree):
    with mp.workdps(90):
        nodes, weights = mp.gauss_quadrature(5, "legendre")
        actual = mp.fsum(w * u**degree for u, w in zip(nodes, weights))
        expected = mp.mpf(1 + (-1) ** degree) / (degree + 1)
        assert abs(actual - expected) < mp.mpf("1e-85")


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_all_exact_scoped_spatial_matching_difference_residuals(name):
    value = audit.residuals()[name]
    values = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(v == 0 for v in values)


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_all_unsupported_spatial_matching_difference_scope_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_exact_spatial_matching_difference_audit_counts():
    assert len(audit.residuals()) == 53
    assert audit.scalar_entry_count() == 238
    assert len(audit.gates()) == 43
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 133


def test_original_frontier_not_closed_by_spatial_pole_identity():
    assert audit.frontier() == audit.previous.frontier()
    assert len(audit.matching()) == len(audit.previous.matching()) + 1
    assert audit.gates()["finite_evanescent_matching_not_claimed_from_pole"]
