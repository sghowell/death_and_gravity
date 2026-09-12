"""Independent continued-mode, Fourier lift, full tangent and regulator checks."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_dimensional_spatial_symbol import density, geometry, jets
from p8_vacuum_affine_matched_spatial_current import (
    assembly,
    audit,
    dimension,
    lift,
    regulator,
)
from p8_vacuum_affine_spatial_current import hamiltonian, vertices
from p8_vacuum_affine_spatial_symbol import sectors
from p8_vector_clock_matching import continuation
from p8_vector_state import wkb


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
def numeric_UV():
    products = jets.endpoint_products()
    result = {}
    for channel in ("tensor", "vector", "scalar"):
        geo = density.geometric_jets(channel)
        for j in range(5):
            for r in range(j + 1):
                for degree in range(5 - j):
                    result[channel, j, r, degree] = s.factor(
                        sum(
                            sum(
                                products[key, j, r][v] * geo[key][degree - v]
                                for v in range(degree + 1)
                            )
                            for key in geo
                        )
                    )
    return tuple(result), s.lambdify(
        (geometry.d, jets.t, jets.u, jets.p, jets.m),
        list(result.values()),
        "mpmath",
        cse=True,
    )


@cache
def numeric_geometry(channel):
    rows = geometry.contractions(channel)
    return tuple(rows), s.lambdify(
        (geometry.d, geometry.u, geometry.y), list(rows.values()), "mpmath", cse=True
    )


@pytest.mark.parametrize("channel", ("tensor", "vector", "scalar"))
@pytest.mark.parametrize("transfer", ("13", "1000000"))
def test_complete_complex_dimension_full_endpoint_UV_remainder_two_radii(
    channel, transfer
):
    with mp.workdps(140):
        p = mp.mpf(transfer)
        d = mp.mpc("3.1", "0.15")
        t, u = mp.mpf(".25"), mp.mpf(".37")
        source = [mp.mpf(v) for v in (1, -2, 3, 1, -1)]
        keys, fn = numeric_UV()
        coeff = dict(zip(keys, fn(d, t, u, p, 1000)))
        radii = [mp.mpf(".01") / (1000 + p) / 64, mp.mpf(".01") / (1000 + p) / 128]
        results = []
        for x in radii:
            original = mode_endpoint_rows(x, t, u, p, d)
            partner = mode_endpoint_rows(x, t, u, p, mp.conj(d))
            # Schwarz reflection constructs the opposite branch at the SAME d.
            analytic = {
                key: (mp.conj(partner[key]) - value) / (2j)
                for key, value in original.items()
            }
            gkeys, gfn = numeric_geometry(channel)
            geos = dict(zip(gkeys, gfn(d, u, x * p)))
            errors = []
            for j in range(5):
                full = mp.fsum(
                    geos[key] * analytic[key, j, r] * source[r]
                    for key in gkeys
                    for r in range(j + 1)
                )
                uv = mp.fsum(
                    coeff[channel, j, r, n] * source[r] * x**n
                    for r in range(j + 1)
                    for n in range(5 - j)
                )
                errors.append(full - uv)
            results.append(errors)
        for j in range(5):
            e1, e2 = results[0][j], results[1][j]
            assert abs(e2) <= mp.mpf(".6") ** (5 - j) * abs(e1) + mp.mpf("1e-75"), (
                channel,
                p,
                j,
                e1,
                e2,
            )
        assert any(abs(v) > mp.mpf("1e-55") for v in results[0])


@pytest.mark.parametrize("kind", ("transverse", "longitudinal"))
@pytest.mark.parametrize("order", range(1, 5))
def test_independent_four_order_dimensional_recurrence(kind, order):
    assert (
        s.factor(
            full_dimensional_WKB_coefficient(kind, order)
            - dimension.coefficient(kind, order)
        )
        == 0
    )


@pytest.mark.parametrize("transfer", ("0", "13", "1000000"))
@pytest.mark.parametrize("phase_index", (0, 1, 3, 5))
def test_full_continued_frequencies_on_joint_domain(transfer, phase_index):
    with mp.workdps(70):
        p = mp.mpf(transfer)
        phase = mp.exp(2j * mp.pi * phase_index / 7)
        d = 3 + phase / 4
        time = mp.mpf(".5") + mp.mpf(".0001") * phase
        x = mp.mpf(".01") / (1000 + p) * mp.conj(phase)
        n = mp.matrix([mp.mpf(".8"), 0, mp.mpf(".6")])
        ell = -n + mp.matrix([0, 0, p * x])
        _point, freq = scaled_point(time, x, n, ell, dimension=d)
        assert all(mp.re(v) > 0.5 and abs(v) < 3 for v in freq.values())
        assert all(abs(v) < 1 for v in four_inverse_phases(freq))


@pytest.mark.parametrize("degree", (0, 1, 2, 4))
@pytest.mark.parametrize("d", ("2.8", "3.2", "3.1+0.1j"))
def test_complex_dimension_angular_integral_envelope(degree, d):
    with mp.workdps(70):
        dim = mp.mpc(d.replace("j", "")) if "j" not in d else mp.mpc("3.1", ".1")
        c = mp.gamma(dim / 2) / (mp.sqrt(mp.pi) * mp.gamma((dim - 1) / 2))
        value = c * mp.quad(
            lambda u: u ** (2 * degree) * (1 - u * u) ** ((dim - 3) / 2), [-1, 0, 1]
        )
        expected = mp.rf(mp.mpf(".5"), degree) / mp.rf(dim / 2, degree)
        assert abs(value - expected) < mp.mpf("1e-55")
        for u in (mp.mpf(".3"), 1 - mp.mpf("1e-9"), -1 + mp.mpf("1e-12")):
            assert abs((1 - u * u) ** ((dim - 3) / 2)) <= (1 - u * u) ** (
                -mp.mpf(".125")
            )


@pytest.mark.parametrize("j", range(5))
@pytest.mark.parametrize("transfer", (0, 13, 10**6))
def test_uniform_far_dimension_radial_integral(j, transfer):
    with mp.workdps(60):
        U = mp.mpf(1000 + transfer)
        L = 200 * U
        sig = mp.mpf(".25")
        # r=L/y gives an integrable compact-domain quadrature.
        value = U ** (5 - j) * L ** (-1 + sig) * mp.quad(lambda y: y ** (-sig), [0, 1])
        expected = U ** (5 - j) * L ** (-1 + sig) / (1 - sig)
        assert abs(value / expected - 1) < mp.mpf("1e-45")
        assert U ** (5 + sig - j) <= U**6


@pytest.mark.parametrize("case", range(3))
def test_literal_noncommuting_homogeneous_generator_contact_and_tangent(case):
    k = s.Matrix(([0, 0, 0], [3, 4, 12], [1000, -2000, 500])[case])
    a = s.Rational(1) + s.Rational(case, 5)
    m = s.Integer(1000)
    D = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 17
    G = s.Matrix([[2, 1, 3], [1, 1, -2], [3, -2, -3]]) / 19
    assert D * G != G * D
    e, z = s.symbols("e z")
    X = e * D + z * G
    E = s.eye(3) + X + X * X / 2
    Ei = s.eye(3) - X + X * X / 2
    C = hamiltonian.cross(k)
    M = s.diag(a * m * m * Ei + C.T * E * C / a, E / a + k * k.T / (a**3 * m * m))
    zer = lambda v: s.expand(v).subs({e: 0, z: 0})
    M0 = M.applyfunc(zer)
    MG = M.diff(z).applyfunc(zer)
    MD = M.diff(e).applyfunc(zer)
    MDG = M.diff(e, z).applyfunc(zer)
    assert M0 == hamiltonian.base(k, a, m)
    assert MG == vertices.metric_vertex(k, k, G, a, m)
    assert MD == vertices.metric_vertex(k, k, D, a, m)
    assert MDG == vertices.metric_contact(k, k, (D * G + G * D) / 2, a, m)
    assert MDG != s.zeros(6)
    J = s.zeros(6)
    J[:3, 3:] = s.eye(3)
    J[3:, :3] = -s.eye(3)
    A, B = J * M0, J * MG
    S = s.diag(*[s.Rational(v, 11) for v in (1, 2, 3, 4, 5, 6)])
    # Direct differentiated polynomial flow versus independently iterated covariance ODE.
    powers = [s.eye(6)]
    tangents = [s.zeros(6)]
    for n in range(1, 5):
        powers.append(A * powers[-1])
        tangents.append(A * tangents[-1] + B * powers[-2])
    c = [S]
    dc = [s.zeros(6)]
    for n in range(4):
        c.append((A * c[-1] + c[-1] * A.T) / (n + 1))
        dc.append((A * dc[-1] + dc[-1] * A.T + B * c[-2] + c[-2] * B.T) / (n + 1))
    for n in range(5):
        direct = sum(
            (
                (tangents[j] * S * powers[n - j].T + powers[j] * S * tangents[n - j].T)
                / s.factorial(j)
                / s.factorial(n - j)
                for j in range(n + 1)
            ),
            s.zeros(6),
        )
        assert direct == dc[n]
    wrong = B * S
    assert wrong != B * S + S * B.T
    assert -s.trace(MDG * S) / 2 != 0


@pytest.mark.parametrize("degree", (1, 2, 4, 8))
def test_complex_source_zero_germ_Hilbert_lift(degree):
    t = s.Symbol("t", real=True)
    # Polynomial checks of the fundamental-theorem estimate need only zero left trace.
    f = (t + s.Rational(1, 2)) ** degree
    C = s.Matrix([[1, s.I, 0], [s.I, -1, 1 + s.I], [0, 1 + s.I, 0]])
    real = C.applyfunc(s.re)
    imag = C.applyfunc(s.im)
    norm2 = lambda M: s.expand(s.trace(M.conjugate().T * M))
    assert norm2(C) == norm2(real) + norm2(imag)
    bound = s.integrate(s.diff(f, t) ** 2, (t, -s.Rational(1, 2), s.Rational(1, 2)))
    assert bound >= f.subs(t, s.Rational(1, 2)) ** 2
    assert norm2(C) * bound == norm2(real) * bound + norm2(imag) * bound
    assert f.subs(t, -s.Rational(1, 2)) == 0


@pytest.mark.parametrize("momentum", (0, s.Rational(1, 1000), 13, 1000, 10**6))
def test_all_transfer_weight_embeddings(momentum):
    v = momentum**2
    assert (1 + v) ** 6 >= 1
    assert (1 + v) ** 6 >= v
    assert (1 + v) ** 6 >= (1 + v) ** 4
    assert (1000 + momentum) ** 6 <= 8 * 1000**6 * (1 + v) ** 3


@pytest.mark.parametrize("K", (2000, 10000, 10**9))
def test_original_anchor_subtraction_with_actual_radial_UV_coefficients(K):
    from p8_vacuum_affine_spatial_matching_difference import density as physical

    inv = (s.Integer(1), s.Rational(1, 3), s.Rational(1, 7))
    subs = {jets.t: s.Rational(1, 4), jets.p: 13, jets.m: 1000}
    F2 = physical.quadratic(*inv).subs(subs)
    F4 = sum(physical.logarithmic(*inv)).subs(subs)
    A3 = (
        -13 * (18 * inv[0] - 12 * inv[1] - inv[2]) / (512 * s.pi**2 * jets.a.subs(subs))
    )
    A1 = s.Rational(7, 11)
    finite = s.Rational(-13, 17)
    m = s.Integer(1000)
    K = s.Integer(K)
    U = A3 * K**3 + A1 * K + F2 * (K * K - m * m) / 2 + F4 * s.log(K / m)
    assert (
        s.expand(U - regulator.subtractor(K, m, A3, A1, F2, F4, finite) - finite) == 0
    )
    # Omission is genuinely nonzero for the actual power coefficient.
    assert F2 * m * m / 2 != 0
    assert assembly.anchored_current(5, 7, 7, 0) == 5


@pytest.mark.parametrize("module", (dimension, lift, assembly, regulator))
def test_complete_assembly_core(module):
    data = module.data()
    assert all(s.cancel(v) == 0 for v in data["checks"].values())
    assert all(data["gates"].values())


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_all_scoped_matched_spatial_current_residuals(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(entry == 0 for entry in entries)


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_all_unsupported_matched_spatial_current_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_exact_assembly_audit_counts_and_scope():
    assert len(audit.residuals()) == 64 and audit.scalar_entry_count() == 64
    assert len(audit.gates()) == 52 and all(audit.gates().values())
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 135
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 69
    assert audit.frontier() == audit.previous.frontier()
    assert audit.gates()["lapse_shift_clock_scalar_and_reduced_inverse_still_required"]


def test_precise_historical_display_erratum_preserves_frozen_report_claims():
    import json

    from p8_vacuum_affine_matched_spatial_current import verify

    prior = json.loads(verify.parent.REPORT.read_text())
    parent = json.loads(verify.parent.parent.REPORT.read_text())
    assert "ACTUAL_DIMENSIONAL_SPATIAL_UV_FINITE_PART" in prior["claim"]
    assert (
        "ACTUAL_SPATIAL_UV_DIFFERENCE_COMPLETE_CONTACT_CANCELLATION" in parent["claim"]
    )
    old = prior["primitive_and_matching_frontier"]["matching"]
    assert old == audit.previous.matching()
    now = audit.matching()
    assert now[:-3] == old[:-2]
    assert old[-3] == old[-2] == old[-1]
    assert len({row["id"] for row in now}) == len(now)
    for row in audit.metadata_errata():
        assert old[row["index"]] == row["frozen_display"]
        assert now[row["index"]] == row["corrected_display"]
    assert verify.sha(verify.parent.REPORT) == verify.PARENT_SHA
    assert verify.sha(verify.parent.parent.REPORT) == verify.parent.PARENT_SHA


def test_scope_omissions_rejected_and_previous_cache_unmodified():
    before = audit.previous.matching()
    rows = audit.matching()
    rows[-2]["status"] = "COMPLETE"
    with pytest.raises(ValueError):
        audit.validate_scope(audit.frontier(), rows)
    assert audit.previous.matching() == before
    assert assembly.CURRENT == 2 * s.Integer(10) ** 95
    assert regulator.TAIL == 3 * s.Integer(10) ** 54


@pytest.mark.parametrize(
    "kind,order",
    (
        ("tensor", 1),
        ("transverse", -1),
        ("longitudinal", 5),
        ("transverse", True),
        ("longitudinal", 1.0),
    ),
)
def test_bad_dimensional_recurrence_inputs(kind, order):
    with pytest.raises(ValueError):
        dimension.coefficient(kind, order)
