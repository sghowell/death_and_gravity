"""Independent actual mixed endpoint coefficients, finite moments and original tails."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_cd_metric_noise import stress
from p8_vacuum_affine_finite_endpoint_polynomial import (
    audit,
    degrees,
    limit,
    radials,
    tail,
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


def complete_endpoint_rows(K, parameter, samples=16):
    k = mp.mpf(K) * mp.matrix([mp.mpf(3) / 5, 0, mp.mpf(4) / 5])
    nu = mp.sqrt(1000**2 + (k.T * k)[0] / (mp.mpf(25) / 16) ** 2)
    direction = mp.matrix([mp.mpf(2) / 3, mp.mpf(1) / 3, mp.mpf(2) / 3])
    P = mp.mpf("1e-6") * nu * parameter * direction
    ell, _, n, _, e1, e2 = joint_frame(k, P)
    nk = k / mp.norm(k)
    ek1, ek2 = real_frame(nk)
    center = mp.mpf("0.25")
    radius = mp.mpf("1e-5")
    D = mp.diag([1, -2, 1])
    G = D.copy()
    G[0, 1] = G[1, 0] = mp.mpf(1) / 3
    MD = contracted_spatial_matrix(D)
    MG = contracted_spatial_matrix(G)
    detector_k = physical_modes(center, k, nk, ek1, ek2, 1)
    detector_l = physical_modes(center, ell, n, e1, e2, 1)
    detector = [
        (detector_k[a][0].T * MD * detector_l[b][0])[0]
        for a in range(3)
        for b in range(3)
    ]
    phases = [mp.exp(2j * mp.pi * l / samples) for l in range(samples)]
    gs = []
    amplitudes = []
    for phase in phases:
        u = center + radius * phase
        first = physical_modes(u, k, nk, ek1, ek2)
        second = physical_modes(u, ell, n, e1, e2)
        gs.append(
            [1 / (first[a][1] + second[b][1]) for a in range(3) for b in range(3)]
        )
        amplitudes.append(
            [(first[a][0].T * MG * second[b][0])[0] for a in range(3) for b in range(3)]
        )
    result = mp.zeros(5, 5)
    for pair in range(9):
        g = [
            mp.fsum(gs[l][pair] * phases[l] ** (-j) for l in range(samples))
            / samples
            / radius**j
            for j in range(5)
        ]
        amp = [
            mp.fsum(amplitudes[l][pair] * phases[l] ** (-j) for l in range(samples))
            / samples
            / radius**j
            for j in range(5)
        ]
        for source_jet in range(5):
            current = [
                amp[l - source_jet] / mp.factorial(source_jet)
                if l >= source_jet
                else mp.mpc(0)
                for l in range(5)
            ]
            for order in range(5):
                result[order, source_jet] += (
                    (-1j) * (1j) ** order * detector[pair] * g[0] * current[0]
                )
                current = [
                    (l + 1) * mp.fsum(g[j] * current[l + 1 - j] for j in range(l + 2))
                    for l in range(len(current) - 1)
                ]
    return result, nu, mp.norm(D) * mp.norm(G)


@cache
def actual_mixed_coefficients(K, points=12):
    with mp.workdps(135):
        radius = mp.mpf("0.5")
        phases = [mp.exp(2j * mp.pi * l / points) for l in range(points)]
        values = [complete_endpoint_rows(K, radius * phase) for phase in phases]
        nu = values[0][1]
        norms = values[0][2]
        rows = {}
        for j, n in degrees.finite_cells():
            row = mp.matrix(
                [
                    mp.fsum(
                        values[l][0][j, r] * phases[l] ** (-n) for l in range(points)
                    )
                    / points
                    / radius**n
                    / (mp.mpf("1e-6") * nu) ** n
                    for r in range(j + 1)
                ]
            )
            rows[j, n] = row
        return rows, nu, norms


def exact_mp(value):
    numerator, denominator = s.fraction(value)
    return mp.mpf(int(numerator)) / int(denominator)


@pytest.mark.parametrize("K", ("1000", "1e6", "1e12"))
@pytest.mark.parametrize("j,n", degrees.finite_cells())
def test_all_ten_actual_spatial_time_coefficients_and_full_time_jet_rows(K, j, n):
    with mp.workdps(130):
        rows, nu, norms = actual_mixed_coefficients(K)
        refined, _, _ = actual_mixed_coefficients(K, 16)
        row = rows[j, n]
        reference = refined[j, n]
        scale = max(mp.mpf("1e-300"), mp.norm(reference))
        assert mp.norm(row - reference) < mp.mpf("1e-45") * scale
        assert (
            0
            < mp.norm(reference)
            < exact_mp(degrees.row_majorant(j, n)) * nu ** (1 - j - n) * norms
        )
        for r, value in enumerate(reference):
            assert (
                abs(value)
                < exact_mp(degrees.jet_majorant(j, n, r)) * nu ** (1 - j - n) * norms
            )


@pytest.mark.parametrize("order", range(4, 8))
def test_independent_complete_massive_radial_moments(order):
    with mp.workdps(85):
        A = mp.mpf(25) / 16
        m = mp.mpf(1000)
        integral = (
            A**3
            * m ** (3 - order)
            / (2 * mp.pi**2)
            * mp.quad(
                lambda theta: mp.sin(theta) ** 2 * mp.cos(theta) ** (order - 4),
                [0, mp.pi / 2],
            )
        )
        expected = mp.mpf(str(s.N(radials.moment(order), 100)))
        assert abs(integral - expected) < mp.mpf("1e-70") * expected
        assert integral < exact_mp(radials.upper(order))


@pytest.mark.parametrize("module", (degrees, radials, tail, limit))
def test_all_exact_core_packets(module):
    d = module.data()
    for value in d["checks"].values():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.cancel(x) == 0 for x in entries)
    assert all(d["gates"].values())


def removed_fraction(k, P, K):
    if k > K:
        return mp.mpf(1)
    if P == 0:
        return mp.mpf(0)
    c = (k * k + P * P - K * K) / (2 * k * P)
    return (min(mp.mpf(1), max(mp.mpf(-1), c)) + 1) / 2


@cache
def removed_moments(K, fraction):
    with mp.workdps(100):
        K = mp.mpf(K)
        P = mp.mpf(fraction) * K
        A = mp.mpf(25) / 16
        m = mp.mpf(1000)
        radii = [m] + [x for x in sorted({abs(K - P), K, K + P}) if x > m]
        limits = [mp.atan(x / (A * m)) for x in radii] + [mp.pi / 2]
        out = {}
        for order in range(4, 8):

            def integrand(theta, order=order):
                if theta == mp.pi / 2:
                    fraction = mp.mpf(1)
                else:
                    fraction = removed_fraction(A * m * mp.tan(theta), P, K)
                return mp.sin(theta) ** 2 * mp.cos(theta) ** (order - 4) * fraction

            out[order] = (
                A**3 * m ** (3 - order) / (2 * mp.pi**2) * mp.quad(integrand, limits)
            )
        return out, P, K


@pytest.mark.parametrize("K", ("1000", "1e6"))
@pytest.mark.parametrize("fraction", ("0", "0.25", "0.5", "1", "2", "100"))
@pytest.mark.parametrize("j,n", degrees.finite_cells())
def test_each_finite_cell_original_removed_union_and_same_spatial_norm(
    K, fraction, j, n
):
    with mp.workdps(90):
        moments, P, K = removed_moments(K, fraction)
        actual = exact_mp(degrees.row_majorant(j, n)) * P**n * moments[j + n - 1]
        if P <= K / 2:
            expected = exact_mp(tail.low_transfer_cell(j, n)) * P**n / K
        else:
            expected = 2 * exact_mp(radials.finite_cell_bound(j, n)) * P ** (n + 1) / K
        assert 0 <= actual <= expected
        coefficient = tail.low_transfer_cell(j, n) + 2 * radials.finite_cell_bound(j, n)
        assert actual <= exact_mp(coefficient) * (1 + P * P) ** 3 / K


@pytest.mark.parametrize("order", range(4, 8))
@pytest.mark.parametrize("K", ("1000", "1e6", "1e12"))
def test_independent_half_band_radial_tail_and_mass_to_one_over_K(order, K):
    with mp.workdps(90):
        A = mp.mpf(25) / 16
        m = mp.mpf(1000)
        K = mp.mpf(K)
        lower = mp.atan(K / (2 * A * m))
        value = (
            A**3
            * m ** (3 - order)
            / (2 * mp.pi**2)
            * mp.quad(
                lambda theta: mp.sin(theta) ** 2 * mp.cos(theta) ** (order - 4),
                [lower, mp.pi / 2],
            )
        )
        power = (
            A**order
            * 2 ** (order - 3)
            * K ** (3 - order)
            / (2 * mp.pi**2 * (order - 3))
        )
        displayed = (
            A**order * 2 ** (order - 3) * m ** (4 - order) / (18 * (order - 3) * K)
        )
        assert 0 < value <= power < displayed


def test_borderline_radial_majorant_is_not_a_current_divergence_claim():
    k, A, m = s.symbols("k A m", positive=True)
    majorant = k * k / (m * m + k * k / (A * A)) ** s.Rational(3, 2)
    assert s.limit(k * majorant, k, s.oo) == A**3
    assert s.integrate(1 / k, (k, m, s.oo)) == s.oo
    assert all(j + n <= 4 for j, n in degrees.ultraviolet_cells())
    # This concerns only the positive majorant. Cancellations in the actual
    # complete current and its contact have not been evaluated here.


BAD_ORDER = (
    True,
    False,
    4.0,
    s.Float(4),
    "4",
    None,
    s.oo,
    -s.oo,
    s.zoo,
    s.I,
    s.nan,
    s.Symbol("n"),
    3,
    8,
)


@pytest.mark.parametrize("value", BAD_ORDER)
def test_radial_order_validation_precedes_prewarmed_cache(value):
    for order in range(4, 8):
        radials.moment(order)
    with pytest.raises(ValueError):
        radials.moment(value)
    with pytest.raises(ValueError):
        radials.upper(value)


@pytest.mark.parametrize(
    "value", (True, False, 1.0, s.Float(1), "1", None, s.oo, s.I, s.Symbol("j"), -1, 5)
)
@pytest.mark.parametrize("position", (0, 1))
def test_exact_endpoint_and_spatial_order_guards(value, position):
    arguments = [1, 4]
    arguments[position] = value
    with pytest.raises(ValueError):
        degrees.row_majorant(*arguments)
    with pytest.raises(ValueError):
        radials.finite_cell_bound(*arguments)
    with pytest.raises(ValueError):
        tail.low_transfer_cell(*arguments)


@pytest.mark.parametrize("j,n", degrees.ultraviolet_cells())
def test_every_candidate_UV_cell_rejected_by_finite_only_apis(j, n):
    assert degrees.row_majorant(j, n) > 0
    with pytest.raises(ValueError):
        radials.finite_cell_bound(j, n)
    with pytest.raises(ValueError):
        tail.low_transfer_cell(j, n)


@pytest.mark.parametrize(
    "value", (True, False, 1.0, s.Float(1), "1", None, s.oo, -1, 5)
)
def test_source_time_jet_order_guards(value):
    with pytest.raises(ValueError):
        degrees.jet_majorant(4, 4, value)


def test_all_finite_cells_vanish_at_exact_zero_transfer():
    assert all(n >= 1 for j, n in degrees.finite_cells())
    assert (
        sum(
            radials.finite_cell_bound(j, n) * s.Integer(0) ** n
            for j, n in degrees.finite_cells()
        )
        == 0
    )


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_every_exact_residual(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(s.cancel(x) == 0 for x in entries)


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_each_invalid_scope_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_complete_counts_and_true_gates():
    assert len(audit.residuals()) == 29 and audit.scalar_entry_count() == 37
    assert len(audit.gates()) == 39 and all(
        value is True for value in audit.gates().values()
    )
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 127


def test_original_primitive_and_matching_statuses_retained():
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert "NOT_REMAINING_UV_CONTACT_MATCHING_INVERSE_OR_V_G_B" in audit.ITEM["status"]
    audit.validate_scope(audit.frontier(), audit.matching())
