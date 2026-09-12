"""Independent joint complex domains, frames and complete physical readouts."""

from functools import cache
from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_cd_metric_noise import stress
from p8_vacuum_affine_reference_spatial_analyticity import (
    audit,
    boundary,
    domain,
    frame,
    readouts,
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


def configuration(K, center, outer=False):
    A = mp.mpf(25) / 16
    k = mp.mpf(K) * mp.matrix([mp.mpf(3) / 5, 0, mp.mpf(4) / 5])
    nu = mp.sqrt(1000**2 + (k.T * k)[0] / A**2)
    direction = mp.matrix(
        [mp.mpc("0.4", "0.2"), mp.mpc("-0.3", "0.1"), mp.mpc("0.5", "-0.4")]
    )
    direction /= mp.norm(direction)
    p = mp.mpf("1e-6") * nu * direction
    radius = mp.mpf("0.00009") if outer else mp.mpf("0.00004")
    u = mp.mpf(center) + radius * mp.exp(mp.mpc(0, "0.7"))
    return k, nu, p, u


@pytest.mark.parametrize("K", ("1000", "1e6", "1e40"))
@pytest.mark.parametrize("center", ("-0.5", "0", "0.5"))
def test_independent_complete_joint_complex_frequency_domain(K, center):
    with mp.workdps(90):
        k, nu, p, u = configuration(K, center, outer=True)
        ell, _, _, _, _, _ = joint_frame(k, p)
        a0 = (1 + mp.mpf(center) ** 2) ** 2
        baseline = 1000**2 + (k.T * k)[0] / a0**2
        ew = mp.mpf(
            str(
                s.N(
                    domain.constants()["complete_relative_squared_frequency_defect"], 95
                )
            )
        )
        for momentum in (k, ell):
            for kind in ("transverse", "longitudinal"):
                _, _, omega, z, W, _ = W_data(kind, u, momentum)
                assert abs(omega**2 - baseline) / baseline <= ew
                assert mp.re(omega) > mp.mpf("0.99") * nu and abs(omega) < 2 * nu
                assert abs(z) < mp.mpf("1.01")
                assert abs(W / omega - 1) < mp.mpf("0.001")
                assert mp.re(W) > nu / 2 and abs(W) < 3 * nu


@pytest.mark.parametrize("K", ("1000", "1e6", "1e40"))
def test_independent_complex_polarization_frame_and_true_Euclidean_norms(K):
    with mp.workdps(90):
        k, nu, p, _ = configuration(K, "0")
        ell, rho, n, R, e1, e2 = joint_frame(k, p)
        n0 = -k / mp.norm(k)
        assert mp.norm(R.T * R - mp.eye(3)) < mp.mpf("1e-80")
        assert mp.norm(R * n0 - n) < mp.mpf("1e-80")
        assert abs(mp.det(R) - 1) < mp.mpf("1e-80")
        assert abs((n.T * n)[0] - 1) < mp.mpf("1e-80")
        assert abs((ell.T * e1)[0]) < mp.mpf("1e-75") * nu
        assert abs((ell.T * e2)[0]) < mp.mpf("1e-75") * nu
        assert mp.norm(n - n0) < mp.mpf("1e-5")
        assert mp.norm(R - mp.eye(3)) < mp.mpf("1e-4")
        assert mp.norm(e1) < 2 and mp.norm(e2) < 2 and mp.norm(n) < 2
        assert abs(rho) > mp.mpf("0.99") * mp.norm(k)


def test_complex_bilinear_unit_does_not_mean_Euclidean_unit():
    with mp.workdps(70):
        t = mp.mpf("0.2")
        n = mp.matrix([1j * mp.sinh(t), 0, mp.cosh(t)])
        assert abs((n.T * n)[0] - 1) < mp.mpf("1e-60")
        assert mp.norm(n) > 1


@pytest.mark.parametrize("K", ("1000", "1e6", "1e40"))
@pytest.mark.parametrize("center", ("-0.5", "0", "0.5"))
@pytest.mark.parametrize("phase_sign", (-1, 1))
def test_independent_full_constrained_readouts_for_both_Schwarz_signs(
    K, center, phase_sign
):
    with mp.workdps(90):
        k, nu, p, u = configuration(K, center)
        ell, _, n, _, e1, e2 = joint_frame(k, p)
        for vector, W in physical_modes(u, ell, n, e1, e2, phase_sign):
            assert max(abs(x) for x in vector) < 1000 * mp.sqrt(nu)
            assert mp.norm(vector) < 4000 * mp.sqrt(nu)
            assert mp.re(W) > nu / 2


@pytest.mark.parametrize("K", ("1000", "1e6", "1e40"))
def test_independent_Schwarz_continuation_preserves_full_field_order(K):
    with mp.workdps(85):
        k, nu, p, u = configuration(K, "0.5")
        ell, _, n, _, e1, e2 = joint_frame(k, p)
        positive = physical_modes(u, ell, n, e1, e2, 1)
        ellc, _, nc, _, e1c, e2c = joint_frame(k, conjugate(p))
        negative = physical_modes(mp.conj(u), ellc, nc, e1c, e2c, -1)
        for (v, _), (w, _) in zip(positive, negative):
            assert mp.norm(v - conjugate(w)) < mp.mpf("1e-70") * mp.sqrt(nu)


@cache
def complete_stress_matrices():
    def real_number(x):
        n, d = s.fraction(x)
        return mp.mpf(int(n)) / int(d)

    return [
        mp.matrix([[real_number(M[i, j]) for j in range(10)] for i in range(10)])
        for M in stress.data()["quadratic_matrices"].values()
    ]


@pytest.mark.parametrize("K", ("1000", "1e6", "1e40"))
@pytest.mark.parametrize("center", ("-0.5", "0.5"))
def test_all_nine_full_stress_pair_tensors_and_inverse_phases(K, center):
    with mp.workdps(80):
        k, nu, p, u = configuration(K, center)
        ell, _, n, _, e1, e2 = joint_frame(k, p)
        nk = k / mp.norm(k)
        ek1, ek2 = real_frame(nk)
        first = physical_modes(u, k, nk, ek1, ek2)
        second = physical_modes(u, ell, n, e1, e2)
        matrices = complete_stress_matrices()
        for (v, W), (w, V) in product(first, second):
            pair = mp.matrix([(v.T * M * w)[0] for M in matrices])
            assert mp.norm(pair) < mp.mpf("1e9") * nu
            assert abs(1 / (W + V)) < 1 / nu


@pytest.mark.parametrize("module", (domain, frame, readouts, boundary))
def test_all_exact_scientific_packets(module):
    d = module.data()
    for value in d["checks"].values():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.cancel(x) == 0 for x in entries)
    assert all(d["gates"].values())


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


def complete_leading_endpoint(K, parameter):
    k = mp.mpf(K) * mp.matrix([mp.mpf(3) / 5, 0, mp.mpf(4) / 5])
    nu = mp.sqrt(1000**2 + (k.T * k)[0] / (mp.mpf(25) / 16) ** 2)
    direction = mp.matrix([mp.mpf(2) / 3, mp.mpf(1) / 3, mp.mpf(2) / 3])
    p = mp.mpf("1e-6") * nu * parameter * direction
    ell, _, n, _, e1, e2 = joint_frame(k, p)
    nk = k / mp.norm(k)
    ek1, ek2 = real_frame(nk)
    u = mp.mpf("0.25")
    D = mp.diag([1, -2, 1])
    G = D.copy()
    G[0, 1] = G[1, 0] = mp.mpf(1) / 3
    MD = contracted_spatial_matrix(D)
    MG = contracted_spatial_matrix(G)
    plus_k = physical_modes(u, k, nk, ek1, ek2, 1)
    plus_l = physical_modes(u, ell, n, e1, e2, 1)
    minus_k = physical_modes(u, k, nk, ek1, ek2, -1)
    minus_l = physical_modes(u, ell, n, e1, e2, -1)
    result = mp.mpc(0)
    for first in range(3):
        for second in range(3):
            vk, Wk = minus_k[first]
            vl, Wl = minus_l[second]
            detector = (plus_k[first][0].T * MD * plus_l[second][0])[0]
            source = (vk.T * MG * vl)[0]
            result += -1j * detector * source / (Wk + Wl)
    return result, nu, mp.norm(D) * mp.norm(G)


@pytest.mark.parametrize("K", ("1000", "1e6", "1e40"))
def test_actual_full_nine_pair_endpoint_spatial_Taylor_remainder(K):
    with mp.workdps(105):
        function = lambda parameter: complete_leading_endpoint(K, parameter)[0]
        coefficients = mp.taylor(function, mp.mpf(0), 4)
        errors = []
        for parameter in (mp.mpf("0.25"), mp.mpf("0.5")):
            value, nu, norms = complete_leading_endpoint(K, parameter)
            error = abs(value - mp.polyval(list(reversed(coefficients)), parameter))
            P = mp.mpf("1e-6") * nu * parameter
            assert 0 < error < mp.mpf("4e57") * P**5 / nu**4 * norms
            errors.append(error)
        # The actual full endpoint has a nonzero fifth coefficient; this is
        # stronger than merely comparing zero against the conservative bound.
        assert 30 < errors[1] / errors[0] < 34


@pytest.mark.parametrize("K", ("1000", "1e6", "1e40"))
def test_Schwarz_continuation_is_not_naive_in_place_conjugation(K):
    with mp.workdps(80):
        value, nu, _ = complete_leading_endpoint(K, mp.mpc("0.3", "0.2"))
        reflected, _, _ = complete_leading_endpoint(K, mp.mpc("0.3", "-0.2"))
        # -i times a real-axis-real analytic contraction has the odd Schwarz
        # relation. Conjugating at the same complex argument breaks it.
        assert abs(value + mp.conj(reflected)) < mp.mpf("1e-65") * nu
        assert abs(value + mp.conj(value)) > mp.mpf("1e-20") * nu


@pytest.mark.parametrize("K", ("1000", "1e6", "1e12"))
def test_independent_fourth_weight_radial_integrals_and_half_band_tail(K):
    with mp.workdps(75):
        A = mp.mpf(25) / 16
        m = mp.mpf(1000)
        K = mp.mpf(K)
        prefactor = A**3 / (2 * mp.pi**2 * m)
        full = prefactor * mp.quad(lambda theta: mp.sin(theta) ** 2, [0, mp.pi / 2])
        assert mp.almosteq(full, A**3 / (8 * mp.pi * m))
        lower = mp.atan(K / (2 * A * m))
        tail = prefactor * mp.quad(lambda theta: mp.sin(theta) ** 2, [lower, mp.pi / 2])
        assert 0 < tail < A**4 / (9 * K)
        assert mp.mpf("4e57") * full < mp.mpf("1e54")
        assert mp.mpf("4e57") * tail < mp.mpf("1e58") / K


@pytest.mark.parametrize("K", ("1000", "1e6"))
@pytest.mark.parametrize("transfer_fraction", ("0", "0.25", "1"))
def test_far_removed_two_leg_geometry_and_actual_angular_radial_integral(
    K, transfer_fraction
):
    with mp.workdps(75):
        A = mp.mpf(25) / 16
        m = mp.mpf(1000)
        delta = mp.mpf("1e-6")
        K = mp.mpf(K)
        P = delta * K * mp.mpf(transfer_fraction)
        far_min = m
        if 2 * P / delta > m:
            far_min = max(m, A * mp.sqrt((2 * P / delta) ** 2 - m * m))
        lower = mp.atan(far_min / (A * m))
        prefactor = A**3 / (2 * mp.pi**2 * m)

        def fraction(theta):
            if theta == mp.pi / 2:
                return mp.mpf(1)
            k = A * m * mp.tan(theta)
            if k > K:
                return mp.mpf(1)
            if P == 0:
                return mp.mpf(0)
            c = (k * k + P * P - K * K) / (2 * k * P)
            return (min(mp.mpf(1), max(mp.mpf(-1), c)) + 1) / 2

        splits = [lower]
        splits.extend(
            mp.atan(k / (A * m)) for k in sorted({K - P, K, K + P}) if k > far_min
        )
        splits.append(mp.pi / 2)
        actual = prefactor * mp.quad(
            lambda theta: mp.sin(theta) ** 2 * fraction(theta), splits
        )
        half_band = prefactor * mp.quad(
            lambda theta: mp.sin(theta) ** 2, [mp.atan(K / (2 * A * m)), mp.pi / 2]
        )
        assert 0 < actual <= half_band
        assert mp.mpf("4e57") * actual < mp.mpf("1e58") / K


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
    assert len(audit.residuals()) == 24 and audit.scalar_entry_count() == 53
    assert len(audit.gates()) == 46 and all(
        value is True for value in audit.gates().values()
    )
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 124


def test_original_primitive_and_matching_statuses_retained():
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert "NOT_NEAR_REGION_LOCAL_MATCHING_INVERSE_OR_V_G_B" in audit.ITEM["status"]
    audit.validate_scope(audit.frontier(), audit.matching())
