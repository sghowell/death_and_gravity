"""Independent actual endpoint time jets, complete near radials and band geometry."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_cd_metric_noise import stress
from p8_vacuum_affine_reference_spatial_remainder import audit, limit, near, real, tail
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
        numerator, denominator = s.fraction(x)
        return mp.mpf(int(numerator)) / int(denominator)

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


def real_modes(u, momentum):
    norm = mp.norm(momentum)
    n = momentum / norm if norm else mp.matrix([0, 0, 1])
    e1, e2 = real_frame(n)
    return physical_modes(u, momentum, n, e1, e2)


def actual_endpoint_jet_row(K, P, center="0.25", samples=16):
    k = mp.mpf(K) * mp.matrix([mp.mpf(3) / 5, 0, mp.mpf(4) / 5])
    transfer = mp.mpf(P) * mp.matrix([mp.mpf(2) / 3, mp.mpf(1) / 3, mp.mpf(2) / 3])
    ell = -k + transfer
    nu = mp.sqrt(1000**2 + (k.T * k)[0] / (mp.mpf(25) / 16) ** 2)
    D = mp.diag([1, -2, 1])
    G = D.copy()
    G[0, 1] = G[1, 0] = mp.mpf(1) / 3
    MD = contracted_spatial_matrix(D)
    MG = contracted_spatial_matrix(G)
    center = mp.mpf(center)
    radius = mp.mpf("1e-5")
    basek = real_modes(center, k)
    basel = real_modes(center, ell)
    detector = [
        mp.conj((basek[a][0].T * MD * basel[b][0])[0])
        for a in range(3)
        for b in range(3)
    ]
    gs = []
    amplitudes = []
    phases = [mp.exp(2j * mp.pi * n / samples) for n in range(samples)]
    for phase in phases:
        u = center + radius * phase
        first = real_modes(u, k)
        second = real_modes(u, ell)
        gs.append(
            [1 / (first[a][1] + second[b][1]) for a in range(3) for b in range(3)]
        )
        amplitudes.append(
            [(first[a][0].T * MG * second[b][0])[0] for a in range(3) for b in range(3)]
        )
    result = [mp.mpc(0) for _ in range(5)]
    for pair in range(9):
        g = [
            mp.fsum(gs[n][pair] * phases[n] ** (-j) for n in range(samples))
            / samples
            / radius**j
            for j in range(5)
        ]
        amplitude = [
            mp.fsum(amplitudes[n][pair] * phases[n] ** (-j) for n in range(samples))
            / samples
            / radius**j
            for j in range(5)
        ]
        for source_jet in range(5):
            current = [
                amplitude[n - source_jet] / mp.factorial(source_jet)
                if n >= source_jet
                else mp.mpc(0)
                for n in range(5)
            ]
            for order in range(5):
                result[source_jet] += (
                    (-1j) * (1j) ** order * detector[pair] * g[0] * current[0]
                )
                current = [
                    (n + 1) * mp.fsum(g[j] * current[n + 1 - j] for j in range(n + 2))
                    for n in range(len(current) - 1)
                ]
    return mp.matrix(result), nu, mp.norm(D) * mp.norm(G)


@pytest.mark.parametrize(
    "K,P", (("0", "0"), ("500", "1234"), ("1000", "1e7"), ("1e6", "0.001"))
)
def test_actual_all_five_endpoint_time_jet_rows_at_all_real_momenta(K, P):
    with mp.workdps(90):
        row, nu, norms = actual_endpoint_jet_row(K, P)
        assert mp.norm(row) < mp.mpf("4e27") * nu * norms
        assert all(abs(value) > 0 for value in row)
        # Independently changing the Cauchy quadrature resolves the same full
        # time-jet coefficients, including noncommuting source/detector tensors.
        refined, _, _ = actual_endpoint_jet_row(K, P, samples=20)
        assert mp.norm(row - refined) < mp.mpf("1e-45") * max(
            mp.mpf(1), mp.norm(refined)
        )


@pytest.mark.parametrize(
    "K,P", (("1000", "1"), ("1000", "1e6"), ("1e6", "1"), ("1e6", "1e6"))
)
def test_actual_full_leading_endpoint_outside_small_spatial_Cauchy_ball(K, P):
    with mp.workdps(105):
        function = lambda parameter: complete_leading_endpoint(K, parameter)[0]
        coefficients = mp.taylor(function, mp.mpf(0), 4)
        _, nu, norms = complete_leading_endpoint(K, 0)
        parameter = mp.mpf(P) / (mp.mpf("1e-6") * nu)
        value, _, _ = complete_leading_endpoint(K, parameter)
        polynomial = mp.polyval(list(reversed(coefficients)), parameter)
        majorant = (
            mp.mpf("2e27") * nu * (2 + sum(parameter**j for j in range(5))) * norms
        )
        assert 0 < abs(value - polynomial) < majorant


def near_upper(P):
    A = mp.mpf(25) / 16
    m = mp.mpf(1000)
    d = mp.mpf("1e-6")
    if 2 * P / d <= m:
        return mp.mpf(0)
    return A * mp.sqrt((2 * P / d) ** 2 - m * m)


def radial_power(order, lower, upper):
    if upper <= lower:
        return mp.mpf(0)
    A = mp.mpf(25) / 16
    m = mp.mpf(1000)
    start = mp.asinh(lower / (A * m))
    stop = mp.asinh(upper / (A * m))
    return (
        A**3
        * m ** (4 - order)
        / (2 * mp.pi**2)
        * mp.quad(lambda x: mp.sinh(x) ** 2 * mp.cosh(x) ** (2 - order), [start, stop])
    )


def numeric_constant(name):
    return mp.mpf(str(s.N(near.constants()[name], 110)))


@pytest.mark.parametrize("P", ("0.001", "1", "1e6", "1e12"))
@pytest.mark.parametrize("order", range(5))
def test_complete_near_radial_bounds_for_each_Taylor_order(P, order):
    with mp.workdps(80):
        A = mp.mpf(25) / 16
        m = mp.mpf(1000)
        d = mp.mpf("1e-6")
        P = mp.mpf(P)
        L = 2 * A * P / d
        actual = radial_power(order, m, near_upper(P))
        bounds = [
            L**4 / (4 * mp.pi**2),
            L**3 / (6 * mp.pi**2),
            A * L**2 / (4 * mp.pi**2),
            A * A * L / (2 * mp.pi**2),
            A**3 * mp.log(L / m) / (2 * mp.pi**2),
        ]
        assert 0 < actual < bounds[order]


@pytest.mark.parametrize("P", ("0", "1e-6", "0.0005", "0.001", "1", "1e6", "1e12"))
def test_combined_raw_and_all_Taylor_near_integrals(P):
    with mp.workdps(85):
        P = mp.mpf(P)
        upper = near_upper(P)
        d = mp.mpf("1e-6")
        E = mp.mpf("2e27")
        actual = 2 * E * radial_power(0, mp.mpf(1000), upper)
        actual += E * sum(
            (P / d) ** j * radial_power(j, mp.mpf(1000), upper) for j in range(5)
        )
        bound = (
            numeric_constant("complete_near_P4_upper") * P**4
            + numeric_constant("complete_near_P5_upper") * P**5
        )
        assert 0 <= actual <= bound
        assert actual < mp.mpf("2e54") * (1 + P * P) ** 3


@pytest.mark.parametrize("module", (real, near, tail, limit))
def test_exact_core_packets(module):
    d = module.data()
    assert all(s.cancel(value) == 0 for value in d["checks"].values())
    assert all(d["gates"].values())


def removed_fraction(k, P, K):
    if k > K:
        return mp.mpf(1)
    if P == 0:
        return mp.mpf(0)
    if k == 0:
        return mp.mpf(1) if P > K else mp.mpf(0)
    c = (k * k + P * P - K * K) / (2 * k * P)
    return (min(mp.mpf(1), max(mp.mpf(-1), c)) + 1) / 2


def removed_radial(order, lower, upper, P, K):
    if upper <= lower:
        return mp.mpf(0)
    A = mp.mpf(25) / 16
    m = mp.mpf(1000)
    breaks = (
        [lower]
        + [x for x in sorted({abs(K - P), K, K + P}) if lower < x < upper]
        + [upper]
    )
    hyperbolic = [mp.asinh(x / (A * m)) for x in breaks]
    return (
        A**3
        * m ** (4 - order)
        / (2 * mp.pi**2)
        * mp.quad(
            lambda x: (
                mp.sinh(x) ** 2
                * mp.cosh(x) ** (2 - order)
                * removed_fraction(A * m * mp.sinh(x), P, K)
            ),
            hyperbolic,
        )
    )


@pytest.mark.parametrize("K", ("1000", "1e6"))
@pytest.mark.parametrize("fraction", ("0.25", "1", "2"))
def test_actual_removed_near_two_leg_union_and_same_six_derivative_tail(K, fraction):
    with mp.workdps(85):
        K = mp.mpf(K)
        d = mp.mpf("1e-6")
        P = d * K * mp.mpf(fraction)
        E = mp.mpf("2e27")
        upper = near_upper(P)
        m = mp.mpf(1000)
        actual = 2 * E * removed_radial(0, m, upper, P, K)
        actual += E * sum(
            (P / d) ** j * removed_radial(j, m, upper, P, K) for j in range(5)
        )
        C4 = numeric_constant("complete_near_P4_upper")
        C5 = numeric_constant("complete_near_P5_upper")
        B = mp.mpf(3125001)
        assert 0 <= actual <= B * P / K * (C4 * P**4 + C5 * P**5)
        assert actual <= B * (C4 + C5) / K * (1 + P * P) ** 3
        if fraction in ("1", "2"):
            assert actual > 0


@pytest.mark.parametrize("K", ("1000", "2000", "1e6"))
@pytest.mark.parametrize("P", ("0", "1", "1000", "1e6"))
def test_actual_unexpanded_low_band_removed_union(K, P):
    with mp.workdps(80):
        K = mp.mpf(K)
        P = mp.mpf(P)
        m = mp.mpf(1000)
        low = mp.mpf(
            str(s.N(real.constants()["unexpanded_low_band_integral_upper"], 90))
        )
        actual = mp.mpf("4e27") * removed_radial(0, mp.mpf(0), m, P, K)
        assert 0 <= actual <= low
        assert actual <= low * (m + P) / K
        assert actual <= low * (m + 1) / K * (1 + P * P) ** 3


@pytest.mark.parametrize("P", ("0", "0.001", "0.5", "1", "2", "1e12"))
def test_same_spatial_weight_controls_every_tail_power(P):
    with mp.workdps(70):
        P = mp.mpf(P)
        for order in (0, 1, 4, 5, 6):
            assert P**order <= (1 + P * P) ** 3


@pytest.mark.parametrize("fraction", ("0", "0.1", "0.5", "1", "2", "2.1"))
def test_polynomial_integrand_does_not_make_moving_intersection_polynomial(fraction):
    with mp.workdps(80):
        K = mp.mpf(1000)
        P = K * mp.mpf(fraction)
        breaks = (
            [mp.mpf(0)] + [x for x in sorted({abs(K - P), K + P}) if 0 < x < K] + [K]
        )
        actual = (
            4
            * mp.pi
            * mp.quad(lambda k: k * k * (1 - removed_fraction(k, P, K)), breaks)
        )
        exact = mp.pi * (4 * K + P) * (2 * K - P) ** 2 / 12 if P <= 2 * K else mp.mpf(0)
        assert abs(actual - exact) < mp.mpf("1e-65") * K**3


def test_opposite_sharp_band_derivatives_and_divergent_surface_artifact():
    K, d = s.symbols("K d", positive=True)
    volume = s.pi * (4 * K + d) * (2 * K - d) ** 2 / 12
    right = s.diff(volume, d).subs(d, 0)
    left = -right
    assert s.factor(right - left) == -2 * s.pi * K * K
    assert s.limit((volume - 4 * s.pi * K**3 / 3) / K**2, K, s.oo) == -s.pi * d
    # This constant-integrand geometric counterexample is not an assertion
    # that the corresponding coefficient survives in the actual full current.


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
    assert len(audit.residuals()) == 29 and audit.scalar_entry_count() == 29
    assert len(audit.gates()) == 38 and all(
        value is True for value in audit.gates().values()
    )
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 125


def test_original_primitive_and_matching_statuses_retained():
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert "NOT_LOCAL_MATCHING_FULL_INVERSE_OR_V_G_B" in audit.ITEM["status"]
    audit.validate_scope(audit.frontier(), audit.matching())
