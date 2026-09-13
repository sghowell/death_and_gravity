"""Independent full scalar evolution, literal variations, geometry, state tails and contacts."""

from fractions import Fraction as F
from functools import cache

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_heavy_curved_state import state
from p8_vacuum_affine_heavy_homogeneous_response import (
    audit,
    estimates,
    evolution,
    renormalization,
)
from scipy.integrate import solve_ivp
from scipy.linalg import expm, expm_frechet


@pytest.mark.parametrize("name,value", tuple(audit.residuals().items()))
def test_every_exact_identity_and_full_matrix(name, value):
    assert (
        all(entry == 0 for entry in value)
        if isinstance(value, s.MatrixBase)
        else value == 0
    ), name


@pytest.mark.parametrize("name,value", tuple(audit.gates().items()))
def test_every_written_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_or_claim_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_entire_covariance_cross_terms_and_pure_normalization():
    data = evolution.data()
    cov = data["exact_pure_graph_covariance"]
    assert cov[0, 1] != 0 and cov[0, 1] == cov[1, 0]
    assert s.factor(cov.det() - s.Rational(1, 4)) == 0
    assert data["physical_cross_covariance_not_declared_even"] != 0
    assert data["complete_covariance_response_with_second_contact"] != 0


def test_full_error_forcing_retains_first_error_square():
    row = evolution.data()["full_second_error_forcing"]
    error1 = next(x for x in row.free_symbols if str(x) == "error1")
    squeeze0 = next(x for x in row.free_symbols if str(x) == "squeeze0")
    assert s.expand(row).coeff(error1, 2) == -2 * squeeze0


def test_entire_mixed_jet_and_reference_majorant_manifest():
    data = estimates.data()
    rows = data["all_finite_Riccati_coefficient_derivative_majorants"]
    assert len(rows) == 195
    assert set(rows) == {
        (n, j, a) for n in range(1, 11) for j in range(12 - n) for a in range(3)
    }
    assert all(value > 0 for value in rows.values())
    assert len(data["complete_mixed_frequency_jets"]) == 39
    assert all(x > 0 for x in data["whole_finite_residual_parameter_bounds"])
    assert data["whole_initial_state_comparison_and_W6_coefficient_tail"] == 10**90
    assert estimates.INITIAL / estimates.K + estimates.C[0] > 0
    assert F(99, 100) * estimates.n > 10**196 > estimates.K**2


def test_every_parameter_bound_uses_whole_initial_and_residual_errors():
    assert estimates.raw_E0 == 1000 * (
        F(estimates.INITIAL, estimates.K) + estimates.C[0]
    )
    assert estimates.raw_E1 > 3 * estimates.E0 * 12
    assert estimates.raw_E2 > 3 * estimates.E1 * 24
    assert all(
        value < limit
        for value, limit in zip(
            (estimates.raw_E0, estimates.raw_E1, estimates.raw_E2),
            (10**40, 10**42, 10**44),
        )
    )
    assert all(value < 10**50 for value in estimates.readout)
    assert all(2 * value < 10**51 for value in estimates.TAIL)


def test_full_covariance_frechet_majorants_on_whole_complex_graph_ball():
    radius = F(1, 10)
    den = 1 - radius**2
    g = (
        1 / den,
        2 * radius / den**2,
        2 * (2 * radius) ** 2 / den**3 + 2 / den**2,
        6 * (2 * radius) ** 3 / den**4 + 6 * (2 * radius) * 2 / den**3,
    )
    assert all(value < bound for value, bound in zip(g, (2, 1, 10, 16)))
    assert ((1 + radius) ** 2 / 2 + radius) / den < 2
    assert 4 * 2 + 2 * 1 < 16
    assert 2 * 2 + 2 * 4 * 1 + 2 * 10 < 128
    assert 3 * 2 * 1 + 3 * 4 * 10 + 2 * 16 < 512


def test_mass_independent_marker_circle_and_full_even_tail():
    assert estimates.CONTOUR == 10**6
    assert all(x < bound for x, bound in zip(estimates.contour, (10, 100, 1000)))
    assert estimates.R[0] < F(1, 100)
    assert estimates.TAIL[2] == estimates.readout[2] + 2 * 1000 * 10**36
    assert estimates.TAIL[2] > 0


def test_complete_radial_integral_and_same_original_frequency_tail():
    # k=4sqrt(n)u and then u=tan(theta) give integral sin²(theta)cos(theta)=1/3.
    theta = s.Symbol("theta", real=True)
    integral = s.integrate(s.sin(theta) ** 2 * s.cos(theta), (theta, 0, s.pi / 2))
    assert integral == s.Rational(1, 3)
    coefficient = (
        64 * integral / (2 * s.pi**2) * s.Rational(100, 99) ** s.Rational(5, 2)
    )
    assert float(coefficient) < 2
    assert F(32, 9) * F(100, 99) ** 2 < 4
    assert all(F(200, 99) * value < 10**51 for value in estimates.TAIL)
    assert estimates.data()[
        "same_frequency_cutoff_tail_coefficient_over_K_squared"
    ] < s.Rational(1, 10**748)


def test_full_nonzero_covariant_boundary_and_unchanged_finite_action():
    data = renormalization.data()
    assert data["complete_general_dimensional_boundary_difference"] != 0
    assert data["all_scalar_radial_poles"] == (1, -s.Rational(1, 3), 1)
    ell = next(
        x
        for x in data[
            "entire_finite_scalar_local_action_per_physical_volume"
        ].free_symbols
        if str(x) == "ell"
    )
    assert data["all_scalar_radial_finite_coefficients"] == (
        s.Rational(3, 2) - ell,
        (ell - 1) / 3,
        -ell,
    )
    assert data["complete_heat_current_two_parameter_derivative_bound"] > 0
    assert data["entire_unchanged_fixed_profile_derivative_bound"] > 0
    assert data[
        "entire_normalized_homogeneous_current_and_response_bound"
    ] < s.Rational(1, 10**350)
    assert data["full_Euler_current_jet_derivative_term_count"] == 1088743692


def test_fixed_profile_and_its_volume_contacts_are_not_zero_alone():
    epsilon = s.Symbol("epsilon", real=True)
    pressure = s.Symbol("fixed_pressure", real=True)
    td, tg = s.symbols("trace_D trace_G", real=True)
    current = -s.exp(epsilon * tg / 2) * pressure * td / 2
    assert s.diff(current, epsilon).subs(epsilon, 0) == -pressure * td * tg / 4
    assert s.diff(current, epsilon, 2).subs(epsilon, 0) == -pressure * td * tg**2 / 8
    projected = s.Symbol("projected_pressure", real=True)
    assert pressure - projected != 0


@pytest.mark.parametrize("j", range(13))
def test_admitted_history_jets_are_each_checked(j):
    assert audit.require_jet_bound(j, s.Rational(1, 100)) == (j, s.Rational(1, 100))


def test_complete_scope_and_matching_frontier():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 100
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 297
    assert audit.require_domain(0, 0) == (0, 0)
    assert audit.require_graph("source_time", 13) == ("source_time", 13)
    assert "nonzero-transfer" in audit.observable()["remaining"]
    assert audit.require_cutoff_squared(4 * state.MASS2) == 4 * state.MASS2


# Independent actual-mass W6 versus order-ten scalar graph coefficient tail.
def add(a, b):
    return [x + y for x, y in zip(a, b)]


def scale(a, c):
    return [c * x for x in a]


def mul(a, b):
    length = min(len(a), len(b))
    return [sum(a[k] * b[j - k] for k in range(j + 1)) for j in range(length)]


def inverse(a):
    out = [1 / a[0]]
    for j in range(1, len(a)):
        out.append(-sum(a[k] * out[j - k] for k in range(1, j + 1)) / a[0])
    return out


def root(a):
    out = [mp.sqrt(a[0])]
    for j in range(1, len(a)):
        out.append(
            (a[j] - sum(out[k] * out[j - k] for k in range(1, j))) / (2 * out[0])
        )
    return out


def derivative(a):
    return [(j + 1) * a[j + 1] for j in range(len(a) - 1)]


def constant(c, length):
    return [mp.mpf(c)] + [mp.mpf(0)] * (length - 1)


def graph_tail(precision, ratio):
    with mp.workdps(precision):
        length = 27
        time = mp.mpf("-0.5")
        n = mp.mpf(10) ** 200 / 512 + 2
        coord = [time, mp.mpf(1)] + [mp.mpf(0)] * (length - 2)
        v = add(constant(1, length), mul(coord, coord))
        H = scale(mul(coord, inverse(v)), 4)
        inva = mul(inverse(v), inverse(v))
        square = add(constant(n, length), scale(mul(inva, inva), n * ratio * ratio))
        omega = root(square)
        lam = mul(derivative(omega), inverse(omega))
        squeeze = scale(add(scale(H, 3), lam), mp.mpf("0.5"))
        phase = scale(inverse(omega), 1 / (2 * mp.j))
        terms = [None, scale(mul(phase, squeeze), -1)]
        for order in range(2, 11):
            value = derivative(terms[order - 1])
            for left in range(1, order - 1):
                value = add(
                    value, mul(squeeze, mul(terms[left], terms[order - 1 - left]))
                )
            terms.append(mul(phase, value))
        reference = sum(term[0] for term in terms[1:])
        U = add(scale(derivative(H), mp.mpf("1.5")), scale(mul(H, H), mp.mpf("2.25")))
        W = omega
        for _ in range(6):
            rate = mul(derivative(W), inverse(W))
            W = root(
                add(
                    square,
                    add(
                        scale(U, -1),
                        add(
                            scale(derivative(rate), -mp.mpf("0.5")),
                            scale(mul(rate, rate), mp.mpf("0.25")),
                        ),
                    ),
                )
            )
        d = mp.mpf("1.5") * H[0] + W[1] / (2 * W[0])
        full = (omega[0] - W[0] + mp.j * d) / (omega[0] + W[0] - mp.j * d)
        Omega = mp.sqrt(n * (1 + mp.mpf(ratio) ** 2 / 16))
        scaled = (full - reference) * Omega**11
        assert abs(scaled) < mp.mpf(10) ** 88
        assert abs(full) < mp.mpf("0.01")
        return +scaled


@pytest.mark.parametrize("ratio", (0, 1, 4))
def test_full_W6_to_order_ten_initial_graph_at_actual_mass(ratio):
    first = graph_tail(1300, ratio)
    second = graph_tail(1500, ratio)
    with mp.workdps(180):
        assert abs(first - second) < mp.mpf("1e-150") * max(1, abs(second))
        assert abs(second) > mp.mpf("1e-30")


# Independent literal Christoffels and Cartesian angular/radial moment checks.
def average(value, coordinates):
    result = s.S.Zero
    d = len(coordinates)
    for powers, coefficient in s.Poly(s.expand(value), *coordinates).terms():
        if any(power % 2 for power in powers):
            continue
        half = sum(powers) // 2
        result += (
            coefficient
            * s.prod(s.factorial2(power - 1) for power in powers)
            / s.prod(d + 2 * j for j in range(half))
        )
    return s.factor(result)


def fixture(d, seed):
    K = s.Matrix(
        d,
        d,
        lambda i, j: (
            s.Rational((i + j + seed) % 5 - 2, 3) if i != j else s.Rational(i + seed, 2)
        ),
    )
    A = s.Matrix(
        d,
        d,
        lambda i, j: (
            s.Rational((i + j + 2 * seed + 1) % 7 - 3, 5)
            if i != j
            else s.Rational(2 - i - seed, 3)
        ),
    )
    return K, A


def curvature_check(K, A):
    d = len(K[:, 0])
    count = d + 1
    sign = [-1] + [1] * d
    second = A + 2 * K * K

    def gamma(r, mu, nu):
        if r == 0 and mu and nu:
            return K[mu - 1, nu - 1]
        if r and bool(mu) != bool(nu):
            return K[r - 1, (mu or nu) - 1]
        return s.S.Zero

    def gamma_t(r, mu, nu):
        if r == 0 and mu and nu:
            return second[mu - 1, nu - 1]
        if r and bool(mu) != bool(nu):
            return A[r - 1, (mu or nu) - 1]
        return s.S.Zero

    R = {}
    for r in range(count):
        for sigma in range(count):
            for mu in range(count):
                for nu in range(count):
                    R[r, sigma, mu, nu] = (
                        (gamma_t(r, nu, sigma) if mu == 0 else 0)
                        - (gamma_t(r, mu, sigma) if nu == 0 else 0)
                        + sum(
                            gamma(r, mu, k) * gamma(k, nu, sigma)
                            - gamma(r, nu, k) * gamma(k, mu, sigma)
                            for k in range(count)
                        )
                    )
    ric = s.Matrix(count, count, lambda i, j: sum(R[r, i, r, j] for r in range(count)))
    scalar = sum(sign[i] * ric[i, i] for i in range(count))
    ric2 = sum(
        sign[i] * sign[j] * ric[i, j] ** 2 for i in range(count) for j in range(count)
    )
    riem2 = sum(
        sign[r] * sign[sigma] * sign[mu] * sign[nu] * R[r, sigma, mu, nu] ** 2
        for r in range(count)
        for sigma in range(count)
        for mu in range(count)
        for nu in range(count)
    )
    theta = s.trace(K)
    X = s.trace(A)
    k2 = s.trace(K * K)
    expectedR = 2 * X + theta**2 + k2
    expectedRic = (X + k2) ** 2 + s.trace((A + theta * K) ** 2)
    expectedRiem = 4 * s.trace((A + K * K) ** 2) + 2 * (k2 * k2 - s.trace(K**4))
    assert s.factor(scalar - expectedR) == 0
    assert s.factor(ric2 - expectedRic) == 0
    assert s.factor(riem2 - expectedRiem) == 0
    return 3


def angular_check(K, A):
    from p8_vacuum_affine_heavy_homogeneous_response import renormalization as action

    d = K.rows
    n = s.Matrix(s.symbols("n0:" + str(d), real=True))
    z = s.Symbol("z", real=True)
    theta = s.trace(K)
    X = s.trace(A)
    kn = (n.T * K * n)[0]
    M = -A + 2 * K * K + theta * K
    squeeze = (theta - z * kn) / 2
    tt = (X + z * (n.T * M * n)[0] - 3 * z * z * kn * kn) / 2
    poly = s.Poly(s.expand(tt * tt + squeeze**4), z)
    actual = s.factor(
        sum(
            average(value, tuple(n))
            * s.rf(s.Rational(d, 2), j)
            / s.rf(s.Rational(3, 2), j)
            for (j,), value in poly.terms()
        )
    )
    subs = {
        action.theta: theta,
        action.X: X,
        action.k2: s.trace(K * K),
        action.k3: s.trace(K**3),
        action.k4: s.trace(K**4),
        action.A2: s.trace(A * A),
        action.AK: s.trace(A * K),
        action.AK2: s.trace(A * K * K),
    }
    assert s.factor(actual - action.whole.subs(subs, simultaneous=True)) == 0
    return 1


@pytest.mark.parametrize("dimension", (2, 3, 4, 5))
@pytest.mark.parametrize("seed", (1, 2))
def test_independent_literal_Christoffels_and_Cartesian_ad4(dimension, seed):
    K, A = fixture(dimension, seed)
    assert K * A != A * K
    assert curvature_check(K, A) == 3
    assert angular_check(K, A) == 1


# Independent prepared scalar covariance IVP, full noncommuting vertices and contacts.
J = np.array([[0.0, 1.0], [-1.0, 0.0]])
G0 = np.array([[0.7, 0.2, -0.1], [0.2, -0.3, 0.15], [-0.1, 0.15, 0.4]])
G1 = np.array([[-0.2, 0.12, 0.3], [0.12, 0.4, -0.1], [0.3, -0.1, -0.1]])
D = np.array([[0.6, -0.2, 0.1], [-0.2, 0.1, 0.25], [0.1, 0.25, -0.4]])


def source(time):
    x = (time + 0.05) / 0.3
    if abs(x) >= 1:
        return np.zeros((3, 3))
    b = np.exp(1 - 1 / (1 - x * x))
    return b * (G0 + x * G1)


def fields(time, k, n, epsilon):
    a = (1 + time * time) ** 2
    G = source(time)
    Q = epsilon * G
    v = a**3 * np.exp(np.trace(Q) / 2)
    inverse_metric = expm(-Q)
    omega2 = n + k @ inverse_metric @ k / a**2
    first = expm_frechet(-Q, -D, compute_expm=False)
    qD = k @ first @ k / a**2
    tauD = np.trace(D)
    M = np.diag([v * omega2, 1 / v])
    MD = np.diag([v * (tauD * omega2 / 2 + qD), -tauD / (2 * v)])
    return M, MD


def full_jets(time, k, n):
    a = (1 + time * time) ** 2
    v = a**3
    G = source(time)
    tg = np.trace(G)
    td = np.trace(D)
    omega2 = n + k @ k / a**2
    qg = -k @ G @ k / a**2
    qgg = k @ G @ G @ k / a**2
    qd = -k @ D @ k / a**2
    qdg = k @ (G @ D + D @ G) @ k / (2 * a**2)
    qdgg = -k @ (G @ G @ D + G @ D @ G + D @ G @ G) @ k / (3 * a**2)
    mm = (
        np.diag([v * omega2, 1 / v]),
        np.diag([v * (tg * omega2 / 2 + qg), -tg / (2 * v)]),
        np.diag([v * (tg * tg * omega2 / 4 + tg * qg + qgg), tg * tg / (4 * v)]),
    )
    base = td * omega2 / 2 + qd
    first = td * qg / 2 + qdg
    second = td * qgg / 2 + qdgg
    md = (
        np.diag([v * base, -td / (2 * v)]),
        np.diag([v * (tg * base / 2 + first), td * tg / (4 * v)]),
        np.diag(
            [v * (tg * tg * base / 4 + tg * first + second), -td * tg * tg / (8 * v)]
        ),
    )
    return mm, md


def full_ivp_diagnostics():
    errors = []
    assert np.linalg.norm(G0 @ G1 - G1 @ G0) > 0.1
    assert np.linalg.norm(D @ G0 - G0 @ D) > 0.1
    for k in (
        np.array([0.0, 0.0, 0.0]),
        np.array([0.7, 1.3, -0.4]),
        np.array([2.0, -1.0, 3.0]),
    ):
        n = 9.0
        start = -0.5
        end = 0.1
        a = (1 + start * start) ** 2
        omega = np.sqrt(n + k @ k / a**2)
        scale = np.diag([1 / np.sqrt(a**3 * omega), np.sqrt(a**3 * omega)])
        inverse = np.linalg.inv(scale)
        r = 0.06 + 0.08j
        den = 1 - abs(r) ** 2
        Cinit = np.array(
            [
                [(1 + 2 * r.real + abs(r) ** 2) / (2 * den), -r.imag / den],
                [-r.imag / den, (1 - 2 * r.real + abs(r) ** 2) / (2 * den)],
            ]
        )

        def transformed(time, epsilon, k=k, n=n, inverse=inverse, scale=scale):
            M, MD = fields(time, k, n, epsilon)
            return inverse @ J @ M @ scale, scale.T @ MD @ scale

        def exact_current(
            epsilon, transformed=transformed, start=start, end=end, Cinit=Cinit
        ):
            def fun(time, value):
                A, _ = transformed(time, epsilon)
                C = value.reshape(2, 2)
                return (A @ C + C @ A.T).ravel()

            result = solve_ivp(
                fun,
                (start, end),
                Cinit.ravel(),
                method="DOP853",
                rtol=5e-14,
                atol=2e-15,
            )
            assert result.success
            _, MD = transformed(end, epsilon)
            return -0.5 * np.trace(MD @ result.y[:, -1].reshape(2, 2))

        def tangent(time, value, k=k, n=n, inverse=inverse, scale=scale):
            mm, _ = full_jets(time, k, n)
            aa = [inverse @ J @ M @ scale for M in mm]
            cc = value.reshape(3, 2, 2)
            zero = aa[0] @ cc[0] + cc[0] @ aa[0].T
            one = aa[0] @ cc[1] + cc[1] @ aa[0].T + aa[1] @ cc[0] + cc[0] @ aa[1].T
            two = (
                aa[0] @ cc[2]
                + cc[2] @ aa[0].T
                + 2 * (aa[1] @ cc[1] + cc[1] @ aa[1].T)
                + aa[2] @ cc[0]
                + cc[0] @ aa[2].T
            )
            return np.stack([zero, one, two]).ravel()

        initial = np.stack([Cinit, np.zeros((2, 2)), np.zeros((2, 2))])
        result = solve_ivp(
            tangent,
            (start, end),
            initial.ravel(),
            method="DOP853",
            rtol=5e-14,
            atol=2e-15,
        )
        assert result.success
        cc = result.y[:, -1].reshape(3, 2, 2)
        _, md = full_jets(end, k, n)
        md = [scale.T @ M @ scale for M in md]
        answer1 = -0.5 * np.trace(md[0] @ cc[1] + md[1] @ cc[0])
        answer2 = -0.5 * np.trace(md[0] @ cc[2] + 2 * md[1] @ cc[1] + md[2] @ cc[0])
        contact1 = -0.5 * np.trace(md[1] @ cc[0])
        contact2 = -0.5 * np.trace(2 * md[1] @ cc[1] + md[2] @ cc[0])
        h = 0.004
        j0 = exact_current(0)
        jp = exact_current(h)
        jm = exact_current(-h)
        jph = exact_current(h / 2)
        jmh = exact_current(-h / 2)
        fd1 = (4 * (jph - jmh) / h - (jp - jm) / (2 * h)) / 3
        fd2 = (4 * (jph - 2 * j0 + jmh) / (h / 2) ** 2 - (jp - 2 * j0 + jm) / h**2) / 3
        assert abs(fd1 - answer1) < 3e-8
        assert abs(fd2 - answer2) < 3e-6
        assert abs(contact1) > 1e-4 and abs(contact2) > 1e-4
        errors.append((abs(fd1 - answer1), abs(fd2 - answer2)))
        print(
            "PREPARED_WHOLE_SCALAR_CURRENT_IVP",
            k.tolist(),
            "response",
            answer1,
            answer2,
            "contacts",
            contact1,
            contact2,
            "errors",
            errors[-1],
            flush=True,
        )
    print("FULL_PREPARED_CURRENT_TWO_DERIVATIVES_AND_CONTACTS_C0", flush=True)
    return errors


def test_complete_prepared_current_first_second_derivatives_and_nonzero_contacts():
    assert len(full_ivp_diagnostics()) == 3


# Independent literal scalar Euler variation versus the complete adiabatic covariance.


@cache
def variational_diagnostics():
    w = (s.Symbol("w0", positive=True),) + s.symbols("w1:5", real=True)
    eta = s.symbols("eta0:5", real=True)

    def dt(value):
        return s.expand(
            sum(
                s.diff(value, row[j]) * row[j + 1] for row in (w, eta) for j in range(4)
            )
        )

    lam = w[1] / w[0]
    squeeze = (lam - eta[1]) / 2
    temp = s.factor(dt(squeeze) - lam * squeeze)
    lagrangians = {
        0: -w[0] / 2,
        2: squeeze**2 / (4 * w[0]),
        4: (temp**2 + squeeze**4) / (16 * w[0] ** 3),
    }
    r = [s.S.Zero, s.I * squeeze / (2 * w[0])]
    for order in range(2, 5):
        value = dt(r[order - 1]) + squeeze * sum(
            r[j] * r[order - 1 - j] for j in range(1, order - 1)
        )
        r.append(s.cancel(value / (2 * s.I * w[0])))
    sharp = [s.conjugate(value) for value in r]

    def product(left, right):
        return [
            s.expand(sum(left[j] * right[n - j] for j in range(n + 1)))
            for n in range(5)
        ]

    quadratic = product(r, sharp)
    denominator = [s.S.One] + [-quadratic[j] for j in range(1, 5)]
    inverse = [s.S.One]
    for j in range(1, 5):
        inverse.append(
            s.expand(-sum(denominator[k] * inverse[j - k] for k in range(1, j + 1)))
        )
    one = [s.S.One] + [s.S.Zero] * 4
    plus = [one[j] + r[j] for j in range(5)]
    plussharp = [one[j] + sharp[j] for j in range(5)]
    minus = [one[j] - r[j] for j in range(5)]
    minussharp = [one[j] - sharp[j] for j in range(5)]
    QQ = [value / 2 for value in product(product(plus, plussharp), inverse)]
    PP = [value / 2 for value in product(product(minus, minussharp), inverse)]
    result = {}
    for order, L in lagrangians.items():

        def euler(row, L=L):
            terms = []
            for j in range(3):
                value = s.diff(L, row[j])
                for _ in range(j):
                    value = -dt(value)
                terms.append(value)
            return s.factor(sum(terms))

        result[order] = (
            s.factor(euler(eta) - w[0] * (QQ[order] - PP[order]) / 2),
            s.factor(euler(w) + QQ[order]),
        )
        assert result[order] == (0, 0)
    assert any(QQ[j] != 0 for j in (2, 4))
    print(
        "LITERAL_FULL_SCALAR_ACTION_COVARIANCE_VARIATIONS_C0",
        len(result) * 2,
        flush=True,
    )
    return result


@pytest.mark.parametrize("order", (0, 2, 4))
def test_independent_literal_scalar_Euler_covariance_variation(order):
    assert variational_diagnostics()[order] == (0, 0)
