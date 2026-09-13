"""Independent full scalar charts, mode covariance, profiles, graph and cutoff controls."""

from fractions import Fraction as F
from math import comb, factorial

import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_heavy_adm_clock_response import (
    audit,
    estimates,
    geometry,
    profile,
)
from p8_vacuum_affine_heavy_curved_state import state
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm


@pytest.mark.parametrize("name,value", tuple(audit.residuals().items()))
def test_every_full_exact_identity(name, value):
    assert (
        all(entry == 0 for entry in value)
        if isinstance(value, s.MatrixBase)
        else value == 0
    ), name


@pytest.mark.parametrize("name,value", tuple(audit.gates().items()))
def test_every_written_full_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_or_claim_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_full_counts_and_unchanged_frontier():
    assert (
        len(audit.residuals()),
        audit.scalar_entry_count(),
        len(audit.gates()),
        len(audit.controls()),
        audit.rejected_inputs(),
    ) == (65, 200, 50, 9, 172)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 102
    assert audit.matching()[:-1] == audit.previous.matching()


def metric_value(n, beta, Q, scale):
    h = scale**2 * expm(Q)
    out = np.zeros((4, 4))
    out[0, 0] = (1 + n) ** 2 - beta @ h @ beta
    out[0, 1:] = -h @ beta
    out[1:, 0] = -h @ beta
    out[1:, 1:] = -h
    return out


def hamiltonian_value(n, beta, Q, scale):
    out = np.zeros((5, 5))
    N = 1 + n
    trace = np.trace(Q)
    out[0, 0] = N * np.exp(-trace / 2)
    out[0, 1:4] = scale * beta
    out[1:4, 0] = scale * beta
    out[1:4, 1:4] = N * np.exp(trace / 2) * expm(-Q)
    out[4, 4] = N * np.exp(trace / 2)
    return out


def to_symbolic(direction):
    n, beta, Q = direction
    return s.Rational(str(n)), s.Matrix(beta), s.Matrix(Q)


@pytest.mark.parametrize("seed", range(12))
def test_literal_entire_matrix_exponential_metric_and_scalar_Hamiltonian_contacts(seed):
    rng = np.random.default_rng(7301 + seed)
    QD = rng.normal(size=(3, 3))
    QD = (QD + QD.T) / 5
    QG = rng.normal(size=(3, 3))
    QG = (QG + QG.T) / 5
    assert np.linalg.norm(QD @ QG - QG @ QD) > 1e-3
    D = (float(rng.normal() / 3), rng.normal(size=3) / 3, QD)
    G = (float(rng.normal() / 3), rng.normal(size=3) / 3, QG)
    t = (seed - 5.5) / 12
    scale = (1 + t * t) ** 2
    eps = 1e-4
    for literal, first, second in (
        (metric_value, geometry.first, geometry.second),
        (hamiltonian_value, geometry.hamiltonian_first, geometry.hamiltonian_second),
    ):

        def value(e, f, literal=literal):
            return literal(
                e * D[0] + f * G[0], e * D[1] + f * G[1], e * D[2] + f * G[2], scale
            )

        d1 = (value(eps, 0) - value(-eps, 0)) / (2 * eps)
        d2 = (
            value(eps, eps) - value(eps, -eps) - value(-eps, eps) + value(-eps, -eps)
        ) / (4 * eps * eps)
        expected1 = np.array(first(to_symbolic(D)).subs(geometry.t, t), dtype=float)
        expected2 = np.array(
            second(to_symbolic(D), to_symbolic(G)).subs(geometry.t, t), dtype=float
        )
        assert np.linalg.norm(d1 - expected1) < 3e-7 * (1 + np.linalg.norm(expected1))
        assert np.linalg.norm(d2 - expected2) < 3e-6 * (1 + np.linalg.norm(expected2))
        assert np.linalg.norm(expected2) > 1e-3


def test_complete_shift_chart_contact_cancels_nonzero_density_Lie_term():
    zero = s.zeros(3)
    D = (s.S.Zero, s.Matrix([1, 0, 0]), zero)
    xi = s.Matrix([0, geometry.t, 0, 0])
    E = geometry.density(s.S.One, -s.S.One)
    lie = geometry.pair(geometry.first(D), geometry.lie_density(E, xi))
    contact = geometry.pair(E, geometry.second(D, geometry.gauge(xi)))
    assert s.cancel(lie + contact) == 0
    assert lie.subs(geometry.t, 0) == -1 and contact.subs(geometry.t, 0) == 1
    assert geometry.hamiltonian_second(D, D) == s.zeros(5)
    assert geometry.second(D, D)[0, 0] != 0


def test_full_wrong_detector_endpoint_flux_is_nonzero():
    E = geometry.density(s.S.One, -s.S.One)
    G = (s.S.One, s.zeros(3, 1), s.zeros(3))
    h = geometry.first(G)
    inv = geometry.metric.inv()
    dE = -(geometry.a**3) * (s.trace(inv * h) * inv / 2 - inv * h * inv) / 2
    left = geometry.primitive(s.S.One, "retarded_source")
    right = geometry.primitive(s.S.One, "advanced_detector")
    wrong = geometry.flux(E, h, dE, s.Matrix([left, 0, 0, 0]))[0]
    good = geometry.flux(E, h, dE, s.Matrix([right, 0, 0, 0]))[0]
    assert wrong.subs(geometry.t, s.Rational(1, 2)) != 0
    assert good.subs(geometry.t, s.Rational(1, 2)) == 0
    assert s.cancel(wrong + geometry.a**3 * left) == 0


def test_full_nonlinear_clock_contact_is_not_deleted_at_finite_projection():
    data = profile.data()
    assert data["complete_Gaussian_nonlinear_clock_contact"] != 0
    assert data["complete_profile_nonlinear_clock_contact"] != 0
    assert data["complete_original_projected_mean_clock_contact"] != 0
    assert (
        s.expand(
            data["complete_Gaussian_nonlinear_clock_contact"]
            + data["complete_profile_nonlinear_clock_contact"]
        )
        == 0
    )
    assert data["whole_physical_profile_one_current"] != 0


@pytest.mark.parametrize("j", range(14))
@pytest.mark.parametrize("b", (6, 7, 28))
def test_all_independent_Fraction_Cauchy_Leibniz_rows(j, b):
    actual = sum(F(comb(j, k) * b * factorial(k) * 4**k) for k in range(j + 1))
    assert actual == F(str(estimates.leibniz(j, b)))


def test_independent_complete_positive_ADM_clock_and_mean_tail_sums():
    def br(j, b):
        return sum(comb(j, k) * b * factorial(k) * 4**k for k in range(j + 1))

    assert 4 * (3 + 4 * br(13, 7) + 2 * br(12, 28)) == 62408287126207901212
    assert 4 * (3 + 4 * br(13, 6)) == 51511602072425569260
    data = estimates.constants()
    kappa = 10**800
    B = 10**400
    tailB = 10**310
    C = F(str(data["complete_actual_spatial_input_normalized"]))
    tail = F(str(data["complete_original_spatial_projection_tail_normalized"]))
    adm = 100 * 10**20 * (C + F(8 * B, kappa)) + F((28520 + 64) * B, kappa)
    clock = 4 * 10**20 * adm
    tailadm = 100 * 10**20 * tail + F(28520 * tailB, kappa)
    tailclock = 4 * 10**20 * tailadm + F(24 * tailB, kappa)
    assert (
        adm
        == F(str(data["complete_actual_ADM_plus_profile_normalized_bound"]))
        < F(1, 10**327)
    )
    assert (
        clock
        == F(str(data["complete_actual_common_clock_normalized_bound"]))
        < F(1, 10**300)
    )
    assert (
        tailclock
        == F(
            str(data["complete_Ward_completed_common_clock_approximant_tail_numerator"])
        )
        < F(1, 10**200)
    )
    assert data["complete_projected_mean_nonlinear_clock_tail_numerator"] > 0
    assert F(16, 9) * 10**301 < 10**310


@pytest.mark.parametrize("cutoff", (2.0, 3.0, 10.0, 100.0, 1000.0, 100000.0))
def test_full_original_frequency_mean_tail_without_a_new_radial_cutoff(cutoff):
    # Mass-scaled Omega=1/x makes the infinite original-frequency tail finite.
    actual = (
        32
        / np.pi**2
        * quad(
            lambda x: x * np.sqrt(1 - x * x), 0, 1 / cutoff, epsabs=1e-25, epsrel=1e-12
        )[0]
    )
    assert actual < 16 / (9 * cutoff * cutoff)


@pytest.mark.parametrize("q", (F(0), F(1, 10**20), F(1), F(1000), F(10**30)))
def test_source_and_detector_spatial_weights_without_detector_time_derivatives(q):
    assert q * q <= (1 + q) ** 2
    assert q * q * (1 + q) ** 6 <= (1 + q) ** 8


@pytest.mark.parametrize("p", (0, s.Rational(1, 10**50), 1, 10**99, 10**500))
def test_all_transfers_and_distinct_full_or_projected_mean_choices(p):
    assert audit.require_domain(0, p) == (0, p)
    assert audit.require_mean("actual_full_reference_mean")
    assert audit.require_mean("original_projected_mean_with_nonlinear_clock_tail")
    assert audit.require_projection(
        "original_two_leg_spatial_one_leg_mean_Ward_completion"
    )
    assert audit.require_cutoff_squared(4 * state.MASS2) == 4 * state.MASS2


def mode_geometry(t):
    a = (1 + t * t) ** 2
    H = 4 * t / (1 + t * t)
    return a, a**3, H


def mode_gauge_profile(t):
    x = t + 0.4
    if x <= 0:
        return 0.0, 0.0, np.zeros(3), np.zeros(3)
    w = np.exp(-1 / x)
    wp = w / x**2
    eta = w * (1 + t / 7)
    etap = wp * (1 + t / 7) + w / 7
    c = np.array([0.2 + 0.1 * t, -0.3 + 0.05 * t * t, 0.4 - 0.06 * t])
    cp = np.array([0.1, 0.1 * t, -0.06])
    return eta, etap, w * c, wp * c + w * cp


def mode_free(t, k, m2):
    a, v, _H = mode_geometry(t)
    return np.array([[0, 1 / v], [-v * (m2 + np.dot(k, k) / a**2), 0]], dtype=complex)


def mode_perturbation(t, k, p, m2):
    a, v, H = mode_geometry(t)
    l = k + p
    eta, n, chi, chip = mode_gauge_profile(t)
    beta = chip - 1j * p * eta / a**2
    Q = 2 * H * eta * np.eye(3) + 1j * (np.outer(p, chi) + np.outer(chi, p))
    trace = np.trace(Q)
    return np.array(
        [
            [1j * np.dot(k, beta), (n - trace / 2) / v],
            [
                -v * ((m2 + np.dot(k, l) / a**2) * (n + trace / 2) - l @ Q @ k / a**2),
                1j * np.dot(l, beta),
            ],
        ],
        dtype=complex,
    )


def mode_initial(k, m2, squeeze):
    a, v, H = mode_geometry(-0.5)
    omega = np.sqrt(m2 + np.dot(k, k) / a**2)
    lam = -H * (np.dot(k, k) / a**2) / (omega * omega)
    phi = 1 / np.sqrt(2 * v * omega)
    pi = (-0.5 * (3 * H + lam) - 1j * omega) * v * phi
    value = np.array([phi, pi])
    alpha = np.sqrt(1 + abs(squeeze) ** 2)
    value = alpha * value + squeeze * value.conjugate()
    assert (
        abs(value[0] * value[1].conjugate() - value[0].conjugate() * value[1] - 1j)
        < 1e-12
    )
    return value


def test_independent_entire_scalar_mode_and_Hermitian_Gram_gauge_tangents():
    largest = 0.0
    largest_gram = 0.0
    wrong = 0.0
    cases = 0
    for k0 in ((0.0, 0.0, 0.0), (0.7, 1.3, -0.4), (2.0, -1.0, 3.0)):
        for p0 in ((0.0, 0.0, 0.0), (0.8, -0.5, 1.1), (-2.0, 0.3, 0.7)):
            for sign in (-1.0, 1.0):
                k = np.array(k0)
                p = sign * np.array(p0)
                l = k + p
                m2 = 9.0
                z0 = mode_initial(k, m2, 0.1 + 0.08j)

                def rhs(t, value, k=k, l=l, p=p, m2=m2):
                    z = value[:2]
                    dz = value[2:]
                    return np.concatenate(
                        (
                            mode_free(t, k, m2) @ z,
                            mode_free(t, l, m2) @ dz
                            + mode_perturbation(t, k, p, m2) @ z,
                        )
                    )

                times = np.linspace(-0.5, 0.5, 61)
                solution = solve_ivp(
                    rhs,
                    (-0.5, 0.5),
                    np.concatenate((z0, np.zeros(2, dtype=complex))),
                    method="DOP853",
                    t_eval=times,
                    rtol=3e-12,
                    atol=3e-13,
                )
                assert solution.success
                for idx, t in enumerate(times):
                    a, _v, _H = mode_geometry(t)
                    eta, _etap, chi, _chip = mode_gauge_profile(t)
                    z = solution.y[:2, idx]
                    actual = solution.y[2:, idx]
                    zdot = mode_free(t, k, m2) @ z
                    # Direct scalar pullback and canonical density variation; not an evolved B matrix.
                    wanted = np.array(
                        [
                            eta * zdot[0] + 1j * np.dot(k, chi) * z[0],
                            eta * zdot[1]
                            + 1j * np.dot(l, chi) * z[1]
                            - a * np.dot(p, k) * eta * z[0],
                        ]
                    )
                    ratio = np.linalg.norm(actual - wanted) / (
                        1 + np.linalg.norm(wanted)
                    )
                    largest = max(largest, ratio)
                    zfull = np.concatenate((z, np.zeros(2, dtype=complex)))
                    afull = np.concatenate((np.zeros(2, dtype=complex), actual))
                    wfull = np.concatenate((np.zeros(2, dtype=complex), wanted))
                    actual_gram = np.outer(afull, zfull.conjugate()) + np.outer(
                        zfull, afull.conjugate()
                    )
                    wanted_gram = np.outer(wfull, zfull.conjugate()) + np.outer(
                        zfull, wfull.conjugate()
                    )
                    gram_ratio = np.linalg.norm(actual_gram - wanted_gram) / (
                        1 + np.linalg.norm(wanted_gram)
                    )
                    largest_gram = max(largest_gram, gram_ratio)
                    missing = np.array(
                        [wanted[0], eta * zdot[1] + 1j * np.dot(l, chi) * z[1]]
                    )
                    wrong = max(wrong, np.linalg.norm(missing - wanted))
                    assert ratio < 2e-10, (k, p, t, ratio)
                    assert gram_ratio < 2e-10, (k, p, t, gram_ratio)
                assert np.linalg.norm(solution.y[2:, -1]) > 1e-4
                cases += 1
    assert wrong > 1e-3
    print(
        "INDEPENDENT_COMPLETE_SCALAR_MODE_AND_HERMITIAN_GRAM_GAUGE_TANGENTS_C0",
        cases,
        flush=True,
    )
    print(
        "MAX_MODE_RELATIVE_ERROR",
        largest,
        "MAX_FULL_GRAM_RELATIVE_ERROR",
        largest_gram,
        "MISSING_DENSITY_GRADIENT_ERROR",
        wrong,
        flush=True,
    )
    assert cases == 18
