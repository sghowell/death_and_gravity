"""Independent full-state, action-variation, auxiliary-jet and scope checks."""

import math
from functools import cache

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from numpy.polynomial.legendre import leggauss
from p8_vacuum_affine_heavy_curved_state import audit, estimates, quantum, state
from scipy.integrate import solve_ivp


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_exact_residual(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_every_full_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_all_unsupported_inputs_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_exact_complete_counts_and_unchanged_frontier():
    assert (
        len(audit.residuals()),
        audit.scalar_entry_count(),
        len(audit.gates()),
        len(audit.controls()),
        audit.rejected_inputs(),
    ) == (57, 57, 74, 9, 343)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 96
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.validate_scope(audit.frontier(), audit.matching())


@pytest.mark.parametrize(
    "u,p,order", [(-1, 0, 0), (1, 10**300, 5), (0, s.Rational(1, 3), 3)]
)
def test_full_closed_reference_estimate_scope(u, p, order):
    assert audit.require_estimate(u, p, order) == (u, p, order)


@pytest.mark.parametrize(
    "u,X,order", [(0, 1, 0), (-1, s.Rational(7, 8), 4), (1, s.Rational(9, 8), 4)]
)
def test_closed_new_clock_four_jet_tube(u, X, order):
    assert audit.require_tube(u, X, order) == (u, X, order)


def test_actual_mass_state_and_strict_energy_cone():
    assert audit.require_state("H8A420-SLE-PRE-T0", state.MASS2) == (
        "H8A420-SLE-PRE-T0",
        state.MASS2,
    )
    assert audit.require_cone(3, 1, 2) == (3, 1, 2)
    assert state.MASS2 == s.Rational(10**200, 512) + 2
    assert 10**196 < state.MASS2 < 10**198
    assert state.KAPPA == 10**800


@pytest.mark.parametrize(
    "stage",
    [
        "specified_exact_heavy_SLE_Hadamard_state",
        "full_dimensional_scalar_counteraction",
        "complete_all_momentum_reference_stress",
        "five_reference_time_derivatives",
        "fixed_QG2_physical_mean_profile",
        "new_combined_clock_four_jet_budget",
        "unchanged_first_loop_vacuum_scattering_jets",
        "complete_flat_Gaussian_one_loop_constant",
    ],
)
def test_precisely_allowed_stage(stage):
    assert audit.require_stage(stage) == stage


@pytest.mark.parametrize("c1,c2", [(3, 1 + 2 * s.I), (1, 0), (7, -2 + 3 * s.I)])
def test_independent_phase_regular_minimizer_and_basis_invariance(c1, c2):
    with mp.workdps(80):

        def choose(energy, pair):
            gap = mp.sqrt(energy**2 - abs(pair) ** 2)
            return (
                mp.sqrt((energy + gap) / (2 * gap)),
                -pair / mp.sqrt(2 * gap * (energy + gap)),
                gap,
            )

        def transform(energy, pair, A, B):
            return (
                (abs(A) ** 2 + abs(B) ** 2) * energy + 2 * mp.re(A * mp.conj(B) * pair),
                A * A * pair + B * B * mp.conj(pair) + 2 * A * B * energy,
            )

        energy = mp.mpf(c1)
        pair = mp.mpc(str(s.re(c2)), str(s.im(c2)))
        a, b, gap = choose(energy, pair)
        exact = state.sle_coefficients(
            s.Integer(c1), s.sympify(c2), s.conjugate(s.sympify(c2))
        )
        assert abs(a - mp.mpf(str(s.N(exact[0], 80)))) < mp.mpf("1e-70")
        assert abs(
            b - mp.mpc(str(s.N(s.re(exact[1]), 80)), str(s.N(s.im(exact[1]), 80)))
        ) < mp.mpf("1e-70")
        e, q = transform(energy, pair, a, b)
        assert abs(abs(a) ** 2 - abs(b) ** 2 - 1) < mp.mpf("1e-70")
        assert abs(e - gap) < mp.mpf("1e-70") and abs(q) < mp.mpf("1e-70")
        A = mp.cosh(mp.mpf(".73"))
        B = mp.exp(mp.j * mp.mpf(".41")) * mp.sinh(mp.mpf(".73"))
        aa, bb, gg = choose(*transform(energy, pair, A, B))
        changed = (aa * A + bb * mp.conj(B), aa * B + bb * mp.conj(A))
        for i in range(2):
            for j in range(2):
                assert abs(
                    (a, b)[i] * mp.conj((a, b)[j]) - changed[i] * mp.conj(changed[j])
                ) < mp.mpf("1e-68")
        assert abs(gg - gap) < mp.mpf("1e-68")
        if s.im(c2) != 0:
            assert abs(transform(energy, pair, a, mp.conj(b))[1]) > mp.mpf(".1")
        else:
            assert b == 0


def test_sampling_is_fixed_before_cauchy_and_not_instantaneous():
    assert state.SAMPLING_INTERVAL == (-s.Rational(3, 4), -s.Rational(1, 2))
    assert state.T0 == state.SAMPLING_INTERVAL[1]
    with mp.workdps(60):
        normal = mp.quad(lambda x: mp.exp(-2 / (1 - x * x)), [-1, 0, 1]) / 8
        assert normal > mp.mpf(1) / 216
        assert normal < mp.mpf(1) / 4
    assert "Do not re-minimize" in state.data()["state_family"]
    assert "NOT the WKB comparison mode" in state.data()["WKB_comparison"]


@pytest.mark.parametrize("momentum", (0, 2, 10))
def test_exact_ODE_state_selection_quadrature_phase_and_basis_convergence(momentum):
    first = reference_sle(25, momentum, 160)
    second = reference_sle(25, momentum, 256)
    assert np.max(abs(first - second)) < 2e-11


@pytest.mark.parametrize("ratio", (0, 1, 4))
def test_full_actual_hierarchy_WKB_residual_is_resolved_not_underflow(ratio):
    low = actual_diagnostic(1800, ratio, "0.333333333333333333333333333333")
    high = actual_diagnostic(2000, ratio, "0.333333333333333333333333333333")
    with mp.workdps(200):
        assert abs(low - high) < mp.mpf("1e-150") * max(1, abs(high))
        assert abs(high) > 1
        assert abs(high) < mp.mpf(10) ** 60


@pytest.mark.parametrize("ratio", (0, 1, 4))
def test_independent_direct_physical_mode_energy_volume_factor(ratio):
    with mp.workdps(100):
        t = mp.mpf(1) / 3
        it, om2, H, inv_a = full_iterates(mp.mpf(10000), mp.mpf(ratio), t)
        W = it[6]
        a = (1 + t * t) ** 2
        mode = a ** (-mp.mpf(3) / 2) / mp.sqrt(2 * W[0])
        velocity = (-3 * H[0] / 2 - W[1] / (2 * W[0]) - mp.j * W[0]) * mode
        energy = (abs(velocity) ** 2 + om2[0] * abs(mode) ** 2) / 2
        pressure = (
            abs(velocity) ** 2 - (mp.mpf(10000) + (om2[0] - 10000) / 3) * abs(mode) ** 2
        ) / 2
        rate = 3 * H[0] / 2 + W[1] / (2 * W[0])
        direct = inv_a[0] ** 3 * (W[0] ** 2 + om2[0] + rate**2) / (4 * W[0])
        direct_p = (
            inv_a[0] ** 3
            * (W[0] ** 2 - mp.mpf(10000) - (om2[0] - 10000) / 3 + rate**2)
            / (4 * W[0])
        )
        assert abs(energy - direct) < mp.mpf("1e-90")
        assert abs(pressure - direct_p) < mp.mpf("1e-90")
        wrong_volume = (
            inv_a[0] ** mp.mpf("1.5") * (W[0] ** 2 + om2[0] + rate**2) / (4 * W[0])
        )
        assert abs(energy - wrong_volume) > 1
        assert abs(
            a**3 * (mode * mp.conj(velocity) - mp.conj(mode) * velocity) - mp.j
        ) < mp.mpf("1e-90")


@pytest.mark.parametrize("dimension", (3, 4, 5))
@pytest.mark.parametrize("label", ("constant", "R", "R2", "Ric2", "Riem2"))
def test_independent_full_lapse_and_scale_Euler_variation(dimension, label):
    assert independent_euler_stress(dimension, label) == (0, 0)


def test_independent_complete_adiabatic_readout_not_just_Riccati_coefficients():
    tau, w = s.symbols("tau omega", positive=True)
    p2, p4, b2, r0, z, D = s.symbols("p2 p4 b2 rate0 z D", real=True)
    full_w = w * (1 + p2 * tau / w**2 + p4 * tau**2 / w**4)
    full_rate = r0 + b2 * tau / (2 * w**2)
    rho = (full_w**2 + w * w + tau * full_rate**2) / (4 * full_w)
    pressure = (full_w**2 - (1 - 2 * z / D) * w * w + tau * full_rate**2) / (4 * full_w)
    ad = quantum.adiabatic()
    sub = {
        p2: ad["P2"],
        p4: ad["P4"],
        b2: quantum.dt(ad["P2"]) + 2 * quantum.H * quantum.z * ad["P2"],
        r0: (quantum.D * quantum.H - quantum.H * quantum.z) / 2,
        z: quantum.z,
        D: quantum.D,
    }
    for name, readout in (("rho", rho), ("pressure", pressure)):
        series = s.series(readout, tau, 0, 3).removeO().expand()
        for order in range(3):
            coefficient = s.cancel(series.coeff(tau, order) * w ** (2 * order - 1))
            assert (
                s.expand(coefficient.subs(sub, simultaneous=True) - ad[name][order])
                == 0
            )


@pytest.mark.parametrize("order", range(3))
def test_complete_subtraction_mode_Ward_identity_in_full_dimension(order):
    ad = quantum.adiabatic()
    rho, pressure = ad["rho"][order], ad["pressure"][order]
    lam = -quantum.H * quantum.z
    residual = (
        quantum.dt(rho) + (1 - 2 * order) * lam * rho + quantum.D * quantum.H * pressure
    )
    assert s.expand(residual) == 0


def test_complete_dimensional_finite_terms_and_false_component_shortcut():
    actual, component_only, checks = quantum.local()
    assert all(value == 0 for value in checks.values())
    H, Hd, Hdd, Hddd = quantum.h[:4]
    ell = quantum.ell
    wanted = {
        0: {"rho": ell - s.Rational(3, 2), "pressure": s.Rational(3, 2) - ell},
        1: {
            "rho": -2 * H * H * (ell - 1),
            "pressure": 2 * (3 * H * H + 2 * Hd) * (ell - 1) / 3,
        },
        2: {
            "rho": ell * (6 * H * H * Hd + 2 * H * Hdd - Hd * Hd),
            "pressure": -ell
            * (18 * H * H * Hd + 12 * H * Hdd + 9 * Hd * Hd + 2 * Hddd)
            / 3,
        },
    }
    for order in range(3):
        for label in ("rho", "pressure"):
            assert s.expand(actual[order][label] - wanted[order][label]) == 0
            assert s.denom(actual[order][label]).has(H) is False
    assert s.expand(actual[2]["rho"] - component_only[2]["rho"]) != 0
    assert quantum.radial(quantum.adiabatic()["rho"][0], 0)[0] == -1


def test_full_state_error_and_exact_radial_bound_survive_all_derivatives():
    data = estimates.data()
    for key in (
        "complete_exact_mode_stress_difference_coefficients",
        "complete_SLE_stress_difference_coefficients",
    ):
        assert set(data[key]) == set(range(6))
        assert 0 < max(data[key].values()) < 10**300
    assert set(data["complete_sampling_derivative_bounds"]) == set(range(11))
    assert (
        data["complete_state_dependent_and_subtraction_stress_upper"]
        == s.Rational(10**310, 1) / state.MASS2
    )
    assert (
        0
        < data["complete_normalized_heavy_stress_all_five_jets_upper"]
        < s.Rational(1, 10**400)
    )
    # The tighter time margin is indispensable for the advertised 10^20 constant.
    assert math.factorial(5) * 2048**5 < 10**20 < math.factorial(5) * 4096**5
    m = s.Symbol("mass", positive=True)
    p = s.Symbol("p", nonnegative=True)
    primitive = 64 * p**3 / (3 * m**2 * (16 * m**2 + p**2) ** s.Rational(3, 2))
    assert (
        s.simplify(
            s.diff(primitive, p) - p**2 * (m * m + p * p / 16) ** (-s.Rational(5, 2))
        )
        == 0
    )
    assert s.limit(primitive, p, s.oo) == 64 / (3 * m * m)


def test_whole_fixed_profile_vacuum_constants_and_clock_Ward():
    profile = quantum.fixed_profile()
    X, u = state.X, state.TIME
    rho, P = profile["rho"], profile["pressure"]
    F = profile["DeltaF"]
    assert (
        s.cancel(F.subs(X, 0) + (profile["p_v_H"] + profile["p_v_Phi"]) / state.KAPPA)
        == 0
    )
    assert s.diff(F.subs(X, 0), u) == 0
    A = F.subs(X, 1)
    B = s.diff(F, X).subs(X, 1)
    assert s.cancel(A + P / state.KAPPA) == 0
    assert s.cancel(2 * B - A + rho / state.KAPPA) == 0
    ward = s.diff(rho, u) + 3 * state.HUBBLE * (rho + P)
    assert (
        s.cancel(
            s.diff(A, u) - 2 * s.diff(B, u) - 6 * state.HUBBLE * B - ward / state.KAPPA
        )
        == 0
    )
    assert profile["p_v_Phi"] == s.Rational(3, 128) / s.pi**2
    assert profile["p_v_H"].subs(state.MASS2, state.MASS2) != 0
    denominator = X**1024 + (1 - X) ** 1024
    assert profile["T"] == X**1024 / denominator
    assert s.factor(denominator * (1 - profile["T"]) - (1 - X) ** 1024) == 0
    assert 2 * 1024 - 2 > 4
    assert s.Rational(3, 10**400) + s.Rational(1, 10**2700) < s.Rational(1, 10**399)


def test_classical_heavy_source_still_has_full_eighth_order_clock_zero():
    X = state.germs.X
    h = state.germs.parent.switch(X, state.germs.parent.A)
    assert h == (1 - X) ** 8 * s.exp(-state.germs.parent.A * X**2)
    assert h.subs(X, 1) == 0
    assert s.diff(h, X).subs(X, 1) == 0
    assert s.diff(h, X, 2).subs(X, 1) == 0


def test_actual_reference_and_new_Hessian_scope_not_erased():
    data = quantum.data()
    assert "CHANGE the classical light Hessian" in data["domain_and_response_boundary"]
    assert "not claimed real analytic" in data["fixed_locality"]
    assert (
        "not adjustable independent functions"
        in data["full_reference_stress_definition"]
    )
    assert "not a finite-gravity quantum-limit" in data["flat_vacuum_matching"]
    assert "original V/G/B/P8 remain open" in audit.observable()["remaining"]
    assert audit.frontier() == audit.previous.frontier()


def reference_sle(mass2, momentum, quadrature_order):
    def scale(t):
        return (1 + t * t) ** 2

    def hubble(t):
        return 4 * t / (1 + t * t)

    def omega2(t):
        return mass2 + momentum**2 / scale(t) ** 2

    start = -1.0
    om = np.sqrt(omega2(start))
    lam = -hubble(start) * (momentum**2 / scale(start) ** 2) / omega2(start)
    mode = 1 / np.sqrt(2 * scale(start) ** 3 * om)
    y0 = np.array(
        [mode, (-1.5 * hubble(start) - lam / 2 - 1j * om) * mode], dtype=complex
    )

    def rhs(t, y):
        return np.array([y[1], -3 * hubble(t) * y[1] - omega2(t) * y[0]])

    solution = solve_ivp(
        rhs,
        (-1, 1),
        y0,
        method="DOP853",
        rtol=1e-12,
        atol=1e-14,
        max_step=0.01,
        dense_output=True,
    )
    assert solution.success
    nodes, weights = leggauss(quadrature_order)
    times = -0.625 + nodes / 8
    sampling = np.exp(-2 / (1 - nodes * nodes))
    measure = weights * sampling / 8
    measure /= np.sum(measure)
    modes = solution.sol(times)

    def moments(values):
        S, Sd = values
        c1 = np.sum(measure * (abs(Sd) ** 2 + omega2(times) * abs(S) ** 2)) / 2
        c2 = np.sum(measure * (Sd**2 + omega2(times) * S**2)) / 2
        return c1, c2

    def choose(c1, c2):
        delta = np.sqrt(c1 * c1 - abs(c2) ** 2)
        return (
            np.sqrt((c1 + delta) / (2 * delta)),
            -c2 / np.sqrt(2 * delta * (c1 + delta)),
            delta,
        )

    c1, c2 = moments(modes)
    assert c1 > abs(c2)
    alpha, beta, delta = choose(c1, c2)
    selected = alpha * modes + beta * np.conj(modes)
    final1, final2 = moments(selected)
    assert abs(final1 - delta) < 1e-12
    assert abs(final2) < 1e-12
    A = np.cosh(0.73)
    B = np.exp(0.41j) * np.sinh(0.73)
    transformed = A * modes + B * np.conj(modes)
    aa, bb, dd = choose(*moments(transformed))
    assert abs(dd - delta) < 1e-11
    target_times = np.linspace(-1, 1, 13)
    exact = solution.sol(target_times)
    selected_exact = alpha * exact + beta * np.conj(exact)
    selected_transformed = aa * (A * exact + B * np.conj(exact)) + bb * np.conj(
        A * exact + B * np.conj(exact)
    )
    for i in range(13):
        first = selected_exact[:, i]
        second = selected_transformed[:, i]
        assert (
            np.max(
                abs(np.outer(first, np.conj(first)) - np.outer(second, np.conj(second)))
            )
            < 1e-10
        )
        wronskian = scale(target_times[i]) ** 3 * (
            first[0] * np.conj(first[1]) - np.conj(first[0]) * first[1]
        )
        assert abs(wronskian - 1j) < 1e-10
    return np.array([c1, c2, delta, beta], dtype=complex)


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
    return [c] + [mp.mpf(0)] * (length - 1)


def full_iterates(mass2, ratio, time, epsilon=1):
    length = 17
    coordinate = [time, mp.mpf(1)] + [mp.mpf(0)] * (length - 2)
    v = add(constant(1, length), mul(coordinate, coordinate))
    inv_v = inverse(v)
    inv_a = mul(inv_v, inv_v)
    omega2 = add(
        constant(mass2, length), scale(mul(inv_a, inv_a), mass2 * ratio * ratio)
    )
    H = scale(mul(coordinate, inv_v), 4)
    U = add(scale(derivative(H), mp.mpf(3) / 2), scale(mul(H, H), mp.mpf(9) / 4))
    current = root(omega2)
    iterates = [current]
    for _ in range(7):
        rate = mul(derivative(current), inverse(current))
        bracket = add(
            scale(U, -1),
            add(
                scale(derivative(rate), -mp.mpf(1) / 2),
                scale(mul(rate, rate), mp.mpf(1) / 4),
            ),
        )
        square = add(omega2, scale(bracket, epsilon * epsilon))
        current = root(square)
        iterates.append(current)
    return iterates, omega2, H, inv_a


def actual_diagnostic(precision, ratio, time):
    with mp.workdps(precision):
        n = mp.mpf(10) ** 200 / 512 + 2
        ts = mp.mpf(time) + mp.j / 2048
        values, _omega2, _H, _inv_a = full_iterates(n, mp.mpf(ratio), ts)
        Omega = mp.sqrt(n * (1 + mp.mpf(ratio) ** 2 / 16))
        defect = (values[7][0] ** 2 - values[6][0] ** 2) * Omega**12
        assert abs(defect) < mp.mpf(10) ** 60
        epsilon = Omega / mp.mpf(10) ** 8 * mp.exp(mp.j * mp.mpf("0.37"))
        aux, aux_omega2, aux_H, aux_inv_a = full_iterates(n, mp.mpf(ratio), ts, epsilon)
        for W in aux:
            assert abs(W[0] / mp.sqrt(aux_omega2[0]) - 1) < mp.mpf("0.01")
        W = aux[6]
        rate = W[1] / W[0]
        d = 3 * aux_H[0] / 2 + rate / 2
        a3inv = aux_inv_a[0] ** 3
        energy = a3inv * (W[0] ** 2 + aux_omega2[0] + epsilon**2 * d * d) / (4 * W[0])
        assert abs(energy / Omega) < 10**4
        return +defect


@cache
def independent_euler_stress(dimension, label):
    a = s.symbols("a0:5", positive=True)
    N = s.symbols("N0:4", positive=True)

    def total(expr):
        return s.expand(
            sum(a[j + 1] * s.diff(expr, a[j]) for j in range(4))
            + sum(N[j + 1] * s.diff(expr, N[j]) for j in range(3))
        )

    H = a[1] / (a[0] * N[0])
    Hd = total(H) / N[0]
    D = s.Integer(dimension)
    R = 2 * D * Hd + D * (D + 1) * H**2
    Ric = D**2 * (Hd + H**2) ** 2 + D * (Hd + D * H**2) ** 2
    Riem = 4 * D * (Hd + H**2) ** 2 + 2 * D * (D - 1) * H**4
    cov = {"constant": s.S.One, "R": R, "R2": R * R, "Ric2": Ric, "Riem2": Riem}[label]
    lag = s.expand(N[0] * a[0] ** D * cov)
    eN = s.diff(lag, N[0]) - total(s.diff(lag, N[1]))
    ea = s.diff(lag, a[0]) - total(s.diff(lag, a[1])) + total(total(s.diff(lag, a[2])))
    rho = s.cancel(-eN / a[0] ** D)
    pressure = s.cancel(ea / (D * N[0] * a[0] ** (D - 1)))
    sample = {
        a[0]: s.Rational(7, 5),
        a[1]: s.Rational(2, 3),
        a[2]: -s.Rational(3, 7),
        a[3]: s.Rational(5, 11),
        a[4]: -s.Rational(7, 13),
        N[0]: 1,
        N[1]: 0,
        N[2]: 0,
        N[3]: 0,
    }
    proper_jets = [H]
    for _ in range(3):
        proper_jets.append(total(proper_jets[-1]) / N[0])
    hsub = dict(zip(quantum.h[:4], [s.cancel(v.subs(sample)) for v in proper_jets]))
    RR, RRic, RRiem = quantum.curvatures()
    target_cov = {
        "constant": s.S.One,
        "R": RR,
        "R2": RR**2,
        "Ric2": RRic,
        "Riem2": RRiem,
    }[label]
    target_rho = quantum.rho_of(target_cov)
    target_p = quantum.pressure_of(target_rho)
    target_sub = {quantum.D: D, **hsub}
    return (
        s.cancel(rho.subs(sample) - target_rho.subs(target_sub)),
        s.cancel(pressure.subs(sample) - target_p.subs(target_sub)),
    )
