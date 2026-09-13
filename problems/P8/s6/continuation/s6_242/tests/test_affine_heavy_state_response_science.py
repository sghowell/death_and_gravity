"""Independent full ADM, ordered quantum, state-profile, regulator and scope checks."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_heavy_state_response import audit, bounds, response, vertices


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_exact_scalar_and_matrix_entry(name):
    value = audit.residuals()[name]
    assert all(e == 0 for e in value) if isinstance(value, s.MatrixBase) else value == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_every_full_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_exact_counts_and_unchanged_frontier():
    assert (
        len(audit.residuals()),
        audit.scalar_entry_count(),
        len(audit.gates()),
        len(audit.controls()),
        audit.rejected_inputs(),
        len(audit.matching()),
    ) == (29, 233, 44, 9, 273, 98)
    assert len(audit.frontier()) == 9
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.validate_scope(audit.frontier(), audit.matching()) is True


@pytest.mark.parametrize(
    "time,momentum", [(-s.Rational(1, 2), 0), (0, 1), (s.Rational(1, 2), 10**300)]
)
def test_all_closed_time_and_momentum_boundaries(time, momentum):
    assert audit.require_domain(time, momentum) == (time, momentum)


@pytest.mark.parametrize(
    "label,order",
    [
        ("full_ADM_feature_spatial_graph", 1),
        ("same_clock_scalar_graph", 2),
        ("external_time_derivatives", 0),
    ],
)
def test_complete_graph_orders(label, order):
    assert audit.require_graph(label, order) == (label, order)


@pytest.mark.parametrize(
    "label,rule",
    [
        ("even", "minus_imaginary_part"),
        ("odd", "plus_i_real_part"),
        ("general", "both_reflected_ordered_products"),
    ],
)
def test_actual_ordered_parity_rules(label, rule):
    assert audit.require_parity(label, rule) == (label, rule)


@pytest.mark.parametrize(
    "value", [4 * audit.state.MASS2, 16 * audit.state.MASS2, 10**400]
)
def test_auxiliary_cutoff_exact_boundary(value):
    assert audit.require_cutoff_squared(value) == value


@pytest.mark.parametrize(
    "stage",
    [
        "full_ADM_Hamiltonian_features",
        "full_second_metric_contact",
        "exact_covariance_and_ordered_Wick_bridge",
        "complete_SLE_minus_exact_comparison",
        "full_all_internal_and_external_momentum_bound",
        "full_two_leg_regulator_tail",
        "unchanged_fixed_reference_profile_split",
        "complete_matched_common_clock_pullback",
    ],
)
def test_supported_stages(stage):
    assert audit.require_stage(stage) == stage


def test_full_isotropic_reference_current_not_profile_only():
    E, L, Z = s.symbols("energy_squared gradient_squared mass_field_squared", real=True)
    covariance = s.diag(E, L / 3, L / 3, L / 3, Z)
    M = vertices.first(vertices.nD, vertices.bD, vertices.D)
    current = -s.trace(M * covariance) / 2
    rho = (E + L + Z) / 2
    pressure = (E - L / 3 - Z) / 2
    wanted = -rho * vertices.nD + pressure * s.trace(vertices.D) / 2
    assert s.expand(current - wanted) == 0
    assert (
        s.expand(current + rho * vertices.nD - pressure * s.trace(vertices.D) / 2) == 0
    )
    assert current.diff(vertices.nD) != 0
    assert current.diff(vertices.D[0, 0]) != 0


def test_literal_profile_full_clock_mixed_variation():
    ed, eg = vertices.ed, vertices.eg
    nd, ng, vd, vg, r, p, delta = s.symbols("nd ng vd vg rho pressure delta", real=True)
    N = 1 + ed * nd + eg * ng
    v = ed * vd + eg * vg
    X = N**-2
    # Exact R clock two-jet from the frozen full parent.
    R = 1 + 2 * delta * (X - 1)
    T = X**audit.state.N / (X**audit.state.N + (1 - X) ** audit.state.N)
    density = (
        N * s.exp(3 * v) * R ** (-s.Rational(3, 4)) * T * (-p - (r + p) * (X - 1) / 2)
    )
    actual = s.diff(density, ed, eg).subs({ed: 0, eg: 0})
    dJ = (21 * delta**2 - 3 * delta) * (-p) / 2 + (1 - 6 * delta) * (-(r + p) / 2)
    wanted = (
        2 * dJ * nd * ng
        + 3 * (r - 3 * delta * p) * (nd * vg + ng * vd)
        - 9 * p * vd * vg
    )
    assert s.expand(actual - wanted) == 0


def test_independent_full_legendre_transform():
    N, volume = s.symbols("N spatial_volume", positive=True)
    velocity, pi, gradient2, potential, bgrad = s.symbols(
        "velocity pi gradient2 potential shift_dot_gradient", real=True
    )
    lag = (
        N
        * volume
        * ((velocity - bgrad) ** 2 / (2 * N * N) - gradient2 / 2 - potential / 2)
    )
    equation = s.diff(lag, velocity)
    solution = N * pi / volume + bgrad
    assert s.expand(equation.subs(velocity, solution) - pi) == 0
    H = (pi * velocity - lag).subs(velocity, solution)
    wanted = N * (pi * pi / volume + volume * (gradient2 + potential)) / 2 + pi * bgrad
    assert s.cancel(H - wanted) == 0
    assert s.diff(wanted, bgrad, 2) == 0
    assert s.diff(wanted, N, 2) == 0
    assert s.diff(wanted, volume, bgrad) == 0


@pytest.mark.parametrize("fixture", [1, 2, 3])
def test_literal_noncommuting_ADM_exponential_derivatives(fixture):
    with mp.workdps(75):
        DD = mp.matrix(
            [
                [mp.mpf(1) / 3, mp.mpf(fixture) / 7, -mp.mpf(1) / 5],
                [mp.mpf(fixture) / 7, -mp.mpf(2) / 9, mp.mpf(1) / 8],
                [-mp.mpf(1) / 5, mp.mpf(1) / 8, mp.mpf(1) / 6],
            ]
        )
        GG = mp.matrix(
            [
                [mp.mpf(1) / 4, -mp.mpf(1) / 6, mp.mpf(fixture) / 11],
                [-mp.mpf(1) / 6, mp.mpf(1) / 5, -mp.mpf(1) / 7],
                [mp.mpf(fixture) / 11, -mp.mpf(1) / 7, -mp.mpf(1) / 3],
            ]
        )
        assert mp.norm(DD * GG - GG * DD) > 0.001
        bd = mp.matrix([mp.mpf(1) / 9, -mp.mpf(1) / 8, mp.mpf(1) / 7])
        bg = mp.matrix([-mp.mpf(1) / 6, mp.mpf(1) / 11, mp.mpf(1) / 13])
        grad = mp.matrix([mp.mpf(1) / 5, -mp.mpf(2) / 7, mp.mpf(3) / 8])
        aa = mp.mpf(5) / 4
        mass = mp.mpf(7) / 3
        field = mp.mpf(2) / 5
        pi = mp.mpf(3) / 7
        nd = mp.mpf(1) / 6
        ng = -mp.mpf(1) / 7

        def full(d, g):
            Q = d * DD + g * GG
            N = 1 + d * nd + g * ng
            beta = d * bd + g * bg
            volume = aa**3 * mp.exp(sum(Q[i, i] for i in range(3)) / 2)
            inv = mp.expm(-Q) / aa**2
            return (
                N
                * (
                    pi * pi / volume
                    + volume * ((grad.T * inv * grad)[0] + mass * mass * field * field)
                )
                / 2
                + pi * (beta.T * grad)[0]
            )

        eps = mp.mpf("1e-15")
        first = (full(eps, 0) - full(-eps, 0)) / (2 * eps)
        mixed = (
            full(eps, eps) - full(eps, -eps) - full(-eps, eps) + full(-eps, -eps)
        ) / (4 * eps * eps)
        td = sum(DD[i, i] for i in range(3))
        tg = sum(GG[i, i] for i in range(3))
        f = mp.matrix([pi / aa**3, *[grad[i] / aa for i in range(3)], mass * field])
        B_D = td * mp.eye(3) / 2 - DD
        B_G = tg * mp.eye(3) / 2 - GG
        MD = mp.zeros(5)
        CDG = mp.zeros(5)
        MD[0, 0] = nd - td / 2
        MD[4, 4] = nd + td / 2
        CDG[0, 0] = td * tg / 4 - (nd * tg + ng * td) / 2
        CDG[4, 4] = td * tg / 4 + (nd * tg + ng * td) / 2
        block = (B_D * B_G + B_G * B_D) / 2 + nd * B_G + ng * B_D
        for i in range(3):
            MD[0, i + 1] = MD[i + 1, 0] = aa * bd[i]
            for j in range(3):
                MD[i + 1, j + 1] = B_D[i, j] + nd * (i == j)
                CDG[i + 1, j + 1] = block[i, j]
        expected_first = aa**3 * (f.T * MD * f)[0] / 2
        expected_mixed = aa**3 * (f.T * CDG * f)[0] / 2
        assert abs(first - expected_first) < mp.mpf("1e-27")
        assert abs(mixed - expected_mixed) < mp.mpf("1e-27")
        assert abs((f.T * CDG * f)[0]) > mp.mpf("1e-5")


@pytest.mark.parametrize("shift", [0, mp.mpf("0.5"), 10, 10000])
def test_full_continuous_momentum_convolution_not_truncated(shift):
    with mp.workdps(55):
        shift = mp.mpf(shift)

        def mean_frequency(z):
            if shift == 0 or z == 0:
                return mp.sqrt(1 + z * z + shift * shift)
            return (
                (1 + (z + shift) ** 2) ** mp.mpf("1.5")
                - (1 + (z - shift) ** 2) ** mp.mpf("1.5")
            ) / (6 * z * shift)

        integral = mp.quad(
            lambda z: z * z * (1 + z * z) ** -5 * mean_frequency(z), [0, 1, 4, mp.inf]
        )
        ratio = integral / (mp.mpf(8) / 105)
        assert ratio <= 1 + shift + mp.mpf("1e-45")
        assert ratio > 0
        if shift == 0:
            assert abs(ratio - 1) < mp.mpf("1e-45")
        if shift == 10000:
            assert ratio > 7000


def test_nonzero_contact_and_exact_same_clock_graph_constants():
    d = bounds.data()
    contact = d["retained_profile_second_clock_contact"]
    vals = {
        sym: (mp.mpf("0.5") if str(sym) == "delta" else 1)
        for sym in contact.free_symbols
    }
    assert abs(float(contact.subs(vals)) - 1.5) < 1e-14
    assert s.expand(contact + d["retained_state_current_second_clock_contact"]) == 0
    assert d["complete_normalized_same_clock_response_bound"] == s.Rational(
        10**103, 1
    ) / (audit.state.MASS2**3 * audit.state.KAPPA)
    assert d["complete_normalized_clock_regulator_tail_coefficient"] == s.Rational(
        10**104, 1
    ) / (audit.state.MASS2 ** s.Rational(5, 2) * audit.state.KAPPA)


def test_reference_modes_are_not_reselected_and_comparison_remains_open():
    d = response.data()
    assert "ZERO initial tangent" in d["full_covariance_tangent"]
    assert "not asserted Hadamard" in d["unchanged_preparation"]
    assert "not set to zero" in d["entire_state_difference"]
    assert "BOTH projected internal legs" in d["regulator_and_limit"]
    assert (
        "No candidate, prescribed sum or finite renormalization changes."
        in d["reference_profile_split"]
    )
    assert "P8" in bounds.data()["not_full_heavy_response"]


def test_full_ordered_lattice_covariance_and_prepared_metric_evolution():
    """Independent finite-lattice Gaussian covariance versus full ordered Fourier current."""
    import numpy as np
    from numpy.testing import assert_allclose
    from scipy.linalg import expm

    L = 7
    mass = np.sqrt(7.0)
    x = 2 * np.pi * np.arange(L) / L
    mom = np.fft.fftfreq(L, d=1 / L).astype(int)
    Fourier = np.exp(1j * x[:, None] * mom[None, :]) / np.sqrt(L)
    derivative = np.real_if_close(Fourier @ np.diag(1j * mom) @ Fourier.conj().T).real
    eye = np.eye(L)
    J = np.block([[np.zeros((L, L)), eye], [-eye, np.zeros((L, L))]])
    M0 = np.block(
        [
            [mass**2 * eye - derivative @ derivative, np.zeros((L, L))],
            [np.zeros((L, L)), eye],
        ]
    )
    omega = np.sqrt(mass**2 + mom**2)

    def modes(time, target):
        b0 = 0.15 * np.exp(1j * (0.4 + 0.03 * mom * mom))
        aa0 = np.sqrt(1 + abs(b0) ** 2)
        f = (
            aa0 * np.exp(-1j * omega * time) + b0 * np.exp(1j * omega * time)
        ) / np.sqrt(2 * omega)
        fd = (
            -1j * omega * aa0 * np.exp(-1j * omega * time)
            + 1j * omega * b0 * np.exp(1j * omega * time)
        ) / np.sqrt(2 * omega)
        if target:
            b1 = 0.07 * np.exp(1j * (0.8 + 0.04 * mom * mom))
            aa1 = np.sqrt(1 + abs(b1) ** 2)
            f, fd = aa1 * f + b1 * f.conj(), aa1 * fd + b1 * fd.conj()
        return f, fd

    def covariance(time, target):
        f, fd = modes(time, target)

        def lift(value):
            return np.real_if_close(Fourier @ np.diag(value) @ Fourier.conj().T).real

        return np.block(
            [
                [lift(abs(f) ** 2), lift(np.real(f * fd.conj()))],
                [lift(np.real(f * fd.conj())), lift(abs(fd) ** 2)],
            ]
        )

    def feature(kind):
        if kind == "N":
            return np.eye(5)
        if kind == "Q":
            return np.diag([-0.5, -0.5, 0.5, 0.5, 0.5])
        if kind == "beta":
            out = np.zeros((5, 5))
            out[0, 1] = out[1, 0] = 1
            return out
        raise ValueError(kind)

    def ham(kind, weight):
        w = np.diag(weight)
        z = np.zeros((L, L), dtype=complex)
        if kind == "N":
            return np.block([[derivative.T @ w @ derivative + mass**2 * w, z], [z, w]])
        if kind == "Q":
            return np.block(
                [
                    [-0.5 * derivative.T @ w @ derivative + 0.5 * mass**2 * w, z],
                    [z, -0.5 * w],
                ]
            )
        if kind == "beta":
            return np.block([[z, derivative.T @ w], [w @ derivative, z]])
        raise ValueError(kind)

    def contact_feature(d, g):
        if "beta" in (d, g) or (d, g) == ("N", "N"):
            return np.zeros((5, 5))
        if (d, g) in (("N", "Q"), ("Q", "N")):
            return feature("Q")
        if (d, g) == ("Q", "Q"):
            return np.eye(5) / 4
        raise ValueError((d, g))

    def contact_ham(d, g, weight):
        c = contact_feature(d, g)
        w = np.diag(weight)
        z = np.zeros((L, L), dtype=complex)
        return np.block(
            [
                [c[1, 1] * derivative.T @ w @ derivative + c[4, 4] * mass**2 * w, z],
                [z, c[0, 0] * w],
            ]
        )

    def readout(time, target):
        f, fd = modes(time, target)
        return np.stack([fd, 1j * mom * f, np.zeros(L), np.zeros(L), mass * f], axis=1)

    def ordered(d, g, P, t, s, target):
        ft, fs = readout(t, target), readout(s, target)
        md, mg = feature(d), feature(g)
        total = 0j
        wrong = 0j
        for ip, p in enumerate(mom):
            iq = int((P - p) % L)
            im = int((-p) % L)
            jm = int((p - P) % L)
            Vd = ft[ip] @ md @ ft[iq]
            Vg = fs[ip] @ mg @ fs[iq]
            reverse = complex(ft[im] @ md @ ft[jm]).conjugate() * (fs[im] @ mg @ fs[jm])
            product = Vd * complex(Vg).conjugate()
            total += 0.5j * (product - reverse)
            wrong += -np.imag(product)
        return total, wrong

    t = 0.37
    s = -0.21
    prop = expm(J @ M0 * (t - s))
    maximum = 0.0
    odd_errors = []
    for P in (0, 1, 2):
        for d in ("N", "Q", "beta"):
            for g in ("N", "Q", "beta"):
                wd = np.exp(-1j * P * x)
                wg = np.exp(1j * P * x)
                md, mg = ham(d, wd), ham(g, wg)
                for target in (False, True):
                    cs = covariance(s, target)
                    force = J @ mg @ cs + cs @ (J @ mg).T
                    tangent = -np.trace(md @ prop @ force @ prop.T) / 2
                    fourier, wrong = ordered(d, g, P, t, s, target)
                    error = abs(tangent - fourier)
                    maximum = max(maximum, error)
                    assert_allclose(
                        tangent,
                        fourier,
                        rtol=1e-11,
                        atol=2e-11,
                        err_msg=str((P, d, g, target)),
                    )
                    ct = (
                        -np.trace(contact_ham(d, g, wd * wg) @ covariance(t, target))
                        / 2
                    )
                    ff = readout(t, target)
                    mode_contact = (
                        -sum(z.conj() @ contact_feature(d, g) @ z for z in ff) / 2
                    )
                    assert_allclose(ct, mode_contact, rtol=1e-11, atol=2e-11)
                    if P and ((d == "beta") != (g == "beta")):
                        odd_errors.append(abs(fourier - wrong))
    print(
        "FULL_ORDERED_FOURIER_AND_CANONICAL_COVARIANCE",
        maximum,
        "WRONG_ODD_PARITY_ERROR",
        max(odd_errors),
        flush=True,
    )
    assert max(odd_errors) > 0.01

    # A separately evolved real prepared metric IVP, with nonzero final contact.
    from numpy.polynomial.legendre import leggauss
    from scipy.integrate import solve_ivp

    def sampling(t):
        xi = (t + 0.1) / 0.3
        return np.exp(-1 / (1 - xi * xi)) if abs(xi) < 1 else 0.0

    start = -0.5
    stop = 0.1
    P = 1
    md = ham("Q", np.cos(P * x)).real
    mg = ham("N", np.cos(P * x)).real
    cdg = contact_ham("Q", "N", np.cos(P * x) ** 2).real
    nodes, weights = leggauss(96)
    sampling_start = -0.4
    times = (stop - sampling_start) * nodes / 2 + (stop + sampling_start) / 2
    weights = (stop - sampling_start) * weights / 2
    max_fd = 0.0
    responses = []
    for target in (False, True):
        Cinitial = covariance(start, target)
        exact_memory = 0.0
        for ss, ww in zip(times, weights):
            us = expm(J @ M0 * (stop - ss))
            cs = covariance(ss, target)
            force = J @ mg @ cs + cs @ (J @ mg).T
            exact_memory += ww * sampling(ss) * (-np.trace(md @ us @ force @ us.T) / 2)
        contact = -sampling(stop) * np.trace(cdg @ covariance(stop, target)) / 2
        exact = exact_memory + contact

        def integrate(eps, Cinitial=Cinitial):
            def rhs(tt, flat):
                cc = flat.reshape(2 * L, 2 * L)
                aa = J @ (M0 + eps * sampling(tt) * mg)
                return (aa @ cc + cc @ aa.T).ravel()

            result = solve_ivp(
                rhs,
                (start, stop),
                Cinitial.ravel(),
                method="DOP853",
                rtol=2e-13,
                atol=2e-14,
            )
            assert result.success
            final = result.y[:, -1].reshape(2 * L, 2 * L)
            return -np.trace((md + eps * sampling(stop) * cdg) @ final) / 2

        eps = 1e-4
        finite = (integrate(eps) - integrate(-eps)) / (2 * eps)
        max_fd = max(max_fd, abs(finite - exact))
        assert_allclose(finite, exact, rtol=2e-7, atol=2e-8)
        assert abs(contact) > 0.01
        responses.append((exact, contact))
    print(
        "PREPARED_FULL_METRIC_IVP",
        max_fd,
        "FULL_STATE_DIFFERENCE",
        responses[1][0] - responses[0][0],
        "RETAINED_CONTACT_DIFFERENCE",
        responses[1][1] - responses[0][1],
        flush=True,
    )
    assert abs(responses[1][0] - responses[0][0]) > 0.0001
    assert abs(responses[1][1] - responses[0][1]) > 0.0001


def test_finite_cutoff_clock_contact_is_retained_until_the_full_limit():
    d = bounds.data()
    contact = d["retained_nonzero_finite_cutoff_clock_contact"]
    values = {}
    for sym in contact.free_symbols:
        name = str(sym)
        values[sym] = (
            s.Rational(1, 2)
            if name == "delta"
            else (0 if name == "fixed_state_selection_pressure_projected_K" else 1)
        )
    assert contact.subs(values) == s.Rational(3, 2)
    values[
        next(
            sym
            for sym in values
            if str(sym) == "fixed_state_selection_pressure_projected_K"
        )
    ] = 1
    assert contact.subs(values) == 0
    assert (
        "nonzero second-clock mean-tail contact is retained" in d["full_regulator_tail"]
    )
