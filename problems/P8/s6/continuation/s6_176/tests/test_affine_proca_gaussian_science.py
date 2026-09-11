"""Independent constrained-action, metric, covariance and source controls."""

import copy

import pytest
import sympy as s
from p8_constant_proca import quantum as ordinary
from p8_proca_stress import estimates, subtraction
from p8_vacuum_affine_proca_gaussian import (
    audit,
    bridge,
    gaussian,
    state,
    stress,
    verify,
)
from p8_vector_hadamard import cutoffs, preparation
from p8_vector_state import wkb


@pytest.mark.parametrize("name,value", tuple(audit.residuals().items()))
def test_every_exact_residual(name, value):
    if isinstance(value, s.MatrixBase):
        assert value == s.zeros(*value.shape), name
    else:
        assert s.simplify(value) == 0, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_strict_rejected_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("kind", ("transverse", "longitudinal"))
@pytest.mark.parametrize(
    "lapse,scale,momentum",
    (
        (s.Rational(1, 2), 1, 0),
        (1, 1, 1),
        (2, s.Rational(3, 2), 1000),
        (s.Rational(3, 2), 2, 10**6),
    ),
)
def test_uneliminated_action_metric_variation(kind, lapse, scale, momentum):
    N, a, k, m = s.symbols("N a k m", positive=True)
    A, A0, velocity = s.symbols("A A0 velocity", real=True)
    if kind == "transverse":
        L = a * velocity**2 / (2 * N) - N * (k * k / a + a * m * m) * A * A / 2
        solution = 0
    else:
        L = (
            a * (velocity - k * A0) ** 2 / (2 * N)
            + a**3 * m * m * A0 * A0 / (2 * N)
            - N * a * m * m * A * A / 2
        )
        solution = s.solve(s.diff(L, A0), A0)[0]
    # Metric derivatives precede BOTH the temporal substitution and canonical map.
    rho = s.factor((-s.diff(L, N) / a**3).subs(A0, solution))
    pressure = s.factor((s.diff(L, a) / (3 * N * a * a)).subs(A0, solution))
    q = k * k / a**2
    g = s.sqrt(a if kind == "transverse" else a * m * m / (q + m * m))
    v, p = s.symbols("v p", real=True)
    mapping = {A: v / g, velocity: N * p / g}
    rho = s.factor(rho.subs(mapping))
    pressure = s.factor(pressure.subs(mapping))
    expected_rho = (p * p + (q + m * m) * v * v) / (2 * a**3)
    expected_p = (
        (p * p + (q - m * m) * v * v) / (6 * a**3)
        if kind == "transverse"
        else ((m * m + 3 * q) / (m * m + q) * p * p - (m * m + q) * v * v) / (6 * a**3)
    )
    assert s.factor(rho - expected_rho) == 0
    assert s.factor(pressure - expected_p) == 0
    values = {
        N: lapse,
        a: scale,
        k: momentum,
        m: 1000,
        v: s.Rational(2, 3),
        p: -s.Rational(4, 5),
    }
    assert rho.subs(values) > 0
    # Treating v,p as fixed too early omits physical metric dependence.
    prematurely_mapped = L.subs(A0, solution).subs(mapping)
    wrong = -s.diff(prematurely_mapped, N) / a**3
    assert s.factor(wrong - rho) != 0


@pytest.mark.parametrize("kind", ("transverse", "longitudinal"))
@pytest.mark.parametrize("k_value", (0, 1, 1000, 10**6))
def test_original_coordinate_equation_to_canonical_frequency(kind, k_value):
    t = wkb.u
    a = (1 + t * t) ** 2
    m = s.Integer(1000)
    q = s.Integer(k_value) ** 2 / a**2
    g2 = a if kind == "transverse" else a * m * m / (q + m * m)
    g = s.sqrt(g2)
    f = s.Function("v")(t)
    A = f / g
    # Coordinate Euler equation of the literal reduced action at N=1.
    coordinate = s.diff(g2 * s.diff(A, t), t) + g2 * (q + m * m) * A
    normalized = s.simplify(coordinate / g)
    U = s.simplify(s.diff(g, t, 2) / g)
    expected = s.diff(f, t, 2) + (q + m * m - U) * f
    assert s.simplify(normalized - expected) == 0
    z = q / (q + m * m)
    assert s.factor(U - wkb.frequency(kind)["U"].subs(wkb.z, z)) == 0


@pytest.mark.parametrize(
    "omega,d,mean",
    (
        (1, 0, (0, 0)),
        (2, -3, (1, 2)),
        (1000, s.Rational(7, 3), (-2, 1)),
        (10**6, -1000, (s.Rational(1, 3), -5)),
    ),
)
def test_literal_positive_coherent_Gram_and_connected_CCR(omega, d, mean):
    v = 1 / s.sqrt(2 * omega)
    p = (-d - s.I * omega) * v
    phase = s.Matrix([v, p])
    C = phase * phase.conjugate().T
    mean = s.Matrix(mean)
    full = C + mean * mean.T
    Omega = s.Matrix([[0, 1], [-1, 0]])
    assert C - C.T == s.I * Omega
    assert full - full.T == s.I * Omega
    assert full - mean * mean.T == C
    for test in (s.Matrix([1, 0]), s.Matrix([1, s.I]), s.Matrix([2 - s.I, 3 + s.I])):
        assert s.simplify((test.conjugate().T * full * test)[0]) >= 0
    CW = C / (bridge.KAPPA * bridge.ZETA)
    assert (CW - CW.T) * bridge.KAPPA * bridge.ZETA == s.I * Omega
    assert CW != C / bridge.KAPPA


def test_literal_tensor_action_normalization_and_missing_zeta_control():
    eta = s.diag(1, -1, -1, -1)
    A = s.Matrix([s.Rational(1, 2), -2, s.Rational(3, 5), s.Rational(7, 11)])
    S = s.Matrix([s.Rational(2, 7), s.Rational(1, 3), -1, s.Rational(4, 9)])
    F = s.Matrix(4, 4, lambda i, j: s.Rational(i - j, 13))
    scale = s.sqrt(bridge.KAPPA * bridge.ZETA)
    W, FW = A / scale, F / scale
    curl = lambda f: sum(
        eta[i, i] * eta[j, j] * f[i, j] ** 2 for i in range(4) for j in range(4)
    )
    old = (
        -bridge.KAPPA * bridge.ZETA * curl(FW) / 4
        + bridge.KAPPA * ((W - S).T * eta * (W - S))[0] / 2
    )
    canonical = (
        -curl(F) / 4
        + bridge.MASS**2 * (A.T * eta * A)[0] / 2
        - s.sqrt(bridge.KAPPA / bridge.ZETA) * (A.T * eta * S)[0]
        + bridge.KAPPA * (S.T * eta * S)[0] / 2
    )
    assert old == canonical
    assert bridge.KAPPA * bridge.ZETA / bridge.KAPPA != 1
    assert 1 / bridge.ZETA == 10**6


@pytest.mark.parametrize("kind", ("transverse", "longitudinal"))
@pytest.mark.parametrize("nu", (1000, 3000, 4000))
def test_exact_selected_Cauchy_preparation_not_finite_order_relabel(kind, nu):
    initial = preparation.initial_data(kind, 1000, nu)
    assert initial["m0_tau"] == 1000 and initial["initial_u"] == -s.Rational(1, 2)
    assert initial["all_order_frequency"] > 0
    assert initial["all_order_half_log_rate"] == initial[
        "all_order_frequency_slope"
    ] / (2 * initial["all_order_frequency"])
    if nu == 1000:
        assert initial["active_higher_derivative_orders"] == []
        assert initial["all_order_frequency"] == initial["frozen_frequency"]
    else:
        assert initial["active_higher_derivative_orders"] == [6]
        assert initial["all_order_frequency"] != initial["frozen_frequency"]
    assert cutoffs.threshold(kind, 3, 1000) == 2000
    assert cutoffs.threshold(kind, 4, 1000) >= 4000
    assert "S6.55" in state.data()["reference_state"]
    assert "finite-order comparison" in state.data()["reference_state"]


@pytest.mark.parametrize("order", range(6))
def test_all_physical_stress_bounds_and_decompositions(order):
    packet = stress.data()
    row = packet["bounds_by_time_derivative_order"][order]
    old = estimates.physical_bounds(s.sqrt(bridge.KAPPA))
    assert (
        row["normalized_energy_decomposition"]
        == old["new_energy_derivative_bounds"][order]
    )
    assert (
        row["normalized_pressure_decomposition"]
        == old["identical_pressure_derivative_bounds"][order]
    )
    for name in ("energy", "pressure"):
        parts = row[f"normalized_{name}_decomposition"]
        assert parts["total"] == sum(v for key, v in parts.items() if key != "total")
        value = row[f"absolute_{name}_derivative_upper"]
        assert value == bridge.KAPPA * parts["total"] < 10**30
        assert value / bridge.KAPPA < s.Rational(1, 10**770)
        assert value != bridge.KAPPA * bridge.ZETA * parts["total"]


@pytest.mark.parametrize("order", range(3))
def test_independent_covariant_density_lapse_and_scale_variation(order):
    a, ad, add, addd, a4 = s.symbols("a ad add addd a4", positive=True)
    H, Hd, Hdd, Hddd = s.symbols("H Hd Hdd Hddd", real=True)
    N, Nd = s.symbols("N Nd", real=True)
    acc = (Hd + H * H) / N**2 - H * Nd / N**3
    exp2 = H * H / N**2
    R = 6 * acc + 6 * exp2
    Ricci2 = 9 * acc**2 + 3 * (acc + 2 * exp2) ** 2
    Riemann2 = 12 * acc**2 + 12 * exp2**2
    scal = R * R / 72 - Ricci2 / 180 + Riemann2 / 180
    inv = (s.Rational(5, 2), s.Rational(5, 3) * R, -4 * scal)[order]
    first = s.diff(N * inv, N).subs({N: 1, Nd: 0})
    rate = s.diff(N * inv, Nd).subs({N: 1, Nd: 0})
    dt = lambda f: s.diff(f, H) * Hd + s.diff(f, Hd) * Hdd + s.diff(f, Hdd) * Hddd
    rho = s.factor(-first + dt(rate) + 3 * H * rate)
    # Independent spatial Euler variation, including both scale derivatives.
    inv_scale = inv.subs(
        {N: 1, Nd: 0, H: ad / a, Hd: add / a - ad**2 / a**2}, simultaneous=True
    )
    density = s.factor(a**3 * inv_scale)
    variables = (a, ad, add, addd)
    rates = (ad, add, addd, a4)
    D = lambda f: sum(s.diff(f, v) * r for v, r in zip(variables, rates))
    euler = s.diff(density, a) - D(s.diff(density, ad)) + D(D(s.diff(density, add)))
    pressure = s.factor(euler / (3 * a * a))
    mapping = {
        ad: a * H,
        add: a * (Hd + H * H),
        addd: a * (Hdd + 3 * H * Hd + H**3),
        a4: a * (Hddd + 4 * H * Hdd + 3 * Hd * Hd + 6 * H * H * Hd + H**4),
    }
    pressure = s.factor(pressure.subs(mapping, simultaneous=True))
    u = wkb.u
    actual_H = 4 * u / (1 + u * u)
    point = {
        H: actual_H,
        Hd: s.diff(actual_H, u),
        Hdd: s.diff(actual_H, u, 2),
        Hddd: s.diff(actual_H, u, 3),
    }
    row = stress.local()["actual_local_coefficients"][order]
    assert s.factor(rho.subs(point) - row["energy"]) == 0
    assert s.factor(pressure.subs(point) - row["pressure"]) == 0
    assert s.factor(dt(rho) + 3 * H * (rho + pressure)) == 0


def test_nonzero_spatial_clock_source_and_cubic_force_not_discarded():
    eps, x = s.symbols("epsilon x", real=True)
    f = s.Function("compact_profile")(x)
    # At the bounce use u=t, N=1+epsilon*f and
    # a=a_CD*exp(epsilon*t*f) near t=0, smoothly cut off away from it.
    # Actual R_X(0,1)=1, H=R_u=uHu=0, Box u/X=3epsilon*f.
    X = (1 + eps * f) ** -2
    source = 3 * eps * f * (X - 1)
    second = s.diff(source, eps, 2).subs(eps, 0) / 2
    assert s.factor(second + 6 * f * f) == 0
    assert s.diff(second, x) != 0
    # The source is temporal but spatially varying, hence generally not closed.
    assert s.diff(second, x) == -12 * f * s.diff(f, x)
    data = gaussian.data()
    assert data["nonzero_second_source_operator_contact"] != 0
    assert "third order" in data["clock_bridge"]
    assert "NOT added" in data["no_profile_transfer"]


def test_constant_mass_counterterms_not_old_clock_only_value_substitution():
    N, alpha, beta, G = s.symbols("N alpha beta G", real=True)
    old_deviation = alpha * (N - 1) + beta * (N - 1) ** 2 / 2
    assert s.diff(old_deviation * G, N).subs(N, 1) == alpha * G
    assert s.diff(s.S.Zero * G, N) == 0
    new = stress.prescription()
    assert "no X denominator" in new["constant_mass_domain"]
    assert new["checks"]["ordinary_pole_has_no_clock_norm"] == 0
    assert ordinary.boundary()["actual_old_first_mass_lapse_jets"]["temporal"] != 0


def test_old_low_reference_failure_control_still_present():
    failure = subtraction.reference_tail("longitudinal", 4, 2)
    assert not failure["integrable_tail_certified"]
    assert failure["low_coefficient_residuals"][4].subs({wkb.u: 0, wkb.z: 1}) == 23808


@pytest.mark.parametrize("t", (-s.Rational(1, 2), 0, s.Rational(1, 2)))
def test_exact_allowed_closed_stress_interval(t):
    assert stress.require_scope(t, 5, bridge.ZETA, bridge.KAPPA) == (t, 5)


def test_exact_payloads_state_family_and_open_frontier():
    for d in (
        bridge.data(),
        bridge.clock(),
        state.data(),
        gaussian.data(),
        stress.prescription(),
        stress.local(),
        stress.data(),
    ):
        verify.serialize(verify.payload(d))
    assert all(v is True for v in audit.gates().values())
    assert "affine Proca algebra" in state.data()["family"]
    rows = copy.deepcopy(audit.matching())
    rows[-1]["status"] = "COMPLETE"
    with pytest.raises(ValueError):
        audit.validate_scope(audit.frontier(), rows)
    with pytest.raises(ValueError):
        verify.serialize(s.Float("0.1"))
