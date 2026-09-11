"""Independent action, chart, continuum selection and scope controls."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_affine_vacuum_domain import family
from p8_vacuum_affine_metric_response import audit, bounds, bridge, chart, limits
from p8_vector_hadamard import preparation

PACKETS = {
    "bridge": bridge.data,
    "homogeneous": bridge.homogeneous,
    "chart": chart.data,
    "envelopes": chart.envelopes,
    "bounds": bounds.data,
    "limits": limits.data,
}


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_exact_residual(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(x == 0 for x in entries), name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_unsupported_scope_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize("name", list(PACKETS))
def test_each_quantitative_packet(name):
    assert all(PACKETS[name]().get("gates", {}).values())


@pytest.mark.parametrize("sector", ("T", "L"))
@pytest.mark.parametrize("momentum", (0, 1, 1000, 10**6))
def test_literal_temporal_constraint_and_both_metric_sources(sector, momentum):
    N, a, k, m, A, A0, V = s.symbols("N a k m A A0 V", positive=True)
    if sector == "T":
        L = a * V * V / (2 * N) - N * (k * k / a + a * m * m) * A * A / 2
    else:
        L = (
            a * (V - k * A0) ** 2 / (2 * N)
            + a**3 * m * m * A0 * A0 / (2 * N)
            - N * a * m * m * A * A / 2
        )
        A0_solution = s.solve(s.diff(L, A0), A0)[0]
        L = s.factor(L.subs(A0, A0_solution))
    g2 = s.diff(L, V, 2)
    w2 = s.factor(-s.diff(L, A, 2) / g2)
    n, v = s.symbols("n v", real=True)
    delta = lambda f: n * s.diff(f, N) + v * a * s.diff(f, a)
    r = s.factor((delta(w2) / (2 * w2)).subs(N, 1))
    rg = s.factor((delta(g2) / (2 * g2)).subs(N, 1))
    z = k * k / (k * k + a * a * m * m)
    expected_r = n - z * v
    expected_rg = (v - n) / 2 if sector == "T" else ((1 + 2 * z) * v - n) / 2
    assert s.factor(r - expected_r) == 0
    assert s.factor(rg - expected_rg) == 0
    assert g2.subs({N: 1, a: s.Rational(25, 16), k: momentum, m: 1000}) > 0
    alpha = s.symbols("old_mass_lapse_jet", real=True)
    assert s.diff(m * m + alpha * (N - 1), N).subs(N, 1) == alpha


@cache
def full_chart_functions():
    literal_R = s.lambdify((family.u, family.X), family.data()["R"], "mpmath")

    def omega(time, lapse):
        return -mp.log(literal_R(time, lapse**-2)) / 4

    return omega


@pytest.mark.parametrize(
    "time",
    (-s.Rational(1, 2), -s.Rational(1, 4), 0, s.Rational(1, 3), s.Rational(1, 2)),
)
def test_literal_full_analytic_chart_first_two_derivatives(time):
    with mp.workdps(80):
        t = mp.mpf(str(s.N(time, 85)))
        literal = lambda lapse: full_chart_functions()(t, lapse)
        val = literal(mp.mpf(1))
        w1 = mp.diff(literal, mp.mpf(1))
        w2 = mp.diff(literal, mp.mpf(1), 2)
        h = (1 + t * t) ** 3
        assert abs(val) < mp.mpf("1e-70")
        assert abs(w1 - 1 / (2 * h)) < mp.mpf("1e-70")
        assert abs(w2 + 3 / (2 * h) - 1 / (h * h)) < mp.mpf("1e-70")
        for lapse in (mp.mpf("0.99"), mp.mpf("1.01")):
            actual = full_chart_functions()(t, lapse)
            assert mp.isfinite(actual)


@pytest.mark.parametrize(
    "time", (-s.Rational(1, 2), 0, s.Rational(1, 4), s.Rational(1, 2))
)
def test_generic_nonstationary_density_second_chain_rule(time):
    N, V, Vhat = s.symbols("N V Vhat", real=True)
    n, vhat = s.symbols("n vhat", real=True)
    h = (1 + time * time) ** 3
    omega = -s.log(1 + (N**-2 - 1) / h) / 4
    F = 2 + 3 * (N - 1) + 5 * V + 7 * (N - 1) ** 2 + 11 * (N - 1) * V + 13 * V * V
    rho = -s.exp(-3 * V) * s.diff(F, N)
    P = s.exp(-3 * V) * s.diff(F, V) / (3 * N)
    w1 = s.diff(omega, N).subs(N, 1)
    qphysical = s.Matrix([n, vhat + w1 * n])
    stresses = s.Matrix([rho, P])
    delta = stresses.jacobian((N, V)).subs({N: 1, V: 0}) * qphysical
    p = chart.data()
    replacements = {chart.u: time}
    values = {
        "rho": rho.subs({N: 1, V: 0}),
        "pressure": P.subs({N: 1, V: 0}),
        "delta_rho": delta[0],
        "delta_pressure": delta[1],
        "n": n,
        "vhat": vhat,
    }
    replacements.update(
        {
            x: values[str(x)]
            for x in p["clock_retarded_response"].free_symbols
            if str(x) in values
        }
    )
    predicted = p["clock_retarded_response"].subs(replacements, simultaneous=True)
    pulled = F.subs(V, Vhat + omega)
    exact = s.hessian(pulled, (N, Vhat)).subs({N: 1, Vhat: 0}) * s.Matrix([n, vhat])
    assert (exact - predicted).applyfunc(s.factor) == s.zeros(2, 1)
    missing = (
        p["nonzero_second_metric_map_contact"]
        .subs(
            {
                x: values[str(x)]
                for x in p["nonzero_second_metric_map_contact"].free_symbols
                if str(x) in values
            }
        )
        .subs(chart.u, time)
    )
    assert missing[0, 0] != 0


@pytest.mark.parametrize("order", range(11))
def test_leibniz_chart_norm_on_whole_declared_interval(order):
    u = chart.u
    f = s.Function("n")(u)
    w1 = chart.data()["omega_N"]
    actual = s.diff(w1 * f, u, order)
    expected = sum(
        s.binomial(order, j) * s.diff(w1, u, order - j) * s.diff(f, u, j)
        for j in range(order + 1)
    )
    assert s.expand(actual - expected) == 0
    bound = chart.envelopes()["coefficient_derivative_envelopes"][order][
        "absolute_upper"
    ]
    for j in range(-8, 9):
        assert abs(s.diff(w1, u, order).subs(u, s.Rational(j, 16))) <= bound


@pytest.mark.parametrize(
    "r", (s.Rational(64, 125), s.Rational(2, 3), s.Rational(3, 4), s.S.One)
)
def test_sharp_pointwise_chart_contact_envelopes(r):
    w1 = r / 2
    w2 = r * r - 3 * r / 2
    assert 0 < w1 <= s.Rational(1, 2)
    assert abs(w2) <= s.Rational(9, 16)
    assert abs(w1 + w2) <= s.Rational(1, 4)
    for rho, P in ((1, 1), (2, 3), (-2, 3), (2, -3)):
        C = s.Matrix(
            [
                [3 * P * (w1 + 3 * w1 * w1 + w2) - 3 * rho * w1, -3 * rho + 9 * w1 * P],
                [3 * P + 9 * w1 * P, 9 * P],
            ]
        )
        assert sum(abs(x) for x in C[0, :]) <= s.Rational(9, 2) * abs(rho) + s.Rational(
            15, 2
        ) * abs(P)
        assert sum(abs(x) for x in C[1, :]) <= s.Rational(33, 2) * abs(P)


@pytest.mark.parametrize("kind", ("transverse", "longitudinal"))
@pytest.mark.parametrize("nu", (1000, 3000))
def test_same_all_order_initial_covariance_not_reset(kind, nu):
    d = preparation.initial_data(kind, 1000, nu)
    assert d["initial_u"] == -s.Rational(1, 2)
    assert d["all_order_frequency"] > 0
    if nu == 3000:
        assert d["active_higher_derivative_orders"] == [6]
        assert d["all_order_frequency"] != d["frozen_frequency"]
    assert "no independent state variation" in bridge.data()["state"]


def test_only_pinned_physical_vector_response_not_old_total_profile():
    report = bounds.read_pinned(bounds.old_response.REPORT, bounds.RESPONSE_SHA)
    d = bounds.data()
    vector = report["new_complete_physical_vector_response_bounds"]
    old_total = report[
        "new_complete_background_cancelled_response_bounds_in_both_charts"
    ]
    for name, value in d[
        "actual_conditional_physical_stress_C10_to_C0_over_kappa"
    ].items():
        assert value == s.Rational(
            vector["complete_metric_response_C10_to_C0_component_bounds"][name]
        )
        assert value < s.Rational(
            old_total["background_cancelled_physical_C10_to_C0_component_bounds"][name]
        )
    assert d["clock_current_C10_to_C0_over_kappa"] != s.Rational(
        old_total["joint_background_cancelled_clock_C10_to_C0_upper"]
    )
    with pytest.raises(ValueError, match="changed"):
        bounds.read_pinned(bounds.old_response.REPORT, "0" * 64)


def test_nonzero_cosmological_contact_and_no_false_tail_verdict():
    d = limits.data()
    assert d["nonzero_second_map_cosmological_contact_at_bounce"] != 0
    assert d["uncancelled_cosmological_density_clock_Hessian_at_bounce"] != s.zeros(2)
    assert "not a total-state asymptotic estimate" in d["finite_time_boundary"]
    assert "does NOT prove" in d["prepared_norm_diagnostic"]


def test_no_mean_response_versus_noise_or_coupled_inverse_promotion():
    assert "no coupled inverse" in audit.observable()["not_inferred"]
    assert "no scalar profile is added" in audit.observable()["parent"]
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 35
    assert len(audit.gates()) == 27
    assert audit.scalar_entry_count() == 74
