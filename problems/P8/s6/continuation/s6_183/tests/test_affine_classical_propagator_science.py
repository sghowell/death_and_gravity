"""Independent Hamiltonian, continuum-envelope, derivative and ODE fixtures."""

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_classical_propagator import audit, perturbation
from p8_vacuum_affine_classical_propagator import hamiltonian as ham
from p8_vacuum_affine_classical_propagator import jet_bounds as jets
from p8_vacuum_affine_classical_propagator import real_bounds as real
from scipy.integrate import solve_ivp


def mpq(value):
    value = s.Rational(value)
    return mp.mpf(int(value.p)) / int(value.q)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_all_exact_residuals(name):
    value = audit.residuals()[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    ), name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_unsupported_scope_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("name", list(audit.packets()))
def test_every_continuous_packet(name):
    assert all(audit.packets()[name]["gates"].values())


@pytest.mark.parametrize("t", [s.Rational(j, 20) for j in range(-10, 11)])
def test_actual_rational_coefficients_and_weighted_norms(t):
    d = ham.data()
    A = d["weighted_phase_matrix"].subs(ham.u, t)
    B = d["weighted_force_matrix"].subs(ham.u, t)
    row = lambda M: max(sum(abs(M[i, j]) for j in range(M.cols)) for i in range(M.rows))
    assert row(A) < 7 and row(B) <= 48
    assert ham.J.subs(ham.u, t) > s.Rational(1, 40)
    assert abs(ham.theta.subs(ham.u, t)) <= abs(ham.H.subs(ham.u, t))
    nrow = d["actual_lapse_state_row"].subs(ham.u, t)
    grows = d["actual_lapse_force_row"].subs(ham.u, t)
    assert row(nrow) <= s.Rational(7, 2)
    assert row(grows) < 30


@pytest.mark.parametrize(
    "t", [-s.Rational(1, 2), -s.Rational(1, 4), 0, s.Rational(1, 4), s.Rational(1, 2)]
)
def test_independent_literal_complex_disc_fixtures(t):
    with mp.workdps(85):
        tt = mpq(t)
        pp = s.lambdify(ham.u, ham.P, "mpmath")
        for i in range(16):
            z = tt + mpq(jets.RADIUS) * mp.exp(2j * mp.pi * i / 16)
            D = 1 + z * z
            delta = 1 / (2 * D**3)
            ell = 1 / (10 * D**6)
            theta = 4 * z / D - z / D**4
            L = 3 * ell**2 * (3 * delta - 1)
            J = pp(z) / (800 * D**18)
            H = 4 * z / D
            A = [
                [theta * L / (2 * J), (theta**2 / (2 * J) - mp.mpf(1) / 6) / 10],
                [10 * (-(L**2) / (2 * J) - 9 * ell**2), -3 * H - theta * L / (2 * J)],
            ]
            B = [
                [theta / (2 * J), delta * theta / (2 * J)],
                [-10 * L / (2 * J), -10 * (1 + delta * L / (2 * J))],
            ]
            assert abs(pp(z) - pp(tt)) < mpq(
                jets.data()["polynomial_disc_variation_upper"]
            )
            assert abs(pp(z)) > 1000 and abs(1 / (2 * J)) < 50
            assert max(sum(abs(x) for x in row) for row in A) < 200000
            assert max(sum(abs(x) for x in row) for row in B) < 30000


@pytest.mark.parametrize("order", range(11))
def test_independent_high_precision_coefficient_derivatives(order):
    # Literal scalar formulas, differentiated without the symbolic report's
    # Cauchy/recurrence construction. Fixtures do not replace that proof.
    with mp.workdps(100):
        for tt in (mp.mpf("-0.5"), mp.mpf(0), mp.mpf("0.5")):
            f = lambda t: 1 / (2 * (1 + t * t) ** 3)
            actual = abs(mp.diff(f, tt, order))
            assert (
                actual < mp.mpf("1.5") * mp.factorial(order) / mpq(jets.RADIUS) ** order
            )
            assert jets.data()["physical_lapse_derivative_upper"][order] < jets.C10
            assert jets.data()["physical_log_scale_derivative_upper"][order] < jets.C10


@pytest.mark.parametrize("order", range(10))
def test_scaled_induction_by_independent_binomial_sum(order):
    r = jets.RADIUS
    M, B = s.Integer(200000), s.Integer(30000)
    ordinary = sum(
        s.binomial(order, k)
        * M
        * s.factorial(k)
        / r**k
        * real.Y0
        * s.factorial(order - k)
        / r ** (order - k)
        for k in range(order + 1)
    )
    ordinary += sum(
        s.binomial(order, k) * B * s.factorial(k) / r**k for k in range(order + 1)
    )
    scaled = ordinary * r ** (order + 1) / s.factorial(order + 1)
    direct = (
        r
        / (order + 1)
        * (
            M * (order + 1) * real.Y0
            + B * sum(r**i / s.factorial(i) for i in range(order + 1))
        )
    )
    assert s.factor(scaled - direct) == 0
    assert scaled < real.Y0


@pytest.mark.parametrize(
    "t,n,v,p,gn,gv", [(s.Rational(1, 4), 2, -3, 5, 7, -11), (0, -1, 2, -3, 4, 5)]
)
def test_literal_original_charge_Routhian_and_source_map(t, n, v, p, gn, gv):
    q = ham.classical.quadratic()
    theta, J, w, ell, delta = [
        expr.subs(ham.u, t) for expr in (q["theta"], q["J"], q["w"], ham.ell, ham.delta)
    ]
    vd, sd, nn, vv = s.symbols("vd sd n v")
    raw = (
        -3 * vd**2
        + (J + w * w / 2 - 3 * theta**2) * nn**2
        + 6 * theta * nn * vd
        + sd**2 / 2
        + w * nn * sd
        + 3 * ell * vv * sd
    )
    rate = -w * nn - 3 * ell * vv
    effective = s.expand(raw.subs(sd, rate))
    expected = (
        -3 * vd**2
        + 6 * theta * nn * vd
        + (J - 3 * theta**2) * nn**2
        - 3 * ell * w * nn * vv
        - s.Rational(9, 2) * ell**2 * vv**2
    )
    # Routh reduction is done in the un-integrated form with cyclic sigma;
    # subtracting its zero fixed perturbative momentum is legitimate here.
    assert s.factor(effective - expected) == 0
    force = gn * nn + gv * (vv + delta * nn)
    lapse = (theta * p + 3 * ell * w * v + gn + delta * gv) / (2 * J)
    assert (
        s.factor(
            s.diff(expected - force, nn).subs(
                {nn: lapse, vv: v, vd: theta * lapse - p / 6}
            )
        )
        == 0
    )
    wrong = (theta * p + 3 * ell * w * v + gn) / (2 * J)
    assert (
        s.factor(
            s.diff(expected - force, nn).subs(
                {nn: wrong, vv: v, vd: theta * wrong - p / 6}
            )
        )
        != 0
    )


def force(t):
    x = (t + 0.4) / 0.3
    b = float(np.exp(4 - 1 / (x * (1 - x)))) if 0 < x < 1 else 0.0
    return np.array([b, -0.4 * b])


@pytest.fixture(scope="module")
def solved_reference():
    A = s.lambdify(ham.u, ham.data()["actual_phase_matrix"], "numpy")
    B = s.lambdify(ham.u, ham.data()["actual_physical_force_matrix"], "numpy")

    def rhs(t, y):
        return A(t) @ y + B(t) @ force(t)

    sol = solve_ivp(
        rhs,
        (-0.5, 0.5),
        [0.0, 0.0],
        method="DOP853",
        rtol=2e-11,
        atol=1e-13,
        max_step=0.01,
        dense_output=True,
    )
    assert sol.success
    return sol, rhs


@pytest.mark.parametrize("t", [-0.35, -0.25, -0.15, 0.0, 0.25, 0.45])
def test_forced_ODE_against_original_Euler_currents(solved_reference, t):
    sol, _ = solved_reference
    v, p = sol.sol(t)
    step = 2e-5
    derivative = (sol.sol(t + step) - sol.sol(t - step)) / (2 * step)
    theta, J, c, ell, delta, H = [
        float(expr.subs(ham.u, t))
        for expr in (ham.theta, ham.J, ham.L, ham.ell, ham.delta, ham.H)
    ]
    gn, gv = force(t)
    n = (theta * p + c * v + gn + delta * gv) / (2 * J)
    En = (
        6 * theta * derivative[0] + 2 * (J - 3 * theta**2) * n - c * v - gn - delta * gv
    )
    Ev = -c * n - 9 * ell**2 * v - derivative[1] - 3 * H * p - gv
    assert abs(En) < 2e-6 and abs(Ev) < 2e-6


def test_prepared_causal_response_not_zero_at_final_endpoint(solved_reference):
    sol, _ = solved_reference
    assert max(abs(sol.sol(-0.45))) < 1e-15
    assert max(abs(sol.sol(0.5))) > 1e-5
    for t in np.linspace(-0.5, 0.5, 101):
        v, p = sol.sol(t)
        theta, J, c, delta = [
            float(expr.subs(ham.u, t)) for expr in (ham.theta, ham.J, ham.L, ham.delta)
        ]
        gn, gv = force(t)
        n = (theta * p + c * v + gn + delta * gv) / (2 * J)
        assert max(abs(n), abs(v + delta * n)) < 53000
    assert "need not vanish" in perturbation.data()["future_extension"]


def test_first_coefficient_is_not_a_finite_coupling_error():
    d = perturbation.data()
    assert real.C0 * perturbation.Q_UPPER * jets.C10 < perturbation.FIRST_UPPER
    assert "No differentiable exact solution branch" in d["not_an_error_bound"]
    assert "BOTH" in d["marker_scope"]
    assert (
        s.diff(
            d["actual_complete_stationary_physical_response"][1],
            s.Symbol("delta_p", real=True),
        )
        == 3
    )
