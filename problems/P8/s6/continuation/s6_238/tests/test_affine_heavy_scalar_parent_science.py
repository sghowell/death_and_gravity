"""Independent full-potential, metric-variation, domain and clock checks."""

import itertools

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_heavy_scalar_loop_coefficients import audit as previous
from p8_vacuum_affine_heavy_scalar_parent import audit, clock, family, heavy


def mpq(value):
    value = s.Rational(value)
    return mp.mpf(int(value.p)) / int(value.q)


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_every_exact_entry(name, value):
    assert s.cancel(value) == 0, name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_every_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_literal_entire_new_functions_not_only_a_germ():
    d = family.coefficients()
    h = family.switch()
    # These pin identical complete defining forms, without polynomially expanding
    # the huge rational switch. Generic rational Ia identities are checked separately.
    assert (
        d["R"]
        == family.original.data()["R"] + family.GAMMA * family.K0 * family.X**2 * h
    )
    assert (
        d["F"]
        == family.profile.data()["new_complete_scalar_coefficient"] + d["delta_F"]
    )
    rx = s.diff(d["R"], family.X)
    assert d["A3"] == rx / family.X
    assert d["A4"] == -rx / family.X - 7 * rx**2 / (4 * d["R"])
    assert d["A5"] == rx**2 / (d["R"] * family.X)
    assert (
        d["normalized_heavy_source"]
        == family.G * s.sqrt(family.K0) * family.u**2 * h / 2
    )


def test_full_source_stress_from_independent_lapse_variation():
    lapse = s.Symbol("lapse", positive=True)
    Y, phi, H, hdot, n, g = s.symbols("Y Phi H Hdot n g", real=True)
    hfun = s.Function("h")
    density = (
        hdot**2 / (2 * lapse**2)
        - n * H**2 / 2
        + g * H * phi**2 * hfun(Y / lapse**2) / 2
    )
    rho = -s.diff(lapse * density, lapse).subs(lapse, 1).doit()
    h, hy = s.symbols("h h_Y", real=True)
    rho = rho.subs({hfun(Y): h, s.diff(hfun(Y), Y): hy}, simultaneous=True)
    J = g * phi**2 * h / 2
    H0 = J / n
    delta = g * Y * h / n**2
    actual = s.cancel(rho.subs(H, H0 - delta) - rho.subs(H, H0))
    wanted = g * g * Y**2 * h * h / (2 * n**3) - g * g * Y**2 * phi**2 * h * hy / n**2
    assert s.cancel(actual - wanted) == 0
    # This regression detects the exact missing-Y error found privately.
    wrong = g * g * Y**2 * h * h / (2 * n**3) - g * g * Y * phi**2 * h * hy / n**2
    assert s.cancel(actual - wrong) != 0
    packet = heavy.data()["full_source_metric_energy_difference"]
    replacement = {
        z: {"positive_Y": Y, "Phi": phi, "g": g, "n": n, "h": h, "h_Y": hy}[str(z)]
        for z in packet.free_symbols
    }
    assert s.cancel(packet.subs(replacement, simultaneous=True) - wanted) == 0


@pytest.mark.parametrize(
    "n,g,Y,t,A,k",
    [
        (13, 2, 3, s.Rational(1, 5), 3, 7),
        (29, 1, 2, -2, 11, 13),
        (7, 3, s.Rational(1, 4), 3, 2, 5),
    ],
)
def test_particular_heavy_solution_against_direct_time_derivatives(n, g, Y, t, A, k):
    tt = s.Symbol("t", real=True)
    x = Y / s.Integer(k)
    h = (1 - x) ** 8 * s.exp(-A * x * x)
    phi = s.sqrt(Y) * tt
    J = g * phi**2 * h / 2
    H = J / n - g * Y * h / n**2
    assert s.cancel(s.diff(H, tt, 2) + n * H - J) == 0
    assert s.diff(H - J / n, tt) == 0
    assert s.diff(H, tt).subs(tt, t) != 0


@pytest.mark.parametrize(
    "x", [s.Rational(1, 100), s.Rational(1, 2), s.Rational(7, 8), s.Rational(9, 8)]
)
@pytest.mark.parametrize("A", [2, 17, 1000])
def test_direct_source_metric_stress_numerical_fixture(x, A):
    with mp.workdps(70):
        k, n, g, phi = mp.mpf(7), mp.mpf(29), mp.mpf("0.125"), mp.mpf(3)
        Y = k * mpq(x)
        hh = lambda yy: (1 - yy / k) ** 8 * mp.exp(-A * (yy / k) ** 2)
        h, hY = hh(Y), mp.diff(hh, Y)
        J = g * phi * phi * h / 2
        H0 = J / n
        Hp = H0 - g * Y * h / n**2
        hd = g * phi * mp.sqrt(Y) * h / n
        # Vary the metric before substituting either field configuration.
        rho = lambda H: (
            -mp.diff(
                lambda N: (
                    N
                    * (
                        hd**2 / (2 * N * N)
                        - n * H * H / 2
                        + g * H * phi * phi * hh(Y / (N * N)) / 2
                    )
                ),
                1,
            )
        )
        actual = rho(Hp) - rho(H0)
        expected = (
            g * g * Y * Y * h * h / (2 * n**3)
            - g * g * Y * Y * phi * phi * h * hY / n**2
        )
        scale = max(abs(expected), mp.mpf("1e-2000"))
        assert abs(actual - expected) < mp.mpf("1e-55") * scale


def test_actual_unrestricted_small_offset_large_metric_stress_counterexample():
    with mp.workdps(300):
        k, A, n, g = map(mpq, (family.K0, family.LOCALIZER, family.MASS2, family.G))
        x = 1 / mp.sqrt(A)
        h = (1 - x) ** 8 * mp.exp(-A * x * x)
        hx = -8 * (1 - x) ** 7 * mp.exp(-A * x * x) - 2 * A * x * h
        Y = k * x
        actual_relative = g * g * Y * h * h / (2 * n**3) - g * g * k * x * h * hx / n**2
        proxy = g * g * k / (mp.sqrt(A) * n**3)
        assert proxy < mp.mpf("1e-8")
        assert actual_relative > mp.mpf("1e395")
        assert actual_relative > g * g * k / (9 * n * n)
    with pytest.raises(ValueError, match="corridor"):
        audit.require_corridor(s.Integer(10) ** 400, s.Rational(1, 10**210))


@pytest.mark.parametrize(
    "phi,x",
    [
        (0, 0),
        (10**190, 0),
        (-(10**190), s.Rational(6, 5)),
        (10**400, s.Rational(1, 2)),
        (-(10**400), 1),
    ],
)
def test_corridor_admissible_anchors(phi, x):
    assert audit.require_corridor(phi, x) == (s.Rational(phi), s.Rational(x))


def test_actual_complete_corridor_bounds_not_only_offset():
    d = heavy.data()
    eta = d["source_response_proxy_squared_upper"]
    low = d["source_metric_stress_low_chart_upper"]
    high = d["source_metric_stress_high_chart_upper"]
    assert 0 < low < s.Rational(1, 10**21)
    assert 0 < high < s.Rational(1, 10**2500)
    assert eta / 2 + max(low, high) < s.Rational(1, 10**8)
    assert d["source_metric_stress_corridor_upper"] == max(low, high)


@pytest.mark.parametrize("X", [s.Rational(9, 8), s.Rational(7, 6), s.Rational(6, 5)])
def test_complete_upper_edge_margin_independent_stable_evaluation(X):
    with mp.workdps(100):
        x = mpq(X)
        z = ((x - 1) / x) ** 1024
        one_minus_B = z / (1 + z) * (1 - 1024 * x * x * mp.exp(-1024 * x * x))
        margin = x - 1
        gap = (x - 1) * one_minus_B
        assert gap > mpq(clock.data()["upper_edge_original_margin"])
        log_delta = (
            mp.log(1024)
            + 2 * mp.log(x)
            + 8 * mp.log(abs(1 - x))
            - mpq(family.LOCALIZER) * x * x
        )
        assert log_delta < mp.log(gap / 2)
        assert margin > 0


@pytest.mark.parametrize("u", [-3, -1, 0, s.Rational(1, 3), 2])
@pytest.mark.parametrize(
    "X", [-s.Rational(1, 8192), 0, s.Rational(1, 2), 1, s.Rational(9, 8)]
)
def test_entire_original_R_and_positive_delta_inside_actual_range(u, X):
    with mp.workdps(90):
        uu, x = mpq(u), mpq(X)
        T = x**1024 / (x**1024 + (1 - x) ** 1024)
        bump = 1024 * x * x * mp.exp(-1024 * x * x)
        B = T + (1 - T) * bump
        R = 1 + B * (x - 1) / (1 + uu * uu) ** 3
        assert R > mpq(s.Rational(1, 2) + s.Rational(1, 8192))
        assert R < mpq(s.Rational(6, 5))
        # The exact additive global upper is a rational, not an underflowed sample.
        assert mpq(clock.data()["global_R_addition_bound"]) < mp.mpf("1e-414")


@pytest.mark.parametrize("j", range(9))
def test_switch_clock_jets_by_independent_shifted_polynomial(j):
    z, A = s.symbols("z A", real=True)
    exponential = s.series(s.exp(-A * (1 + z) ** 2), z, 0, max(1, 9 - j)).removeO()
    whole = s.expand(z**8 * exponential)
    got = whole.coeff(z, j) * s.factorial(j)
    expected = 0 if j < 8 else s.factorial(8) * s.exp(-A)
    assert s.cancel(got - expected) == 0


@pytest.mark.parametrize("j", range(5))
@pytest.mark.parametrize("x", [s.Rational(7, 8), s.Rational(15, 16), s.Rational(9, 8)])
def test_full_switch_derivative_bound_without_underflow(j, x):
    with mp.workdps(100):
        xx = mpq(x)
        for A in (mp.mpf(2), mp.mpf(100), mp.mpf("1e6")):
            f = lambda z, parameter=A: (1 - z) ** 8 * mp.exp(-parameter * z * z)
            got = mp.diff(f, xx, j)
            if got:
                assert (
                    mp.log(abs(got))
                    < mp.log(mp.factorial(j) * 16**j) - mpq(clock.COMPLEX_EXPONENT) * A
                )


@pytest.mark.parametrize(
    "u0,x0", list(itertools.product([-1, 0, 1], [s.Rational(7, 8), s.Rational(9, 8)]))
)
def test_full_complex_R_neighborhood_and_switch_factor(u0, x0):
    with mp.workdps(100):
        for angle in (0, mp.pi / 3, mp.pi / 2, mp.pi):
            u = mpq(u0) + mp.exp(mp.j * angle) / 16
            x = mpq(x0) + mp.exp(mp.j * (angle + mp.pi / 4)) / 16
            q = (1 - x) / x
            T = 1 / (1 + q**1024)
            bump = 1024 * x * x * mp.exp(-1024 * x * x)
            B = T + (1 - T) * bump
            R = 1 + B * (x - 1) / (1 + u * u) ** 3
            assert abs(R - 1) < mp.mpf(17) / 64
            assert abs(R) > mp.mpf("0.5")
            assert mp.re(x * x) >= mpq(clock.COMPLEX_EXPONENT)
            assert abs(1 - x) < 1


def test_every_mixed_jet_majorant_is_actual_and_complete():
    d = clock.data()
    expected = {f"{i}_{j}" for i in range(5) for j in range(5 - i)}
    assert set(d["all_mixed_derivative_factors"]) == expected
    assert d["complete_actual_four_jet_majorant"] < s.Rational(1, 10**2700)
    assert max(d["all_mixed_derivative_factors"].values()) == 25165824


@pytest.mark.parametrize(
    "x", [-s.Rational(1, 8192), s.Rational(1, 1000), s.Rational(1, 2), 1]
)
def test_full_regular_q_parameter_integral_and_ODE_on_both_sides(x):
    with mp.workdps(65):
        u = mp.mpf("0.3")
        a = lambda uu: (1 + uu * uu) / 10
        au = u / 5
        R = lambda xx: 1 + a(u) * xx * xx
        RX = lambda xx: 2 * a(u) * xx
        RU = lambda xx: au * xx * xx
        q = lambda xx: (
            3
            * xx
            * R(xx) ** mp.mpf("0.75")
            * mp.quad(
                lambda t: (
                    mp.sqrt(t) * RX(t * xx) * RU(t * xx) * R(t * xx) ** mp.mpf("-1.75")
                ),
                [0, 1],
            )
            / 4
        )
        xx = mpq(x)
        residual = (
            mp.diff(q, xx)
            + (1 / (2 * xx) - 3 * RX(xx) / (4 * R(xx))) * q(xx)
            - 3 * RX(xx) * RU(xx) / (4 * R(xx))
        )
        assert abs(residual) < mp.mpf("1e-55")
        if abs(xx) < mp.mpf("0.002"):
            assert abs(q(xx) / xx**4 - a(u) * au / 3) < mp.mpf("1e-8")


@pytest.mark.parametrize(
    "vector", [(1, 0, 0, 0), (0, s.sqrt(2) / 128, 0, 0), (1, 1, 0, 0)]
)
def test_full_covariant_metric_inverse_for_timelike_spacelike_and_null(vector):
    eta = s.diag(1, -1, -1, -1)
    v = s.Matrix(vector)
    X = (v.T * eta * v)[0]
    R = 1 + X * X / 10
    C = 1 / s.sqrt(R)
    D = X / (10 * s.sqrt(R) * (s.sqrt(R) + 1))
    physical = C * eta + D * v * v.T
    assert s.simplify(physical.det() + C**3) == 0
    assert s.simplify((v.T * physical.inv() * v)[0] - X) == 0
    inverse = physical / C - X * v * v.T / (10 * (s.sqrt(R) + 1))
    assert all(s.simplify(e) == 0 for e in inverse - eta)


def test_unmodified_original_frontier_and_new_separate_record():
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert len(audit.matching()) == 94
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 282
    assert "do not automatically transfer" in audit.observable()["quantum_boundary"]
    assert "No new heavy state" in audit.observable()["quantum_boundary"]


@pytest.mark.parametrize("phi", [s.Rational(1, 10), 1, 10**100, 10**400])
@pytest.mark.parametrize(
    "valley_factor", [0, s.Rational(1, 2), 1, s.Rational(5, 4), -1]
)
def test_full_actual_potential_against_literal_negative_function(phi, valley_factor):
    with mp.workdps(700):
        pp = mpq(phi)
        k, n, g, C, gamma = map(
            mpq, (family.K0, family.MASS2, family.G, family.CONTACT, family.GAMMA)
        )
        u = pp / mp.sqrt(k)
        tree = -(
            8800 * u**14
            + 53600 * u**12
            + 136800 * u**10
            + 180050 * u**8
            + 123200 * u**6
            + 38099 * u**4
            + 4573 * u**2
            + 224
        ) / (200 * (1 + u * u) ** 8)
        vac = -u * u / 2 + (-mp.mpf(1024) / 3 + mp.mpf(28) / 25) * u**4
        old = -k * (mp.exp(-(u**4)) * vac - mp.expm1(-(u**4)) * tree)
        H = mpq(valley_factor) * g * pp * pp / (2 * n)
        V = (
            old
            - gamma * pp**4 / 3
            - C * pp**4 / 24
            + n * H * H / 2
            - g * H * pp * pp / 2
        )
        bound = H * H / 6 + mpq(heavy.data()["actual_quartic_q"]) * pp**4 / 2
        assert V > bound
        beta = mpq(heavy.data()["actual_coercive_quartic_beta"])
        square = (
            old
            + H * H / 6
            + (n - mp.mpf(1) / 3)
            * (H - g * pp * pp / (2 * (n - mp.mpf(1) / 3))) ** 2
            / 2
            + beta * pp**4
        )
        assert abs(V - square) < mp.mpf("1e-450") * max(abs(V), 1)
