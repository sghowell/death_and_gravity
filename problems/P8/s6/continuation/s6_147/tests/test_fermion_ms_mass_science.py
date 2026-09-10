"""Independent dimensional, mass-integral, finite-part and scope checks."""

import copy

import mpmath as mp
import pytest
import sympy as sp
from p8_vacuum_fermion_ms_mass import (
    audit,
    calibration,
    correction,
)
from p8_vacuum_fermion_ms_slopes import audit as previous


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_exact_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_invalid_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def gamma_matrices():
    sigma = (
        mp.matrix([[0, 1], [1, 0]]),
        mp.matrix([[0, -1j], [1j, 0]]),
        mp.matrix([[1, 0], [0, -1]]),
    )
    g = []
    for x in sigma:
        g.append(
            mp.matrix(
                [
                    [0, 0, -1j * x[0, 0], -1j * x[0, 1]],
                    [0, 0, -1j * x[1, 0], -1j * x[1, 1]],
                    [1j * x[0, 0], 1j * x[0, 1], 0, 0],
                    [1j * x[1, 0], 1j * x[1, 1], 0, 0],
                ]
            )
        )
    g.append(mp.matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]))
    g.append(g[0] * g[1] * g[2] * g[3])
    return g


@pytest.mark.parametrize("dim", (3, 4, 5))
@pytest.mark.parametrize("sector", ("scalar", "gauge"))
@pytest.mark.parametrize("config", (0, 1))
def test_independent_matrix_vacuum_derivative_and_three_words(dim, sector, config):
    with mp.workdps(45):
        gam = gamma_matrices()[:dim]
        k = [mp.mpf(j + 1 + config) / 5 for j in range(dim)]
        l = [mp.mpf((-1) ** j * (j + 2)) / (7 + config) for j in range(dim)]
        trace = lambda a: sum(a[j, j] for j in range(4))

        def S(m, v):
            return (
                m * mp.eye(4) - 1j * sum((g * x for g, x in zip(gam, v)), mp.zeros(4))
            ) / (m * m + sum(x * x for x in v))

        def vacuum(m):
            sk, sl = S(m, k), S(m, l)
            if sector == "scalar":
                return trace(sk * sl) / 2
            return -sum(trace(sk * g * sl * g) for g in gam) / 2

        m = mp.mpf("1.3")
        sk, sl = S(m, k), S(m, l)
        if sector == "scalar":
            words = trace(sk**3 * sl) + trace(sk**2 * sl**2) + trace(sk * sl**3)
        else:
            words = -sum(
                trace(sk**3 * g * sl * g)
                + trace(sk**2 * g * sl**2 * g)
                + trace(sk * g * sl**3 * g)
                for g in gam
            )
        assert abs(mp.diff(vacuum, m, 2) - words) < mp.mpf("1e-38")


@pytest.mark.parametrize(
    "sector,pole2,pole1,finite", (("scalar", 36, -48, -56), ("gauge", -72, 24, 40))
)
def test_independent_complex_MS_mass_Laurent_extraction(sector, pole2, pole1, finite):
    with mp.workdps(70):
        values = []
        for j in range(32):
            e = mp.mpf(".0001") * mp.exp(2j * mp.pi * j / 32)
            A = mp.exp(mp.euler * e) * mp.gamma(1 + e)
            if sector == "scalar":
                V = 4 / ((1 - e) * (2 * e - 1)) + 1 / (1 - e) ** 2
                ct = 2 * (6 - 3 * e) * (4 - 2 / (e - 1))
            else:
                V = -4 / ((1 - e) * (2 * e - 1)) + (2 - 2 * e) / (1 - e) ** 2
                ct = -4 * (6 - 3 * e) * (4 - 2 / (e - 1))
            raw = (4 - 4 * e) * (3 - 4 * e) * V
            values.append((A * A * raw + A * ct - pole2) / e**2 - pole1 / e)
        assert abs(sum(values) / 32 - finite) < mp.mpf("1e-50")


def test_independent_fixed_mu_full_field_mass_and_vertex_variation():
    h, z, eta, nu, e, m, mu = sp.symbols("h z eta nu e m mu", positive=True)
    changed_m = m * (1 + h * eta) / (1 + h * z)
    original = m**2 * (mu / m) ** (2 * e)
    changed = (
        (1 + h * nu) ** 2
        / (1 + h * z) ** 2
        * changed_m**2
        * (mu / changed_m) ** (2 * e)
    )
    expected = 2 * (nu - z) * original + (eta - z) * m * sp.diff(original, m)
    assert sp.simplify(sp.diff(changed, h).subs(h, 0) - expected) == 0
    assert sp.simplify(m * sp.diff(original, m) - (2 - 2 * e) * original) == 0
    assert sp.simplify(m * sp.diff(original, m) - 2 * original) != 0
    # Setting mu=m first loses a finite epsilon times pole term.
    assert sp.diff(original.subs(mu, m), m) == 2 * m


@pytest.mark.parametrize("e", (".17", ".37", ".6"))
@pytest.mark.parametrize("v,b", (("7", ".2"), ("20", "1"), ("3", "2")))
def test_general_mass_triangle_by_independent_radial_integral(e, v, b):
    with mp.workdps(60):
        e, v, b = map(mp.mpf, (e, v, b))
        radial = mp.quad(
            lambda x: x ** (1 - e) / ((x + v) * (x + b) ** 2), [0, b, v, mp.inf]
        ) / mp.gamma(2 - e)
        bracket = b ** (-e) / (v - b) - (v ** (1 - e) - b ** (1 - e)) / (
            (1 - e) * (v - b) ** 2
        )
        exact = mp.gamma(e) * bracket
        assert abs(radial / exact - 1) < mp.mpf("1e-30")


@pytest.mark.parametrize("v,b", (("7", ".2"), ("20", "1"), ("3", "2"), ("1", ".0001")))
def test_general_mass_triangle_finite_limit(v, b):
    with mp.workdps(60):
        v, b = map(mp.mpf, (v, b))
        radial = mp.quad(lambda x: x / ((x + v) * (x + b) ** 2), [0, b, v, mp.inf])
        exact = (v * mp.log(v / b) - v + b) / (v - b) ** 2
        assert abs(radial / exact - 1) < mp.mpf("1e-45")


@pytest.mark.parametrize("x", (".001", ".2", "1", "10", "1e8"))
def test_complete_scalar_propagator_mass_transfer(x):
    with mp.workdps(55):
        x = mp.mpf(x)
        integral = mp.quad(lambda b: 1 / (x + b) ** 2, [0, mp.mpf(".001"), 1])
        assert abs((1 / (x + 1) - 1 / x) + integral) < mp.mpf("1e-43")
        assert abs(integral * x * (x + 1) - 1) < mp.mpf("1e-40")


@pytest.mark.parametrize("e", (".2", ".6", "-.2"))
def test_quartic_counterterm_mass_integral_not_scaleless_at_b_one(e):
    with mp.workdps(70):
        e = mp.mpf(e)
        integral = mp.gamma(e) * mp.quad(lambda b: b ** (-e), [0, mp.mpf(".001"), 1])
        tad = mp.gamma(e) / (e - 1)
        assert abs(integral / tad + 1) < mp.mpf("1e-25")
        assert tad != 0


def hard_B(e):
    return (
        mp.sqrt(mp.pi)
        * mp.exp(2 * mp.euler * e)
        * 4 ** (-e)
        * (1 - 2 * e / 3)
        * (1 - 4 * e)
        * mp.gamma(1 + e)
        * mp.gamma(1 + 2 * e)
        / ((1 - e) * (1 + 2 * e) * mp.gamma(mp.mpf(".5") + e))
    )


@pytest.mark.parametrize("ell", ("0", "3", "12"))
def test_regulator_first_mass_integral_finite_part(ell):
    with mp.workdps(70):
        ell = mp.mpf(ell)
        values = []
        for j in range(32):
            e = mp.mpf(".0001") * mp.exp(2j * mp.pi * j / 32)
            A = mp.exp(mp.euler * e) * mp.gamma(1 + e)
            normal = -mp.exp(ell * e) * (
                (12 - 32 * e + 16 * e * e) * A * A - 12 * A
            ) / (1 - e) + 6 * hard_B(e)
            values.append((normal - 6) / e**2 + 2 / e)
        assert abs(sum(values) / 32 - (78 + 32 * ell)) < mp.mpf("1e-48")


@pytest.mark.parametrize("e", (".13", ".23"))
@pytest.mark.parametrize("ell", ("0", "7"))
def test_nonzero_regulator_b_integral_before_finite_extraction(e, ell):
    with mp.workdps(65):
        e, ell = map(mp.mpf, (e, ell))
        A = mp.exp(mp.euler * e) * mp.gamma(1 + e)
        soft1 = mp.exp(ell * e) * ((12 - 32 * e + 16 * e * e) * A * A - 12 * A)
        integrand = lambda b: -(soft1 * b ** (-e) - 6 * hard_B(e))
        numerical = mp.quad(integrand, [0, mp.mpf(".001"), 1])
        exact = -soft1 / (1 - e) + 6 * hard_B(e)
        assert abs(numerical - exact) < mp.mpf("1e-40")


def finite_ratio_delta(r):
    # z=1-u^2 removes the square-root endpoint. h(k) is evaluated
    # as k(2-k)/(1-k)^2, avoiding cancellation at tiny r.
    def integrand(u):
        z = 1 - u * u
        k = r * z
        return (
            -24
            * r
            * (1 - 2 * u * u)
            / ((1 - k) ** 2)
            * ((2 - k) * (-mp.log(k) - 1) + 1)
        )

    return mp.quad(integrand, [0, mp.mpf(".5"), 1])


@pytest.mark.parametrize("r", ("0.0625", ".01", ".0001", "1e-20", "1e-100"))
def test_whole_finite_ratio_remainder_integral_against_majorant(r):
    with mp.workdps(65):
        r = mp.mpf(r)
        delta = finite_ratio_delta(r)
        assert abs(delta) < r * (300 + 100 * mp.log(1 / r))


@pytest.mark.parametrize("T", ("16", "100"))
def test_integrated_finite_ratio_remainder_against_full_b_majorant(T):
    with mp.workdps(35):
        T = mp.mpf(T)
        result = mp.quad(lambda b: finite_ratio_delta(b / T), [0, mp.mpf(".1"), 1])
        assert abs(result) < (175 + 50 * mp.log(T)) / T


@pytest.mark.parametrize("T", ("16", "100", "1e400"))
def test_independent_b_log_integral_constant(T):
    with mp.workdps(55):
        T = mp.mpf(T)
        result = mp.quad(
            lambda b: b * (300 + 100 * mp.log(T / b)), [0, mp.mpf(".01"), 1]
        )
        assert abs(result - (175 + 50 * mp.log(T))) < mp.mpf("1e-45")


@pytest.mark.parametrize(
    "m,Y,Q",
    ((720, 0, 144), (721, sp.Rational(1, 3), 100), (10000, sp.Rational(1, 10**9), 144)),
)
def test_exact_enclosure_scaling_and_dyadic_domains(m, Y, Q):
    b = correction.enclosure(m, Y, Q)
    assert 2 ** (b["dyadic_m_squared"] - 1) < m * m <= 2 ** b["dyadic_m_squared"]
    assert 2 ** (b["dyadic_threshold"] - 1) < 4 * m * m <= 2 ** b["dyadic_threshold"]
    assert b["complete_scalar_mass_difference_absolute_upper"] >= 0
    doubled = correction.enclosure(m, 2 * Y, Q)
    assert (
        doubled["complete_scalar_mass_difference_absolute_upper"]
        == 4 * b["complete_scalar_mass_difference_absolute_upper"]
    )


def test_actual_reference_and_on_shell_decomposition():
    d = calibration.data()
    assert all(v is True for v in d["bounds"].values())
    assert d["actual_on_shell_MS_mass_reference_absolute_upper"] < sp.Rational(
        1, 10**10
    )
    assert (
        d["actual_on_shell_MS_mass_reference_absolute_upper"]
        > d["actual_zero_momentum_MS_mass_reference_absolute_upper"]
        > 0
    )
    assert (
        d["actual_zero_momentum_MS_mass_reference_absolute_upper"]
        > d["massless_gauge_reference_absolute_upper"]
    )
    assert d["actual_zero_soft_MS_slope_absolute_upper"] > 0
    assert d["previous_soft_tail_upper"] > 0


@pytest.mark.parametrize("degree", (2, 3, 4, 7))
def test_on_shell_mass_keeps_both_affine_coefficients(degree):
    x, f0, f1, c = sp.symbols("s f0 f1 c")
    f = f0 + x * f1 + c * x**degree
    tail = f - f.subs(x, 0) - x * sp.diff(f, x).subs(x, 0)
    assert sp.expand(f.subs(x, 1) - f0 - f1 - tail.subs(x, 1)) == 0
    assert sp.diff(tail, x).subs(x, 0) == 0
    assert tail.subs(x, 0) == 0


def test_exact_two_row_frontier_change_and_no_parent_mutation():
    old = copy.deepcopy(previous.frontier())
    new = audit.frontier()
    changed = [a["id"] for a, b in zip(old, new) if a != b]
    assert changed == list(audit.TARGETS)
    assert previous.frontier() == old
    assert len([r for r in new if r["status"] == "UNEVALUATED"]) == 2
    assert all(r["status"] != "COMPLETE" for r in new)
    old_older = next(r for r in old if r["id"] == "scalar_Phi2_W1_F0")
    assert next(r for r in new if r["id"] == old_older["id"]) == old_older


def test_research_boundary_controls():
    d = audit.controls()
    assert len(d) == 9
    assert d["finite_reference_not_a_new_adjustable_parameter"] is True
    assert d["two_vacuum_rows_still_open"] is True
    assert d["original_P8_not_closed"] is True
