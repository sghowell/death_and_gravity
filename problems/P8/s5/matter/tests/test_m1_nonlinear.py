import pytest
import sympy as sp
from p8_m1 import njets, nonlinear, regressions, series
from p8_s5.lapse_series import generic_checks


def test_full_chart_legendre_and_prior_linear_bridges():
    assert all(value == 0 for value in nonlinear.audit_checks().values())


def test_independent_literal_bounce_all_forty_jets():
    assert all(value == 0 for value in regressions.bounce_checks().values())


def test_boundary_cannot_be_dropped_at_the_bounce():
    result = regressions.boundary_negative_control()
    assert result["at_bounce"] == "18"
    omitted, full = nonlinear.lapse_jets(omit_boundary=True), nonlinear.lapse_jets()
    # These innocent-looking bridges would not detect the missing term.
    for key in ("A", "B", "C", "L"):
        assert sp.cancel(omitted[key][1]-full[key][1]) == 0
    assert sp.cancel(omitted["A"][2]-full["A"][2]) != 0


def test_jet_algebra_is_ordinary_taylor_not_derivative_arithmetic():
    n = sp.Symbol("n")
    assert njets.derivative_coefficients(njets.Npower(3)) == (1, 3, 6, 6, 0)
    for p in (sp.Rational(-7, 4), sp.Rational(3, 4), -2, 3):
        result = njets.power((1, 2, 3, 4, 5), p)
        direct = sp.series((1+2*n+3*n**2+4*n**3+5*n**4)**p, n, 0, 5).removeO()
        assert sp.expand(direct-sum(v*n**k for k, v in enumerate(result))) == 0
    with pytest.raises(ValueError):
        njets.power((2, 0, 0, 0, 0), -1)


def test_generic_stationary_series_including_second_order_lapse():
    assert all(value == 0 for value in generic_checks().values())
    symbols, result = series.symbolic_series()
    assert [len(series.coefficients(expr)) for expr in result["hamiltonian"]] == [1, 3, 8, 16, 30]
    assert [len(series.coefficients(expr)) for expr in result["lapse"]] == [3, 8, 16]
    for name, orders in result.items():
        start = 1 if name == "lapse" else 0
        for k, expr in enumerate(orders, start):
            for powers in series.coefficients(expr):
                assert sum(int(p)*w for p, w in zip(powers.split(","), nonlinear.WEIGHTS)) == k
    # The generic quartic correction from n2 must survive, including matter.
    n2 = result["lapse"][1]
    assert n2.has(nonlinear.eta, nonlinear.z)
    correction = -symbols["A"][2]*n2**2/2
    assert correction != 0 and correction.has(nonlinear.eta, nonlinear.z)


def test_compact_scaling_and_uniform_coefficient_bounds():
    assert all(value == 0 for value in series.compact_checks().values())
    report = series.coefficient_report()
    assert sum(len(row) for rows in report.values() for row in rows) == 85
    for rows in report.values():
        for row in rows:
            for entry in row.values():
                assert abs(sp.Rational(entry["at_bounce"])) <= sp.Rational(entry["uniform_absolute_upper"])


@pytest.mark.parametrize("bad", [1/nonlinear.u, 1/(nonlinear.u**2-1), nonlinear.u, nonlinear.u**2])
def test_invalid_compact_rationals_rejected(bad):
    with pytest.raises(ValueError):
        series.compact_even(bad)


def test_unexpected_stationary_denominator_is_not_certified():
    symbols, _ = series.symbolic_series()
    with pytest.raises(ValueError):
        series.coefficient_bound(1/(symbols["A"][2]+1), {})


def test_compact_bound_at_non_bounce_points():
    compact = series.compact_jets()
    actual = nonlinear.lapse_jets()
    scale_degree = {"A": 2, "B": 1, "C": 0, "L": 1, "D": 0, "E": 0, "M": 0, "Z": 0}
    for u, x, y in ((sp.Rational(3, 4), sp.Rational(3, 5), sp.Rational(4, 5)),
                    (sp.Rational(-4, 3), sp.Rational(-4, 5), sp.Rational(3, 5))):
        for key in actual:
            for a, c in zip(actual[key], compact[key]):
                assert sp.cancel(a.subs(nonlinear.u, u)/y**scale_degree[key]-c.subs({series.x: x, series.y: y})) == 0
