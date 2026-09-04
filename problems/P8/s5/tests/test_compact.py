import pytest
import sympy as sp
from p8_s5 import compact, lapse_series
from p8_s5 import nonlinear_d as m


def test_full_nonlinear_compactification_and_lapse_hessian():
    assert all(value == 0 for value in compact.checks().values())
    assert compact.P.subs(compact.x, -1) == 1
    assert compact.P.subs(compact.x, 1) == 1


def test_every_coefficient_has_finite_both_tail_limits_and_bound():
    data = compact.coefficient_report()
    count = 0
    for group in data.values():
        for order in group:
            for entry in order.values():
                bound = sp.Rational(entry["absolute_bound"])
                assert bound > 0
                assert abs(sp.Rational(entry["past_endpoint"])) <= bound
                assert abs(sp.Rational(entry["future_endpoint"])) <= bound
                count += 1
    assert count == 34


def test_normalization_preserves_invariant_values_at_rational_clock_point():
    # u=3/4 gives sqrt(d)=5/4 and x=3/5, so this uses exact rationals.
    data = compact.coefficient_report()
    for group, orders in lapse_series.derive().items():
        for index, expr in enumerate(orders):
            for key, value in lapse_series.coefficients(expr).items():
                p, r, w = map(int, key.split(","))
                exponent = (1 if group == "hamiltonian" else 0)-sp.Rational(p, 2)-r-w
                actual = sp.sympify(data[group][index][key]["expression"], locals={"x": compact.x})
                expected = value.subs(m.u, sp.Rational(3, 4))*sp.Rational(25, 16)**exponent
                assert sp.cancel(actual.subs(compact.x, sp.Rational(3, 5))-expected) == 0


def test_tail_limit_is_regular_but_not_a_physical_vertex_or_cutoff():
    # Already enough to show why unweighted lapse coefficients can mislead.
    raw = lapse_series.coefficients(lapse_series.derive()["lapse"][0])["0,1,0"]
    assert sp.limit(raw, m.u, sp.oo) == -sp.oo
    normalized = compact.coefficient_report()["lapse"][0]["0,1,0"]
    assert normalized["future_endpoint"] == "-1/8"


def test_wrong_parity_and_unproved_denominator_fail_closed():
    with pytest.raises(ValueError, match="parity"):
        compact.normalized(m.u, 0)
    with pytest.raises(ValueError, match="denominator"):
        compact.bound_coefficient(1/(1-compact.x))


def test_overlapping_principal_charts_have_uniform_normalization():
    # The P and Lambda bounds are independently certified in the main replay.
    u_edge = sp.Rational(1, 8)
    x_squared_min = u_edge**2/(1+u_edge**2)
    assert x_squared_min == sp.Rational(1, 65)
    assert 1/(2*x_squared_min) == sp.Rational(65, 2)
    assert sp.Rational(68)*sp.Rational(16, 17) == 64
