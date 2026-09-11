"""Independent retarded integration, Cauchy recurrence and spatial controls."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_retarded_reference_boundary import (
    audit,
    boundary,
    bounds,
    ibp,
    majorants,
)


@cache
def literal_reference():
    t = s.Symbol("t", real=True)
    g = 1 / (3 + 2 * t)
    b = (1 + s.I * t + t * t) * t**7 * (s.Rational(1, 4) - t) ** 7
    rows = [b]
    for _ in range(6):
        rows.append(s.cancel(s.diff(g * rows[-1], t)))
    return t, g, rows, tuple(s.lambdify(t, v, "mpmath") for v in rows)


@pytest.mark.parametrize("order", range(1, 7))
@pytest.mark.parametrize("at", (s.Rational(1, 16), s.Rational(1, 8), s.Rational(3, 16)))
def test_literal_variable_phase_retarded_integral_and_every_boundary(order, at):
    with mp.workdps(100):
        _t, _g, _rows, fn = literal_reference()
        t = mp.mpf(str(at.p)) / int(at.q)
        theta = lambda x: 3 * x + x * x
        direct = mp.exp(-1j * theta(t)) * mp.quad(
            lambda x: mp.exp(1j * theta(x)) * fn[0](x), [0, t]
        )
        terms = [-1j * (1j**j) * fn[j](t) / (3 + 2 * t) for j in range(order)]
        remainder = (
            (1j**order)
            * mp.exp(-1j * theta(t))
            * mp.quad(lambda x: mp.exp(1j * theta(x)) * fn[order](x), [0, t])
        )
        assert abs(direct - sum(terms) - remainder) < mp.mpf("1e-85")
        assert abs(direct) > mp.mpf("1e-22")
        if order == 6:
            assert abs(sum(terms[:5]) + remainder - direct) > mp.mpf("1e-20")


@pytest.mark.parametrize("order", range(1, 7))
def test_literal_common_preparation_eliminates_all_initial_terms(order):
    t, g, rows, _fn = literal_reference()
    assert all(s.simplify((g * rows[j]).subs(t, 0)) == 0 for j in range(order))
    # A source not vanishing near preparation must keep its lower boundary.
    assert (-s.I * g * (1 + t)).subs(t, 0) == -s.I / 3


@pytest.mark.parametrize("order", range(7))
def test_independent_rational_derivatives_reconstruct_every_Cauchy_coefficient(order):
    t = s.Symbol("t", real=True)
    # Their derivatives at t0 equal the factorial inputs to the
    # algebraic recursion; no bounded-disc extremality is asserted.
    h = s.symbols("h0:7")
    a = 1 / (1 - t)
    g = 1 / (1 - t)
    polynomial = sum(h[j] * t**j / s.factorial(j) for j in range(7))
    expression = a * polynomial
    for _ in range(order):
        expression = s.diff(g * expression, t)
    result = s.expand(expression.subs(t, 0))
    assert (
        tuple(result.coeff(h[j]) for j in range(order + 1))
        == majorants.recurrence()[order]
    )


@pytest.mark.parametrize("power", (5, 6))
def test_exact_zero_transfer_radial_bulk_power_boundary(power):
    x = s.Symbol("x", positive=True)
    radial = x * x * (1 + x * x) ** (1 - s.Rational(power, 2)) / 2**power
    value = s.integrate(radial, (x, 0, s.oo))
    if power == 5:
        assert value == s.oo
    else:
        assert s.simplify(value - s.pi / 256) == 0


@pytest.mark.parametrize("K", (1000, 10**6, 10**16))
def test_massive_removed_tail_and_complete_pair_numerators(K):
    with mp.workdps(80):
        A = mp.mpf(25) / 16
        m = mp.mpf(1000)
        R = mp.mpf(K)
        J4tail = (
            A**4
            * mp.quad(lambda y: 1 / (1 + (A * m * y / R) ** 2) ** 2, [0, 1])
            / (2 * mp.pi**2 * R)
        )
        assert 0 < J4tail < A**4 / (2 * mp.pi**2 * R)
        c = bounds.constants()
        coefficient = mp.mpf(
            str(c["sixth_bulk_pair_numerator"] + c["fifth_endpoint_pair_numerator"])
        )
        assert coefficient * J4tail / 2 < mp.mpf("1e52") / R
        assert (
            4 * bounds.tail_bound(K) / s.Integer(10) ** 800
            == 4 * s.Rational(1, 10) ** 748 / K
        )


@pytest.mark.parametrize("degree", range(1, 9))
def test_equal_time_inverse_phase_has_unbounded_spatial_series_degree(degree):
    x = s.Symbol("x", real=True)
    coefficient = s.binomial(s.Rational(1, 2), degree + 1)
    expression = 1 / (1 + s.sqrt(1 + x * x))
    actual = s.series(expression, x, 0, 2 * degree + 1).removeO().coeff(x, 2 * degree)
    assert s.simplify(actual - coefficient) == 0 and actual != 0


def test_complex_reference_cannot_drop_odd_boundary_terms_before_contraction():
    t = s.Symbol("t", real=True)
    g = 1 / (3 + 2 * t)
    amplitude = 1 + s.I * t
    detector = 1 - s.I
    first = s.conjugate(detector) * g * s.diff(g * amplitude, t)
    assert s.simplify(s.im(first).subs(t, s.Rational(1, 8))) != 0
    assert ibp.data()["boundary_coefficients"] == [-s.I, 1, s.I, -1, -s.I, 1]


@pytest.mark.parametrize("name", ("ibp", "majorants", "bounds", "boundary"))
def test_exact_scientific_packet(name):
    p = {"ibp": ibp, "majorants": majorants, "bounds": bounds, "boundary": boundary}[
        name
    ].data()
    assert all(s.cancel(v) == 0 for v in p["checks"].values())
    assert all(p["gates"].values())


def test_complete_known_finite_piece_uses_both_predecessor_norms():
    assert bounds.DISPLAY + 2 * 10**23 < 2 * 10**48
    assert bounds.TAIL + 10**27 < 2 * 10**52
    assert "N61[Gamma]" in boundary.data()["complete_known_piece"]
    assert (
        "does not bound the remaining contact"
        in boundary.data()["complete_known_piece"]
    )


@pytest.mark.parametrize(
    "name,value",
    audit.residuals().items(),
    ids=lambda x: x if isinstance(x, str) else None,
)
def test_all_exact_residuals(name, value):
    assert s.cancel(value) == 0, name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=lambda x: x if isinstance(x, str) else None
)
def test_all_scope_mutations_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_complete_counts_and_unclosed_frontier():
    assert len(audit.residuals()) == 57 and audit.scalar_entry_count() == 57
    assert len(audit.gates()) == 39 and all(audit.gates().values())
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 120
    assert len(audit.frontier()) == 9 and audit.matching()[-1] == audit.ITEM


def test_fixed_prescription_and_projector_boundaries_are_explicit():
    assert "original covariant prescription" in boundary.data()["spatial_warning"]
    assert "one-momentum band" in boundary.data()["separate_projectors"]
    assert (
        "no independent reference Ward identity" in boundary.data()["reference_warning"]
    )
