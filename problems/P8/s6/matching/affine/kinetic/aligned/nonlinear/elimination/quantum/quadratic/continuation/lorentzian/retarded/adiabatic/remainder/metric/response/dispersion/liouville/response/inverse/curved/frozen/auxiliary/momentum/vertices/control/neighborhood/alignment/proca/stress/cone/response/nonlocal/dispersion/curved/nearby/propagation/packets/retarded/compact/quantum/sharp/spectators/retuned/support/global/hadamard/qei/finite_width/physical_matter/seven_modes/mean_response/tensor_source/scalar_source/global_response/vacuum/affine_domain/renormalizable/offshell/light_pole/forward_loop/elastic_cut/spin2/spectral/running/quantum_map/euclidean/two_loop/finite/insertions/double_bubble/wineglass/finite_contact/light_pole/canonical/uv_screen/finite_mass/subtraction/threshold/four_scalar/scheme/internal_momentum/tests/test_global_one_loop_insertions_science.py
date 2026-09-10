"""Independent fraction, parameter-gap, routing and nonclosure regressions."""

from fractions import Fraction as F

import pytest
import sympy as sp
from p8_vacuum_global_one_loop_insertions import (
    audit,
    calibration,
    halfplane,
    routing,
    spacelike,
    window,
)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_named_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_every_explicit_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[v[0] for v in audit.bad_cases()]
)
def test_reject_every_unsupported_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "m,sigma", ((1, 1), (1, 3), (2, 1), (2, 15), (36, 2), (10**200, 1))
)
@pytest.mark.parametrize("Y,Q", ((0, 144), (F(1, 17), 144), (1, 256)))
def test_halfplane_bound_from_independent_fraction_moments(m, sigma, Y, Q):
    m, sigma, Y, Q = map(F, (m, sigma, Y, Q))
    gap = m * m - sigma / 4
    expected = 6 * Y / Q * (F(1, 6) / gap + m * m / (30 * gap * gap))
    d = halfplane.enclosure(m, Y, Q, sigma)
    assert d["parameter_gap_lower"] == sp.Rational(gap)
    assert d["unscaled_OS_insertion_uniform_upper"] == sp.Rational(expected)
    assert expected == 6 * Y * (6 * m * m - 5 * sigma / 4) / (30 * Q * gap * gap)


@pytest.mark.parametrize("x", (F(0), F(1, 4), F(1, 2), F(3, 4), F(1)))
@pytest.mark.parametrize(
    "real,imag", ((-(10**100), 10**200), (-1, 0), (0, 10**12), (1, 0), (1, 10**400))
)
def test_complex_halfplane_parameter_gap_with_unbounded_momentum_examples(
    x, real, imag
):
    m, sigma = F(36), F(1)
    A = x * (1 - x)
    gap = m * m - sigma / 4
    reD = m * m - A * real
    imD = -A * imag
    assert reD >= gap > 0
    assert reD * reD + imD * imD >= gap * gap
    assert 1 - 4 * A == (2 * x - 1) ** 2 >= 0
    assert reD - gap == A * (sigma - real) + sigma * (F(1, 4) - A)


@pytest.mark.parametrize(
    "m,x,t",
    (
        (2, F(1, 4), 0),
        (2, F(1, 2), 1),
        (36, F(1, 4), 100),
        (36, F(1, 2), 10**4),
        (10**200, F(1, 2), 0),
    ),
)
def test_spacelike_sign_with_independent_exact_log_series(m, x, t):
    m, x, t = map(F, (m, x, t))
    A = x * (1 - x)
    K = 4 * m * m - 1
    d = t + 1
    b = A / (m * m - A)
    h = d * b
    ratio = h / (2 + h)
    # Exact positive atanh lower bound, sufficient for these independent points.
    lo = 2 * sum(ratio ** (2 * j + 1) / F(2 * j + 1) for j in range(32))
    assert (K + d) * lo - K * h > 0
    assert 0 < b <= 1 / K
    assert h > 0 and lo > 0


@pytest.mark.parametrize(
    "m,Y,Q,t", ((2, 0, 144, 0), (2, F(1, 3), 144, 0), (36, F(1, 10**8), 256, 10**10))
)
def test_global_envelope_prefactors_independently(m, Y, Q, t):
    m, Y, Q, t = map(F, (m, Y, Q, t))
    d, K = t + 1, 4 * m * m - 1
    result = spacelike.envelope(m, Y, Q, t)
    expected = sp.Rational(12 * Y / Q) * sp.log(1 + sp.Rational(d / K))
    assert result["insertion_decay_envelope"] == expected / sp.Rational(d)
    assert result["minus_fermion_remainder_upper"] == expected * sp.Rational(d)
    assert result["coarse_uniform_insertion_upper"] == sp.Rational(12 * Y / (Q * K))


@pytest.mark.parametrize(
    "a,b", ((F(1, 2), F(0)), (F(3, 4), F(0)), (F(3, 4), F(1, 8)), (F(3, 4), F(-1, 8)))
)
@pytest.mark.parametrize(
    "q0,q1", ((0, 0), (1, -2), (-(10**200), 10**100), (10**400, 0))
)
@pytest.mark.parametrize("sign", (-1, 1))
def test_symmetric_routing_independent_complex_fraction_arithmetic(a, b, q0, q1, sign):
    real_s, imag_s = 4 * (a * a - b * b), 8 * a * b
    assert (real_s - 2) ** 2 + imag_s**2 <= 1
    # q0+sign*iE has real part q0-sign*b and imaginary part sign*a.
    real_line = a * a - (F(q0) - sign * b) ** 2 - F(q1) ** 2
    assert real_line <= a * a <= F(3, 4) < 1


@pytest.mark.parametrize(
    "center,radius,sigma", ((2, 1, 1), (0, 0, 1), (2, 2, 1), (2, 3, 1), (10, 2, 3))
)
def test_routing_interface_does_not_force_a_pass(center, radius, sigma):
    d = routing.channel_margin(center, radius, sigma)
    expected = F(center + radius, 4)
    assert d["line_invariant_real_part_upper"] == sp.Rational(expected)
    assert d["halfplane_margin"] == sp.Rational(F(sigma) - expected)
    assert d["selected_symmetric_bubble_routing_covered"] == (expected <= sigma)


@pytest.mark.parametrize(
    "M,x,t",
    (
        (3, F(0), 0),
        (3, F(1, 4), 1),
        (36, F(1, 2), 100),
        (10**197, F(3, 4), 10**800),
        (36, F(1), 0),
    ),
)
def test_scalar_parameter_gap_on_global_spacelike_domain(M, x, t):
    M, x, t = map(F, (M, x, t))
    A = x * (1 - x)
    Delta1 = x * M + (1 - x) ** 2
    b = A / Delta1
    d = t + 1
    assert Delta1 >= 1
    assert Delta1 + A * d == Delta1 * (1 + d * b)
    assert b >= 0


@pytest.mark.parametrize("m", (2, 36, 10**200))
def test_reference_energy_log_window_with_independent_fractions(m):
    m = F(m)
    energy = m * m
    K = 4 * m * m - 1
    argument = 1 + (energy * energy + 1) / K
    assert 1 < argument < m * m
    assert m * m - argument == (3 * m**4 - 5 * m * m) / K


@pytest.mark.parametrize(
    "r,C,ell,k",
    (
        (0, 0, 0, 1),
        (F(1, 10), F(1, 10), 1, 1),
        (0, 2, 1, 1),
        (1, 0, 0, 1),
        (F(1, 10**208), F(1, 10**207), 1000, F(999, 1000)),
    ),
)
def test_finite_window_defect_from_fraction_arithmetic(r, C, ell, k):
    r, C, ell, k = map(F, (r, C, ell, k))
    epsilon = (r + C * ell) / k
    d = window.defect_bound(r, C, ell, k)
    assert d["fractional_inverse_defect_upper"] == sp.Rational(epsilon)
    assert d["normalized_inverse_lower"] == sp.Rational(1 - epsilon)
    assert d["strictly_positive_on_declared_log_window"] == (epsilon < 1)


def test_full_inverse_does_not_drop_scalar_sector_or_fermion_sign():
    r, fp, S, Ferm, t = sp.symbols("r fp S Ferm t")
    k = 1 + r - fp
    gamma = 1 + t - (S + Ferm) / k
    assert sp.factor((1 - gamma / (1 + t)) - (S + Ferm) / (k * (1 + t))) == 0
    wrong = 1 + t + (Ferm - S) / k
    assert sp.factor(gamma - wrong) == -2 * Ferm / k


def test_actual_bounds_and_named_window():
    d = calibration.data()
    assert d["reference_energy_not_a_derived_cutoff"] == 10**400
    assert d["Euclidean_momentum_squared_window"] == [0, 10**800]
    assert 0 < d["canonical_halfplane_insertion_upper"] < sp.Rational(1, 10**606)
    w = d["complete_one_loop_reference_window_inverse_enclosure"]
    assert 0 < w["fractional_inverse_defect_upper"] < sp.Rational(1, 10**202)
    assert w["normalized_inverse_lower"] > 1 - sp.Rational(1, 10**202)


def test_valid_but_uncovered_or_inconclusive_domains():
    assert (
        routing.channel_margin(2, 3)["selected_symmetric_bubble_routing_covered"]
        is False
    )
    assert (
        window.defect_bound(0, 2, 1, 1)["strictly_positive_on_declared_log_window"]
        is False
    )
    assert (
        window.defect_bound(1, 0, 0, 1)["strictly_positive_on_declared_log_window"]
        is False
    )


def test_scope_keeps_the_outer_loop_and_cutoff_open():
    assert "not a full two-loop bound" in halfplane.data()["scope"]
    assert "not arbitrary graph" in routing.channel_margin(2, 1)["scope"]
    assert "not a derived Wilsonian cutoff" in window.defect_bound(0, 0, 0, 1)["scope"]
    assert "reflection positivity" in window.data()["scope"]
    assert audit.controls()["original_P8_not_closed"] is True


def test_exact_counts():
    assert len(audit.residuals()) == 51
    assert len(audit.gates()) == 31
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 196
