"""Independent tensor, full-D pole, physical sign and retained-scope tests."""

import pytest
import sympy as s
from p8_affine import verify as base
from p8_vacuum_affine_massive_dimensional_cut_completion import dimensional
from p8_vacuum_affine_massive_gravity_pole_completion import (
    audit,
    graphs,
    poles,
    soft,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not value.atoms(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_every_unsupported_input(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("D", (4, 5, 6))
@pytest.mark.parametrize("choice", (1, 2, 3))
def test_independent_literal_box_currents_and_exact_product(D, choice):
    E = s.Integer(choice + 3)
    Q = s.Integer(choice)
    W = s.Integer(1)
    eta = s.diag(1, *([-1] * (D - 1)))
    p = s.Matrix([E, -Q, W, *([0] * (D - 3))])
    r = s.Matrix([E, Q, -W, *([0] * (D - 3))])
    q = s.Matrix([0, 2 * Q, *([0] * (D - 2))])
    k = s.Matrix([s.Rational((i + 1) * choice, 7) for i in range(D)])
    dot = lambda x, y: (x.T * eta * y)[0]
    mu = dot(p, p)
    ss = dot(p + r, p + r)
    assert mu > 0
    T = lambda x, y: x * y.T + y * x.T - eta * (dot(x, y) - mu)
    H = lambda x, y: (
        s.trace(eta * x * eta * y) - s.trace(eta * x) * s.trace(eta * y) / (D - 2)
    )
    actual1 = H(T(p, p - k), T(r, r + k))
    actual2 = H(T(p + q, p - k), T(r - q, r + k))
    d1 = dot(k, k)
    d2 = dot(p - k, p - k) - mu
    d3 = dot(k + q, k + q)
    d4 = dot(r + k, r + k) - mu
    v = poles.eikonal(ss, mu, D)
    z = ss - 2 * mu
    assert actual1 == v + z * (d1 - d2 - d4)
    assert actual2 == v + z * (d3 - d2 - d4)
    assert (
        s.factor(
            actual1 * actual2
            - v * v
            - v * z * (d1 + d3 - 2 * d2 - 2 * d4)
            - z * z * (d1 * d3 - (d1 + d3) * (d2 + d4) + (d2 + d4) ** 2)
        )
        == 0
    )


@pytest.mark.parametrize("D", (4, s.Rational(9, 2), 5))
@pytest.mark.parametrize("mu", (s.Integer(1), s.Rational(3, 2)))
def test_independent_frozen_master_and_complete_crossed_principal_parts(D, mu):
    D = s.sympify(D)
    a, b = s.symbols("test_channel test_other", real=True)
    e = (D - 4) / 2
    actual = dimensional.master_coefficients(a, b, mu, e)
    breg = s.cancel(a * a * actual["Bmm"])
    creg = s.cancel(a * actual["C0mumu"])
    b2 = s.factor(breg.subs(a, 0))
    b1 = s.factor(s.diff(breg, a).subs(a, 0))
    c1 = s.factor(creg.subs(a, 0))
    r = s.cancel(a * a * poles.completion(a, b, mu, D))
    assert s.factor(r.subs(a, 0) - b2) == 0
    expected = (
        b1
        - e * b2 / (6 * mu)
        + c1 / (2 * mu)
        + 2 * mu * poles.residue_ratio(D) * poles.eikonal(b, mu, D)
    )
    assert s.factor(s.diff(r, a).subs(a, 0) - expected) == 0


def test_evanescent_tree_cannot_be_dropped_from_finite_lsz_pole():
    a, b, mu, e = s.symbols("a b mu e", positive=True)
    expr = (
        2
        * mu
        * (-1 / e)
        * poles.residue_ratio(4 + 2 * e)
        * (poles.newton_shape(a, b, mu, 4 + 2 * e) - poles.newton_shape(a, b, mu, 4))
    )
    finite = s.factor(s.limit(expr, e, 0))
    assert s.factor(finite + 12 * poles.completion_basis(a, b, mu)[0]) == 0
    assert finite != 0


@pytest.mark.parametrize("i", range(4))
def test_complete_finite_raw_coefficients_not_only_D4(i):
    e, L = s.symbols("epsilon finite_log", real=True)
    exact = poles.completion_coefficients(4 + 2 * e)[i]
    finite = s.series((-1 / e + L) * exact, e, 0, 1).removeO()
    r0 = poles.data()["whole_D4_meromorphic_coefficients"][i]
    r1 = poles.data()["whole_EP_linear_meromorphic_coefficients"][i]
    assert s.factor(finite - (-r0 / e + L * r0 - r1)) == 0
    assert r1 != 0


@pytest.mark.parametrize("mu", (s.Integer(1), s.Integer(2), s.Rational(3, 2)))
@pytest.mark.parametrize("Q", (s.Integer(1), s.Integer(7)))
@pytest.mark.parametrize(
    "fraction", (s.Rational(1, 4), s.Rational(1, 2), s.Rational(3, 4))
)
@pytest.mark.parametrize("x", (s.Integer(0), s.Rational(1, 3), s.Integer(1)))
def test_independent_exact_positive_density_and_global_majorants(mu, Q, fraction, x):
    tau = Q * fraction
    P = 1 + 6 * x * x + x**4
    value = s.factor(soft.attenuation_density(Q, tau, mu, x))
    lower = 8 * mu * mu * tau * (Q - tau) * P / (4 * mu + Q) ** 3
    upper = 8 * mu * mu * tau * (Q - tau) * P / (4 * mu) ** 3
    assert 0 < lower <= value <= upper
    lo, hi = soft.attenuation_bounds(Q, tau, mu)
    assert 0 < lo <= hi


@pytest.mark.parametrize("endpoint", (0, 1))
def test_both_physical_zero_attenuation_endpoints(endpoint):
    assert soft.attenuation_bounds(7, endpoint * 7, 2) == (0, 0)
    assert soft.attenuation_density(7, endpoint * 7, 2, s.Rational(2, 3)) == 0


@pytest.mark.parametrize(
    "bad", (True, False, 1.0, s.Float(1), s.I, s.oo, None, "1", -1, 2)
)
def test_physical_transfer_domain_rejects_inexact_and_outside(bad):
    with pytest.raises((TypeError, ValueError)):
        soft.attenuation_bounds(1, bad)


@pytest.mark.parametrize("bad", (0, -1, True, 1.0, s.Float(1), None))
def test_positive_mass_and_excess_domains_reject(bad):
    with pytest.raises((TypeError, ValueError)):
        soft.attenuation_bounds(bad, 0, 1)
    with pytest.raises((TypeError, ValueError)):
        soft.attenuation_bounds(1, 0, bad)


@pytest.mark.parametrize("beta", (s.Rational(1, 3), s.Rational(3, 5), s.Rational(4, 5)))
def test_literal_massive_forward_phase_not_zero_after_endpoint_charge(beta):
    mu = s.Integer(1)
    energy = 4 * mu / (1 - beta * beta)
    v = energy**2 - 4 * mu * energy + 2 * mu * mu
    ms = (-s.atanh(beta) + s.I * s.pi / 2) / (energy * beta)
    crossed = s.atanh(beta) / (energy * beta)
    full = s.simplify(2 * (v * ms + mu / 2 + v * crossed) - mu)
    assert full == s.I * s.pi * v / (energy * beta)
    assert full != 0
    assert s.re(full) == 0


def test_entire_minimal_topology_inventory_and_bound():
    rows = graphs.topology_inventory()
    assert len(rows) == 7
    for M, G, extra, ends, Iphi, Ih in rows:
        assert M + G <= 4 and M >= 2 and ends >= M
        assert 3 * G + extra + ends == 2 * Ih
        assert Iphi == M - 2 and Ih == G + 2
        assert Iphi + Ih - (M + G) + 1 == 1


def test_whole_original_frontier_and_matching_not_overwritten():
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 144
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())


def test_parent_packet_cache_warmth_cannot_change_new_residuals():
    before = dict(audit.residuals())
    audit.previous.packets()
    assert audit.residuals() == before
    assert source.data()["checks"] is not source.previous.data()["checks"]
    base.certify_residuals(before)


def test_three_finite_anchors_are_not_assigned():
    text = audit.observable()["not_established"]
    assert "finite values or signs" in text.lower()
    assert "Regge" in text
    assert "P8" in text
    assert "light_Newton_and_regular_local_anchors_not_chosen" in poles.data()["gates"]


def test_integer_math_apis_preserve_exact_rationals():
    values = (
        poles.eikonal(6, 1, 4),
        poles.tree_numerator(6, -1, -1, 1, 4),
        poles.residue_ratio(4),
        *poles.completion_coefficients(4),
        *poles.completion_basis(6, -1, 1),
        poles.completion(6, -1, 1, 4),
        poles.newton_shape(6, -1, 1, 4),
        soft.density(-1, 1, 0),
        soft.attenuation_density(4, 1, 1, 0),
        soft.coulomb_eta(8, 1, 1),
    )
    assert all(
        isinstance(value, s.Basic) and not value.has(s.Float) for value in values
    )
