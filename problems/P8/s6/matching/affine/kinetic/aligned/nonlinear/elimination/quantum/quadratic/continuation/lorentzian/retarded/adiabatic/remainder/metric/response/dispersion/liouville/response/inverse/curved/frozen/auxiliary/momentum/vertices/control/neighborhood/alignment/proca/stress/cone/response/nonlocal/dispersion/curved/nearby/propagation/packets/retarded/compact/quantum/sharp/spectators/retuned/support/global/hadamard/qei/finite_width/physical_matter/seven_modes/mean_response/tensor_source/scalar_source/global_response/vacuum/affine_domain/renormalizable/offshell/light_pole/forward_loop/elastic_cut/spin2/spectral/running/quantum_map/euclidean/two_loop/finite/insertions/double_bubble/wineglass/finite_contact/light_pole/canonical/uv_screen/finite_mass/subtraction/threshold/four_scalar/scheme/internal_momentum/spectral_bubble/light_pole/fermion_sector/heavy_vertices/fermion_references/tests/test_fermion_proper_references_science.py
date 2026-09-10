"""Independent Dirac, finite-regulator, moment and reference tests."""

import itertools

import mpmath as mp
import pytest
import sympy as sp
from p8_vacuum_fermion_proper_references import (
    anchors,
    audit,
    bounds,
    conversion,
    dirac,
)


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_every_exact_identity_and_matrix_entry(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, sp.MatrixBase) else [value]
    assert all(entry == 0 for entry in entries)


@pytest.mark.parametrize("name", tuple(audit.gates()))
def test_every_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=lambda v: v if isinstance(v, str) else None
)
def test_all_unsupported_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "m,q", tuple(itertools.product((2, 3), ((0, 0, 0, 0), (1, 2, 3, 4), (-2, 1, 0, 3))))
)
def test_exact_resolvent_mass_derivative_is_negative_square(m, q):
    gammas = dirac.gamma_matrices()
    mass = sp.Symbol("independent_mass", positive=True)
    slash = sum((a * b for a, b in zip(q, gammas)), sp.zeros(4))
    S = (mass * sp.eye(4) - sp.I * slash) / (mass**2 + sum(v * v for v in q))
    assert (S.diff(mass) + S * S).subs(mass, m).applyfunc(sp.simplify) == sp.zeros(4)


@pytest.mark.parametrize(
    "rvalue",
    (
        sp.Rational(1, 4),
        sp.Rational(1, 16),
        sp.Rational(1, 2),
        sp.Rational(9, 10),
        sp.Rational(1, 1000000),
    ),
)
def test_closed_moments_against_independent_quadrature(rvalue):
    data = anchors.data()
    with mp.workdps(50):
        r = mp.mpf(int(rvalue.p)) / int(rvalue.q)
        values = (
            mp.quad(lambda x: mp.log(r + (1 - r) * x), [0, 1]),
            mp.quad(lambda x: (1 - x) * mp.log(r + (1 - r) * x), [0, 1]),
            mp.quad(lambda x: x / (r + (1 - r) * x), [0, 1]),
        )
        for key, target in zip(("J0_closed", "J1_closed", "R_closed"), values):
            expr = data[key]
            actual = expr.subs(
                next(s for s in expr.free_symbols if str(s) == "r"), rvalue
            ).evalf(50)
            assert mp.almosteq(mp.mpf(str(actual)), target)
        assert -1 <= values[0] <= 0 and -mp.mpf(3) / 4 <= values[1] <= 0
        assert mp.mpf(1) / 2 <= values[2] <= 1


@pytest.mark.parametrize(
    "key,at_zero,at_one",
    (
        ("J0_closed", -1, 0),
        ("J1_closed", -sp.Rational(3, 4), 0),
        ("R_closed", 1, sp.Rational(1, 2)),
    ),
)
def test_mass_ratio_endpoint_limits(key, at_zero, at_one):
    expr = anchors.data()[key]
    r = next(s for s in expr.free_symbols if str(s) == "r")
    assert sp.limit(expr, r, 0, dir="+") == at_zero
    assert sp.limit(expr, r, 1) == at_one


def test_early_four_dimensional_gauge_numerator_loses_finite_constants():
    eps = sp.Symbol("eps")
    correct_mass = sp.expand((4 - 2 * eps) * (1 / eps + 1) - 4 / eps).coeff(eps, 0)
    wrong_mass = sp.expand(4 * (1 / eps + 1) - 4 / eps).coeff(eps, 0)
    correct_kin = sp.expand(
        (2 - 2 * eps) * (1 / (2 * eps) + sp.Rational(3, 4)) - 1 / eps
    ).coeff(eps, 0)
    assert correct_mass == 2 and wrong_mass == 4 and correct_kin == sp.Rational(1, 2)


def test_scale_must_remain_fixed_during_scalar_background_derivative():
    m, mu = sp.symbols("m mu", positive=True)
    mass = m * (2 - 4 * sp.log(m * m / mu**2))
    assert sp.diff(mass, m).subs(mu, m) == -6
    assert sp.diff(mass.subs(mu, m), m) == 2


@pytest.mark.parametrize(
    "j0,j1,R,w",
    tuple(
        itertools.product(
            (-1, 0),
            (-sp.Rational(3, 4), 0),
            (sp.Rational(1, 2), 1),
            (-sp.Rational(1, 2), sp.Rational(1, 2)),
        )
    ),
)
def test_selected_canonical_bounds_at_independent_anchor_box_corners(j0, j1, R, w):
    Y, a, Q = sp.Rational(1, 100), sp.Rational(1, 200), sp.Integer(144)
    Cf = sp.Rational(4, 3)
    z = (-Y * j1 + a * Cf / 2) / Q
    eta = (Y * j0 + 2 * a * Cf) / Q
    ups = (Y * (j0 + 2 * R) - 6 * a * Cf) / Q
    result = bounds.enclosure(2, Y, a, sp.Rational(1, 2), Q)
    assert (
        abs((1 + eta) / (1 + z) - 1)
        <= result["selected_mass_ratio_minus_one_absolute_upper"]
    )
    value = (1 + ups) / ((1 + z) * sp.sqrt(1 + w))
    assert abs(value - 1) <= result["selected_Yukawa_ratio_minus_one_absolute_upper"]
    assert value > 0


def test_large_valid_enclosure_is_not_falsely_positive():
    d = bounds.enclosure(2, 1000, 1000, 0, 1)
    assert d["positive_selected_mass_reference"] is False
    assert d["positive_selected_Yukawa_reference"] is False


def test_inert_scalar_vertex_zero_but_gauge_mass_not_zero():
    data = conversion.data()
    assert data["checks"]["inert_flavor_scalar_Yukawa_vertex_is_zero"] == 0
    mass = data["finite_mass_relative_eta"]
    Y = next(s for s in mass.free_symbols if str(s) == "Y")
    assert mass.subs(Y, 0) != 0


def test_scope_and_counts():
    assert len(audit.residuals()) == 76
    assert audit.scalar_entry_count() == 357
    assert len(audit.gates()) == 30
    assert audit.rejected_inputs() == 68
    assert len(audit.controls()) == 9
    assert "off-shell" in anchors.data()["scope"].lower()
    assert "not" in conversion.data()["scope"].lower()
