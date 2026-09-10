"""Independent Gaussian, noncommuting, color and exact-frontier regressions."""

import copy
from fractions import Fraction as F

import pytest
import sympy as sp
from p8_vacuum_two_loop_fermion_ledger import (
    audit,
    boundary,
    counterterms,
    families,
    gauge,
    legendre,
)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_named_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_every_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[v[0] for v in audit.bad_cases()]
)
def test_every_unsupported_degree_or_frontier(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def pairings(items):
    if not items:
        return [()]
    return [
        ((items[0], items[j]),) + rest
        for j in range(1, len(items))
        for rest in pairings(items[1:j] + items[j + 1 :])
    ]


@pytest.mark.parametrize("power", (0, 2, 4, 6, 8))
@pytest.mark.parametrize("a", (F(1), F(3, 2), F(7)))
def test_gaussian_moments_from_independent_Wick_pairings(power, a):
    count = len(pairings(tuple(range(power))))
    expected = F(count) / a ** (power // 2)
    assert legendre.gaussian_moment(power, sp.Rational(a)) == sp.Rational(expected)


def test_exact_zero_dimensional_fermion_determinant_with_source():
    # A pair of opposite Yukawas gives det=m^2-y^2*x^2.
    h, phi = sp.symbols("h phi")
    a, m, y = sp.symbols("a m y", positive=True)
    delta = m * m - y * y * phi * phi
    determinant_gaussian_average = delta - h * y * y / a
    Fdet = -sp.log(delta)
    W2 = sp.diff(sp.log(determinant_gaussian_average), h).subs(h, 0)
    Gamma2 = -W2 + sp.diff(Fdet, phi) ** 2 / (2 * a)
    assert sp.factor(Gamma2 - sp.diff(Fdet, phi, 2) / (2 * a)) == 0
    # Omitting the mean-background shift would retain the wrong result.
    assert sp.factor(-W2 - sp.diff(Fdet, phi, 2) / (2 * a)) != 0


@pytest.mark.parametrize("size", (2, 3))
@pytest.mark.parametrize("degree", (0, 2, 4))
def test_noncommuting_background_coefficients_from_exact_matrix_inverse(size, degree):
    z = sp.symbols("z")
    D = sp.Matrix(
        [[sp.Integer(3 if i == j else 1) for j in range(size)] for i in range(size)]
    )
    W = sp.Matrix([[sp.Integer(i + j + 1) for j in range(size)] for i in range(size)])
    F0 = sp.diag(*range(1, size + 1))
    F2 = sp.Matrix(
        [[sp.Integer(2 if i == j else -1) for j in range(size)] for i in range(size)]
    )
    F4 = sp.eye(size)
    exact = sp.trace((D.inv() + z * z * W).inv() * (F0 + z * z * F2 + z**4 * F4)) / 2
    coefficients = {
        0: sp.trace(D * F0) / 2,
        2: sp.trace(D * F2 - D * W * D * F0) / 2,
        4: sp.trace(D * F4 - D * W * D * F2 + D * W * D * W * D * F0) / 2,
    }
    actual = sp.diff(exact, z, degree).subs(z, 0) / sp.factorial(degree)
    assert sp.factor(actual - coefficients[degree]) == 0


def test_illegal_commuting_of_operator_factors_changes_the_trace():
    D = sp.diag(1, 2)
    W = sp.Matrix([[0, 1], [1, 0]])
    F0 = sp.eye(2)
    right = sp.trace(D * W * D * W * D * F0) / 2
    wrong = sp.trace(D**3 * W**2 * F0) / 2
    assert right == 3 and wrong == sp.Rational(9, 2)


@pytest.mark.parametrize("degree,count", ((0, 2), (2, 3), (4, 4)))
def test_every_background_degree_has_all_scalar_and_gauge_rows(degree, count):
    rows = boundary.rows_at_degree(degree)
    assert len(rows) == count
    assert len([r for r in rows if r["contracted_boson"] == "gauge"]) == 1
    scalar = [r for r in rows if r["contracted_boson"] == "Phi"]
    assert {r["tree_W_insertions"] for r in scalar} == set(range(degree // 2 + 1))
    for r in scalar:
        assert (
            r["fermion_Hessian_background_degree"] + 2 * r["tree_W_insertions"]
            == degree
        )
        assert r["scalar_Y_power"] == r["fermion_Hessian_background_degree"] // 2 + 1
        assert r["trace_coefficient"] == sp.Rational((-1) ** r["tree_W_insertions"], 2)
    g = next(r for r in rows if r["contracted_boson"] == "gauge")
    assert g["gauge_a_power"] == 1 and g["scalar_Y_power"] == degree // 2


@pytest.mark.parametrize("index", range(8))
def test_literal_generators_are_traceless_Hermitian(index):
    T = gauge.data()["fundamental_generators"][index]
    assert T == T.conjugate().T
    assert sp.trace(T) == 0
    assert sp.trace(T * T) == sp.Rational(1, 2)


def test_active_gauge_factor_is_eight_not_six():
    d = gauge.data()
    total = sum(sp.trace(T * T) for T in d["fundamental_generators"])
    assert total == 4
    assert 2 * total == d["gauge_active_color_flavor_factor"] == 8
    assert 14 * total == d["gauge_all_flavor_vacuum_color_factor"] == 56
    assert d["scalar_active_color_flavor_factor"] == 6
    assert (
        d["gauge_active_color_flavor_factor"] != d["scalar_active_color_flavor_factor"]
    )


@pytest.mark.parametrize("degree", (1, 2, 3, 4, 5, 6))
def test_opposite_Yukawa_flavor_projection_independently(degree):
    y = sp.symbols("y")
    all_flavors = [y, -y] + [sp.Integer(0)] * 12
    expected = 0 if degree % 2 else 2 * y**degree
    assert sp.expand(sum(v**degree for v in all_flavors) - expected) == 0


@pytest.mark.parametrize(
    "a0,a1,a2,p1,p2,f0",
    ((2, 3, 5, 7, 11, 13), (1, -2, 3, 4, 5, -6), (0, 1, 1, 1, 1, 0)),
)
def test_regulator_sensitive_finite_products(a0, a1, a2, p1, p2, f0):
    eps = sp.symbols("eps")
    expression = (a0 + eps * a1 + eps * eps * a2) * (p2 / eps**2 + p1 / eps + f0)
    finite = sp.expand(expression).coeff(eps, 0)
    assert finite == a0 * f0 + a1 * p1 + a2 * p2
    early_projection = sp.expand(a0 * (p2 / eps**2 + p1 / eps + f0)).coeff(eps, 0)
    assert finite - early_projection == a1 * p1 + a2 * p2
    assert finite != early_projection


def test_assigned_inner_counterterm_is_included_exactly_once():
    Fraw, CT, rest = sp.symbols("Fraw CT rest", commutative=False)
    raw = Fraw + CT + rest
    paired = (Fraw + CT) + rest
    wrong = paired + CT
    assert sp.expand(raw - paired) == 0
    assert sp.expand(wrong - raw) == CT
    assert sp.expand(wrong - raw) != 0


def test_second_order_external_field_and_cubic_products():
    h, G, G1, G2, z1, z2, A0, A1, A2 = sp.symbols("h G G1 G2 z1 z2 A0 A1 A2")
    square = sp.expand((G + h * G1 + h * h * G2) ** 2).coeff(h, 2)
    assert square == G1 * G1 + 2 * G * G2
    amp = (A0 + h * A1 + h * h * A2) / (1 + h * z1 + h * h * z2) ** 2
    second = sp.diff(amp, h, 2).subs(h, 0) / 2
    assert sp.factor(second - A2 + 2 * z1 * A1 - (3 * z1 * z1 - 2 * z2) * A0) == 0


@pytest.mark.parametrize("index", range(9))
def test_no_primitive_row_can_be_silently_deleted(index):
    rows = boundary.ownership_frontier()
    del rows[index]
    with pytest.raises(ValueError, match="differs"):
        boundary.validate_frontier(rows)


def test_falsely_complete_frontier_is_rejected():
    rows = copy.deepcopy(boundary.ownership_frontier())
    for row in rows:
        row["status"] = "BOUNDED_PAIRED_SECTOR"
    with pytest.raises(ValueError, match="differs"):
        boundary.validate_frontier(rows)


def test_exact_frontier_retains_one_full_one_partial_and_seven_unevaluated():
    rows = boundary.ownership_frontier()
    assert boundary.validate_frontier(rows) is True
    assert sum(r["status"] == "BOUNDED_PAIRED_SECTOR" for r in rows) == 1
    assert sum(r["status"] == "LOCAL_QUARTIC_SUBSET_ONLY" for r in rows) == 1
    assert sum(r["status"] == "UNEVALUATED" for r in rows) == 7
    assert len({r["id"] for r in rows}) == 9
    assert len(boundary.data()["separate_uncomputed_counterterm_tasks"]) == 3


def test_returned_frontier_does_not_mutate_the_catalog():
    rows = boundary.ownership_frontier()
    rows[0]["status"] = "changed outside"
    assert boundary.ownership_frontier()[0]["status"] == "UNEVALUATED"
    assert len(families.catalog()) == 9


def test_scope_excludes_integral_bounds_and_finite_matching():
    assert "not evaluation or an error bound" in legendre.data()["scope"]
    assert "not completed integrals" in families.data()["scope"]
    assert "not supplied" in counterterms.data()["scope"]
    assert "Finite-gravity" in gauge.data()["scope"]
    assert "seven primitive rows" in boundary.data()["nonclosure"]
    assert audit.controls()["original_P8_not_closed"] is True


def test_exact_counts():
    assert len(audit.residuals()) == 141
    assert len(audit.gates()) == 28
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 39
