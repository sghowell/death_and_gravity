"""Independent finite-index, simplex, spectral-tail and scope tests."""

import copy
import itertools

import mpmath as mp
import pytest
import sympy as sp
from p8_vacuum_fermion_heavy_vertices import (
    audit,
    bounds,
    calibration,
    ownership,
    parametric,
    renormalization,
)


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_every_exact_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", tuple(audit.gates()))
def test_every_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args",
    audit.bad_cases(),
    ids=lambda value: value if isinstance(value, str) else None,
)
def test_every_unsupported_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "v,M,s",
    tuple(
        itertools.product(
            (16, 64),
            (24, 100),
            (sp.Integer(2), 2 + sp.Rational(3, 5) + sp.I * sp.Rational(4, 5), 2 - sp.I),
        )
    ),
)
def test_full_complex_parameter_gap_on_independent_simplex_grid(v, M, s):
    for a0, a1, a2 in itertools.product(range(5), repeat=3):
        if a0 + a1 + a2 > 4:
            continue
        a0, a1, a2 = (sp.Rational(a, 4) for a in (a0, a1, a2))
        a3 = 1 - a0 - a1 - a2
        H = a2 + a3
        for m0, m1 in ((v, 1), (1, v)):
            Delta = (
                m0 * a0
                + m1 * a1
                + M * H
                - s * a0 * a1
                - (4 - s) * a2 * a3
                - (a0 + a1) * H
            )
            marked = a0 if m0 == v else a1
            lower = sp.Rational(1, 4) + (M - 1) * H + (v - 1) * marked
            assert sp.re(Delta) >= lower > 0


@pytest.mark.parametrize("v", (16, 17, 64, 256, 10000, 1000000))
def test_direct_simplex_integrals_below_their_analytic_majorants(v):
    with mp.workdps(50):
        vv = mp.mpf(v)
        tri = mp.quad(lambda x: (1 - x) / (mp.mpf(1) / 4 + (vv - 1) * x), [0, 1])
        box = mp.quad(
            lambda x: (1 - x) ** 2 / (2 * (mp.mpf(1) / 4 + (vv - 1) * x) ** 2), [0, 1]
        )
        assert 0 < tri < 2 * mp.log(4 * vv) / vv
        assert 0 < box < 4 / vv


@pytest.mark.parametrize("T", (16, 64, 10000))
def test_infinite_spectral_tail_by_independent_compactification(T):
    with mp.workdps(50):
        tt = mp.mpf(T)
        transformed = mp.quad(lambda x: mp.log(4 * tt / x) / tt, [0, 1])
        exact = (mp.log(4 * tt) + 1) / tt
        assert mp.almosteq(transformed, exact)
        tri_bound = 8 * transformed
        assert mp.almosteq(tri_bound, 8 * (mp.log(4 * tt) + 1) / tt)


@pytest.mark.parametrize("n", (1, 2, 6, 32, 100, 1333))
@pytest.mark.parametrize("offset", (-1, 0, 1))
def test_dyadic_integer_bound_at_rational_boundaries(n, offset):
    value = sp.Integer(2) ** n + sp.Rational(offset, 3)
    got = bounds.dyadic_log_upper(value)
    expected = n + (offset > 0)
    assert got == expected
    assert 2**got >= value
    assert got == 0 or 2 ** (got - 1) < value


@pytest.mark.parametrize("value", (sp.Rational(3, 2), 2, 3, 16, sp.Rational(257, 8)))
def test_dyadic_bound_dominates_high_precision_log(value):
    n = bounds.dyadic_log_upper(value)
    with mp.workdps(50):
        x = (
            mp.mpf(int(value.p)) / int(value.q)
            if isinstance(value, sp.Rational)
            else mp.mpf(value)
        )
        assert mp.log(x) < int(n)


def test_counterterm_cannot_be_replaced_by_a_pure_quartic_contact():
    d = renormalization.data()
    value = d["finite_nonlocal_reference_change_b2"]
    sub = {
        symbol: {
            "quartic_L": 2,
            "cubic_squared_g": 1,
            "M": 24,
            "finite_change_in_paired_bubble_reference": 1,
        }[str(symbol)]
        for symbol in value.free_symbols
    }
    assert value.subs(sub) != 0


def test_full_vertex_contains_all_nine_terms():
    terms = ownership.data()["vertex_product_terms"]
    assert len(terms) == 9
    assert [kind for kind, _ in terms].count("triangle") == 4
    assert [kind for kind, _ in terms].count("box") == 4
    assert len({str(term) for _, term in terms}) == 9


def test_noncommuting_covariances_must_not_be_collapsed():
    D = sp.Matrix([[2, 1], [1, 3]])
    P = sp.Matrix([[5, 2], [2, 1]])
    W = sp.Matrix([[1, 2], [2, 4]])
    assert sp.trace(P * W * D * W) != sp.trace(P * D * W**2)


@pytest.mark.parametrize("index", range(9))
def test_no_primitive_id_may_be_changed(index):
    rows = audit.frontier()
    rows[index]["id"] = "not_the_derived_family"
    with pytest.raises(ValueError):
        audit.validate_frontier(rows)


@pytest.mark.parametrize("index", range(9))
def test_no_primitive_row_may_be_removed(index):
    rows = audit.frontier()
    del rows[index]
    with pytest.raises(ValueError):
        audit.validate_frontier(rows)


def test_exactly_one_frontier_row_changes_and_old_source_is_unchanged():
    old = audit.parent.ownership_frontier()
    new = audit.frontier()
    changed = [(a, b) for a, b in zip(old, new) if a != b]
    assert len(changed) == 1
    assert changed[0][0]["id"] == "scalar_Phi4_W2_F0"
    assert changed[0][0]["status"] == "LOCAL_QUARTIC_SUBSET_ONLY"
    assert changed[0][1]["status"] == "BOUNDED_IN_DECLARED_PARENT_SUBTRACTION"
    assert sum(r["status"] == "UNEVALUATED" for r in new) == 7
    copy_of_old = copy.deepcopy(old)
    new[0]["status"] = "arbitrary mutation"
    assert audit.parent.ownership_frontier() == copy_of_old


def test_complete_family_is_a_sum_not_an_extra_copy_of_the_subset():
    result = bounds.enclosure(2, 1, 2, 3, 24, 144)
    assert result["full_family_disc_amplitude_and_b2_upper"] == sum(
        result[k]
        for k in (
            "all_channel_bubble_amplitude_upper",
            "all_channel_triangle_amplitude_upper",
            "all_channel_box_amplitude_upper",
        )
    )


def test_actual_relative_bound_and_exact_log_cover():
    d = calibration.data()
    r = d["literal_full_family_enclosure"]
    assert r["dyadic_exponent"] == 1333
    assert r["log_4T_plus_one_upper"] == 1334
    assert 0 < d["full_family_upper_relative_to_tree"] < sp.Rational(1, 10**22)


def test_absence_of_yukawa_coupling_removes_every_term():
    result = bounds.enclosure(2, 0, 1, 1, 24, 144)
    for key in (
        "all_channel_bubble_amplitude_upper",
        "all_channel_triangle_amplitude_upper",
        "all_channel_box_amplitude_upper",
        "full_family_disc_amplitude_and_b2_upper",
    ):
        assert result[key] == 0


def test_absence_of_cubic_coupling_removes_heavy_vertex_products():
    result = bounds.enclosure(2, 1, 1, 0, 24, 144)
    assert result["all_channel_triangle_amplitude_upper"] == 0
    assert result["all_channel_box_amplitude_upper"] == 0
    assert result["all_channel_bubble_amplitude_upper"] > 0


def test_scope_leaves_finite_matching_and_other_primitives_open():
    assert "unevaluated" in renormalization.data()["canonical_conversion"]
    assert "not the complete" in calibration.data()["scope"]
    assert "No sign" in bounds.data()["scope"]
    assert "1/Q" in parametric.data()["scope"]


def test_expected_check_counts():
    assert len(audit.residuals()) == 60
    assert len(audit.gates()) == 29
    assert audit.rejected_inputs() == 106
    assert len(audit.controls()) == 9
