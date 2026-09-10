"""Independent finite-sector, elementary integral and actual-bound checks."""

from fractions import Fraction
from math import factorial

import pytest
import sympy as sp
from p8_vacuum_two_loop_finite import (
    audit,
    calibration,
    ordered,
    parametric,
    sectors,
    selection,
)


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_selection_integral_or_calibration_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_continuous_integral_domain_or_scope_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[row[0] for row in calibration.bad_cases()],
)
def test_reject_unsupported_or_subtraction_dependent_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_control_or_original_scope(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "kind,choices",
    selection.cases(),
    ids=[
        kind + "_" + "".join(map(str, choices)) for kind, choices in selection.cases()
    ],
)
def test_every_finite_graph_sector_counts_and_positive_log_polynomial(kind, choices):
    d = sectors.data(kind, choices)
    h = d["heavy_edge_count"]
    assert d["edge_order_sector_count"] == factorial(4 + h)
    assert d["minimum_small_tail_exponent"] > 0
    assert d["maximum_log_power"] == 2
    assert all(
        isinstance(c, Fraction) and c > 0 for c in d["box_log_polynomial_coefficients"]
    )
    assert len(d["complete_sector_inventory_sha256"]) == 64
    p = calibration.polynomial(kind, choices, 600)
    assert p["scaled_anisotropic_box_integral_upper"] == sum(
        sp.Rational(c.numerator, c.denominator) * 600**j
        for j, c in enumerate(d["box_log_polynomial_coefficients"])
    )


@pytest.mark.parametrize(
    "kind,choices",
    selection.cases(),
    ids=[
        kind + "_" + "".join(map(str, choices)) for kind, choices in selection.cases()
    ],
)
def test_independent_sector_piece_formula_in_both_extreme_orders(kind, choices):
    graph = selection.require_finite(kind, choices)
    size = len(graph["edges"])
    for order in (tuple(range(size)), tuple(reversed(range(size)))):
        d = sectors.single(kind, choices, order)
        p = d["monomial_powers"]
        assert p.count(-2) == 2
        assert sum(x + 1 for x in p) == size - 4
        for piece in d["pieces"]:
            k = piece["large_light_count"]
            assert all(e < 4 for e in order[:k])
            prefix = piece["prefix_exponents"]
            tails = piece["small_tail_exponents"]
            assert all(a > 0 for a in tails)
            assert all(a <= 0 for a in prefix)
            denom = factorial(prefix.count(0))
            for a in tails:
                denom *= a
            for a in prefix:
                if a < 0:
                    denom *= -a
            assert piece["rational_coefficient"] == Fraction(1, denom)
            # Independent nested integral for this actual small pattern.
            small = ordered.small_data(p[k:])
            b = sp.Symbol("positive_small_upper", positive=True)
            assert small["exact_small_ordered_integral"].subs(b, 1) == sp.Rational(
                1, sp.prod(tails)
            )


@pytest.mark.parametrize("powers", ordered.LARGE_PATTERNS)
def test_large_ordered_integral_anchor_and_named_upper(powers):
    d = ordered.large_data(powers)
    B = d["B"]
    assert d["exact_ordered_integral"].subs(B, 1) == (1 if not powers else 0)
    assert d["log_power"] <= 2
    assert all(v == 0 for v in d["checks"].values())
    # A diagnostic evaluation, not the continuous positivity proof.
    assert (
        sp.N((d["positive_log_upper"] - d["exact_ordered_integral"]).subs(B, 64), 50)
        >= 0
    )


@pytest.mark.parametrize("powers", ordered.small_patterns())
def test_all_supported_small_ordered_integrals(powers):
    d = ordered.small_data(powers)
    assert sp.simplify(d["exact_small_ordered_integral"] - d["closed_formula"]) == 0
    assert all(a > 0 for a in d["positive_tail_exponents"])


@pytest.mark.parametrize("h", (1, 2, 3))
def test_two_loop_layer_factor_and_coupling_scaling(h):
    d = parametric.data()["per_heavy_edge_count"][h]
    assert d["edge_count"] == 4 + h
    assert d["positive_layer_factor"] == factorial(h + 2)
    assert d["fourfold_parameter_scaling"] == 4 ** (h + 2)
    assert (
        d["rational_coupling_and_loop_prefactor"]
        == sp.Rational(17 * factorial(h + 2) * 16, 41472) * sp.Rational(2, 3) ** h
    )


@pytest.mark.parametrize("ell", (0, Fraction(1, 2), 1, 600))
def test_logarithm_polynomial_conditional_bound_api(ell):
    d = calibration.polynomial("double_bubble", (0, 0, 2), ell)
    assert d["log_upper"] == sp.Rational(ell)
    assert (
        d["scaled_anisotropic_box_integral_upper"]
        == 38 + 28 * sp.Rational(ell) + 4 * sp.Rational(ell) ** 2
    )


def test_cached_validation_does_not_accept_equal_inexact_values():
    calibration.polynomial("double_bubble", (0, 0, 2), 1)
    for bad in (True, 1.0, sp.Float(1)):
        with pytest.raises((TypeError, ValueError)):
            calibration.polynomial("double_bubble", (0, 0, 2), bad)
    sectors.single("double_bubble", (0, 0, 2), (0, 1, 2, 3, 4))
    with pytest.raises((TypeError, ValueError)):
        sectors.single("double_bubble", (0, 0, 2), (False, 1, 2, 3, 4))


def test_actual_rational_error_bound_and_nonzero_log_squared_term():
    d = calibration.data()
    E = d["actual_finite_subsector_b2_absolute_upper"]
    assert d["aggregate_log_polynomial_coefficients"] == (
        sp.Rational(13688587, 17496),
        sp.Rational(515831, 972),
        sp.Rational(59789, 972),
    )
    assert d["aggregate_constant_at_log_upper"] == sp.Rational(393017383387, 17496)
    assert 0 < E < sp.Rational(3, 10**607)
    assert 0 < E / d["actual_tree_b2"] < sp.Rational(1, 10**7)
    assert d["actual_heavy_box_scale"] < 10**200
    assert d["elementary_positive_exp_three_partial_sum"] > 10


def test_selection_and_integration_counts():
    d = selection.data()
    assert d["finite_refinement_count"] == 88
    assert d["excluded_subtraction_dependent_count"] == 104
    assert sum(r["edge_order_sector_count"] for r in sectors.summaries()) == 258480
    assert sum(r["heavy_cutoff_piece_count"] for r in sectors.summaries()) == 526320


def test_integrated_subset_is_not_a_complete_two_loop_amplitude():
    assert (
        "not the signed correction or the complete two-loop coefficient"
        in calibration.data()["scope"]
    )
    assert "Counterterm diagrams" in selection.data()["scope"]
    assert audit.controls()["complete_two_loop_and_original_P8_not_closed"]


def test_exact_counts():
    assert len(audit.residuals()) == 171
    assert len(audit.gates()) == 32
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 124
