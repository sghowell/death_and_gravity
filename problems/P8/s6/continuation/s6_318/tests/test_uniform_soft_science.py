"""All-order structural checks, exact majorants and failed unweighted controls."""

import pytest
import sympy as s
from p8_vacuum_affine_uniform_soft_current_bound import (
    audit,
    bounds,
    geometry,
    grading,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_all_exact_residuals(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_all_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_original_scope_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("n", range(1, 11))
def test_all_finite_majorant_calibrations(n):
    C = bounds.coefficient(n)
    assert C.is_Integer and C > 0
    assert C <= bounds.envelope(n)
    assert C == bounds.majorant_data()["whole_first_ten_majorant_coefficients"][n]


@pytest.mark.parametrize(
    "invalid",
    (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        0,
        -1,
        s.Rational(3, 2),
        s.oo,
        s.I,
        s.nan,
        s.Symbol("unknown"),
        [],
        {},
    ),
)
def test_positive_multiplicity_validation_precedes_cache(invalid):
    assert bounds.coefficient(1) == 1
    for call in (bounds.coefficient, bounds.envelope):
        with pytest.raises((TypeError, ValueError)):
            call(invalid)


def test_arbitrary_Rosen_family_is_not_assumed_to_solve_the_uu_equation():
    record = geometry.rosen()
    assert record["whole_vanishing_action_contraction"] == 0
    assert record["whole_Ricci_uu"] != 0
    assert all(
        value == 0
        for key, value in record["checks"].items()
        if "Ricci_" in key and not key.endswith("uu_Ricci")
    )


@pytest.mark.parametrize("valence", (3, 4, 5))
@pytest.mark.parametrize("kind,minimum", (("TT", 2), ("TL", 1), ("LL", 0)))
def test_low_valence_parity_and_required_angular_orders(valence, kind, minimum):
    row = grading.data()["whole_independent_graded_vertices"][f"EH{valence}_{kind}"]
    assert row["nonzero"] and row["minimum_degree"] == minimum


def test_plain_Frobenius_control_without_longitudinal_weights_is_insufficient():
    assert grading.data()["whole_wrong_unweighted_child_defect"] == 144


def test_exact_barrier_proves_more_than_the_first_ten_samples():
    record = bounds.majorant_data()["whole_exact_cap_radius_barrier"]
    cap, radius = record["cap"], record["radius"]
    assert radius == cap / 2
    assert 0 < 416 * cap < 1
    assert (
        record["margin"]
        == cap - radius - 2048 * ((1 - 416 * cap) ** -2 - 1 - 832 * cap)
        > 0
    )


@pytest.mark.parametrize(
    "label,expected", (("0000", s.Rational(2509, 225)), ("0101", -s.Rational(858, 25)))
)
def test_complete_four_ray_weighted_limits_at_both_collinear_ends(label, expected):
    row = bounds.four_ray_calibration()["whole_two_independent_26_tree_families"][label]
    assert row["endpoint_values"]["0"] == row["endpoint_values"]["oo"]
    assert row["endpoint_values"]["0"][-1] == expected
    assert all(value.is_finite is True for value in row["endpoint_values"]["0"])


def test_explicit_energy_and_canonical_scaling():
    first = (s.Rational(1, 64), s.Rational(1, 32))
    assert bounds.energy_majorant(first, 16) == bounds.coefficient(2) * s.Rational(9, 8)
    assert bounds.energy_majorant(
        tuple(2 * w for w in first), 16
    ) == bounds.energy_majorant(first, 16)
    triple = (s.Rational(1, 64), s.Rational(1, 32), s.Rational(3, 64))
    assert bounds.energy_majorant(triple, 16) == s.Rational(9, 4) * bounds.coefficient(
        3
    )
    assert bounds.energy_majorant(
        triple, source.KAPPA
    ) * source.KAPPA == 36 * bounds.coefficient(3)


@pytest.mark.parametrize(
    "energies",
    (
        [],
        (),
        None,
        "1",
        (True,),
        (1.0,),
        (s.Float(1),),
        (0,),
        (-1,),
        (s.Rational(1, 4),),
        (s.oo,),
        ([1],),
    ),
)
def test_energy_domain_and_exactness_rejections(energies):
    with pytest.raises((TypeError, ValueError)):
        bounds.energy_majorant(energies)


@pytest.mark.parametrize(
    "invalid", (True, False, 1.0, s.Float(1), "1", None, 0, -1, s.oo, s.I)
)
def test_coupling_rejections(invalid):
    with pytest.raises((TypeError, ValueError)):
        bounds.energy_majorant((s.Rational(1, 64),), invalid)


def test_singleton_majorant_and_first_nonlinear_counts():
    assert bounds.energy_majorant((s.Rational(1, 64),)) == 1
    assert bounds.coefficient(2) == 2126512128
    assert bounds.coefficient(3) == 13566165030109446144


def test_all_original_frontiers_and_nonclosure_remain():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 174
    assert len(audit.frontier()) == 9 and len(audit.qualifications()) == 6
    scope = audit.observable()["not_established"]
    assert all(
        word in scope for word in ("hard", "inclusive", "Regge", "Energy-independent")
    )
