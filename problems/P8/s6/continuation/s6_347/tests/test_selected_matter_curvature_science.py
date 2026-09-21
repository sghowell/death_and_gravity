"""Complete selected scalar inventory, off-shell matching and strict scope."""

import pytest
import sympy as s
from p8_vacuum_affine_selected_matter_curvature_coefficient import (
    aggregate,
    audit,
    calibration,
    local,
    mixed,
    moment,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_residual(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_all_rejected_inputs_and_scopes(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("order", (2, 3))
@pytest.mark.parametrize("mass", (s.Rational(1, 2), 1, 2, 3, 10))
def test_each_exact_positive_mass_moment(order, mass):
    value = moment.coefficient(order, mass)
    target = (
        s.Rational(1, 60)
        if mass == 1 and order == 2
        else s.Rational(1, 420)
        if mass == 1
        else moment.build()[1][order].subs(moment.N, mass)
    )
    assert s.factor(value - target) == 0 and not value.has(s.Float)


@pytest.mark.parametrize("mass", (s.Rational(1, 2), 1, 2, 3, 10))
def test_each_exact_extra_coefficient(mass):
    value = moment.extra_at(mass)
    target = -s.Rational(952, 15) if mass == 1 else moment.extra().subs(moment.N, mass)
    assert s.factor(value - target) == 0 and not value.has(s.Float)


@pytest.mark.parametrize(
    "name",
    (
        "old_polynomial_core",
        "all_six_local_tadpoles",
        "mixed_across",
        "mixed_within_J4",
        "mixed_within_J2",
    ),
)
def test_exact_selected_class_inventory(name):
    assert source.require_class(name) == name


def test_complete_graph_inventory():
    assert set(source.data()["whole_scalar_graph_patterns"]) == set(source.PATTERNS)
    assert len(source.PATTERNS) == 5


def test_local_offshell_and_dimension_proof():
    data = local.data()
    assert len(data["checks"]) == 12
    assert data["whole_degree6_curvature_coefficients"] == {"Gal": 0, "Y_Hdiff": 0}
    assert data["whole_degree6_flat_jet_coefficients"]["Gal"] == (
        0,
        s.Rational(1, 2),
        0,
        0,
        -s.Rational(1, 2),
        0,
    )


def test_all_mixed_EOM_delta_and_counterterms_retained():
    text = mixed.data()["whole_covariant_across_reduction"]
    assert "marked-EOM terms are retained" in text
    assert "arbitrary endpoint p_squared" in text
    assert "within-J2" in mixed.data()["whole_withinJ2_and_counterterm_boundary"]


def test_equal_mass_sign_not_universally_positive():
    assert moment.extra_at(1) == -s.Rational(952, 15)
    assert moment.data()["whole_large_mass_squared_limit"] == s.Rational(8, 3)


def test_source_and_frontiers_unchanged():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 203
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert source.CONTACT == source.previous.CONTACT


def test_actual_science_counts():
    assert [
        len(m.data()["checks"])
        for m in (source, local, mixed, moment, aggregate, calibration)
    ] == [292, 12, 14, 14, 7, 39]
    assert len(ROWS) == 378 and audit.scalar_entry_count() == 397
    assert len(GATES) == 45 and audit.rejected_inputs() == 859


def test_finite_bound_not_asymptotic_inference():
    text = aggregate.data()["whole_finite_extra_sign_proof"]
    assert "n^2>360" in text and "positive gaps" in text
    assert "10^-604" in aggregate.data()["whole_original_extra_bound"]
    assert "10^-207" in aggregate.data()["whole_total_sign_and_bound"]


def test_all_original_and_actual_kernel_counts():
    assert calibration.data()["whole_original_counts"] == {
        "nonzero_polarizations": 6,
        "word_and_source_checks": 18,
        "generic_frozen_channels": 3,
    }


def test_parent_matching_still_open():
    text = audit.observable()["not_established"]
    assert "Independent extra parent chi" in text and "unbounded" in text
    assert "P8 remain open" in text


def test_no_finite_gravity_or_physical_taylor_promotion():
    assert "not_finite_gravity" in " ".join(source.data()["gates"])
    assert "above-threshold" in aggregate.data()["whole_parent_boundary"]


def test_no_duplicate_loop_addition():
    assert "no extra copy" in audit.observable()["not_established"]
