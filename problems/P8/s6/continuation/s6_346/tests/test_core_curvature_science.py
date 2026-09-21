"""Whole common-basis conversion, aggregate and finite-domain scope tests."""

import pytest
import sympy as s
from p8_vacuum_affine_core_curvature_coefficient import (
    audit,
    bounds,
    calibration,
    conversion,
    core,
    moments,
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


@pytest.mark.parametrize("powers", conversion.POWERS)
def test_all_exact_word_powers(powers):
    assert conversion.require_powers(powers) == powers


@pytest.mark.parametrize("order", range(4))
@pytest.mark.parametrize("mass", (s.Rational(1, 2), 1, 2, 3, 10))
def test_each_exact_positive_mass_moment(order, mass):
    value = moments.coefficient(order, mass)
    assert not value.has(s.Float)
    target = (
        s.integrate(moments.weights()[order], (moments.Z, 0, 1))
        if mass == 1
        else moments.triangle_moments()[order].subs(moments.N, mass)
    )
    assert s.factor(value - target) == 0


def test_full_conversion_table_not_inferred_from_samples():
    checks, tau, coefs, K = conversion.build()
    assert tuple(tau.values()) == conversion.EXPECTED
    assert len(tau) == len(coefs) == 10 and len(K) == 4 and len(checks) == 45
    assert all(len(v) == 6 for v in coefs.values())


def test_entire_bubble_normalization():
    assert conversion.data()["whole_bubble_conversion_coefficient"] == -4
    assert core.data()["checks"]["complete_A_squared_B_degree6_coefficient"] == 0


def test_all_moments_and_nonzero_equal_mass_values():
    assert len(moments.data()["whole_triangle_integrals"]) == 4
    assert all(moments.coefficient(j, 1) != 0 for j in range(4))


def test_source_and_frontiers_not_retuned():
    n, g = source.HEAVY_MASS2, source.CUBIC
    assert s.factor(source.CONTACT - g * g * core.contact_ratio(n)) == 0
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 202
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6


def test_actual_science_module_counts():
    assert [
        len(m.data()["checks"])
        for m in (source, conversion, moments, core, bounds, calibration)
    ] == [287, 45, 27, 12, 4, 72]
    assert len(ROWS) == 447 and audit.scalar_entry_count() == 465
    assert len(GATES) == 47 and audit.rejected_inputs() == 795


def test_entire_generic_sum():
    assert s.factor(core.generic_coefficient() - sum(core.sectors())) == 0
    assert (
        s.factor(
            core.original_coefficient()
            - core.generic_coefficient().subs(core.C, core.contact_ratio(core.N))
        )
        == 0
    )


def test_mass_limit_is_separate_from_finite_bound():
    assert core.data()["whole_original_heavy_mass_limit"] == -s.Rational(23, 105)
    packet = bounds.data()
    assert isinstance(packet["whole_exact_tail_budget"], s.Rational)
    assert packet["whole_positive_exponential_partial_sum"] > 10
    assert "not" in packet["whole_original_parameter_and_Born_bound"]


def test_all_original_diagnostic_counts():
    assert calibration.data()["whole_original_counts"] == {
        "word_checks": 60,
        "nonzero_curvature_polarizations": 6,
        "frozen_bubble_checks": 6,
    }


def test_full_matter_and_parent_matching_not_claimed():
    text = audit.observable()["not_established"]
    assert "not the entire known matter coefficient" in text
    assert "Independent extra parent chi" in text and "P8 remain open" in text


def test_no_physical_taylor_replacement_or_double_counting():
    text = core.data()["whole_physical_coefficient"]
    assert "no extra copy" in text
    assert "above-threshold" in bounds.data()["whole_original_parameter_and_Born_bound"]


def test_complete_source_scheme_explicit():
    assert "No finite derivative subtraction" in source.data()["whole_prescription"]
    assert (
        "New local tadpole and mixed-source" in source.data()["whole_sector_boundary"]
    )
