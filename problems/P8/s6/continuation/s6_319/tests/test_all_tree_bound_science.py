"""Exact core inventories, all-valence bounds and original full amplitude."""

import pytest
import sympy as s
from p8_vacuum_affine_uniform_all_tree_bound import (
    audit,
    bounds,
    core,
    source,
    vertices,
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


@pytest.mark.parametrize("n", range(6))
def test_independent_core_counts_are_not_full_counts(n):
    value = core.inventory(n)
    assert value == core.data()["whole_core_inventories"][n]
    if n >= 2:
        assert value.subs({core.C: 1, core.g: 1}) < core.full.topology_closed(n).subs(
            {core.C: 1, core.g: 1}
        )


@pytest.mark.parametrize("n", range(1, 7))
@pytest.mark.parametrize("sector", ("matter", "gravity"))
def test_composed_coefficients_under_all_order_barrier(n, sector):
    value = bounds.sector_coefficient(n, sector)
    radius = bounds.barrier()["radius"]
    cap = 3 if sector == "matter" else 2 * 10**11
    assert value > 0 and value < cap * s.factorial(n) / radius**n


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
def test_public_multiplicity_validation_precedes_cache(invalid):
    bounds.sector_coefficient(1, "matter")
    for call, args in (
        (bounds.coefficient_envelope, (invalid,)),
        (bounds.sector_coefficient, (invalid, "matter")),
    ):
        with pytest.raises((TypeError, ValueError)):
            call(*args)


@pytest.mark.parametrize(
    "invalid",
    (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
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
def test_core_count_validation_precedes_cache(invalid):
    core.inventory(1)
    with pytest.raises((TypeError, ValueError)):
        core.inventory(invalid)


@pytest.mark.parametrize(
    "sector",
    (True, False, 1.0, None, "inclusive", "UV_complete", "zero_matching", [], {}),
)
def test_unsupported_core_sector_rejected(sector):
    with pytest.raises((TypeError, ValueError)):
        bounds.sector_coefficient(1, sector)


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
        bounds.amplitude_upper(energies)


@pytest.mark.parametrize("n", range(1, 7))
def test_soft_energy_scaling_retains_one_pole_per_leaf(n):
    ws = (s.Rational(1, 16 * n),) * n
    assert bounds.amplitude_upper(
        tuple(w / 2 for w in ws)
    ) == 2**n * bounds.amplitude_upper(ws)
    assert bounds.amplitude_upper(ws) * s.prod(ws) == bounds.coefficient_envelope(n)


def test_exact_composition_barrier():
    row = bounds.barrier()
    cap, radius = row["cap"], row["radius"]
    assert radius == cap / 2
    assert row["margin"] == cap - radius - row["soft_nonlinearity"] > 0
    assert row["core_values"]["V"] < s.Rational(1, 1000)
    assert row["core_values"]["matter"] < 3
    assert row["core_values"]["gravity"] < 2 * 10**11


def test_complete_original_5116_not_a_leading_soft_proxy():
    row = bounds.original_calibration()
    assert row["whole_exact_calibration"]["graph_count"] == 5116
    assert all(bool(v) for v in row["gates"].values())
    actual = row["whole_exact_calibration"]["scaled_unit_polarization_amplitude_square"]
    assert 0 < actual < bounds.coefficient_envelope(3) ** 2


@pytest.mark.parametrize("n,exponent", ((3, 755), (4, 1138)))
def test_original_numeric_envelopes(n, exponent):
    assert bounds.coefficient_envelope(n) < s.Rational(1, 10**exponent)


def test_original_source_parameter_record_unchanged():
    assert source.original_parameters() == source.previous.original_parameters()
    assert source.KAPPA == 10**800


def test_trace_reversal_norm_needs_no_conservation_shortcut():
    row = vertices.data()
    assert row["checks"]["trace_reversal_exact_Frobenius_isometry"] == s.zeros(16)
    assert row["gates"]["trace_reversal_norm_does_not_require_current_conservation"]


def test_bare_bound_remains_IR_divergent():
    text = bounds.data()["whole_infrared_boundary"]
    assert all(word in text for word in ("dwi/wi", "infrared", "real-virtual", "Regge"))


def test_original_frontiers_and_nonclosure():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 175
    assert len(audit.frontier()) == 9 and len(audit.qualifications()) == 6
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    text = audit.observable()["not_established"]
    assert all(
        word in text
        for word in ("inclusive", "Regge", "Energy-independent", "summability")
    )
