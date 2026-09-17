"""Exact complex continuation, relative domains and soft-face exclusions."""

import pytest
import sympy as s
from p8_vacuum_affine_relative_energy_complex_tube import (
    audit,
    bounds,
    geometry,
    source,
)
from p8_vacuum_affine_relative_energy_complex_tube import continuation as c

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


@pytest.mark.parametrize("n", range(1, 9))
def test_complex_current_coefficients_and_envelope(n):
    value = bounds.pure_coefficient(n)
    assert value.is_Integer and 0 < value <= bounds.pure_envelope(n)
    assert value == bounds.data()["whole_first_eight_complex_coefficients"][n]


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
def test_exact_positive_order_validation_precedes_cache(invalid):
    assert bounds.pure_coefficient(1) == 1
    for call in (
        bounds.pure_coefficient,
        bounds.pure_envelope,
        bounds.scaled_coefficient,
    ):
        with pytest.raises((TypeError, ValueError)):
            call(invalid)


@pytest.mark.parametrize(
    "invalid",
    (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        s.zoo,
        s.nan,
        s.Symbol("unknown"),
        s.sqrt(2),
        [],
        {},
    ),
)
def test_exact_gaussian_rejections(invalid):
    with pytest.raises((TypeError, ValueError)):
        c.exact_gaussian(invalid)


@pytest.mark.parametrize(
    "value",
    (0, 1, -1, s.I, -s.I, s.Rational(2, 3) + s.I / 5, (3 + s.I) / (7 - 2 * s.I)),
)
def test_exact_gaussian_acceptance_and_no_floats(value):
    result = c.exact_gaussian(value)
    assert c.canonical(result - value) == 0
    assert not result.has(s.Float)
    assert all(v.is_Rational for v in result.as_real_imag())


def test_complex_projection_temporal_idempotent_and_gauge_annihilating():
    Q = s.Matrix([1 + s.I / 100, 1, 0, 0])
    A = s.diag(1, s.I, 2, 3)
    made = c.project_complex(A, Q)
    assert made[0, :] == s.zeros(1, 4)
    assert c.project_complex(made, Q) == made
    q = c.ETA * Q
    xi = s.Matrix([s.I / 7, 1, -1, 2])
    assert c.project_complex(q * xi.T + xi * q.T, Q) == s.zeros(4)


def test_frozen_real_projection_guard_was_not_weakened():
    from p8_vacuum_affine_temporal_tree_reorganization import chart

    with pytest.raises((TypeError, ValueError)):
        chart.project(s.eye(4), s.Matrix([1 + s.I, 1, 0, 0]))


@pytest.mark.parametrize(
    "Q", (s.zeros(4, 1), s.ones(3, 1), s.eye(4), [True, 1, 0, 0], [1.0, 1, 0, 0])
)
def test_complex_projection_energy_and_shape_rejections(Q):
    with pytest.raises((TypeError, ValueError)):
        c.project_complex(s.eye(4), Q)


@pytest.mark.parametrize(
    "A",
    (
        s.ones(3),
        s.ones(4, 1),
        s.Matrix([[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
        [[1.0] * 4] * 4,
    ),
)
def test_complex_projection_tensor_rejections(A):
    with pytest.raises((TypeError, ValueError)):
        c.project_complex(A, s.Matrix([1, 1, 0, 0]))


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
def test_center_energy_rejections(energies):
    with pytest.raises((TypeError, ValueError)):
        bounds.bare_upper(energies)


def test_relative_tube_has_independent_complex_energy_shifts():
    ws = (s.Rational(1, 64), s.Rational(1, 32))
    zs = (ws[0] * (1 + s.I * bounds.EPSILON / 2), ws[1] * (1 - bounds.EPSILON / 3))
    assert c.require_relative_tube(ws, zs) == (ws, zs)
    for bad in (
        (0, zs[1]),
        (ws[0] * (1 + 2 * s.I * bounds.EPSILON), zs[1]),
        zs[:1],
        None,
    ):
        with pytest.raises((TypeError, ValueError)):
            c.require_relative_tube(ws, bad)


@pytest.mark.parametrize(
    "alpha",
    (
        None,
        "0",
        (),
        (0, 0),
        (True,),
        (-1,),
        (1.0,),
        (s.Float(1),),
        (s.Rational(1, 2),),
        (s.I,),
        ([],),
    ),
)
def test_Cauchy_multiindex_rejections(alpha):
    with pytest.raises((TypeError, ValueError)):
        bounds.cauchy_majorant((s.Rational(1, 64),), alpha)


@pytest.mark.parametrize("alpha", ((0, 0), (1, 0), (0, 1), (1, 1), (2, 3)))
def test_Cauchy_factorials_and_energy_scaling(alpha):
    ws = (s.Rational(1, 64), s.Rational(1, 32))
    target = bounds.scaled_coefficient(2) * s.prod(
        s.factorial(a) / (bounds.EPSILON * w) ** a for a, w in zip(alpha, ws)
    )
    assert bounds.cauchy_majorant(ws, alpha) == target
    assert (
        bounds.cauchy_majorant(tuple(w / 2 for w in ws), alpha)
        == 2 ** sum(alpha) * target
    )


def test_complete_current_pole_rejects_global_W_disc_not_real_domain():
    row = geometry.pole_data()
    record = row["whole_global_total_energy_tube_obstruction"]
    assert record["TT_residues"] == (s.Rational(47089, 2601), 0)
    assert record["required_energy_shift"] < record["global_total_energy_radius"]
    assert record["pole_energies"][-1] < 0
    assert row["gates"]["negative_energy_pole_not_called_physical_real_radiation"]


def test_genuinely_complex_pure_and_original_complete_sources():
    pure, full = c.pure_calibration(), c.full_calibration()
    assert all(bool(v) for v in pure["gates"].values())
    assert all(bool(v) for v in full["gates"].values())
    assert pure["whole_three_independent_energy_shifts"]["graph_count"] == 4
    assert full["whole_original_complex_5116"]["graph_count"] == 5116


def test_exact_barriers_do_not_reach_zero_energy_faces():
    row = bounds.barrier_data()
    assert row["pure_margin"] > 0 and row["hard_margin"] > 0
    assert row["core_values"]["V"] < s.Rational(1, 1000)
    assert 0 < bounds.EPSILON < 1
    assert "do not reach soft faces" in bounds.data()["whole_nonclosure"]


def test_original_parameters_and_all_frontiers_unchanged():
    assert source.original_parameters() == source.previous.original_parameters()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 176
    assert len(audit.frontier()) == 9 and len(audit.qualifications()) == 6
    scope = audit.observable()["not_established"]
    assert all(
        word in scope for word in ("soft", "inclusive", "Regge", "Energy-independent")
    )
