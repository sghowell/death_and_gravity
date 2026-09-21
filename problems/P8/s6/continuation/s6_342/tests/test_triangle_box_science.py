"""Complete selected labelled radiation, exact contours and strict scope tests."""

import pytest
import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import tensor
from p8_vacuum_affine_minimal_gravity_radiation import vertices as original
from p8_vacuum_affine_triangle_box_radiation import (
    audit,
    bounds,
    calibration,
    functional,
    labels,
    radiation,
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
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_all_rejected_scopes(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("permutation", labels.PERMUTATIONS)
def test_every_distinct_label(permutation):
    assert len(permutation) == 4 and set(permutation) == set(range(4))


@pytest.mark.parametrize("sector", bounds.SECTORS)
def test_both_selected_remainder_interferences(sector):
    value = bounds.finite_interference_upper(
        s.Rational(25, 4), 0, s.Rational(1, 384), sector
    )
    assert 0 < value < s.Rational(1, 10**592)
    assert not value.has(s.Float)


@pytest.mark.parametrize("index", range(3))
def test_hard_kernel_constant_diagnostic_only(index):
    ps, k, _ = original.sample(index)
    q = [ps[0] + k, *ps[1:]]
    constC = lambda *args: s.S.One
    constD = lambda *args: s.S.One
    actual = labels.hard(q, constC, constD)
    expected = (
        2
        * source.CUBIC**2
        * sum(
            labels.endpoint(radiation.dot(q[0] + q[j], q[0] + q[j])) for j in (1, 2, 3)
        )
        + 6 * source.CUBIC**4
    )
    assert s.cancel(s.expand(actual - expected)) == 0


@pytest.mark.parametrize("split", (-1, 3, True, 1.0))
def test_split_line_domain_rejected(split):
    ps, k, _ = original.sample(0)
    with pytest.raises((TypeError, ValueError)):
        radiation.setup(
            [ps[0] + ps[1], ps[2], ps[3]],
            [1, 1, source.HEAVY_MASS2],
            k,
            split,
            s.Rational(1, 2),
            source.HEAVY_MASS2,
        )


@pytest.mark.parametrize("gamma", (-1, 2, True, 0.5))
def test_split_parameter_domain_rejected(gamma):
    ps, k, _ = original.sample(0)
    with pytest.raises((TypeError, ValueError)):
        radiation.setup(
            [ps[0] + ps[1], ps[2], ps[3]],
            [1, 1, source.HEAVY_MASS2],
            k,
            0,
            gamma,
            source.HEAVY_MASS2,
        )


def test_wrong_polarization_rejected_before_TT_formula():
    ps, k, _ = original.sample(0)
    with pytest.raises(ValueError, match="TT"):
        radiation.line_density(
            [ps[0] + ps[1], ps[2], ps[3]],
            [1, 1, source.HEAVY_MASS2],
            k,
            radiation.ETA,
            0,
            s.Rational(1, 5),
            s.Rational(1, 7),
            s.Rational(2, 5),
            s.Rational(2, 3),
            source.HEAVY_MASS2,
        )


def test_one_actual_line_has_full_measure_and_no_float():
    ps, k, _ = original.sample(0)
    eps, norm = tensor.transverse_polarizations(tensor.frame(k), radiation.ETA)[0]
    assert norm > 0
    row = radiation.line_density(
        [ps[0] + ps[1], ps[2], ps[3]],
        [1, 1, source.HEAVY_MASS2],
        k,
        eps,
        1,
        s.Rational(1, 5),
        s.Rational(1, 7),
        s.Rational(2, 5),
        s.Rational(2, 3),
        source.HEAVY_MASS2,
    )
    assert s.im(s.expand(row["denominator"])) < 0
    assert not row["density"].has(s.Float)
    assert s.expand(sum(row["weights"]) - 1) == 0


def test_all_functional_positions_and_noncommuting_kernels():
    data = functional.data()
    assert len(data["whole_triangle_metric_positions"]) == 4
    assert len(data["whole_box_metric_positions"]) == 4
    assert (
        data["whole_diagnostic_matrices"]["light"]
        * data["whole_diagnostic_matrices"]["heavy"]
        != data["whole_diagnostic_matrices"]["heavy"]
        * data["whole_diagnostic_matrices"]["light"]
    )


def test_all_original_counts_and_coincidences():
    data = calibration.data()
    assert data["whole_literal_density_counts"] == {
        "routes": 672,
        "physical_TT": 1344,
        "nonzero_TT": 1284,
        "timelike": 224,
        "crossed": 448,
    }
    assert data["whole_coincident_polarized_channels"] == 4
    assert data["whole_wrong_sign_rejections"] == 8
    assert data["whole_nonzero_original_TT_controls"] == 8


def test_all_original_budget_values_are_exact():
    assert bounds.HARD_BUDGET == 205824000
    assert bounds.CHANGE_BUDGET == 263738720256000
    assert bounds.REMAINDER_BUDGET == 17515528820736000
    assert (
        bounds.unrounded_remainder_ratio() + bounds.known.REMAINDER_RATIO
        < bounds.REMAINDER_RATIO
    )
    assert bounds.REMAINDER_RATIO / s.sqrt(source.KAPPA) == s.Rational(1, 10**196)


def test_original_parameters_and_frontiers_unchanged():
    n, g = source.HEAVY_MASS2, source.CUBIC
    assert source.CONTACT == -g * g * (3 / (n - 2) - 2 / (n - 2) ** 2)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 198
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "internal-graviton" in audit.observable()["not_established"]
    assert "curved" in audit.observable()["not_established"]


def test_all_module_counts_and_counterterm_boundary():
    assert len(functional.data()["checks"]) == 7
    assert len(labels.data()["checks"]) == 3
    assert len(radiation.data()["checks"]) == 5
    assert len(bounds.data()["checks"]) == 9
    assert len(calibration.data()["checks"]) == 42
    assert (
        "not double-counted" in radiation.data()["whole_counterterm_and_UV_completion"]
    )
    assert "not" in bounds.data()["whole_large_hard_boundary"].lower()
