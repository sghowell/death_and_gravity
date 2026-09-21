"""Complete local contractions, radiation guards and actual finite-rate bounds."""

import pytest
import sympy as s
from p8_vacuum_affine_local_tadpole_radiation import (
    audit,
    bounds,
    calibration,
    contractions,
    radiation,
    source,
    vertices,
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


@pytest.mark.parametrize("name", vertices.BASIS)
def test_literal_complete24_quartic_vertex(name):
    gram, _ = vertices.symbolic_grams()
    expected = s.Add(
        *(vertices.quartic_word(name, p, gram) for p in vertices.FOUR_ASSIGNMENTS)
    )
    assert s.expand(vertices.quartic_vertex(name, gram) - expected) == 0


@pytest.mark.parametrize("name", source.NAMES)
def test_whole_finite_evanscence_kept(name):
    dim = contractions.DIM
    basis = radiation.raw_basis_factors(dim)[name]
    written = s.Add(
        *(co * vertices.quartic_polynomials()[key] for key, co in basis.items())
    )
    assert s.expand(written - contractions.raw_quartic_factors()[name]) == 0
    finite = s.Add(
        *(
            co * vertices.quartic_polynomials()[key]
            for key, co in radiation.finite_basis_factors()[name].items()
        )
    )
    assert s.expand(finite - contractions.finite_quartic_factors()[name]) == 0


@pytest.mark.parametrize(
    "ss,z", ((s.Rational(25, 4), 0), (9, s.Rational(9, 10)), (16, s.Rational(-99, 100)))
)
def test_actual_all_resolution_finite_bound(ss, z):
    delta = bounds.born.transfer_gap(ss, z)
    a = delta / 192
    assert bounds.finite_interference_upper(ss, z, 0) == 0
    small = bounds.finite_interference_upper(ss, z, a / 2)
    edge = bounds.finite_interference_upper(ss, z, a)
    large = bounds.finite_interference_upper(ss, z, s.Rational(1, 8))
    assert 0 < small < edge < large
    assert small <= edge / 2
    assert edge < s.Rational(1, 10**1587)
    assert not large.has(s.Float)
    assert (
        s.limit(
            bounds.REMAINDER_RATIO
            * (
                128 * s.Symbol("x") + bounds.TREE_REMAINDER * s.Symbol("x") ** 2 / delta
            ),
            s.Symbol("x"),
            0,
        )
        == 0
    )


@pytest.mark.parametrize("bad", (-1, s.Rational(1, 7), 1, s.Float(0), True, "0", None))
def test_resolution_guard(bad):
    with pytest.raises((TypeError, ValueError)):
        bounds.finite_interference_upper(9, 0, bad)


@pytest.mark.parametrize(
    "ss,z", ((6, 0), (17, 0), (9, 1), (9, -1), (True, 0), (9, s.Float(0)))
)
def test_original_nonforward_domain_guard(ss, z):
    with pytest.raises((TypeError, ValueError)):
        bounds.finite_interference_upper(ss, z, 0)


@pytest.mark.parametrize("bad", (-1, 1, 3, True, 1.0, None))
def test_isotropic_rank_guard(bad):
    with pytest.raises((TypeError, ValueError)):
        contractions.denominator(bad)


@pytest.mark.parametrize("bad", (-1, 4, True, 1.0, None))
def test_exact_external_leg_guard(bad):
    gram, _ = vertices.generic_radiative_data()
    ak = tuple(gram[i, 4] for i in range(4))
    with pytest.raises((TypeError, ValueError)):
        vertices.shifted_vertex("Y2", gram, ak, bad)


def test_original_known_coefficients_not_matching_choices():
    assert source.KAPPA == s.Integer(10) ** 800
    assert source.CUBIC == s.Rational(1, 8192)
    assert source.LAMBDA == s.Rational(1, 10**600)
    assert 0 < bounds.coefficient_budget() < s.Rational(1, 10**1398)
    assert (
        0
        < bounds.remainder_upper() / s.Rational(1, 10**600)
        < bounds.REMAINDER_RATIO / s.sqrt(source.KAPPA)
    )
    assert radiation.four_polynomial(s.Rational(4, 3), s.Rational(4, 3)) == 0


def test_complete_original_calibrations_and_unchanged_scope():
    assert len(calibration.data()["whole_original_recoil_calibrations"]) == 3
    assert len(vertices.SIX_ASSIGNMENTS) == 720 and len(vertices.FOUR_ASSIGNMENTS) == 24
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 193
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "curvature" in audit.observable()["not_established"].lower()
    assert "other" in audit.observable()["not_established"].lower()
