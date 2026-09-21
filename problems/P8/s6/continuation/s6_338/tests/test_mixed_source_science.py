"""Exact source radiation, domain guards and complete finite-rate bounds."""

import pytest
import sympy as s
from p8_vacuum_affine_local_tadpole_radiation import vertices
from p8_vacuum_affine_mixed_source_radiation import (
    audit,
    bounds,
    calibration,
    contractions,
    mixed,
    radiation,
    source,
    triangles,
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


@pytest.mark.parametrize("name", source.CLASSES)
def test_all_original_source_classes(name):
    assert name in ("within_J2", "within_J4", "across_J2_J4")
    assert source.derivative_coefficient() == -4 * source.CUBIC / source.KAPPA


@pytest.mark.parametrize(
    "ss,z", ((s.Rational(25, 4), 0), (9, s.Rational(9, 10)), (16, s.Rational(-99, 100)))
)
@pytest.mark.parametrize("sector", bounds.SECTORS)
def test_original_all_resolution_finite_interference(ss, z, sector):
    delta = bounds.born.transfer_gap(ss, z)
    a = delta / 192
    assert bounds.finite_interference_upper(ss, z, 0, sector) == 0
    small = bounds.finite_interference_upper(ss, z, a / 2, sector)
    edge = bounds.finite_interference_upper(ss, z, a, sector)
    large = bounds.finite_interference_upper(ss, z, s.Rational(1, 8), sector)
    assert 0 < small < edge < large and small <= edge / 2
    assert edge < s.Rational(1, 10 ** (1394 if sector == "mixed_source" else 1393))
    assert not large.has(s.Float)


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


@pytest.mark.parametrize(
    "bad",
    (-1, True, 1.0, None, "all_loops", "unknown_chi_equals_zero", "full_inclusive"),
)
def test_known_sector_guard(bad):
    with pytest.raises((TypeError, ValueError)):
        bounds.finite_interference_upper(9, 0, 0, bad)


@pytest.mark.parametrize("channel", (1, 2, 3))
def test_full_traceful_Ward_before_physical_TT(channel):
    G, _ = vertices.generic_radiative_data()
    a = tuple(G[i, 4] for i in range(4))
    xi = s.symbols("xi_dot_p0:4")
    beta = -sum(xi)
    H = s.Matrix(4, 4, lambda i, j: a[i] * xi[j] + a[j] * xi[i])
    L = (0, channel)
    R = tuple(i for i in range(4) if i not in L)
    full, _, _, _ = radiation.one_channel(
        L, R, G, H, a, s.Symbol("n"), 2 * beta, tuple(beta * v for v in a)
    )
    assert s.factor(full) == 0 and not full.has(s.Float)


def test_actual_tadpole_and_EOM_terms_not_dropped():
    c = contractions.data()
    t = triangles.data()
    assert c["whole_source_tadpole_vertex"] != 0
    assert c["whole_source_metric_tadpole"] != 0
    assert c["whole_source_loop_insertion"] == 0
    assert t["whole_complete_EOM_external"] != 0 and t["whole_complete_EOM_metric"] != 0
    assert (
        s.factor(t["whole_complete_EOM_external"] + t["whole_complete_EOM_metric"]) == 0
    )
    assert mixed.data()["whole_internal_kernel_emission"] != 0


@pytest.mark.parametrize("degree", range(1, 9))
def test_all_operator_insertions_not_single_endpoint(degree):
    pp, rr = s.symbols("p_squared r_squared")
    complete = sum(rr**j * pp ** (degree - 1 - j) for j in range(degree))
    assert s.expand((rr - pp) * complete - (rr**degree - pp**degree)) == 0
    if degree > 1:
        assert s.expand((rr - pp) * pp ** (degree - 1) - (rr**degree - pp**degree)) != 0


def test_original_contact_and_selected_factor_unchanged():
    n, g = source.HEAVY_MASS2, source.CUBIC
    assert source.CONTACT == -g * g * (3 / (n - 2) - 2 / (n - 2) ** 2)
    assert radiation.centered_contact() != source.CONTACT
    assert radiation.source_factor() == 8 * (4 - n) / (16 * s.pi**2 * source.KAPPA)
    assert radiation.four_correction(s.Rational(4, 3), s.Rational(4, 3)) == 0


def test_all_original_calibrations_and_frontiers():
    assert len(calibration.data()["whole_original_recoil_calibrations"]) == 3
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 194
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "curvature" in audit.observable()["not_established"].lower()
    assert "internal-graviton" in audit.observable()["not_established"].lower()
