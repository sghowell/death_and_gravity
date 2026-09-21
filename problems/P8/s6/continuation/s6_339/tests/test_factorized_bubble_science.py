"""Exact factorized bubbles, physical branch guards and finite interference."""

import pytest
import sympy as s
from p8_vacuum_affine_factorized_bubble_radiation import (
    audit,
    bounds,
    branch,
    calibration,
    kernel,
    loops,
    source,
)
from p8_vacuum_affine_local_tadpole_radiation import vertices

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


@pytest.mark.parametrize("endpoint", source.ENDPOINTS)
def test_all_four_ordered_endpoint_types(endpoint):
    assert endpoint in (
        "contact_contact",
        "contact_heavy",
        "heavy_contact",
        "heavy_heavy",
    )
    assert source.CUBIC == s.Rational(1, 8192)


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
    assert edge < s.Rational(1, 10**602) and not large.has(s.Float)


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


@pytest.mark.parametrize(
    "bad", (-13, 1, 4, 5, 17, s.Float(0), True, "0", None, s.I, s.oo, s.Symbol("v"))
)
def test_physical_bubble_branch_guard(bad):
    with pytest.raises((TypeError, ValueError)):
        branch.bubble(bad)


@pytest.mark.parametrize("value", (-12, -1, 0, s.Rational(45, 8), 9, 16))
def test_actual_continuous_equal_invariant_quotient(value):
    assert branch.divided_difference(value, value) == branch.derivative(value)
    assert not branch.divided_difference(value, value).has(s.Float)


@pytest.mark.parametrize("left,right", ((-1, 9), (9, -1), (0, 16), (16, 0)))
def test_channel_segment_cannot_cross_threshold(left, right):
    with pytest.raises(ValueError):
        branch.divided_difference(left, right)


@pytest.mark.parametrize("channel", (1, 2, 3))
def test_general_kernel_traceful_Ward(channel):
    G, _ = vertices.generic_radiative_data()
    a = tuple(G[i, 4] for i in range(4))
    xi = s.symbols("xi_dot_p0:4")
    beta = -sum(xi)
    H = s.Matrix(4, 4, lambda i, j: a[i] * xi[j] + a[j] * xi[i])
    L = (0, channel)
    R = tuple(i for i in range(4) if i not in L)
    v = sum(G[i, j] for i in L for j in L)
    u = sum(G[i, j] for i in R for j in R)
    f0, f1 = s.symbols("F0 F1")
    full = kernel.channel(
        G, H, a, L, R, f0, f1, (f1 - f0) / (u - v), 2 * beta, tuple(beta * q for q in a)
    )
    assert s.factor(full) == 0 and not full.has(s.Float)


@pytest.mark.parametrize("degree", range(1, 9))
def test_all_operator_insertions_not_single_endpoint(degree):
    v, u = s.symbols("v u")
    dd = sum(u**j * v ** (degree - 1 - j) for j in range(degree))
    assert s.expand((u - v) * dd - (u**degree - v**degree)) == 0
    if degree > 1:
        assert s.expand((u - v) * v ** (degree - 1) - (u**degree - v**degree)) != 0


def test_cut_and_symmetric_subtraction_are_retained():
    assert s.im(branch.bubble(9)) == s.sqrt(5) * s.pi / 3
    assert s.im(branch.bubble(-1)) == 0
    assert branch.bubble(0) == 0 and branch.derivative(0) == s.Rational(1, 6)
    assert branch.symmetric_value() == 2 - 2 * s.sqrt(2) * s.atan(1 / s.sqrt(2))
    assert loops.data()["whole_existing_linear_finite_contact"] != 0


def test_one_outer_insertion_is_not_the_full_product():
    Al, Ar, Bl, Br, DA, DB = s.symbols("Al Ar Bl Br DA DB")
    full = DA * Br * Ar + Al * DB * Ar + Al * Bl * DA
    assert s.expand(full - (DA * Br * Ar + Al * DB * Ar)) != 0
    assert s.expand(full - Al * DB * Ar) != 0


def test_original_contact_and_bound_exponent_unchanged():
    n, g = source.HEAVY_MASS2, source.CUBIC
    assert source.CONTACT == -g * g * (3 / (n - 2) - 2 / (n - 2) ** 2)
    assert bounds.REMAINDER_RATIO == s.Integer(10) ** 194
    assert bounds.REMAINDER_RATIO / s.sqrt(source.KAPPA) == s.Rational(1, 10**206)
    assert bounds.HARD_RATIO == s.Integer(10) ** 190 > 1


def test_all_original_calibrations_and_frontiers():
    data = calibration.data()
    assert len(data["whole_exact_original_absorptive_calibrations"]) == 4
    assert data["whole_actual_equal_invariant_cases"] == 4
    assert len(calibration.component_checks()) == 30
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 195
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "curvature" in audit.observable()["not_established"].lower()
    assert "internal-graviton" in audit.observable()["not_established"].lower()
