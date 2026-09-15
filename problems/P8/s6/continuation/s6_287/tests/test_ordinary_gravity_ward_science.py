"""Independent whole-metric, parameter-domain, raw-sign and scope tests."""

import pytest
import sympy as s
from p8_affine import verify as base
from p8_vacuum_affine_massive_graviton_cut_isolated_replay import amplitude as EH
from p8_vacuum_affine_massive_matter_graviton_vertex import ward as earlier_ward
from p8_vacuum_affine_ordinary_gravity_ward_bridge import (
    audit,
    continuity,
    gauge,
    soft,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not value.atoms(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_every_unsupported_input(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("D", (4, 5, 6))
@pytest.mark.parametrize("direction", ("trace", "off_diagonal", "mixed"))
def test_independent_full_metric_inverse_and_determinant_variation(D, direction):
    eta = s.diag(1, *([-1] * (D - 1)))
    p = s.Matrix([D + 2, 1, -1, *([0] * (D - 3))])
    k = s.Matrix([1, 1, 2, *([0] * (D - 3))])
    r = p - k
    mu = (p.T * eta * p)[0]
    dr = (r.T * eta * r)[0] - mu
    assert dr != 0
    T = p * r.T + r * p.T - eta * ((p.T * eta * r)[0] - mu)
    H = eta * T * eta - eta * s.trace(eta * T) / (D - 2)
    B = s.zeros(D)
    if direction == "trace":
        B = eta
    elif direction == "off_diagonal":
        B[0, 2] = B[2, 0] = 1
    else:
        B[0, 0] = 2
        B[1, 1] = -1
        B[0, 1] = B[1, 0] = 3
    eps = s.Symbol("metric_variation", real=True)
    g = eta + eps * B
    gi = g.inv()
    kcov = eta * k
    F = H * gi * kcov - kcov * s.trace(gi * H) / 2
    density = s.sqrt((-1) ** (D - 1) * g.det()) * (F.T * gi * F)[0]
    exact = s.simplify(s.diff(density, eps).subs(eps, 0))
    dinv = -eta * B * eta
    dF = H * dinv * kcov - kcov * s.trace(dinv * H) / 2
    quotient = dr * (s.trace(eta * B) * mu / 2 - (p.T * B * p)[0]) - 2 * (p.T * dF)[0]
    assert s.simplify(exact / dr - quotient) == 0


@pytest.mark.parametrize(
    "tau", (s.Rational(1), s.Rational(1, 4), s.Rational(1, 16), s.Rational(1, 64))
)
@pytest.mark.parametrize(
    "z,v",
    (
        (s.Rational(1, 4), s.Rational(1, 3)),
        (s.Rational(1, 2), s.Rational(2, 3)),
        (s.Rational(3, 4), s.Rational(1, 2)),
    ),
)
@pytest.mark.parametrize("ep", (0, s.Rational(1, 4)))
def test_independent_parameter_domain_and_z_squared_difference(tau, z, v, ep):
    h = (1 - v * v) / 4
    A = z * z + tau * (1 - z) ** 2 * h
    A0 = z * z
    assert 0 < A0 <= A <= 1
    if z <= s.Rational(1, 2):
        assert A >= z * z + tau * h / 4
    else:
        assert z / A <= 1 / z
    ratio = A0 / A
    # Prove the quarter-power monotonic comparison using positive integer powers.
    if ep == s.Rational(1, 4):
        assert ratio**3 >= ratio**4
    exact = s.cancel(z * z * (1 / A0 - 1 / A))
    assert exact == tau * (1 - z) ** 2 * h / A
    assert exact <= tau / (4 * A)


@pytest.mark.parametrize(
    "tau", (s.Rational(1), s.Rational(1, 4), s.Rational(1, 16), s.Rational(1, 64))
)
def test_all_recorded_majorants_positive_and_exact(tau):
    values = continuity.bounds(continuity.require_transfer(tau))
    assert len(values) == 7
    for value in values.values():
        assert value.is_positive is True
        assert not value.has(s.Float)


@pytest.mark.parametrize(
    "bad", (0, -1, 2, True, False, 1.0, s.Float(1), s.I, s.oo, None, "1")
)
def test_transfer_domain_rejects_inexact_and_outside_inputs(bad):
    with pytest.raises((TypeError, ValueError)):
        continuity.require_transfer(bad)


@pytest.mark.parametrize("mass", (1, 2, s.Rational(3, 2)))
def test_literal_negative_triangle_and_endpoint_IR_sign(mass):
    M0 = s.Rational(1, 4) / mass
    V0 = soft.eikonal(0, mass, 0)
    literal_relative = ((-s.I) ** 3 * s.I**3 * s.I) / (-s.I)
    triangle_residue = -2 * M0
    assert s.simplify(literal_relative * V0 * triangle_residue - mass) == 0
    assert s.simplify((-literal_relative) * V0 * triangle_residue - mass) != 0


@pytest.mark.parametrize("power", (0, 1, 2, 3, 4))
def test_one_massless_triangle_all_tensor_radial_shifts(power):
    ell, ep = s.symbols("ell ep", positive=True)
    primitive = ell ** (power + 2 * ep) / (power + 2 * ep)
    assert s.simplify(s.diff(primitive, ell) - ell ** (power - 1 + 2 * ep)) == 0


def test_independent_full_massless_bubble_Gamma_finite_constant():
    ep, tau = s.symbols("ep tau", positive=True)
    H = (
        s.gamma(1 - ep)
        * (4 * s.pi) ** (-ep)
        * s.gamma(1 + ep) ** 2
        / s.gamma(2 + 2 * ep)
    )
    assert H.subs(ep, 0) == 1
    derivative = s.simplify(s.diff(H, ep).subs(ep, 0))
    assert s.simplify(derivative - (s.EulerGamma - s.log(4 * s.pi) - 2)) == 0
    assert (
        s.simplify(
            -derivative - s.log(tau) - (2 - s.EulerGamma + s.log(4 * s.pi / tau))
        )
        == 0
    )


def test_literal_complete_EH_cubic_soft_homogeneity():
    fields = [
        s.diag(1, 2, -1, 3),
        s.Matrix([[1, 2, 0, 1], [2, -1, 1, 0], [0, 1, 2, 1], [1, 0, 1, 0]]),
        s.diag(2, 1, 3, -2),
    ]
    momenta = [
        s.Matrix([1, 2, -1, 3]),
        s.Matrix([2, -1, 1, 1]),
        s.Matrix([-3, -1, 0, -4]),
    ]
    scale = s.Symbol("soft_scale", real=True)
    full = EH.einstein_cubic(fields, momenta)
    scaled = EH.einstein_cubic(fields, [scale * q for q in momenta])
    assert full != 0
    assert s.expand(scaled - scale * scale * full) == 0


@pytest.mark.parametrize(
    "mu,tau", ((1, s.Rational(1, 4)), (2, s.Rational(1, 16)), (s.Rational(3, 2), 1))
)
def test_direct_Breit_F1_component_is_nonsingular_even_with_singular_F2(mu, tau):
    mu, tau = map(s.sympify, (mu, tau))
    eta = s.diag(1, -1, -1, -1)
    P = s.Matrix([s.sqrt(mu + 1 + tau / 4), 0, 1, 0])
    q = s.Matrix([0, s.sqrt(tau), 0, 0])
    F1 = 1 + s.sqrt(tau)
    F2 = 1 / tau**2
    G = 2 * P * P.T * F1 + (q * q.T + tau * eta) * F2
    assert s.simplify(G[0, 2] / (2 * P[0] * P[2]) - F1) == 0
    assert (2 * P[0] * P[2]) ** 2 >= 4 * (mu + 1)


@pytest.mark.parametrize("coefficient", (-3, -1, 1, 7))
def test_analytic_transverse_slope_freedom_is_not_removed(coefficient):
    eta = s.diag(1, -1, -1, -1)
    P = s.Matrix([3, 0, 1, 0])
    q = s.Matrix([0, 2, 0, 0])
    t = (q.T * eta * q)[0]
    Q = coefficient * earlier_ward.transverse_ricci(P, q)
    assert Q * eta * q == s.zeros(4, 1)
    assert Q[0, 2] / (2 * P[0] * P[2]) == coefficient * t / 2


def test_source_graph_completeness_and_gauge_scope_are_explicit():
    d = source.data()
    assert len(d["complete_pure_GR_proper_graph_inventory"]) == 5
    assert "Both orientations" in d["complete_pure_GR_proper_graph_inventory"]["D"]
    assert "metric reducible" in d["graph_completeness_boundary"]
    assert "no scalar term" in d["background_split"]
    assert "flat harmonic quadratic" in d["specified_ordinary_vertex"]
    assert d["checks"] is not source.previous.data()["checks"]


def test_every_general_D_metric_component_check_is_exact():
    for D in (4, 5, 6):
        rows = gauge.dimension_checks(D)
        assert len(rows) == 5
        assert (
            len(rows["every_symmetric_metric_scalar_denominator_cancellation"])
            == D * (D + 1) // 2
        )
        assert all(value == s.zeros(*value.shape) for value in rows.values())
    assert "not interpolation" in gauge.data()["general_dimension_proof"]


def test_complete_pole_and_finite_cancellation_does_not_exchange_IR_limits():
    d = soft.data()
    assert d["complete_raw_GR_proper_F1_at_zero_minus_one"].has(s.EulerGamma)
    assert "still has IR pole" in d["finite_t_IR_boundary"]
    assert "not a finite physical observable" in d["finite_t_IR_boundary"]
    assert "finite7" in d["ordinary_charge_result"]
    assert "light curvature" in d["matching_boundary"]


def test_boundaries_are_not_silently_completed():
    text = audit.observable()["not_established"]
    for word in ("finite-transfer", "F1 slope", "F2", "S285", "Regge", "V/G/B/P8"):
        assert word in text
    assert "not a finite slope" in continuity.data()["projection_scope"]


def test_parent_warmup_does_not_mutate_packet_contract():
    before = base.serialize(audit.packets())
    audit.previous.packets()
    assert base.serialize(audit.packets()) == before
    assert len(audit.gates()) == 36


def test_all_original_and_historical_frontiers_preserved():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 143 and len(audit.frontier()) == 9
    assert audit.qualifications() == audit.previous.qualifications()
    assert len(audit.qualifications()) == 6
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    assert audit.controls()["rejected_inputs"] == 73
    assert audit.require_parameters(audit.parameters()) == audit.parameters()
