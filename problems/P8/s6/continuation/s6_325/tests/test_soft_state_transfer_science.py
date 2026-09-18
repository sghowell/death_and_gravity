"""Independent current, continuity, regulator, phase and scope checks."""

from itertools import product

import pytest
import sympy as s
from p8_vacuum_affine_radiative_soft_state_transfer import (
    audit,
    continuity,
    matching,
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
def test_every_uniform_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_scope_and_invalid_inputs_rejected(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("radial", (s.S.Zero, s.Rational(3, 5), s.S.One))
def test_massive_current_derivative_at_nontrivial_projected_momentum(radial):
    z = s.Symbol("z", real=True)
    n = s.ImmutableMatrix([s.Rational(3, 5), 0, s.Rational(4, 5)])
    v = s.ImmutableMatrix([s.Rational(1, 5), -s.Rational(2, 7), s.Rational(3, 11)])
    p = s.ImmutableMatrix([s.Rational(2, 5), s.Rational(1, 7), s.Rational(3, 10)])
    light = s.ImmutableMatrix([*(radial * n), s.sqrt(1 - radial**2)])
    P = s.eye(4) - light * light.T
    pz = p + z * v
    p4 = s.ImmutableMatrix([*pz, 0])
    den = s.sqrt(1 + pz.dot(pz)) - radial * pz.dot(n)
    literal = P * p4 * p4.T * P / den
    actual = s.diff(literal, z).subs(z, 0)
    E = s.sqrt(1 + p.dot(p))
    D = E - radial * p.dot(n)
    dD = p.dot(v) / E - radial * v.dot(n)
    q = P * s.ImmutableMatrix([*p, 0])
    dq = P * s.ImmutableMatrix([*v, 0])
    expected = (dq * q.T + q * dq.T) / D - q * q.T * dD / D**2
    assert (actual - expected).applyfunc(s.factor) == s.zeros(4)
    assert s.factor(s.trace(actual * actual)) < 144**2 * v.dot(v)


@pytest.mark.parametrize(
    "R,AB",
    tuple(
        product(
            (s.Rational(1, 8), s.Rational(1, 100), s.Rational(1, 10**100)),
            ((100000, 64000), (60000, 42000)),
        )
    ),
)
def test_minimum_majorant_integral_on_independent_exact_scales(R, AB):
    A, B = map(s.Integer, AB)
    tau = s.Symbol("tau", positive=True)
    cut = (A * R / B) ** 4
    assert 0 < cut < 1
    low = s.integrate(B * tau ** s.Rational(-3, 4), (tau, 0, cut))
    high = s.integrate(A * R / tau, (tau, cut, 1))
    assert (
        s.expand_log(
            low + high - 4 * A * R * (1 + s.log(B / (A * R))), force=True
        ).expand()
        == 0
    )
    assert B / A < 1
    assert 4 * A * R * (1 - s.log(R)) > low + high


def test_original_radial_calibrations_cover_collinear_splitting():
    row = continuity.calibration()
    records = row["whole_exact_recoil_radial_angular_calibrations"]
    assert len(records) == 18
    assert {r[0] for r in records} == {"one", "two", "three_split"}
    assert {r[2] for r in records} == {0, s.Rational(3, 5), 1}
    assert (
        len([k for k in row["checks"] if k.endswith("exact_collinear_splitting")]) == 6
    )
    assert all(row["gates"].values())
    assert "Frobenius" in row["whole_calibration_boundary"]


@pytest.mark.parametrize(
    "R", (0, s.Rational(1, 10**50), s.Rational(1, 128), s.Rational(1, 8))
)
def test_state_continuity_bounds_include_exact_Born_boundary(R):
    assert continuity.index_change_upper(R) == 1400 * R / source.KAPPA
    finite = continuity.finite_conversion_change_upper(R)
    uniform = continuity.regulator_quotient_upper(R)
    assert 2 * uniform == 3 * finite
    if R == 0:
        assert finite == uniform == 0
    else:
        assert finite == 10000 * R * (1 - s.log(R)) / source.KAPPA


@pytest.mark.parametrize(
    "bad",
    (
        True,
        False,
        None,
        "free",
        "1/64",
        [1],
        (1,),
        {"energy": 1},
        s.Matrix([1]),
        1.0,
        s.Float(1),
        s.I,
        -1,
        s.oo,
        s.Symbol("R"),
        s.Rational(1, 4),
    ),
)
def test_state_energy_guards_precede_logarithms(bad):
    for call in (
        continuity.index_change_upper,
        continuity.finite_conversion_change_upper,
        continuity.regulator_quotient_upper,
    ):
        with pytest.raises((TypeError, ValueError)):
            call(bad)


def test_zero_energy_continuity_is_not_assumed_from_a_constant_bound():
    R = s.Symbol("R", positive=True)
    assert s.limit(R * (1 - s.log(R)), R, 0) == 0
    assert s.diff(R * (1 - s.log(R)), R) == -s.log(R)
    eps = s.Symbol("eps", positive=True)
    assert s.limit(s.integrate(1 / R, (R, eps, s.Rational(1, 8))), eps, 0) == s.oo
    assert s.integrate(1 - s.log(R), (R, 0, s.Rational(1, 8))) == (2 + s.log(8)) / 8


def test_full_trace_and_radial_contributions_enter_both_budgets():
    row = continuity.bounds()
    assert row["whole_finite_and_regulator_numerators"] == (615000, 1045000)
    assert row["whole_current_and_quadratic_difference_caps"] == (866, 50000, 30000)
    assert row["whole_radial_interpolation_caps"] == {
        "H": (100000, 64000, 400000),
        "U": (60000, 42000, 240000),
    }
    assert s.Rational(615000, 72) < 10000
    assert s.Rational(1045000, 72) < 15000


def test_same_state_virtual_normalization_and_retained_phase():
    e, change, phi = s.symbols("e change phi", positive=True)
    V = -change / (4 * e) + s.I * phi / e
    real = change / (2 * e)
    assert s.expand(V + s.conjugate(V) + real) == 0
    assert s.limit(e * real, e, 0) == change / 2
    assert s.im(V) == phi / e
    wrong = change / (2 * e) - change / (4 * e)
    assert s.limit(e * wrong, e, 0) == change / 4


def test_regulator_limit_of_an_exact_integrated_nonconstant_state_change():
    e, x, A, D = s.symbols("e x A D", positive=True)
    integrated = (A + 2 * e * D) * x ** (1 + 2 * e) / (2 * e * (1 + 2 * e)) - A * x / (
        2 * e
    )
    expected = x * (D + A * (s.log(x) - 1))
    assert s.simplify(s.limit(integrated, e, 0) - expected) == 0
    # Replacing x-b by x misses the finite endpoint contribution.
    wrong = x * (D + A * s.log(x))
    assert s.simplify(wrong - expected) == A * x


def test_uniform_power_difference_inequality_has_zero_initial_value_and_positive_curvature():
    z = s.Symbol("z", nonnegative=True)
    h = s.exp(-z) - 1 + z
    assert h.subs(z, 0) == s.diff(h, z).subs(z, 0) == 0
    assert s.diff(h, z, 2) == s.exp(-z)
    assert s.exp(-z).is_positive
    assert (
        matching.data()["gates"]["additional_soft_projector_and_phase_full_D"] is True
    )


def test_Bose_calorimetric_finite_term_from_independent_log_Gamma_derivative():
    z = s.Symbol("z")
    logF = -s.EulerGamma * z - s.log(s.gamma(1 + z))
    assert s.diff(logF, z).subs(z, 0) == 0
    assert s.diff(logF, z, 2).subs(z, 0) == -(s.pi**2) / 6
    assert matching.data()["gates"]["Gamma_total_energy_term_nonzero"] is True


@pytest.mark.parametrize(
    "x", (s.Rational(1, 10**50), s.Rational(1, 128), s.Rational(1, 8))
)
def test_original_transfer_bound_and_resolution_limit(x):
    assert (
        matching.transfer_upper(x)
        == 200000 * x * (1 - s.log(x)) / s.Integer(10) ** 1600
    )
    assert matching.transfer_upper(x) > 0


@pytest.mark.parametrize(
    "bad",
    (
        True,
        False,
        None,
        "1",
        1.0,
        s.Float(1),
        s.I,
        0,
        -1,
        s.oo,
        s.Symbol("x"),
        s.Rational(1, 4),
    ),
)
def test_transfer_cutoff_guards(bad):
    with pytest.raises((TypeError, ValueError)):
        matching.transfer_upper(bad)


def test_original_couplings_and_scope_are_preserved():
    pars = source.original_parameters()
    assert (
        pars["heavy"] == s.Integer(10) ** 200 / 512 + 2
        and pars["kappa"] == s.Integer(10) ** 800
    )
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 181
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6 and audit.rejected_inputs() == 83
    assert len(audit.controls()) == 8
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    for word in (
        "hard loop",
        "NNLO",
        "evanescent",
        "all-N",
        "state",
        "Regge",
        "bounce",
    ):
        assert word in audit.observable()["not_established"]
