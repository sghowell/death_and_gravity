"""Exact enclosures and independent full-mode finite-parameter diagnostics."""

from fractions import Fraction
from functools import cache

import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_coupled_principal_obstruction import coupled as c
from p8_vacuum_affine_finite_window_growth import audit, background, window
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate, fraction
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_whole_exact_identity(name):
    value = audit.residuals()[name]
    assert all(
        entry == 0
        for entry in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_domain_and_scope_rejection(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_counts_and_complete_fixed_scope():
    assert (
        len(audit.residuals()),
        audit.scalar_entry_count(),
        len(audit.gates()),
        len(audit.controls()),
        audit.rejected_inputs(),
        len(audit.matching()),
    ) == (62, 297, 109, 9, 277, 112)
    assert all(value is True for value in audit.gates().values())
    assert audit.require_parameters(audit.parameters()) == audit.parameters()
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()


@pytest.mark.parametrize(
    "value", (True, False, 1.0, s.Float(1), "1", None, s.oo, s.nan, s.I, s.Symbol("x"))
)
def test_interval_rejects_nonexact_and_nonfinite_inputs(value):
    with pytest.raises(TypeError):
        fraction(value)


@pytest.mark.parametrize("left,right", ((-3, -1), (-2, 4), (1, 3), (0, 2)))
@pytest.mark.parametrize("power", range(6))
def test_independent_integer_interval_endpoints_and_interior(left, right, power):
    interval = I(left, right) ** power
    for index in range(9):
        x = Fraction(left) + Fraction(index, 8) * (right - left)
        assert interval.lo <= x**power <= interval.hi


@pytest.mark.parametrize("exponent", (-6000, -300, -10, 0, 10, 300, 6000))
@pytest.mark.parametrize("order", (2, 3, 4))
def test_adaptive_positive_root_encloses_huge_and_tiny_rationals(exponent, order):
    x = Fraction(7, 3) * Fraction(10) ** exponent
    root = I(x) ** Fraction(1, order)
    assert 0 < root.lo <= root.hi
    assert root.lo**order <= x <= root.hi**order
    reciprocal = I(x) ** Fraction(-1, order)
    assert reciprocal.lo**order <= 1 / x <= reciprocal.hi**order
    assert root.hi / root.lo < Fraction(1001, 1000)


@pytest.mark.parametrize("lo,hi", ((-1, 1), (0, 1), (-1, 0)))
def test_zero_crossing_inverse_and_fractional_roots_rejected(lo, hi):
    with pytest.raises(AssertionError):
        I(lo, hi).inverse()
    with pytest.raises(AssertionError):
        I(lo, hi) ** s.Rational(1, 2)


def test_whole_exact_expression_enclosure_retains_dependency_safely():
    x, y = s.symbols("x y")
    expr = (x - y) ** 2 / (x + y) + s.sqrt(x)
    bound = evaluate(expr, {x: I(1, 2), y: I(3, 4)})
    for xx in (1, s.Rational(3, 2), 2):
        for yy in (3, s.Rational(7, 2), 4):
            exact = expr.subs({x: xx, y: yy})
            assert s.Rational(bound.lo.numerator, bound.lo.denominator) <= exact
            assert exact <= s.Rational(bound.hi.numerator, bound.hi.denominator)
    with pytest.raises(TypeError):
        evaluate(s.sin(x), {x: I(1, 2)})


@pytest.mark.parametrize("slot", (0, 1, 2))
def test_independent_joint_lapse_and_all_three_Euler_equations(slot):
    # An independent exact four-by-four solve, not a selected Euler component.
    D, Z = s.Rational(5, 4), s.Rational(6, 5)
    cross = s.Matrix([s.Rational(slot - 1, 7), s.Rational(2, 9), -s.Rational(3, 11)])
    force = s.Matrix([s.Rational(7, 13), -s.Rational(2, 5), s.Rational(1, 17)])
    rest = s.Rational(2, 19)
    kinetic = s.diag(-6 * D, Z, Z)
    Cnn = s.Rational(8, 5)
    J = Cnn - (cross.T * kinetic.inv() * cross)[0] / 2
    joint = kinetic.row_join(cross).col_join(cross.T.row_join(s.Matrix([[2 * Cnn]])))
    direct = joint.inv() * force.col_join(s.Matrix([-rest]))
    Nd = (-rest - (cross.T * kinetic.inv() * force)[0]) / (2 * J)
    expected = (kinetic.inv() * (force - cross * Nd)).col_join(s.Matrix([Nd]))
    assert (direct - expected).applyfunc(s.factor) == s.zeros(4, 1)
    assert J > 0 and direct[2] != 0


def test_mass_adaptation_keeps_huge_mass_without_raw_mass_error():
    # Full fixed coefficients at a=Z=Y=1, not an integration at huge mass.
    mass = audit.parameters()["mass_squared"]
    P = background.MOMENTUM_MIN
    A = mass + P**2
    core = s.Matrix([[0, 1], [-A, 0]])
    W = s.diag(s.sqrt(A), 1)
    transformed = (W * core * W.inv()).applyfunc(s.simplify)
    assert transformed + transformed.T == s.zeros(2)
    assert core[1, 0] < -(10**197)
    assert transformed[0, 1] ** 2 == A
    assert (2 * background.MOMENTUM_MAX) ** 2 < mass
    assert window.growth_data()["principal_fast_rate_uniform_upper_bound"] ** 2 < mass


def test_all_time_connections_survive_independent_differentiation():
    t = s.Symbol("t", real=True)
    K = s.diag(s.exp(t), s.exp(-t))
    A = s.Matrix([[0, 7], [-7, 0]])
    full = K * A * K.inv() + K.diff(t) * K.inv()
    removed = K * A * K.inv()
    assert s.simplify(full - removed) == s.diag(1, -1)
    x = s.Matrix([s.exp(t) * s.cos(7 * t), -s.exp(-t) * s.sin(7 * t)])
    assert (x.diff(t) - full * x).applyfunc(s.simplify) == s.zeros(2, 1)
    assert (x.diff(t) - removed * x).subs(t, 0) != s.zeros(2, 1)


def test_retained_source_and_mass_are_separate_full_hessian_inputs():
    packet = window.matrix_data()
    H = packet["whole_on_shell_eight_phase_Hessian"]
    mass = s.Symbol("positive_heavy_potential_mass2", positive=True)
    for item in (mass, c.cs[1], c.ws[1], c.ds[1]):
        assert H.has(item)
        assert H.diff(item) != s.zeros(8)
    assert packet["whole_perfect_square_shear_derivative"] != s.zeros(8)
    assert len(packet["whole_normalized_mass_adapted_remainder"]) == 64


@cache
def cone_diagnostic(frequency):
    # Dimensionless, finite parameters only. The actual giant mass is
    # certified by the exact energy identity and enclosures, not this ODE.
    cycle = np.array(
        [[0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0]], dtype=float
    )
    slow = np.array([[0, 0], [-1 / 4, 0]], dtype=float)
    omega = np.array([[0, 1], [-1, 0]], dtype=float)
    core = block_diag(cycle, slow, *[w * omega for w in (frequency, 3, 5, 7, 11)])
    direction = np.r_[np.ones(4) / 2, np.zeros(12)]
    raw = np.arange(1, 257, dtype=float).reshape(16, 16)
    perturbation = raw / np.linalg.norm(raw) / 64

    def rhs(t, y):
        return (core + np.sin(2 * t) * perturbation) @ y

    times = np.linspace(0, 2, 41)
    run = solve_ivp(
        rhs, (0, 2), direction, t_eval=times, method="DOP853", rtol=2e-10, atol=2e-12
    )
    assert run.success
    f = direction @ run.y
    complement = run.y - direction[:, None] * f
    return times, f, np.linalg.norm(complement, axis=0), core


@pytest.mark.parametrize("frequency", (13, 137, 1021))
def test_independent_full_sixteen_phase_nonautonomous_cone(frequency):
    times, f, complement, core = cone_diagnostic(frequency)
    assert core.shape == (16, 16)
    assert np.max(complement - f) < -0.9
    assert np.min(f - np.exp(times / 2)) > -2e-9
    assert np.linalg.norm(core[6:8, 6:8]) > 10
    assert core[6:8, 6:8][1, 0] == -frequency


def test_actual_lifetime_has_quantitative_initial_margin_not_bare_continuity():
    packet = background.continuation_data()
    for margin in packet["persistent_initial_coefficient_margins"].values():
        assert margin > 2 * packet["uniform_coefficient_displacement"]
    assert len(packet["persistent_initial_coefficient_margins"]) == 14
    energy = background.energy_data()
    assert energy["heavy_root_energy_improved_bound"] == s.Rational(4, 10**2560)
    assert window.growth_data()["whole_growth_exponent_lower_bound"] == 5 * 10**31


@pytest.mark.parametrize("factor", (1, s.Rational(3, 2), 2))
def test_finite_band_endpoints_and_interior_are_admitted(factor):
    P = factor * background.MOMENTUM_MIN
    assert audit.require_momentum(P) == P
    assert audit.require_time(background.TIME_LENGTH / 2) == background.TIME_LENGTH / 2
    assert window.growth_data()[
        "physical_phase_conversion_prefactor_lower_bound"
    ] == s.Rational(1, 10**300)
