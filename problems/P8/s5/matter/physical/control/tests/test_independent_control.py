"""Independent Hamiltonian, cosmic-time and finite-q known-answer controls."""

import pytest
import sympy as sp
from p8_m1 import nonlinear
from p8_m1_control import independent
from p8_m1_control import model as compact
from p8_m1_physical import quadratic


def assert_zero(value):
    entries = value if isinstance(value, sp.MatrixBase) else (value,)
    assert all(sp.cancel(entry) == 0 for entry in entries)


def test_generic_coupled_canonical_and_energy_identities():
    for name, residual in independent.generic_checks().items():
        try:
            assert_zero(residual)
        except AssertionError as error:
            raise AssertionError(name) from error


def test_generic_omission_controls_are_nonzero():
    for name, residual in independent.generic_negative_controls().items():
        entries = residual if isinstance(residual, sp.MatrixBase) else (residual,)
        assert any(sp.cancel(entry) != 0 for entry in entries), name


def test_jet_map_matches_direct_time_dependent_matrix_differentiation():
    t = sp.Symbol("test_time", real=True)
    T = sp.Matrix([[1+t*t, 0], [t+t*t, 2+t]])
    inverse = T.inv()
    A = inverse*inverse.T
    B = sp.Matrix([[t*t, 1+t], [2-t, t**3]])
    C = sp.Matrix([[3+t, 2-t*t], [2-t*t, 5+2*t]])
    H = t/(1+t*t)
    Bhat = T*B*inverse+3*H*sp.eye(2)/2+T.diff(t)*inverse
    S, Omega = (Bhat+Bhat.T)/2, -(Bhat-Bhat.T)/2
    W = inverse.T*C*inverse-S*S+S*Omega-Omega*S-S.diff(t)
    at = lambda value: value.subs(t, sp.Rational(1, 3))
    actual = independent.oscillator_from_jets(
        at(A), at(B), at(C), at(B.diff(t)), at(T), at(T.diff(t)),
        at(T.diff(t, 2)), at(H), at(sp.diff(H, t)))
    assert_zero(actual["kinetic_residual"])
    assert_zero(actual["potential"]-at(W))
    assert_zero(actual["connection"]-at(Omega))
    assert_zero(actual["momentum_boundary"]-at(S))


def test_exact_bounce_hamiltonian_formula():
    q = sp.Symbol("q_bounce", positive=True)
    derived = independent.bounce_from_cosmic_time(q)
    expected = independent.bounce_expected(q)
    for name in ("d1", "d2", "chi", "ratio", "potential11", "potential22",
                 "potential12_factor", "connection_factor", "connection_dot_factor"):
        assert_zero(derived[name]-expected[name])
    for residual in derived["vanishing_checks"].values():
        assert_zero(residual)
    determinant = (derived["potential11"]*derived["potential22"]
                   -derived["ratio"]*derived["potential12_factor"]**2)
    assert_zero(determinant-expected["determinant"])


@pytest.mark.parametrize("q", (sp.Integer(8), sp.Integer(1000)))
def test_bounce_closed_formula_matches_full_cosmic_jet_route(q):
    actual = independent.cosmic_time_normalization("gamma", sp.Integer(0), q)
    expected = independent.bounce_expected(q)
    off = sp.sqrt(expected["ratio"])*expected["potential12_factor"]
    assert_zero(actual["potential"]-sp.Matrix([
        [expected["potential11"], off], [off, expected["potential22"]]]))
    assert_zero(actual["connection"])


def test_bounce_finite_q_mass_is_not_the_principal_cone():
    at_eight = independent.bounce_expected(sp.Integer(8))
    assert at_eight["d1"] > 0 and at_eight["d2"] > 0
    assert at_eight["determinant"] < 0
    assert at_eight["potential11"] < 0
    at_thousand = independent.bounce_expected(sp.Integer(1000))
    assert at_thousand["potential11"] > 0 and at_thousand["determinant"] > 0
    assert at_thousand["potential12_factor"] != 0
    # Instantaneous W indefiniteness at q=8 is not a dynamical-instability
    # theorem, and is different from the q=6 velocity-chart singularity.
    with pytest.raises(ValueError, match="q>6"):
        independent.bounce_expected(sp.Integer(6))
    with pytest.raises(ValueError, match="q\\*lambda"):
        independent.cosmic_time_normalization("gamma", sp.Integer(0), sp.Integer(6))
    with pytest.raises(ValueError, match="theta=0"):
        independent.cosmic_time_normalization("unitary", sp.Integer(0), sp.Integer(1000))


def test_bounce_high_q_mass_limits_retain_finite_corrections():
    q = sp.Symbol("q_bounce_limit", positive=True)
    data = independent.bounce_expected(q)
    assert sp.limit(data["potential11"]-q, q, sp.oo) == -sp.Rational(1139301, 59950)
    assert sp.limit(data["potential22"]-q, q, sp.oo) == -sp.Rational(301, 50)
    assert sp.limit(data["ratio"]*data["potential12_factor"], q, sp.oo) == -sp.Rational(199, 5995)


POINTS = (("unitary", sp.Rational(-3, 4)), ("unitary", sp.Rational(3, 4)),
          ("gamma", sp.Rational(-9, 40)), ("gamma", sp.Integer(0)),
          ("gamma", sp.Rational(9, 40)))


@pytest.mark.parametrize("chart,u_value", POINTS)
def test_original_cosmic_time_normalization_and_compact_weights(chart, u_value):
    data = independent.cosmic_time_normalization(chart, u_value, sp.Integer(10000))
    ell, delta = data["ell"], data["delta"]
    values = data["jets"][0]
    mapping = {compact.x: data["x"], compact.z: 1/data["local_q"],
               compact.l: ell*values[quadratic.l]}
    at = lambda expression: sp.cancel(expression.subs(mapping, simultaneous=True))
    matrices = compact.coefficients(chart)
    weights = sp.diag(ell**-delta, 1)
    assert_zero(data["alpha"]-weights*matrices["alpha"].applyfunc(at)*weights)
    assert_zero(data["B"]-weights.inv()*matrices["B"].applyfunc(at)*weights/ell)
    assert_zero(data["C"]-weights*matrices["C"].applyfunc(at)*weights/ell**2)
    assert_zero(data["kinetic_residual"])
    assert_zero(data["T"].T*data["T"]-data["alpha"])
    factors = {key: at(value) if isinstance(value, sp.Expr) else value
               for key, value in compact.factors(chart).items()}
    root_ratio = sp.sqrt(factors["ratio"])
    Tbar = sp.Matrix([[sp.sqrt(factors["d1"]), 0],
                     [factors["chi"]*sp.sqrt(factors["d2"]), sp.sqrt(factors["d2"])]])
    assert_zero(data["T"]-Tbar*weights)
    Fbar = sp.Matrix([[factors["f1"], 0], [root_ratio*factors["h"], factors["f2"]]])
    assert_zero(ell*data["volume_connection"]-Fbar)


@pytest.mark.parametrize("u_value", (sp.Rational(-3, 4), sp.Rational(3, 4)))
def test_nonzero_unitary_connection_sign_from_hamiltonian(u_value):
    data = independent.cosmic_time_normalization("unitary", u_value, sp.Integer(10000))
    expected = independent.unitary_connection_expected(u_value)
    assert_zero(data["connection"][0, 1]-expected)
    assert expected < 0  # Positive Cholesky orientation on both exterior pieces.
    assert_zero(data["connection"]+data["connection"].T)
    assert data["connection"][0, 1] != 0


@pytest.mark.parametrize("u_value", (sp.Rational(-9, 40), sp.Rational(3, 4)))
@pytest.mark.parametrize("dimension", (0, 1, 2, 4))
def test_compact_first_and_second_derivatives_against_original_cosmic_time(u_value, dimension):
    u, d = nonlinear.u, nonlinear.d
    ell, kcom_squared = sp.sqrt(d), sp.Integer(10000)
    polynomial = (compact.x**2+compact.z*compact.l+compact.l**2
                  +compact.x*compact.z**2+compact.l**3)
    mapping = {compact.x: u/ell, compact.z: d**3/kcom_squared,
               compact.l: 1/(10*d**sp.Rational(11, 2))}
    physical = polynomial.subs(mapping, simultaneous=True)/ell**dimension
    first = compact.derivative(polynomial, dimension)
    second = compact.derivative(first, dimension+1)
    for order, value in ((1, first), (2, second)):
        actual = (ell**(dimension+order)*sp.diff(physical, u, order)).subs(u, u_value)
        expected = value.subs(mapping, simultaneous=True).subs(u, u_value)
        assert_zero(actual-expected)


def test_physical_momentum_and_volume_drifts_are_not_local_momentum_drifts():
    u, d = nonlinear.u, nonlinear.d
    ell, kcom_squared = sp.sqrt(d), sp.Symbol("kcom_squared", positive=True)
    x, qlocal, qphysical = u/ell, kcom_squared/d**3, kcom_squared/d**4
    volume_root = d**3
    assert_zero(ell*sp.diff(qlocal, u)+6*x*qlocal)
    assert_zero(ell**3*sp.diff(qphysical, u)+8*x*qlocal)
    assert_zero(ell*sp.diff(volume_root, u)/volume_root-6*x)
    # The tempting fixed-local-q derivative drops a real drift at fixed kcom.
    assert (ell*sp.diff(qlocal, u)).subs(u, sp.Rational(3, 4)) != 0
    assert (ell**3*sp.diff(qphysical, u)+6*x*qlocal).subs(u, sp.Rational(3, 4)) != 0
