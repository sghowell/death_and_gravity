import pytest
import sympy as sp
from p8_m1_control import model as m
from p8_m1_control import oscillator


def zero(value):
    assert m.canonical(value) == 0


def test_exact_parity_field_against_independent_reduction():
    values = (m.l**2-m.L2, (1+m.l)*(1-m.l),
              (2+m.l)/(3-m.l), (1+m.x+m.l)**5/(2+m.z+m.l*m.l)**3)
    for value in values:
        numerator, denominator = sp.fraction(sp.cancel(value))
        n, d = m.reduce_l(numerator), m.reduce_l(denominator)
        # Independent polynomial remainder, with no fraction-field parser.
        a, b = sp.expand(d).coeff(m.l, 0), sp.expand(d).coeff(m.l, 1)
        independent = m.reduce_l(sp.expand(n*(a-b*m.l)))/(a*a-m.L2*b*b)
        assert sp.cancel(m.canonical(value)-independent) == 0
    with pytest.raises(ValueError, match="Unsupported"):
        m.canonical(sp.sin(m.l))


def test_derivative_respects_relation_and_physical_weights():
    zero(m.derivative(m.l*m.l)-m.derivative(m.L2))
    zero(m.derivative(1/m.z)+6*m.x/m.z)
    zero(m.derivative(1/m.z, 2)+8*m.x/m.z)
    value = (m.x+m.l)/(1+m.z)
    expected = (1-m.x*m.x)*sp.diff(value, m.x)+6*m.x*m.z*sp.diff(value, m.z)-11*m.x*m.l*sp.diff(value, m.l)
    zero(m.derivative(value)-expected)
    assert m.no_high_frequency_pole(value)
    with pytest.raises(ValueError, match="pole"):
        m.no_high_frequency_pole(1+m.l/m.z)


@pytest.mark.parametrize("chart", ["unitary", "gamma"])
def test_complete_cholesky_and_velocity_bridge(chart):
    data, f = m.coefficients(chart), m.factors(chart)
    d1, d2, chi = (f[k] for k in ("d1", "d2", "chi"))
    target = sp.Matrix([[d1+chi*chi*d2, chi*d2], [chi*d2, d2]])
    for value in data["alpha"]-target:
        zero(value)
    for value in data["alpha"]*data["B"]+data["beta"]:
        zero(value)


@pytest.mark.parametrize("chart", ["unitary", "gamma"])
def test_exact_luminal_remainder_not_a_sampled_limit(chart):
    data = oscillator.derive(chart)
    for value in data["unwhitened_mass_residual"]:
        assert m.no_high_frequency_pole(value)
    for key in ("mass11", "mass22", "mass12_factor", "connection_factor",
                "covariant_mass11", "covariant_mass22", "covariant_mass12_factor"):
        assert m.no_high_frequency_pole(data[key])
    # Deleting the matched principal q*alpha subtraction reintroduces a pole.
    with pytest.raises(ValueError, match="pole"):
        m.no_high_frequency_pole(data["unwhitened_mass_residual"][0, 0]+m.coefficients(chart)["alpha"][0, 0]/m.z)


@pytest.mark.parametrize("x_value", [-1, 1])
def test_compact_tail_scalar_and_tensor_limits(x_value):
    data = oscillator.physical_matrices("unitary", x_value, 1000)
    assert data["mass"] == -30*sp.eye(2)
    assert data["connection"] == sp.zeros(2)
    assert oscillator.tensor()["mass"].subs(m.x, x_value) == -30


def test_invalid_point_and_velocity_chart_rejected():
    for args in (("gamma", 0, 8), ("gamma", sp.Rational(1, 2), 1000),
                 ("unitary", 0, 1000), ("unitary", 2, 1000)):
        with pytest.raises(ValueError):
            oscillator.physical_matrices(*args)
