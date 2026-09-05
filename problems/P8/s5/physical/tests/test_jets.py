import sympy as sp
from p8_physical import jets as j


def test_labelled_repeated_fields_keep_factorial():
    c = j.Context([(1, 0, 0), (0, 1, 0), (-1, -1, 0)])
    field = sum(c.leg(i) for i in range(3))
    assert (field**3).coefficient(c.full) == 6
    assert (c.leg(0)**2).is_zero()


def test_spatial_derivative_leibniz_and_closure():
    c = j.Context([(1, 2, 3), (-2, 1, 1), (1, -3, -4)])
    a, b, d = (c.leg(i) for i in range(3))
    assert ((a*b).derivative(1)-a.derivative(1)*b-a*b.derivative(1)).is_zero()
    assert (a*b*d).derivative(2).is_zero()


def test_inverse_and_square_root_in_exact_ring():
    c = j.Context([(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, -1, -1)])
    f = 1+sum(c.leg(i) for i in range(4))
    assert (f*f.power(-1)-1).is_zero()
    assert (f.power(sp.Rational(1, 2))**2-f).is_zero()


def test_symbolic_time_rational_coefficients_stay_exact():
    x = sp.Symbol("x", real=True)
    c = j.Context([(1, 0, 0), (-1, 0, 0)], parameter=x)
    value = c.leg(0, 1/(1+x*x))*c.leg(1, 1+x*x)
    assert value.coefficient(3) == 1
