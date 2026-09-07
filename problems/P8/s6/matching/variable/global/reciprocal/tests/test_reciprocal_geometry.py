from fractions import Fraction

import pytest
import sympy as sp
from p8_reciprocal_geometry import identities


def test_exact_source_and_integrating_factor_inputs():
    assert len(identities.checks()) == 8
    assert set(identities.checks().values()) == {0}


def test_center_sign_requires_positive_lapse_separation_without_h_division():
    y0, c0, hp0 = sp.symbols("y0 c0 hp0", positive=True)
    P0 = 2*y0**3*hp0/(c0*(c0-y0))
    assert sp.factor(P0*(c0-y0)-2*y0**3*hp0/c0) == 0
    assert (2*y0**3*hp0/c0).is_positive
    # Exact concrete controls: c>y and c<y have opposite P signs.
    assert P0.subs({y0: 2, c0: 4, hp0: 4}) == 8
    assert P0.subs({y0: 2, c0: 1, hp0: 4}) == -64


def test_independent_differential_function_derivation():
    u = sp.Symbol("u", real=True)
    a, c = sp.Function("a", positive=True)(u), sp.Function("c", positive=True)(u)
    alpha, P = sp.symbols("alpha P", positive=True)
    h, y = sp.diff(a, u)/a, alpha/a**2
    z = h/c
    D, B = h/y-z, sp.diff(h/y, u)
    # Rearrange the actual f equation before substitution, preserving c'.
    source_residual = sp.diff(z, u)-(c-y)*P/(2*y**3)
    result = sp.diff(D, u)+P*D/(2*y*y*z)-B+source_residual
    assert sp.factor(result) == 0
    assert sp.factor(B-(sp.diff(a, u)**2+a*sp.diff(a, u, 2))/alpha) == 0


def test_cd_geometric_bound_has_positive_coefficients_not_only_point_samples():
    x = sp.Symbol("x", nonnegative=True)
    polynomial = sp.Poly(2*(1+x)**2*(1+7*x), x)
    assert polynomial.all_coeffs() == [14, 30, 18, 2]
    assert min(polynomial.all_coeffs()) > 0


def test_independent_affine_and_timelike_geodesic_constraints():
    b, c, momentum = sp.symbols("b c momentum", positive=True)
    null_du = momentum/(b*c)
    null_dx = momentum/b**2
    assert sp.factor(c*c*null_du**2-b*b*null_dx**2) == 0
    time_du = sp.sqrt(1+momentum**2/b**2)/c
    assert sp.factor(c*c*time_du**2-b*b*null_dx**2-1) == 0
    squared_ratio = sp.factor((1/time_du/(b*c/momentum))**2)
    assert squared_ratio == momentum**2/(b*b+momentum**2)
    assert sp.factor(1-squared_ratio) == b*b/(b*b+momentum**2)


def test_affine_bound_is_uniform_over_later_z_values():
    d = identities.derive()
    a, alpha, h, _hp, z, z0, _P = d["variables"]
    excess = sp.Symbol("nonnegative_z_increment", nonnegative=True)
    remainder = (alpha*h*(z0-z)/(a*z*z0)).subs(z, z0+excess)
    assert remainder.is_nonpositive
    assert Fraction(2, 1)/(Fraction(3, 2)*Fraction(5, 4)) == Fraction(16, 15)


def test_exact_tail_calibration_is_not_an_unchecked_numeric_background():
    assert identities.affine_tail_bound(2, Fraction(5, 4), Fraction(3, 2)) == Fraction(16, 15)
    assert identities.affine_tail_bound(2, Fraction(5, 4), Fraction(3, 2), 3) == Fraction(16, 45)


@pytest.mark.parametrize("index", range(4))
@pytest.mark.parametrize("bad", [0, -1, True, 0.01])
def test_invalid_tail_bound_inputs(index, bad):
    values = [2, 1, 1, 1]
    values[index] = bad
    with pytest.raises((TypeError, ValueError)):
        identities.affine_tail_bound(*values)


def test_added_f_null_stress_changes_the_invariant_equation():
    a, alpha, h, z, P, stress = sp.symbols("a alpha h z P f_null_stress", positive=True)
    hp = sp.Symbol("h_prime", real=True)
    y, c = alpha/a**2, h/z
    # With f stress, z' gains c*stress/2. Thus D'+F D=B no longer holds.
    zp = (c-y)*P/(2*y**3)+c*stress/2
    D = h/y-z
    derivative = a*h*sp.diff(D, a)+hp*sp.diff(D, h)+zp*sp.diff(D, z)
    correction = sp.factor(derivative+P*D/(2*y*y*z)-(hp+2*h*h)/y)
    assert correction == -c*stress/2
    assert correction.is_negative
