"""Full Proca, angular, parameter and independent zero-momentum controls."""
import pytest
import sympy as sp
from p8_vector_bubble import independent, pole


def test_all_ordered_tensor_and_pole_identities():
    for module in (pole, independent):
        for value in module.checks().values():
            assert all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_rotational_moment_normalizations():
    l = sp.symbols("l0:4", real=True)
    radial = sp.Symbol("r2", positive=True)
    average = independent.angular_average
    assert average(sum(value**2 for value in l), l, radial) == radial
    assert average(sp.expand(sum(value**2 for value in l)**2), l, radial) == radial**2
    assert average(l[0]**4, l, radial) == radial**2/8
    assert average(l[0]**2*l[1]**2, l, radial) == radial**2/24
    assert average(l[0]*l[1], l, radial) == 0
    with pytest.raises(ValueError, match="degree four"):
        average(l[0]**6, l, radial)


def test_actual_mass_direction_and_all_external_momentum_coefficients():
    d = pole.actual_mass_direction()
    h, t, s = d["h"], d["time2"], d["space2"]
    assert d["alpha"] == 4/(9*h)
    assert d["beta"] == 28/(81*h)
    assert sp.factor(d["constant"]-904*pole.m2**2/(2187*h**2)) == 0
    assert sp.factor(d["second"]-416*pole.m2*(4*t+3*s)/(6561*h**2)) == 0
    assert sp.factor(d["fourth"]-32*(120*t**2+220*t*s+101*s**2)/(98415*h**2)) == 0


@pytest.mark.parametrize("time,space", ((1, 0), (0, 1), (1, 2), (2, 3)))
def test_independent_parameter_integral_in_actual_direction(time, space):
    d = independent.direct_shift()
    actual = sp.integrate(d["explicit_integrand"].subs({d["alpha"]: sp.Rational(4, 9),
                              d["beta"]: sp.Rational(28, 81), d["time"]: time, d["space"]: space}), (pole.x, 0, 1))
    original = pole.actual_mass_direction()
    expected = sum(original[name] for name in ("constant", "second", "fourth"))
    expected = expected.subs({original["h"]: 1, original["time2"]: time**2, original["space2"]: space**2})
    assert sp.factor(actual-expected) == 0


def test_discarding_longitudinal_numerator_loses_all_momentum_poles():
    # The incorrect propagator I/(p²+m²) leaves only m^4 tr(A²)*I0.
    incorrect = pole.m2**2*pole.TA2
    data = pole.parameter_integral()
    assert sp.factor(data["actual"]-incorrect) != 0
    assert data["fourth"] != 0


def test_two_insertion_zero_momentum_is_not_the_full_second_mass_variation():
    from p8_aligned_quantum import potential
    d = pole.actual_mass_direction()
    jets = potential.clock_jets()
    complete_quadratic_C = jets["pole_weight"]["N_second"].subs(jets["h"], d["h"])/2
    linear_mass_only = d["constant"]/pole.m2**2
    assert sp.factor(complete_quadratic_C-linear_mass_only) != 0
