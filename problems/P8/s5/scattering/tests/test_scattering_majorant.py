from fractions import Fraction

import pytest
import sympy as sp
from p8_scattering import majorant as m


def test_positive_series_product_and_analytic_majorants():
    v = m.Series.at(1, 3)+m.Series.at(2, 2)
    u = sp.Symbol('u')
    expected = sp.Poly((3*u+2*u*u)**3, u)
    assert (v**3).coefficients == tuple(Fraction(str(expected.nth(i))) for i in range(5))
    for exponent in (sp.Rational(1, 2), sp.Rational(-1, 2)):
        actual = (1+3*u+2*u*u)**exponent
        coefficients = sp.Poly(actual.series(u, 0, 5).removeO(), u)
        bound = v.analytic(exponent)
        assert all(abs(Fraction(str(coefficients.nth(i)))) <= bound.coefficients[i] for i in range(5))
    with pytest.raises(ValueError, match="vanish"):
        (v+1).analytic(-1)


def test_uniform_domain_and_scale_inequalities():
    out = m.build()
    scale = int(out['sufficient_M_tau'])
    C3 = Fraction(out['cubic_transition_block_bound_numerator'])
    C4 = Fraction(out['quartic_connected_tree_block_bound_numerator'])
    assert C3/scale <= Fraction(1, 1000)
    assert C4/(scale*scale) <= Fraction(1, 1000)
    assert len(out['momentum_correction_majorants']) == 3
    assert all(Fraction(v) > 0 for v in out['momentum_correction_majorants'])
    assert all(-2 <= v['unit_weight'] <= 6 for v in out['stationary_coefficient_majorants'])


def test_exponential_energy_estimate_is_strictly_below_two():
    exponent = Fraction(20, 99)
    assert 1/(1-exponent) < 2
    # Independent first ten Taylor terms and a geometric tail enclose exp(v).
    terms = [exponent**n/Fraction(sp.factorial(n)) for n in range(10)]
    tail = exponent**10/Fraction(sp.factorial(10))/(1-exponent/Fraction(11))
    assert sum(terms)+tail < 2


def test_continuum_fiber_schur_majorant_not_a_delta_function_bound():
    out = m.build()['window']
    N = int(out['momentum_species_measure_upper'])
    D = int(out['Schur_measure_majorant'])
    assert D > 3+3*N
    assert D*D > 3*(3*N)  # rectangular 1-to-2 Schur estimate
    assert 'fixed total momentum' in out['Fourier_measure']
