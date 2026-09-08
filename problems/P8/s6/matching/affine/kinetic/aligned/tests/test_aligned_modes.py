"""Independent Euler equations, time normalization and momentum-chart controls."""
import pytest
import sympy as sp
from p8_affine_aligned import modes as m
from p8_affine_retuned.bounds import units


def test_exact_mode_identities_and_continuous_bounds():
    assert all(value == 0 for value in m.checks().values())
    assert all(value is True for value in m.proof_checks().values())


@pytest.mark.parametrize("sector", ("transverse", "longitudinal"))
def test_canonical_equation_from_original_weighted_Euler_equation(sector):
    data = m.canonical()
    weight = data[sector]["weight"]
    # Substitute X=v/g in (weight*X')'+weight*(q+1/zeta)*X=0.
    # Keep g abstract until after the chain rule, independently of log_rate.
    g, gd, gdd, v, vd, vdd = sp.symbols("g gd gdd v vd vdd", nonzero=True)
    Xd = vd/g-v*gd/g**2
    Xdd = vdd/g-2*vd*gd/g**2-v*gdd/g**2+2*v*gd**2/g**3
    transformed = sp.expand((g**2*Xdd+2*g*gd*Xd+g**2*(data["q"]+1/m.ZETA)*v/g)/g)
    assert sp.factor(transformed-(vdd+(data["q"]+1/m.ZETA-gdd/g)*v)) == 0
    independent_correction = sp.diff(weight, m.u, 2)/(2*weight)-sp.diff(weight, m.u)**2/(4*weight**2)
    assert sp.factor(independent_correction-data[sector]["correction"]) == 0


@pytest.mark.parametrize("time", (-10, -1, 0, sp.Rational(1, 2), 1, 10))
@pytest.mark.parametrize("comoving_squared", (0, sp.Rational(1, 100), 1, 10000))
def test_exact_frequency_floor_including_homogeneous_chart(time, comoving_squared):
    result = m.frequency_bounds(sp.Rational(1, 2000), time, comoving_squared)
    assert result["canonical_frequency_squared"] >= result["uniform_frequency_squared_lower"]
    assert result["uniform_frequency_squared_lower"] == result["q"]+1985
    assert result["stationary_gap_or_nonlinear_cutoff_claim"] is False
    assert (result["chart"] == "three_homogeneous_coordinate_vectors") == (comoving_squared == 0)


def test_holding_physical_momentum_fixed_or_dropping_normalization_is_detected():
    data = m.canonical()
    q0 = sp.Symbol("fixed_physical_q", positive=True)
    wrong_weight = data["a"]**3*m.ZETA*q0/(1+m.ZETA*q0)
    wrong_rate = sp.diff(wrong_weight, m.u)/(2*wrong_weight)
    wrong_U = sp.diff(wrong_rate, m.u)+wrong_rate**2
    assert sp.factor(wrong_U-data["longitudinal"]["correction"]) != 0
    assert data["transverse"]["correction"].subs(m.u, 0) == 2
    assert data["longitudinal"]["correction"].subs({m.u: 0, m.KCOM2: 2000, m.ZETA: sp.Rational(1, 2000)}) == 4


def test_nonunit_physical_normalization():
    conversion = units(3, 2, sp.Rational(1, 1000))
    assert conversion["normalized_zeta"] == sp.Rational(1, 12000)
    result = m.frequency_bounds(conversion["normalized_zeta"], sp.Rational(1, 2), sp.Rational(4375, 256))
    assert result["q"] == 7
    assert result["uniform_frequency_squared_lower"]/4 == 2998
    assert result["uniform_frequency_squared_lower"]/4 == sp.Rational(7, 4)+3000-sp.Rational(15, 4)


@pytest.mark.parametrize("values", ((0, 0, 1), (-1, 0, 1), (sp.Rational(1, 1000), 0, 1),
                                   (sp.Rational(1, 2000), 0, -1), (True, 0, 1),
                                   (0.0005, 0, 1), (sp.Rational(1, 2000), sp.oo, 1)))
def test_invalid_frequency_domains_rejected(values):
    with pytest.raises((TypeError, ValueError)):
        m.frequency_bounds(*values)
