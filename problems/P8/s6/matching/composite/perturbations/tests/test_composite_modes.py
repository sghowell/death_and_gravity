import pytest
import sympy as sp
from p8_composite_modes import jets, model, shift, tensor, vector


@pytest.mark.parametrize("module", [tensor, shift, vector, jets])
def test_exact_literal_and_background_identities(module):
    assert module.checks()
    assert all(sp.simplify(value) == 0 for value in module.checks().values())


def test_scale_restoration_keeps_two_independent_mass_scales():
    assert all(value == 0 for value in jets.scale_checks().values())


@pytest.mark.parametrize("module", [tensor, shift, vector])
def test_genuine_omission_controls(module):
    assert all(sp.simplify(value) != 0 for value in module.negative_controls().values())


def test_source4_8_elimination_does_not_import_unstated4_12_field_redefinition():
    result = vector.source_normalization_checks()
    assert result["4_11_C_matches_literal"] == 0
    assert result["printed_4_12_K_over_literal"] == 2
    assert sp.cancel(result["printed_4_12_U_over_literal"]
                     -2*model.F*model.y**2/(model.G+model.F*model.y**2)) == 0
    assert result["printed_frequency_over_literal"] != 1


def test_positive_shift_requires_stated_source_domain():
    assert all(value.is_positive for value in shift.positive_factors().values())
    d = shift.derive()
    initial = {model.alpha: 1, model.beta: 1, model.y: 1, model.c: 1,
               model.rho: model.M4/2, model.p: model.M4/2,
               **dict(zip(model.betas, (0, 0, 1, 0, 0), strict=True))}
    assert sp.simplify(d["Xi"].subs(initial, simultaneous=True)-6*model.M4) == 0


@pytest.mark.parametrize("acceleration,expected", [(sp.Rational(7, 36), -sp.Rational(1, 108)),
                                                  (4, -sp.Rational(23, 18)),
                                                  (sp.Rational(4, 25), sp.Rational(1, 450)),
                                                  (sp.Rational(1, 6), 0)])
def test_actual_clock_vector_jet_and_scale_guard(acceleration, expected):
    d = jets.derive()
    assert d["physical_speed_quadratic"].subs(d["A"], acceleration) == expected


def test_vector_kinetic_survives_zero_algebraic_mass():
    d = vector.derive()
    assert d["K"].subs(d["mu"], 0).is_positive is True
    assert d["U"].subs(d["mu"], 0) == 0
    assert d["K_T"].subs(d["mu"], 0).is_positive is True


def test_negative_shift_coefficient_is_not_covered_by_positive_spring_theorem():
    d = vector.derive()
    bad = {model.G: 1, model.F: 1, model.y: 1, model.c: 1,
           model.a: 1, model.Ng: 1, d["k"]: 2, d["Xi"]: -1}
    assert d["K"].subs(bad) < 0


def test_two_clocks_must_not_be_interchanged_in_normalization():
    t = sp.Symbol("t", real=True)
    kt = sp.exp(2*t)
    # Physical T=exp(t): K_T=K*Ne=exp(3t)=T^3.
    T = sp.Symbol("T", positive=True)
    physical_frequency = vector.canonical_frequency(T**3, 0, T)
    naively_rescaled = vector.canonical_frequency(kt, 0, t).subs(t, sp.log(T))/T**2
    assert sp.simplify(physical_frequency+sp.Rational(3, 4)/T**2) == 0
    assert sp.simplify(physical_frequency-naively_rescaled) == 1/(4*T**2)
