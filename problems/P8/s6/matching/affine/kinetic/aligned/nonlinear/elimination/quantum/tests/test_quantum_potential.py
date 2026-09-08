"""Dimensional multiplicity, physical mass jets and continuous bounds."""
import pytest
import sympy as sp
from p8_aligned_quantum import bounds, kernel, potential


def test_dimensional_pole_finite_and_scale_identities():
    assert all(value == 0 for value in potential.checks().values())
    assert all(value is True for value in bounds.proof_checks().values())


def test_independent_cutoff_and_three_dimensional_integral_sign_controls():
    assert all(value == 0 for value in potential.independent_sign_controls().values())


def test_missing_dimensional_polarization_term_changes_the_finite_constant():
    data = potential.coefficients()
    # Three four-dimensional scalar determinants give -9/2, whereas
    # the dimensionally continued massive vector gives -5/2.
    wrong = 3*potential.Lscale-sp.Rational(9, 2)
    correct = data["finite_weight"].subs({kernel.a: 1, kernel.b: 1})
    assert sp.factor(correct-wrong) == 2


@pytest.mark.parametrize("mass_factor", (sp.Rational(1, 2), 1, 2, 3))
def test_independent_isotropic_mass_rescaling(mass_factor):
    c = sp.sympify(mass_factor)
    data = potential.coefficients()
    replacement = {kernel.a: c, kernel.b: c}
    assert sp.factor(data["pole_weight"].subs(replacement)-3*c**2) == 0
    expected = 3*c**2*(potential.Lscale+sp.log(c)-sp.Rational(5, 6))
    assert sp.factor(data["finite_weight"].subs(replacement)-expected) == 0


@pytest.mark.parametrize("h", (1, 2, sp.Rational(125, 64)))
def test_independent_clock_N_jets_of_actual_mass_functions(h):
    N = sp.Symbol("N", positive=True)
    p = sp.sqrt((h-1+N**-2)/(4*h))
    masses = kernel.masses()
    replacement = {kernel.a: masses["a"].subs(kernel.geometry.P, p),
                   kernel.b: masses["b"].subs(kernel.geometry.P, p)}
    jets = potential.clock_jets()
    for name in ("pole_weight", "finite_weight"):
        expression = potential.coefficients()[name].subs(potential.Lscale, 0).subs(replacement, simultaneous=True)
        for order, key in ((0, "value"), (1, "N_first"), (2, "N_second")):
            direct = sp.factor(sp.diff(expression, N, order).subs(N, 1))
            assert sp.factor(direct-jets[name][key].subs(jets["h"], h)) == 0


def test_actual_clock_coefficient_second_jets():
    jets = potential.clock_jets()
    h = jets["h"]
    assert sp.factor(jets["pole_weight"]["N_second"]+4*(3645*h-236)/(2187*h**2)) == 0
    assert sp.factor(jets["finite_weight"]["N_second"]-2*(8019*h+4784)/(6561*h**2)) == 0
    assert jets["pole_weight"]["N_first"] != 0
    assert jets["finite_weight"]["N_first"] != 0


def test_explicit_scale_example_keeps_all_normalizations():
    d = bounds.scale_bound(10**12, 1000)
    assert d["m0_over_M"] == sp.Rational(1, 10**9)
    assert d["normalized_zeta"] == sp.Rational(1, 10**6)
    assert d["physical_zeta"] == 10**18
    assert d["local_potential_over_reference_M_squared_over_tau_squared"] == sp.Rational(1, 144*10**12)
    assert d["scale_choice_is_a_cutoff_or_complete_quantum_bound"] is False


@pytest.mark.parametrize("value", (0, -1))
def test_nonpositive_scale_products_are_rejected(value):
    with pytest.raises(ValueError):
        bounds.scale_bound(value, 1)
    with pytest.raises(ValueError):
        bounds.scale_bound(1, value)
