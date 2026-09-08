"""Focused differentiated vector evolution, subtraction and clock tests."""
import sympy as sp
from p8_vector_evolution import estimates, majorants, mixing, ward


def assert_zero(values):
    for value in values.values():
        assert all(entry == 0 for entry in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_actual_coefficient_derivatives_and_inverse_reference_chain_rule():
    assert_zero(majorants.checks())


def test_exact_polynomials_replay_frozen_residual_and_full_subtraction():
    assert_zero(majorants.replay_checks())


def test_phase_removal_products_and_physical_time_derivative():
    assert_zero(mixing.checks())


def test_continuous_mixing_and_integrated_derivative_bounds():
    assert all(value is True for value in mixing.proof_checks().values())
    assert mixing.MIXING_CONSTANT == 676000000


def test_ordinary_conservation_in_dimensional_modes_and_full_adiabatic_orders():
    assert_zero(ward.ordinary_exact_mode_balance())
    assert_zero(ward.ordinary_adiabatic_conservation())
    assert_zero(ward.local_mass_balance())


def test_direct_clock_variation_includes_the_mass_energy_exchange():
    assert_zero(ward.clock_variation())
    assert_zero(ward.checks())


def test_actual_local_clock_source_not_spurious_vacuum_energy_derivative():
    u = ward.wkb.u
    expected = {0: -44*u/(9*(1+u**2)**4),
                1: 128*u*(79*u**2+36)/(243*(1+u**2)**6),
                2: -64*u*(449*u**4+898*u**2+81)/(81*(1+u**2)**8)}
    for n, value in expected.items():
        assert sp.factor(ward.local_coefficients()[n]["negative_clock_source"]-value) == 0


def test_physical_derivative_and_clock_source_example_is_not_full_solution():
    assert all(value is True for value in estimates.proof_checks().values())
    d = estimates.physical_bounds(10**12, 1000)
    assert d["all_order_Hadamard_or_quantum_corrected_solution_claim"] is False
    doubled = estimates.physical_bounds(2*10**12, 1000)
    key = "clock_source_matched_over_reference_scalar_equation"
    assert 4*doubled[key] == d[key]


def test_polynomial_majorant_differentiates_all_independent_factors():
    p, r, si = majorants.p, majorants.r, majorants.si
    expression = -2*p**2+3*p*r-si
    expected = 2*p**2+3*p*r+si
    assert sp.expand(majorants.polynomial_majorant(expression)-expected) == 0
    envelopes = {variable: (sp.Integer(1), sp.Integer(1)) for variable in majorants.VARIABLES}
    assert majorants.value_and_slope(expression, envelopes) == {"value_upper": 6, "derivative_upper": 11}
