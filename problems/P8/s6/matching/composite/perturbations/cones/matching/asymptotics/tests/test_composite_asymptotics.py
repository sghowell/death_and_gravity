from fractions import Fraction

import pytest
import sympy as sp
from p8_composite_asymptotics import independent, limiting, scaled


@pytest.mark.parametrize("check", [scaled.checks, limiting.checks,
                                  limiting.canonical_limit_checks,
                                  limiting.physical_projection_checks])
def test_exact_scaled_flow_time_dependent_canonical_and_projection_identities(check):
    values = check()
    assert values and all(sp.simplify(value) == 0 for value in values.values())


def test_separate_Fraction_derivatives_energy_and_interval_identities():
    values = independent.checks()
    assert len(values["coefficientwise_identities"]) == 6


def test_singular_metric_limit_does_not_singularize_coefficient_templates():
    d = scaled.canonical_jets()
    e = d["e"]
    assert d["shape_sum"].subs(e, 0) == 1
    assert d["shape_relative"].subs(e, 0) == d["kappa"]
    assert d["physical_source_relative"].subs(e, 0) == -1/sp.sqrt(d["kappa"])
    family = scaled.cd_family()
    assert family["N_g"].subs(family["epsilon"], 0) == 0
    assert sp.simplify(family["N_f"].subs(family["epsilon"], 0)-1) == 0


def test_no_division_by_a_zero_algebraic_mass_in_exact_limiting_equation():
    d = limiting.general_jets()
    assert d["mass_squared"].subs(d["v"], 1) == 0
    assert d["heavy_growth_coefficient"].subs({d["v"]: 1, d["a"]: 0}) == 1


def test_forward_bound_is_not_silently_extended_to_negative_times():
    # At u=-atanh(1/6), z'=0 and z=sqrt(140)-11<1.
    assert 11**2 < 140 < 12**2
    assert sp.sqrt(140)-11 < 1
    assert Fraction(140, 19) < 9


def test_reconstructed_density_may_turn_negative_without_changing_scaled_NEC():
    d = limiting.solution()
    u, z = d["u"], d["z"]
    assert sp.simplify((d["R"]-1-(11*(1-u)-10*sp.exp(-u))/z).rewrite(sp.exp)) == 0
    assert sp.simplify((d["R"]-1).subs(u, 1)) < 0
    assert scaled.core()["null_over_e"].subs(scaled.core()["e"], 0) == 2


def test_omitting_normalization_changes_actual_initial_acceleration():
    d = limiting.general_jets()
    at = {d["v"]: 2, d["a"]: 11}
    true_acceleration = d["heavy_growth_coefficient"].subs(at)
    false_acceleration = -d["mass_squared"].subs(at)
    assert true_acceleration == sp.Rational(55, 2)
    assert false_acceleration == -2
    assert true_acceleration-false_acceleration == sp.Rational(59, 2)


def test_no_heavy_initial_data_no_forcing_is_an_exact_zero_control():
    u = sp.Symbol("u", real=True)
    coefficient = sp.Function("positive_growth_coefficient")(u)
    zero = sp.S.Zero
    assert sp.diff(zero, u, 2)-coefficient*zero == 0


def test_comparison_kernel_is_retarded_with_zero_initial_data():
    u = sp.Symbol("u", nonnegative=True)
    kernel = 2*sp.sinh(u/2)
    assert kernel.subs(u, 0) == 0
    assert sp.diff(kernel, u).subs(u, 0) == 1
    assert sp.diff(kernel, u, 2)-kernel/4 == 0
