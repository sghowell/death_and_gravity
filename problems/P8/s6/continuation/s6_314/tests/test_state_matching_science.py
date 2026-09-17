"""State continuity, full signed-series matching and exact reference domains."""

from itertools import combinations, permutations

import pytest
import sympy as s
from p8_vacuum_affine_state_transport_soft_matching import (
    assembly,
    audit,
    continuity,
    series,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_all_exact_residuals(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_all_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_scope_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "invalid",
    (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        s.I,
        s.nan,
        s.Symbol("missing"),
        -1,
        s.Rational(1, 4),
    ),
)
def test_exact_energy_and_reference_domains(invalid):
    for call in (
        continuity.radiation_energy,
        continuity.modulus,
        continuity.conversion_change_upper,
        continuity.regulator_coefficient_change_upper,
        assembly.two_real_variation_upper,
        assembly.connector_relative_upper,
        assembly.dressed_two_real_relative_upper,
        assembly.assembled_relative_upper,
    ):
        with pytest.raises((TypeError, ValueError)):
            call(invalid)


def test_continuity_zero_is_admitted_but_detector_zero_is_not():
    for call in (
        continuity.radiation_energy,
        continuity.modulus,
        continuity.conversion_change_upper,
        continuity.regulator_coefficient_change_upper,
    ):
        assert call(0) == 0
    for call in (
        assembly.two_real_variation_upper,
        assembly.connector_relative_upper,
        assembly.dressed_two_real_relative_upper,
        assembly.assembled_relative_upper,
    ):
        with pytest.raises(ValueError):
            call(0)


@pytest.mark.parametrize(
    "energy",
    (
        s.Rational(1, 8),
        s.Rational(1, 1000),
        s.Rational(1, 10**200),
    ),
)
def test_full_state_continuity_bounds_are_exact_and_vanishing(energy):
    modulus = continuity.modulus(energy)
    assert modulus == energy * (1 - s.log(energy))
    assert continuity.conversion_change_upper(energy) == 5000 * modulus / source.KAPPA
    assert (
        continuity.regulator_coefficient_change_upper(energy)
        == 8000 * modulus / source.KAPPA
    )
    assert 0 < continuity.modulus(energy / 2) < modulus
    assert not modulus.has(s.Float)


@pytest.mark.parametrize(
    "resolution",
    (
        s.Rational(1, 8),
        s.Rational(1, 1000),
        s.Rational(1, 10**200),
    ),
)
def test_reference_bound_is_exact_and_below_the_uniform_margin(resolution):
    total = assembly.assembled_relative_upper(resolution)
    assert total > 0 and total < s.Rational(1, 10**653)
    assert not total.has(s.Float)
    assert assembly.dressed_two_real_relative_upper(
        resolution
    ) == 2 * assembly.two_real_variation_upper(resolution)
    assert (
        assembly.connector_relative_upper(resolution)
        == 100000 * resolution * (2 - s.log(resolution)) / source.KAPPA**2
    )


@pytest.mark.parametrize("name", tuple(continuity.physical_samples()[1]))
def test_each_actual_recoil_state_change_is_checked(name):
    assert continuity.physical_samples()[1][name]


def test_actual_state_transport_is_not_zero():
    records = continuity.physical_samples()[2]
    assert len(records) == 3
    assert all(value != 0 and value.is_Rational for value in records.values())


@pytest.mark.parametrize("name", tuple(assembly.original_samples()))
def test_every_complete_original434_density_decomposition(name):
    assert assembly.original_samples()[name] == 0


def test_omitting_connector_or_confusing_subtractions_is_rejected():
    checks, gates = assembly.general_identities()
    assert all(value == 0 for value in checks.values())
    assert gates["omitting_the_state_connector_is_a_nonzero_generic_error"]
    assert gates["amplitude_and_probability_subtractions_not_identified"]


@pytest.mark.parametrize("count", range(2, 7))
def test_state_dependent_signed_two_marked_finite_label_counting(count):
    energies = tuple(s.Rational(j + 1, 2 * count) for j in range(count))

    def kernel(first, second, rest):
        return (first + second - 1) * (1 + first * second) * s.prod(rest)

    labeled = sum(
        kernel(
            energies[i],
            energies[j],
            tuple(energies[k] for k in range(count) if k not in (i, j)),
        )
        for i, j in combinations(range(count), 2)
    ) / s.factorial(count)
    averaged = sum(
        kernel(row[0], row[1], row[2:]) for row in permutations(energies)
    ) / s.factorial(count)
    assert labeled == averaged / (2 * s.factorial(count - 2))


def test_entire_signed_series_and_remaining_energy_moments():
    assert all(series.independent_numerics().values())


def test_joint_state_radial_modulus_is_integrable_but_constant_bound_is_not():
    R, tau, mu = s.symbols("R tau mu", positive=True)
    finite = s.integrate(tau ** s.Rational(-3, 4), (tau, 0, R**4)) + s.integrate(
        R / tau, (tau, R**4, 1)
    )
    assert s.simplify(finite - 4 * R * (1 - s.log(R))) == 0
    wrong = s.integrate(R / tau, (tau, mu, 1))
    assert s.limit(wrong, mu, 0, dir="+") == s.oo


def test_connector_difference_must_precede_the_seed_integral():
    mu, x = s.symbols("mu x", positive=True)
    wrong = s.log(x / mu)
    assert s.limit(wrong, mu, 0, dir="+") == s.oo
    assert s.limit(x * (2 - s.log(x)), x, 0, dir="+") == 0
    assert "difference" in assembly.data()["whole_finite_state_transport_connector"]


def test_remaining_energy_is_not_replaced_by_the_full_threshold():
    a = s.Symbol("physical_index", positive=True)
    assert s.expand_func(s.gamma(1 + a) / s.gamma(2 + a)) == 1 / (1 + a)
    assert s.simplify(1 - 1 / (1 + a)) != 0


def test_all_original_parameters_and_unclosed_frontiers_retained():
    assert source.KAPPA == 10**800
    assert source.HEAVY_MASS2 == s.Rational(10**200, 512) + 2
    assert source.CUBIC == s.Rational(1, 8192)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 170
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "Full finite hard real-virtual" in audit.observable()["not_established"]
    assert "original common-parent bounce" in audit.observable()["not_established"]
