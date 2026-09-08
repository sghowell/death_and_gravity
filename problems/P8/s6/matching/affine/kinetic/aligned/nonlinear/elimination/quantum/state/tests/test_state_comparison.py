"""Exact WKB residual, physical chain rules and convergent energy envelope."""
import pytest
import sympy as sp
from p8_vector_state import comparison, wkb


def test_all_exact_WKB_and_physical_frequency_identities():
    assert all(value == 0 for value in wkb.checks().values())
    assert all(value == 0 for value in wkb.physical_chain_rule_checks().values())


def test_variation_of_constants_and_continuous_norm_proofs():
    for value in comparison.checks().values():
        assert all(entry == 0 for entry in value) if isinstance(value, sp.MatrixBase) else value == 0
    assert all(value is True for value in comparison.proof_checks().values())


def test_explicit_reference_scale_has_a_finite_evolution_error_not_a_full_quantum_verdict():
    d = comparison.physical_bounds(10**12, 1000)
    assert d["physical_lapse_energy_evolution_difference_over_reference_density"] == sp.Rational(3, 125*10**18)
    assert d["pressure_evolution_difference_over_reference_density"] == sp.Rational(27, 625*10**18)
    assert d["full_renormalized_energy_or_all_order_Hadamard_claim"] is False


def test_second_order_reference_does_not_have_the_fourth_order_residual():
    for kind in ("transverse", "longitudinal"):
        # Omitting P4 leaves 2*P4/omega^2 rather than an O(omega^-4)
        # oscillator residual; the all-momentum estimate cannot reuse C.
        d = wkb.frequency(kind)
        lam = wkb.background()["lambda"]
        leading = -d["P2"]**2-wkb.D(d["B2"])/2+sp.Rational(3, 2)*lam*d["B2"]
        assert sp.factor(leading-2*d["P4"]) == 0
        assert d["P4"] != 0


def test_coefficient_bound_is_continuous_not_a_momentum_grid():
    for kind in ("transverse", "longitudinal"):
        for name in ("P2", "P4", "B2", "B4", "A", "B"):
            d = wkb.box_bound(wkb.frequency(kind)[name])
            assert d["absolute_upper"].is_Rational
            assert d["reconstruction"] == 0


@pytest.mark.parametrize("bad", (True, 1.0, "1000", sp.oo, 0, 999))
def test_invalid_or_outside_mass_floor_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        comparison.physical_bounds(10**12, bad)
