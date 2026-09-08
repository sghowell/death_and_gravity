"""Finite prepared-response identities and continuous source-norm controls."""
import pytest
import sympy as sp
from p8_vector_mass_remainder import (
    evolution,
    proofs,
    source,
    tail,
    tangent,
    verify,
)


def test_all_exact_identities():
    verify.affine.certify_residuals(proofs.residuals())


@pytest.mark.parametrize("sector", ("T", "L"))
def test_varied_residual_and_readout_tail_cancel_low_orders(sector):
    assert all(value == 0 for value in tangent.reference(sector)["varied_residual_low_coefficients"].values())
    assert all(value == 0 for value in tail.reference_tail(sector)["low_tail_numerator_coefficients"].values())


@pytest.mark.parametrize("sector", ("T", "L"))
def test_tenth_source_derivative_control(sector):
    assert tangent.reference(sector)["tenth_source_derivative_fixture"] == -sp.Rational(7, 20736)


@pytest.mark.parametrize("order", (1, 2, 3, 4))
def test_regular_zero_momentum_tangent(order):
    assert source.clean((tangent.varied_P("T", order)-tangent.varied_P("L", order)).subs(source.z, 0)) == 0


def test_source_norm_rejects_uncontrolled_jets_and_nonlinearity():
    for value in (source.n[11], source.n[0]*source.n[1], sp.Float(1)*source.n[0]):
        with pytest.raises(ValueError):
            source.linear_bound(value)


def test_nonzero_initial_mixing_retained():
    data = evolution.constants()["frozen_initial_and_evolved_mixing"]
    assert data["initial_mixing_envelopes"][6] == 1813229
    assert data["evolution_mixing_envelope"] == 5347035781757616


def test_complete_prepared_mass_response_bound():
    data = evolution.bound(10**24, 1000)
    assert 0 < data["nonlocal_subtracted_C10_to_C0_upper_bound"] < sp.Rational(1, 10**43)
    assert 0 < data["complete_mass_response_C10_to_C0_upper_bound"] < sp.Rational(1, 10**39)


def test_native_reference_order_validation_after_warming_cache():
    tangent.coefficient("T", 1)
    source.baseline("T", 1)
    for value in (True, 1.0, sp.Integer(1)):
        with pytest.raises(ValueError):
            tangent.coefficient("T", value)
        with pytest.raises(ValueError):
            source.baseline("T", value)


def test_continuous_proofs_and_rejected_inputs():
    assert all(proofs.checks().values())
    assert verify.controls()["rejected_inputs"] == 129
