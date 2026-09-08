"""Independent constrained retarded mass-response checks."""
import pytest
import sympy as sp
from p8_vector_retarded import bubble, canonical, causal, clock, flat, modes, verify


@pytest.mark.parametrize("module", (canonical, modes, bubble, flat, causal, clock))
def test_exact_science(module):
    verify.affine.certify_residuals(module.checks())


@pytest.mark.parametrize("sector", bubble.CHANNELS)
def test_explicit_polarizations_match_general_momentum_form(sector):
    assert modes.angular(bubble.direct(sector)-bubble.closed(sector)) == 0


def test_contact_omission_is_detected():
    data = flat.controls()
    point = {flat.alpha: 1, flat.beta: 0, flat.m: 1, flat.w: 2}
    assert data["omitted_longitudinal_contact_per_momentum"].subs(point) == sp.Rational(3, 4)
    assert data["omitted_contact_UV_residue_times_pi_squared"].subs(point) == sp.Rational(3, 64)


def test_literal_Kubo_covariance_and_realness():
    data = causal.covariance()
    assert data["Kubo_equals_exact_canonical_covariance"] == 0
    assert data["real_response"] == 0


def test_native_index_validation_after_warming_cache():
    modes.data(0, 0)
    modes.polarization(0, 0)
    bubble.vertex(0)
    flat.coefficient(0)
    for value in (False, 0.0, sp.Integer(0)):
        with pytest.raises(ValueError):
            modes.data(value, 0)
        with pytest.raises(ValueError):
            modes.polarization(0, value)
        with pytest.raises(ValueError):
            bubble.vertex(value)
        with pytest.raises(ValueError):
            flat.coefficient(value)


def test_proof_and_rejection_controls():
    assert all(verify.proof_checks().values())
    assert verify.controls()["rejected_inputs"] == 140
