"""Physical photon counting, thermodynamic units and named-scheme controls."""

import pytest
import sympy as sp
from p8a_maxwell_thermal import state


def test_thermal_physical_field_and_clock_identities():
    assert all(sp.simplify(value) == 0 for value in state.identities().values())
    assert state.thermal_amplitude(hbar=1, b_T=2) == sp.pi**2/240


def test_one_scalar_polarization_is_not_the_photon_coefficient():
    assert state.calibration()["physical_transverse_polarizations"] == 2
    assert state.thermal_amplitude(hbar=1, b_T=1) != sp.pi**2/30


def test_scalar_gamma_is_not_an_input_for_the_photon_finite_scheme():
    with pytest.raises(TypeError):
        state.prescription(gamma=0, cosmological_constant=0)
    for beta, lam in ((1, 0), (0, 1), (-1, 0)):
        with pytest.raises(ValueError):
            state.prescription(beta_m=beta, cosmological_constant=lam)


def test_nonzero_hbar_temperature_and_conformal_weights_cannot_be_dropped():
    # b_T is a length: hbar=2,b_T=1,a=3 has rho_thermal=2*pi²/(15*81).
    data = state.zero_type_d_stress(3, 0, 0, hbar=2, thermal_Q=2*sp.pi**2/15)
    assert data["rho"] == 2*sp.pi**2/1215
    assert data["rho"] != 2*sp.pi**2/135  # Wrong covariant a^-2 as physical a^-4.
    assert data["rho"] != sp.pi**2/(15*2**3*81)  # Wrong b_T=1/(k_B T) with hbar retained.


@pytest.mark.parametrize("bad", [True, 0, -1, 0.1, sp.oo, sp.nan, sp.Symbol("b_T")])
def test_numeric_thermal_parameters_are_exact_positive_finite(bad):
    with pytest.raises((TypeError, ValueError)):
        state.thermal_amplitude(hbar=1, b_T=bad)
