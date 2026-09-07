"""Independent field-strength, SEE, proper-clock and continuous-domain audit."""

from fractions import Fraction as Q

import pytest
import sympy as sp
from p8a_maxwell_thermal import bounds, dynamics, state


def test_physical_plane_wave_and_conformal_density_weights():
    eta = sp.diag(1, -1, -1, -1)
    k = sp.Matrix([3, -1, -2, -2])
    polar = sp.Matrix([0, 2, -1, 0])/sp.sqrt(5)
    assert (k.T*eta*k)[0] == 0
    assert (k.T*eta*polar)[0] == 0
    field = k*polar.T-polar*k.T
    scale = sp.Rational(5, 3)
    inverse = eta/scale**2
    square = sum(field[i, j]*(inverse*field*inverse)[i, j]
                 for i in range(4) for j in range(4))
    tensor = -field*inverse*field.T+scale**2*eta*square/4
    assert square == 0
    assert tensor[0, 0]/scale**2 == 9/scale**4
    # Per polarization: hbar/(2k) times the occupation excess2n
    # multiplies this plane-wave tensor; two physical polarizations give2.
    assert sp.Rational(1, 6)*2*9*2 == 6


def test_hbar_temperature_and_spatial_normalization_are_distinct():
    assert state.thermal_amplitude(hbar=4, b_T=2) == sp.pi**2/60
    hb, a, length = sp.symbols("hbar a b_T", positive=True)
    energy_temperature = hb/(a*length)
    rho = sp.pi**2*energy_temperature**4/(15*hb**3)
    assert sp.simplify(rho-hb*sp.pi**2/(15*(a*length)**4)) == 0
    assert sp.simplify(rho.subs({a: 3*a, length: length/3}, simultaneous=True)-rho) == 0


def test_full_thermal_anomaly_stress_satisfies_density_and_spatial_see():
    y, lam, kappa, tau = sp.symbols("y lambda kappa tau", positive=True)
    a4 = 4*(1-4*lam)/(y*y*(1-lam*y*y))
    hb = 1440*sp.pi**2*tau*tau*lam/(31*kappa)
    thermal = 12*(1-4*lam)/(kappa*tau*tau)
    rate = 2*y*y*(1-lam*y*y)/(1-2*lam*y*y)
    data = state.zero_type_d_stress(a4**sp.Rational(1, 4), -y/tau,
                                  -rate/tau**2, hbar=hb, thermal_Q=thermal)
    assert sp.cancel(3*y*y/tau**2-kappa*data["rho"]) == 0
    assert sp.cancel((-2*rate+3*y*y)/tau**2+kappa*data["pressure"]) == 0
    ricci = 6*(-rate+2*y*y)/tau**2
    assert sp.cancel(ricci-kappa*data["trace_FK"]) == 0
    assert sp.cancel(kappa*data["EED"]*tau**2-3*y*y/(1-2*lam*y*y)) == 0


def test_h4_omission_fails_the_same_00_equation():
    h, b = sp.symbols("H b", positive=True)
    radiation_source_over_three = h*h-b*h**4
    assert sp.expand(h*h-radiation_source_over_three) == b*h**4


def test_clock_scale_factor_and_hubble_jets_are_consistent():
    y, lam = sp.symbols("y lambda", positive=True)
    rate = dynamics.velocity(y, lam)
    a4 = 4*(1-4*lam)/(y*y*(1-lam*y*y))
    assert sp.simplify(sp.diff(dynamics.clock(y, lam), y)+1/rate) == 0
    assert sp.cancel(sp.diff(sp.log(a4), y)*rate/4+y) == 0
    direct = y
    for target in dynamics.jet_functions(y, lam):
        assert sp.cancel(target-direct) == 0
        direct = sp.diff(direct, y)*rate


def test_entire_jet_numerator_intervals_have_positive_bernstein_controls():
    # Independent degree-two Bernstein coefficients on z∈[0,1/4].
    controls = [(Q(3), Q(2), Q(11, 8)),
                (Q(2), Q(11, 8), Q(9, 8)),
                (Q(9), Q(25, 4), Q(35, 8))]
    u = sp.Symbol("u", real=True)
    z = u/4
    targets = [6*z*z-8*z+3, 6*z*z-5*z+2, 14*z*z-22*z+9]
    for coefficients, target, upper in zip(controls, targets, (3, 2, 9), strict=True):
        assert min(coefficients) > 0 and max(coefficients) <= upper
        reconstructed = coefficients[0]*(1-u)**2+2*coefficients[1]*u*(1-u)+coefficients[2]*u*u
        assert sp.expand(reconstructed-target) == 0


def test_continuous_history_budget_is_anchored_and_strict():
    result = bounds.history(Q(1, 10**8))
    coefficients = (Q(16, 25), Q(1728, 25), Q(155136, 25), Q(14770176, 25))
    lam = Q(31, 180*10**8)
    expected = [coefficient*lam for coefficient in coefficients]
    assert list(map(sp.Rational, expected)) == result["C3_error_upper"]
    assert all(error < cap for error, cap in zip(expected, (Q(1, 100), Q(1, 2), 3, 16), strict=True))
    assert result["anchored_at_observer"]
    assert result["actual_observer_H0_times_tau"] == 2
    assert not result["A18_future_bounds_proved_beyond_this_branch_endpoint"]


def test_endpoint_lower_bound_uses_a_regular_subinterval():
    lam = dynamics.dimensionless_coupling(Q(1, 10**8))
    assert lam < sp.Rational(1, 64)
    # For2<=y<=4, D>1/2 and1-lambda*y²<=1, so1/F>1/(4y²).
    y = sp.Symbol("y", positive=True)
    assert sp.integrate(1/(4*y*y), (y, 2, 4)) == sp.Rational(1, 16)
    out = dynamics.endpoint_data(Q(1, 10**8))
    assert out["critical_y"] > 4
    assert out["endpoint_x_strict_lower"] == sp.Rational(1, 16)
    assert out["endpoint_x_strict_upper"] == sp.Rational(1, 4)
    assert out["a_endpoint_fourth"] == 16*lam*(1-4*lam) > 0


@pytest.mark.parametrize("y", [100000, 110000])
def test_critical_or_high_branch_is_rejected_even_when_algebraic_a4_is_positive(y):
    delta = Q(90, 31*100000**2)
    lam = Q(1, 2*100000**2)
    assert 4*(1-4*lam)/(y*y*(1-lam*y*y)) > 0
    with pytest.raises(ValueError, match="low-curvature"):
        dynamics.branch_point(y, delta)


def test_named_scheme_is_not_a_scalar_gamma_or_temperature_parameter():
    assert state.prescription(beta_m=0, cosmological_constant=0)["scalar_gamma_used"] is False
    with pytest.raises(ValueError):
        state.prescription(beta_m=1, cosmological_constant=0)
    with pytest.raises(TypeError):
        state.prescription(gamma=0, cosmological_constant=0)


def test_actual_thermal_branch_does_not_claim_an_all_state_see_or_quantum_endpoint():
    assert not state.calibration()["all_Hadamard_states_solve_this_metric"]
    assert dynamics.calibration()["actual_branch_has_positive_EED"]
    assert not dynamics.calibration()["this_branch_is_a_new_QEI_only_endpoint_argument"]
    assert not dynamics.endpoint_data(Q(1, 10**8))["fundamental_EFT_validity_at_endpoint_proved"]
