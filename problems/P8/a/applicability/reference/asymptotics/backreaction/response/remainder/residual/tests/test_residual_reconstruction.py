import pytest
import sympy as sp
from p8a_remainder import stress as prior_stress
from p8a_residual import independent, reconstruction


def test_residual_generic_trace_conservation_anomaly_and_finite_tensor():
    assert all(sp.simplify(value) == 0 for value in reconstruction.identities().values())


def test_residual_generic_trace_matches_pinned_raw_scheme():
    eta = sp.Symbol("eta", positive=True)
    a = sp.Function("a", positive=True)(eta)
    k = sp.Function("K")(eta)
    hbar = sp.Symbol("hbar", positive=True)
    gamma = sp.Symbol("gamma", real=True)
    wick = -hbar*k/(8*sp.pi**2*a**2)
    assert sp.simplify(reconstruction.trace(a, wick, eta, hbar=hbar, gamma=gamma)
                       -prior_stress.born_trace(a, k, eta, hbar=hbar, gamma=gamma)) == 0


def test_residual_pinned_nonlinear_error_is_exactly_the_linear_reconstruction_difference():
    eta = sp.Symbol("eta", positive=True)
    a = sp.Function("a", positive=True)(eta)
    s, r, j, jr, ell = (sp.Function(name)(eta) for name in ("S_B", "Rmode", "J_B", "J_R", "L"))
    actual = reconstruction.components(a, s+r, j+jr, ell, eta)
    born = reconstruction.components(a, s, j, ell, eta)
    expected = prior_stress.components(a, r, jr, eta)
    for name in actual:
        assert sp.simplify(actual[name]-born[name]-expected[name]/sp.pi**2) == 0


def test_residual_actual_radiation_past_fixes_homogeneous_energy_datum():
    eta, normalization, hbar = sp.symbols("eta A hbar", positive=True)
    values = reconstruction.components(normalization*eta, 0, 0, 0, eta, hbar=hbar)
    denominators = {"density": 960, "pressure": 576, "EED": 320}
    for name, value in values.items():
        assert sp.simplify(value-hbar/(denominators[name]*sp.pi**2*normalization**4*eta**8)) == 0
    assert reconstruction.fix_radiation_constant(normalization*eta, values["density"], values["density"]) == 0
    shifted = reconstruction.components(normalization*eta, 0, 0, 0, eta, radiation_constant=1)
    wrong = sp.simplify(shifted["density"]-reconstruction.components(normalization*eta, 0, 0, 0, eta)["density"])
    assert wrong == 1/(normalization*eta)**4


def test_residual_histories_have_the_exact_derivatives_and_initial_values():
    eta, start = sp.symbols("eta start", positive=True)
    a, s = eta**2, eta**5
    values = reconstruction.histories(a, s, eta, start)
    u, h = reconstruction.potential(a, eta), sp.diff(a, eta)/a
    assert sp.simplify(sp.diff(values["J"], eta)-sp.diff(u, eta)*s) == 0
    assert sp.simplify(sp.diff(values["L"], eta)-h*u**2) == 0
    assert all(value.subs(eta, start).doit() == 0 for value in values.values())


def test_residual_both_histories_and_anomaly_have_nonzero_exclusion_controls():
    controls = reconstruction.controls()
    assert controls.pop("wrong_past_constant_can_still_be_conserved") is True
    assert all(value != 0 for value in controls.values())


def test_residual_generic_gamma_is_not_mistaken_for_a_zero_curvature_term():
    eta = sp.Symbol("eta", positive=True)
    a = eta**2
    base = reconstruction.components(a, 0, 0, 0, eta)
    changed = reconstruction.components(a, 0, 0, 0, eta, gamma=1)
    tensor = reconstruction.curvature_tensor(a, eta)
    assert sp.simplify(tensor["density"]) != 0
    for name in base:
        assert sp.simplify(changed[name]-base[name]-tensor[name]) == 0


@pytest.mark.parametrize("power,s_power", [(2, 5), (2, 7), (3, 5), (4, 6)])
def test_residual_independent_laurent_reconstruction_matches_public_api(power, s_power):
    audit = independent.powerlaw_control(power, s_power)
    eta = sp.Symbol("eta", positive=True)
    s = eta**s_power
    j = sp.Rational(2*power*(power-1), s_power-2)*eta**(s_power-2)
    ell = -sp.Rational(power**3*(power-1)**2, 4)*eta**-4
    values = reconstruction.components(eta**power, s, j, ell, eta)
    assert sp.simplify(values["density"].subs(eta, 1)*sp.pi**2) == sp.Rational(audit["rho_at_eta_one_in_hbar_over_pi2_units"])
    assert sp.simplify(values["pressure"].subs(eta, 1)*sp.pi**2) == sp.Rational(audit["pressure_at_eta_one_in_hbar_over_pi2_units"])
    assert not audit["control_is_a_physical_flat_past_state"]


def test_residual_effective_S_preserves_exact_zero_inputs():
    assert reconstruction.effective_s(0, 0, 0) == 0
    assert not reconstruction.effective_s(1, 2, 3).has(sp.Float)
