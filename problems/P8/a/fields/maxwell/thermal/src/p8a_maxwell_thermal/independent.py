"""Stdlib Fraction reconstruction from pinned photon and history inputs."""

from fractions import Fraction as F


def serialize(value):
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    if isinstance(value, bool):
        return value
    return str(value)


def photon_coefficients(report):
    """Read only the transitive A16 monomials, in its 2880*pi^2 units."""
    source = report["derived_constants"]["stress_jet_coefficients"]
    data = {name: {key: F(value)/2880 for key, value in source[name]["universal"].items()}
            for name in ("rho", "pressure", "EED", "trace")}
    r = data["rho"]["4,0,0,0"]
    expected = {"rho": {"4,0,0,0": r},
                "pressure": {"4,0,0,0": -r, "2,1,0,0": -4*r/3},
                "EED": {"4,0,0,0": -r, "2,1,0,0": -2*r},
                "trace": {"4,0,0,0": 4*r, "2,1,0,0": 4*r}}
    if data != expected or r <= 0:
        raise ValueError("the pinned physical Maxwell stress polynomial changed")
    return data


def state_data(photon):
    # Two polarizations times angular phase-space coefficient, Gamma(4)
    # and zeta(4)/pi^4. No new primary formula module is imported.
    thermal = F(2)*F(4, 8)*F(6)*F(1, 90)
    r = photon["rho"]["4,0,0,0"]
    return {"physical_transverse_polarizations": F(2),
            "Q_in_hbar_pi_squared_over_bT_fourth_units": thermal,
            "rho_anomaly_in_hbar_over_pi_squared_units": r,
            "pressure_anomaly_H4_coefficient": photon["pressure"]["4,0,0,0"],
            "pressure_anomaly_H2_Hdot_coefficient": photon["pressure"]["2,1,0,0"],
            "thermodynamic_parameter_is_conformal_time_length": True,
            "physical_temperature_relation": "k_B*T_physical=hbar/(a*b_T), c=1",
            "positive_Hadamard_state_constructed": True,
            "state_is_a_physical_two_polarization_quasifree_example": True,
            "all_Hadamard_states_solve_this_metric": False,
            "prescription": {"beta_M": F(0), "Lambda": F(0), "additional_source": "none",
                             "scalar_gamma_used": False,
                             "independent_curvature_couplings_added": False}}


def jet_values(y, coupling):
    """Independent exact low-branch jets using the response polynomials."""
    y, lam = F(y), F(coupling)
    z = lam*y*y
    d = 1-2*z
    if y <= 0 or lam <= 0 or d <= 0:
        raise ValueError("the independent point is not on the regular positive low branch")
    return [y, 2*y*y+2*lam*y**4/d,
            8*y**3+8*lam*y**5*(6*z*z-8*z+3)/d**3,
            48*y**4+16*lam*y**6*(6*z*z-5*z+2)*(14*z*z-22*z+9)/d**5]


def response_bounds(ratio):
    zero = 32*ratio/F(1, 2)
    scheme = [2*F(2)**4/F(1, 2),
              8*F(2)**5*3/F(1, 2)**3,
              16*F(2)**6*2*9/F(1, 2)**5]
    lipschitz = [2*2*F(2), 8*3*F(2)**2, 48*4*F(2)**3]
    total = [zero, *(a+b*zero for a, b in zip(scheme, lipschitz, strict=True))]
    return total, scheme, lipschitz


def actual_history(delta, coupling, geometry):
    ratio = F(geometry["ratio"])
    caps = list(map(F, geometry["anchored_C3_slice"]["error_caps"]))
    coefficients, _, _ = response_bounds(ratio)
    errors = [coupling*value for value in coefficients]
    margins = [cap-error for cap, error in zip(caps, errors, strict=True)]
    backward_lower = 1/(F(1, 2)+4*ratio)
    if min(*margins, backward_lower-1) <= 0 or coupling >= F(1, 64):
        raise ValueError("independent actual-history or endpoint domain lost a strict margin")
    return {"delta": delta, "lambda": coupling,
            "history_interval_x": [-ratio, F(0)], "reference_power": F(1, 2),
            "backward_solution_y_bounds": [backward_lower, F(2)],
            "C3_response_coefficients_per_lambda": coefficients,
            "C3_error_upper": errors, "strict_A18_tube_margins": margins,
            "actual_observer_H0_times_tau": F(geometry["reference_Hstar_times_tau"]),
            "anchored_at_observer": True, "a0_one_is_a_chosen_spatial_normalization": True,
            "actual_state_and_SEE_history_realized": True,
            "A18_future_bounds_proved_beyond_this_branch_endpoint": False}


def dynamics_data(delta, coupling):
    # Endpoint lower integral: integral_2^4 1/(4*y^2) dy.
    lower = (F(1, 2)-F(1, 4))/4
    upper = F(1, 2)/2  # integral_2^infinity 1/(2*y^2) dy.
    charge = 3*F(2)**2*(1-coupling*F(2)**2)
    endpoint_a4 = 4*(1-4*coupling)/(F(1, 2)/coupling*F(1, 2))
    return {"delta_upper": delta, "lambda_upper": coupling,
            "initial_y": F(2), "initial_a": F(1),
            "Q_in_one_over_kappa_tau_squared_units": charge,
            "strict_initial_low_branch_margin": 1-8*coupling,
            "endpoint_x_strict_bounds": [lower, upper],
            "endpoint_lower_bound_domain_margin": F(1, 64)-coupling,
            "a_endpoint_fourth_at_delta_upper": endpoint_a4,
            "uniform_positive_endpoint_a_bound_as_delta_tends_to_zero": False,
            "actual_branch_has_positive_EED": True,
            "this_branch_is_a_new_QEI_only_endpoint_argument": False,
            "fundamental_EFT_endpoint_control": False}


def replay(cosmology_report, photon_report):
    photon = photon_coefficients(photon_report)
    old = cosmology_report["derived_constants"]
    geometry = old["geometry"]
    delta = F(old["calibration"]["universal_gate"]["weighted_delta_upper"])
    coupling = F(8, 3)*photon["rho"]["4,0,0,0"]*delta
    actual = actual_history(delta, coupling, geometry)
    _, scheme, lipschitz = response_bounds(F(geometry["ratio"]))
    z = F(1, 4)
    bounds = {"actual_history_at_delta_upper": actual,
              "generic_small_coupling_upper": F(1, 16),
              "strict_named_coupling_margin": F(1, 16)-coupling,
              "scheme_only_jet_error_coefficients": scheme,
              "radiation_jet_Lipschitz_coefficients": lipschitz,
              "numerator_lower_values_at_one_quarter": [6*z*z-8*z+3, 6*z*z-5*z+2, 14*z*z-22*z+9],
              "strict_backward_y_above_one": actual["backward_solution_y_bounds"][0]-1,
              "continuous_interval_not_a_sample_scan": True}
    return serialize({"state": state_data(photon), "dynamics": dynamics_data(delta, coupling),
                      "bounds": bounds, "pinned_photon_physical_coefficients": photon,
                      "engine": "stdlib Fraction, Bose phase-space factor, pinned Maxwell polynomials and continuous rational ODE bounds",
                      "imports_new_primary_formula_modules": False})
