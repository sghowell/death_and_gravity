"""Physical two-polarization thermal state and an explicitly named prescription."""

import sympy as sp
from p8a_maxwell.domain import nonnegative, rational


def prescription(*, beta_m, cosmological_constant):
    """Only this child's beta_M=0,Lambda=0 specialization is admitted.

The thermodynamic b_T is a conformal-time length, not beta_M or an older
scalar gamma. No gamma keyword substitutes for the physical photon scheme.
    """
    beta, lam = map(rational, (beta_m, cosmological_constant))
    if beta != 0 or lam != 0:
        raise ValueError("the actual thermal branch is certified only for beta_M=0 and Lambda=0")
    return {"beta_M": beta, "Lambda": lam, "additional_source": "none",
            "scalar_gamma_used": False, "independent_curvature_couplings_added": False}


def thermal_amplitude(*, hbar, b_T):
    """Exact Q for n(k)=(exp(b_T*k)-1)^-1; hbar>0 and b_T>0."""
    hbar, length = (nonnegative(value, positive=True) for value in (hbar, b_T))
    return hbar*sp.pi**2/(15*length**4)


def zero_type_d_stress(a, hubble, hubble_dot, *, hbar, thermal_Q):
    """Formal physical rho,p,EED and FK trace on an arbitrary smooth metric."""
    a, h, hd, hb, q = map(sp.sympify, (a, hubble, hubble_dot, hbar, thermal_Q))
    r = 31*hb/(480*sp.pi**2)
    rho = q/a**4+r*h**4
    pressure = q/(3*a**4)-r*h**4-4*r*h*h*hd/3
    return {"rho": rho, "pressure": pressure,
            "EED": sp.expand((rho+3*pressure)/2),
            "trace_FK": sp.expand(rho-3*pressure)}


def identities():
    a, scale, length, hb, kb, temp = sp.symbols("a scale b_T hbar k_B temperature", positive=True)
    h, hd, q = sp.symbols("H Hdot Q", real=True)
    data = zero_type_d_stress(a, h, hd, hbar=hb, thermal_Q=q)
    derivative_rho = sp.diff(data["rho"], a)*a*h+sp.diff(data["rho"], h)*hd
    blackbody = hb*sp.pi**2/(15*length**4)
    return {"two_physical_polarization_measure": sp.simplify(2*4*sp.pi/(2*sp.pi)**3-1/sp.pi**2),
            "Bose_integral_Gamma_zeta": sp.simplify(sp.gamma(4)*sp.zeta(4)-sp.pi**4/15),
            "thermal_proper_density_from_energy_temperature": sp.simplify(
                (sp.pi**2*(kb*temp)**4/(15*hb**3)).subs(temp, hb/(kb*a*length))-blackbody/a**4),
            "conformal_coordinate_rescaling": sp.simplify(
                (blackbody/a**4).subs({a: scale*a, length: length/scale}, simultaneous=True)-blackbody/a**4),
            "thermal_plus_anomaly_conservation": sp.simplify(derivative_rho+3*h*(data["rho"]+data["pressure"])),
            "thermal_difference_is_traceless": sp.simplify(data["trace_FK"]-data["trace_FK"].subs(q, 0)),
            "thermal_difference_EED_equals_density": sp.simplify(data["EED"]-data["EED"].subs(q, 0)-q/a**4),
            "nonzero_vacuum_anomaly_retained": sp.simplify(data["trace_FK"].subs(q, 0)-31*hb*(h**4+h*h*hd)/(120*sp.pi**2))}


def calibration():
    return {"physical_transverse_polarizations": sp.Integer(2),
            "Q_in_hbar_pi_squared_over_bT_fourth_units": sp.Rational(1, 15),
            "rho_anomaly_in_hbar_over_pi_squared_units": sp.Rational(31, 480),
            "pressure_anomaly_H4_coefficient": -sp.Rational(31, 480),
            "pressure_anomaly_H2_Hdot_coefficient": -sp.Rational(31, 360),
            "thermodynamic_parameter_is_conformal_time_length": True,
            "physical_temperature_relation": "k_B*T_physical=hbar/(a*b_T), c=1",
            "positive_Hadamard_state_constructed": True,
            "state_is_a_physical_two_polarization_quasifree_example": True,
            "all_Hadamard_states_solve_this_metric": False,
            "prescription": prescription(beta_m=0, cosmological_constant=0)}
