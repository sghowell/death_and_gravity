"""Exact named-prescription photon SEE low branch and its regular open domain."""

import sympy as sp
from p8a_maxwell.domain import nonnegative

from . import state

DELTA_MAX = sp.Rational(1, 10**8)


def dimensionless_coupling(delta):
    delta = nonnegative(delta, positive=True)
    if delta > DELTA_MAX:
        raise ValueError("delta exceeds the named thermal history gate")
    return sp.Rational(31, 180)*delta


def velocity(y, coupling):
    """Formal y'=F(y), y=-tau*H_normal, x=s/tau."""
    y, lam = map(sp.sympify, (y, coupling))
    return 2*y*y*(1-lam*y*y)/(1-2*lam*y*y)


def jet_functions(y, coupling):
    """Formal y derivatives through3, with the exact clock retained."""
    z = sp.Dummy("y", positive=True)
    lam = sp.sympify(coupling)
    f = velocity(z, lam)
    second = sp.factor(f*sp.diff(f, z))
    third = sp.factor(f*sp.diff(second, z))
    return [sp.sympify(y), f.subs(z, y), second.subs(z, y), third.subs(z, y)]


def clock(y, coupling):
    """Formal dimensionless primitive Phi; x=Phi(2)-Phi(y) on the low branch."""
    y, lam = map(sp.sympify, (y, coupling))
    return 1/(2*y)+sp.sqrt(lam)*sp.atanh(sp.sqrt(lam)*y)/2


def branch_point(y, delta):
    """Exact rational interior point, rejecting critical/high/zero branches."""
    y = nonnegative(y, positive=True)
    lam = dimensionless_coupling(delta)
    denominator = 1-2*lam*y*y
    if denominator <= 0:
        raise ValueError("the point is not on the regular low-curvature branch")
    a4 = 4*(1-4*lam)/(y*y*(1-lam*y*y))
    jets = jet_functions(y, lam)
    return {"y": y, "lambda": lam, "branch_denominator": denominator,
            "a_fourth_with_a0_one": a4,
            "normal_Hubble_jets_scaled_by_tau": [-value for value in jets],
            "Ricci_scalar_times_tau_squared": -12*lam*y**4/denominator,
            "kappa_EED_times_tau_squared": 3*(jets[1]-y*y),
            "actual_SEE_point_in_named_thermal_state": True}


def endpoint_data(delta):
    lam = dimensionless_coupling(delta)
    critical = 1/sp.sqrt(2*lam)
    return {"lambda": lam, "critical_y": critical,
            "endpoint_x_exact": clock(2, lam)-clock(critical, lam),
            "endpoint_x_strict_lower": sp.Rational(1, 16),
            "endpoint_x_strict_upper": sp.Rational(1, 4),
            "a_endpoint_fourth": 16*lam*(1-4*lam),
            "endpoint_is_excluded_from_smooth_metric_domain": True,
            "finite_positive_endpoint_scale_factor_for_each_delta": True,
            "Ricci_scalar_diverges_to_minus_infinity": True,
            "fundamental_EFT_validity_at_endpoint_proved": False}


def identities():
    a, hb, kappa, b, tau, y, lam = sp.symbols("a hbar kappa b tau y lambda", positive=True)
    h, hd, q = sp.symbols("H Hdot Q", real=True)
    data = state.zero_type_d_stress(a, h, hd, hbar=hb, thermal_Q=q)
    b_actual = 31*kappa*hb/(1440*sp.pi**2)
    hd_solution = -2*h*h*(1-b*h*h)/(1-2*b*h*h)
    pressure_reduced = (2*hd+3*h*h+kappa*data["pressure"]).subs(q, a**4*(3*h*h/kappa-31*hb*h**4/(480*sp.pi**2)))
    f = velocity(y, lam)
    a4 = 4*(1-4*lam)/(y*y*(1-lam*y*y))
    delta = kappa*hb/(8*sp.pi**2*tau*tau)
    return {"density_SEE_reduction": sp.simplify(
                3*h*h-kappa*data["rho"]-3*(h*h-b_actual*h**4-kappa*q/(3*a**4))),
            "pressure_SEE_reduction": sp.simplify(pressure_reduced-(2*(1-2*b_actual*h*h)*hd+4*h*h*(1-b_actual*h*h))),
            "exact_pressure_solution": sp.factor((2*(1-2*b*h*h)*hd+4*h*h*(1-b*h*h)).subs(hd, hd_solution)),
            "delta_coupling_dictionary": sp.simplify(b_actual/tau**2-sp.Rational(31, 180)*delta),
            "proper_clock_derivative": sp.simplify(sp.diff(clock(y, lam), y)*f+1),
            "density_first_integral": sp.simplify(y*y*(1-lam*y*y)*a4-4*(1-4*lam)),
            "scale_factor_clock": sp.simplify(sp.diff(sp.log(a4), y)*f+4*y),
            "anchored_scale_factor": sp.simplify(a4.subs(y, 2)-1),
            "exact_Ricci_scalar": sp.simplify(6*(-f+2*y*y)+12*lam*y**4/(1-2*lam*y*y)),
            "endpoint_positive_scale_factor": sp.simplify(a4.subs(y, 1/sp.sqrt(2*lam))-16*lam*(1-4*lam)),
            "thermal_amplitude_constraint": sp.simplify(
                (12*(1-4*lam)/(kappa*tau*tau))/3*kappa*tau*tau-4*(1-4*lam))}


def calibration():
    lam = dimensionless_coupling(DELTA_MAX)
    return {"delta_upper": DELTA_MAX, "lambda_upper": lam,
            "initial_y": sp.Integer(2), "initial_a": sp.Integer(1),
            "Q_in_one_over_kappa_tau_squared_units": 12*(1-4*lam),
            "strict_initial_low_branch_margin": 1-8*lam,
            "endpoint_x_strict_bounds": [sp.Rational(1, 16), sp.Rational(1, 4)],
            "endpoint_lower_bound_domain_margin": sp.Rational(1, 64)-lam,
            "a_endpoint_fourth_at_delta_upper": 16*lam*(1-4*lam),
            "uniform_positive_endpoint_a_bound_as_delta_tends_to_zero": False,
            "actual_branch_has_positive_EED": True,
            "this_branch_is_a_new_QEI_only_endpoint_argument": False,
            "fundamental_EFT_endpoint_control": False}
