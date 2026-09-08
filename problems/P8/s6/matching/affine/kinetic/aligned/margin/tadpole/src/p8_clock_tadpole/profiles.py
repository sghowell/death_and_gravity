"""A new finite lower-scalar action with profiles fixed before variation."""
from functools import cache

import sympy as sp
from p8_vector_regularity import estimates
from p8_vector_state import wkb

from . import window

u = wkb.u
x = sp.Symbol("physical_clock_norm", real=True)


@cache
def covariant():
    rho, pressure = sp.Function("selected_energy")(u), sp.Function("selected_pressure")(u)
    P = -pressure+(rho+pressure)*(x+1)/2
    H = sp.Function("original_Hubble")(u)
    Px = sp.diff(P, x)
    energy = (2*x*Px-P).subs(x, -1)
    spatial = P.subs(x, -1)
    clock = sp.diff(P.subs(x, -1), u)+2*sp.diff(Px.subs(x, -1), u)+6*H*Px.subs(x, -1)
    return {"rho": rho, "pressure": pressure, "H": H, "P_on_flat_clock_tube": P,
            "global_P": window.cutoff(x)*P,
            "energy": energy, "spatial_pressure": spatial, "clock_source": sp.simplify(clock),
            "checks": {"selected_energy_cancelled": sp.expand(energy+rho),
                       "selected_pressure_cancelled": sp.expand(spatial+pressure),
                       "selected_Ward_clock_source_cancelled": sp.expand(clock-sp.diff(rho, u)-3*H*(rho+pressure)),
                       "no_added_second_clock_norm_derivative_on_flat_tube": sp.diff(P, x, 2)}}


@cache
def point_chart():
    N, h = sp.symbols("physical_lapse h", positive=True)
    rho, pressure = sp.symbols("fixed_energy_profile fixed_pressure_profile", real=True)
    v = sp.Symbol("hat_log_scale", real=True)
    omega = -sp.log((h-1+N**-2)/h)/4
    inner = -pressure+(rho+pressure)*(1-N**-2)/2
    density = N*sp.exp(3*omega+3*v)*inner
    point = {N: 1, v: 0}
    coeff_rho = (3/h-1)/2
    coeff_pressure = -sp.Rational(1, 2)+sp.Rational(9, 4)/h-sp.Rational(21, 8)/h**2
    lapse_square = coeff_rho*rho+coeff_pressure*pressure
    exact = {"value": -pressure, "N_first": rho-3*pressure/(2*h), "v_first": -3*pressure,
             "N_square_coefficient": lapse_square, "Nv_coefficient": 3*rho-9*pressure/(2*h),
             "v_square_coefficient": -9*pressure/2}
    checks = {"profile_point_chart_value": sp.simplify(density.subs(point)-exact["value"]),
              "profile_point_chart_lapse_force": sp.simplify(sp.diff(density, N).subs(point)-exact["N_first"]),
              "profile_point_chart_scale_force": sp.simplify(sp.diff(density, v).subs(point)-exact["v_first"]),
              "profile_point_chart_lapse_square": sp.simplify(sp.diff(density, N, 2).subs(point)/2-lapse_square),
              "profile_point_chart_mixed_square": sp.simplify(sp.diff(density, N, v).subs(point)-exact["Nv_coefficient"]),
              "profile_point_chart_scale_square": sp.simplify(sp.diff(density, v, 2).subs(point)/2-exact["v_square_coefficient"])}
    y = sp.Symbol("inverse_h", nonnegative=True)
    Q = 21*y**2-18*y+4
    checks.update({"pressure_coefficient_strict_sign_square": sp.expand(Q-21*(y-sp.Rational(3, 7))**2-sp.Rational(1, 7)),
                   "pressure_coefficient_endpoint_envelope": sp.expand(7-Q-(1-y)*(21*y+3))})
    return {"h": h, "rho": rho, "pressure": pressure, "coefficients": exact,
            "lapse_square_energy_coefficient": coeff_rho,
            "lapse_square_pressure_coefficient": coeff_pressure, "checks": checks}


def profile_bounds(planck_time_product, mass_time_product):
    prior = estimates.physical_bounds(planck_time_product, mass_time_product)
    stress = prior["normalized_derivative_bounds"]
    eta = {order: max(stress[name][order]["total"] for name in ("energy", "pressure")) for order in range(6)}
    return {"M_tau": prior["M_tau"], "m0_tau": prior["m0_tau"],
            "fixed_profile_derivative_bounds_on_clock_tube": {
                order: {"P": sp.Rational(11, 10)*upper, "P_x": upper, "P_xx": sp.Integer(0)}
                for order, upper in eta.items()},
            "all_clock_norm_mixed_derivative_bounds_total_order_at_most_five": {
                str(j)+","+str(k): (sp.Rational(5, 4)*window.derivative_bound(k)
                                     +(k*window.derivative_bound(k-1) if k else 0))*eta[j]
                for j in range(6) for k in range(6-j)},
            "lapse_square_coefficient_absolute_upper": sp.Rational(15, 8)*eta[0],
            "mixed_lapse_scale_coefficient_absolute_upper": sp.Rational(15, 2)*eta[0],
            "scale_square_coefficient_absolute_upper": sp.Rational(9, 2)*eta[0],
            "profile_held_fixed_under_subsequent_state_and_metric_variations": True,
            "full_quantum_or_UV_completion_claim": False}
