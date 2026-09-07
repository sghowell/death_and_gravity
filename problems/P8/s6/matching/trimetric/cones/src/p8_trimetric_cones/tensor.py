"""Literal TT auxiliary completion and physical-clock common tensor action.

One real unit TT polarization has tr(E²)=1. S_probe=int sqrt|h| j gamma_u;
j includes any source-strength prefactor. The external source is distinct
from the homogeneous canonical scalar and is treated to quadratic response.
"""

from functools import cache

import sympy as sp


@cache
def derive():
    pg, pf = sp.symbols("P_g P_f", positive=True)
    ge, gf, gu, j = sp.symbols("gamma_e gamma_v gamma_u j", real=True)
    weights = pg+pf
    average = (pg*ge+pf*gf)/weights
    stationary = average+2*j/weights
    algebraic = -(pg*(ge-gu)**2+pf*(gf-gu)**2)/4+j*gu
    reduced = -pg*pf*(ge-gf)**2/(4*weights)+j*average+j**2/weights
    g, f, n, a, big_n, big_a, q, k = sp.symbols("G F n a N A q k", positive=True)
    gedot, gfdot = sp.symbols("gamma_e_dot gamma_v_dot", real=True)
    # This density uses coordinate t; D_T=(N*n)^-1 partial_t, a_h=A*a.
    quadratic = a**3/(8*n)*(g*gedot**2+f*gfdot**2-n**2*k**2/a**2*(g*ge**2+f*gf**2)
                            -2*q*n**2*big_n*big_a**2*((ge-gu)**2+(gf-gu)**2))
    common, relative, cdot, rdot = sp.symbols("gamma delta gamma_dot delta_dot", real=True)
    symmetric = quadratic.subs({f: g, ge: common+relative, gf: common-relative,
                                 gu: common, gedot: cdot+rdot, gfdot: cdot-rdot}, simultaneous=True)
    gt, ft = 2*g*big_n/big_a**3, 2*g/(big_n*big_a)
    return {"P_g": pg, "P_f": pf, "gamma_e": ge, "gamma_v": gf, "gamma_u": gu, "j": j,
            "weight_sum": weights, "source_free_auxiliary_map": average,
            "sourced_auxiliary_map": stationary, "algebraic_density_over_h_volume": algebraic,
            "reduced_density_over_h_volume": reduced, "source_contact": j**2/weights,
            "G": g, "F": f, "n": n, "a": a, "N": big_n, "A": big_a, "q": q, "k": k,
            "gamma_e_dot": gedot, "gamma_v_dot": gfdot,
            "full_symmetric_background_TT_density": quadratic,
            "gamma": common, "delta": relative, "gamma_dot": cdot, "delta_dot": rdot,
            "exchange_even_odd_density": sp.expand(symmetric),
            "G_T_h": gt, "F_T_h": ft, "common_speed_squared": big_a**2/big_n**2,
            "relative_algebraic_mass_h_squared": 2*q*big_a**2/(g*big_n),
            "symmetric_source_contact": big_a*j**2/(2*q),
            "symmetric_sourced_auxiliary_map": common+big_a*j/q}


@cache
def checks():
    d = derive()
    pg, pf, ge, gf, gu, j = (d[key] for key in ("P_g", "P_f", "gamma_e", "gamma_v", "gamma_u", "j"))
    lag = d["algebraic_density_over_h_volume"]
    source = d["sourced_auxiliary_map"]
    # Literal diagonal exponential expansion. E=diag(1,-1,0) has norm²=2;
    # divide its raw quadratic coefficient by 2 to obtain a unit polarization.
    lam = sp.Symbol("lambda", real=True)
    trace_e = 1+2*sp.cosh(lam*(ge-gu)/2)
    trace_f = 1+2*sp.cosh(lam*(gf-gu)/2)
    literal = -2*(pg*trace_e+pf*trace_f)
    literal_q2 = sp.diff(literal, lam, 2).subs(lam, 0)/4
    g, n, a, big_n, big_a, q, k = (d[key] for key in ("G", "n", "a", "N", "A", "q", "k"))
    gamma, delta, cdot, rdot = (d[key] for key in ("gamma", "delta", "gamma_dot", "delta_dot"))
    even_odd_expected = a**3/(8*n)*(2*g*(cdot**2+rdot**2)-2*g*n**2*k**2/a**2*(gamma**2+delta**2)
                                        -4*q*n**2*big_n*big_a**2*delta**2)
    dhgamma, dhdelta = sp.symbols("D_h_gamma D_h_delta", real=True)
    physical = d["exchange_even_odd_density"].subs({cdot: big_n*n*dhgamma, rdot: big_n*n*dhdelta})/(big_n*n)
    physical_expected = (big_a*a)**3/8*(d["G_T_h"]*(dhgamma**2+dhdelta**2)
                         -d["F_T_h"]*k**2/(big_a*a)**2*(gamma**2+delta**2)
                         -d["G_T_h"]*d["relative_algebraic_mass_h_squared"]*delta**2)
    return {
        "literal_diagonal_exponential_unit_TT": sp.expand(literal_q2+pg*(ge-gu)**2/4+pf*(gf-gu)**2/4),
        "sourced_auxiliary_stationarity": sp.cancel(sp.diff(lag, gu).subs(gu, source)),
        "full_auxiliary_completed_square": sp.cancel(lag-d["reduced_density_over_h_volume"]
                                           +(pg+pf)*(gu-source)**2/4),
        "source_response_includes_contact": sp.cancel(sp.diff(d["reduced_density_over_h_volume"], j)-source),
        "symmetric_auxiliary_weight": sp.cancel(d["source_free_auxiliary_map"].subs({pg: q/big_a, pf: q/big_a})-(ge+gf)/2),
        "symmetric_contact_normalization": sp.cancel(d["source_contact"].subs({pg: q/big_a, pf: q/big_a})-d["symmetric_source_contact"]),
        "symmetric_even_odd_action": sp.expand(d["exchange_even_odd_density"]-even_odd_expected),
        "exact_physical_clock_and_volume": sp.cancel(physical-physical_expected),
        "common_tensor_physical_cone": sp.cancel(d["F_T_h"]/d["G_T_h"]-d["common_speed_squared"]),
        "vacuum_relative_mass_dictionary": sp.cancel(d["relative_algebraic_mass_h_squared"].subs({big_a: 1, big_n: 1})-2*q/g),
    }


def controls():
    d = derive()
    point = {d["A"]: 2, d["N"]: sp.Rational(2, 7), d["q"]: 1, d["j"]: 1}
    relative = {d["P_g"]: 1, d["P_f"]: 1, d["gamma_e"]: d["delta"], d["gamma_v"]: -d["delta"]}
    return {"dropping_actual_clock_at_A2_N2over7_misses_squared_speed": d["common_speed_squared"].subs(point)-1,
            "dropping_probe_contact_at_A2_q1_j1": d["symmetric_source_contact"].subs(point),
            "source_free_relative_projection": d["source_free_auxiliary_map"].subs(relative),
            "unequal_Einstein_even_odd_kinetic_cross_coefficient": (d["G"]-d["F"])/4,
            "auxiliary_tensor_denominator_zero_outside_domain": d["weight_sum"].subs({d["P_g"]: -1, d["P_f"]: 1})}
