"""Exhaustive pointwise Bianchi cases, NEC obstruction and exact controls."""

import sympy as sp

from . import background as bg

Y, C, TAU = sp.symbols("y c tau", positive=True)
HDOT, NGNULL, NFNULL = sp.symbols("Hdot_g rho_g_plus_p_g rho_f_plus_p_f", real=True)


def stationary_equations():
    polynomial = bg.interaction_polynomial(Y)
    g_null = -2*bg.MG2*HDOT-NGNULL-(Y-C)*polynomial
    f_null = -2*bg.MF2*HDOT/(C*Y)-NFNULL-(C-Y)*polynomial/(C*Y**3)
    return {"P": polynomial, "dynamic_g_null": g_null, "dynamic_f_null": f_null,
            "dynamic_weighted_null": -2*(bg.MG2+bg.MF2*Y**2)*HDOT-NGNULL-C*Y**3*NFNULL,
            "root_g_null": -2*bg.MG2*HDOT-NGNULL,
            "dynamic_positive_inertia": 2*(bg.MG2+bg.MF2*Y**2),
            "root_positive_inertia": 2*bg.MG2,
            "second_matter_positive_weight": C*Y**3}


def branch_checks():
    data = stationary_equations()
    checks = {"dynamic_weighted_cancellation": sp.factor(data["dynamic_g_null"]+C*Y**3*data["dynamic_f_null"]-data["dynamic_weighted_null"]),
              "root_g_equation_requires_no_dynamic_relation": sp.factor(data["dynamic_g_null"]+(Y-C)*data["P"]-data["root_g_null"]),
              "g_only_dynamic_specialization": sp.factor(data["dynamic_weighted_null"].subs(NFNULL, 0)+2*(bg.MG2+bg.MF2*Y**2)*HDOT+NGNULL)}
    _, beta1, beta2, beta3, _ = bg.BETAS
    checks["simple_root_at_y_one"] = data["P"].subs({beta1: -2, beta2: 1, beta3: 0, Y: 1})
    checks["double_root_at_y_one"] = data["P"].subs({beta1: 1, beta2: -1, beta3: 1, Y: 1})
    checks["double_root_factorization"] = sp.factor(data["P"].subs({beta1: 1, beta2: -1, beta3: 1})-bg.M4*(Y-1)**2)
    checks["identically_zero_polynomial_with_free_cosmological_terms"] = data["P"].subs({beta1: 0, beta2: 0, beta3: 0})
    # The pointwise root formula needs no inference about which neighboring
    # branch is followed, root multiplicity, or isolated branch switching.
    scale = sp.Symbol("a_stationary", positive=True)
    eq = bg.equations()
    point = {bg.a: scale, bg.b: Y*scale, bg.Ng: 1, bg.Nf: C,
             bg.KG: NGNULL, bg.KF: C**2*NFNULL,
             sp.diff(bg.a, bg.t): 0, sp.diff(bg.b, bg.t): 0,
             sp.diff(bg.Ng, bg.t): 0, sp.diff(bg.Nf, bg.t): 0,
             sp.diff(bg.a, bg.t, 2): scale*HDOT,
             sp.diff(bg.b, bg.t, 2): C*scale*HDOT}
    checks["dynamic_g_null_from_full_Euler_equations"] = sp.factor(
        (eq["EL_a"]-eq["EL_Ng"]).subs(point, simultaneous=True)+data["dynamic_g_null"])
    checks["dynamic_f_null_from_full_Euler_equations"] = sp.factor(
        (eq["EL_b"]-eq["EL_Nf"]).subs(point, simultaneous=True)+data["dynamic_f_null"])
    point[sp.diff(bg.b, bg.t, 2)] = sp.Symbol("independent_b_ddot", real=True)
    checks["root_g_null_without_f_acceleration_assumption"] = sp.factor(
        (eq["EL_a"]-eq["EL_Ng"]).subs(point, simultaneous=True).subs(beta1, -2*beta2*Y-beta3*Y**2)+data["root_g_null"])
    # Differentiate Hf=Hg/y only on a neighborhood where P != 0.
    h, y, c = (sp.Function(name)(bg.t) for name in ("H", "y", "c"))
    checks["dynamic_derivative_at_stationarity"] = sp.factor(
        (sp.diff(h/y, bg.t)/c).subs({h: 0, sp.diff(h, bg.t): HDOT}, simultaneous=True)-HDOT/(c*y))
    checks["complete_dynamic_derivative_identity"] = sp.factor(
        sp.diff(h/y, bg.t)/c-sp.diff(h, bg.t)/(c*y)+h*sp.diff(y, bg.t)/(c*y**2))
    return checks


def window_checks():
    u, w = sp.symbols("u w", real=True)
    hubble = 4*u/(TAU*(1+u**2))
    derivative = sp.diff(hubble, u)/TAU
    return {"CD_bounce_derivative": sp.factor(derivative.subs(u, 0)-4/TAU**2),
            "CD_left_window_H": sp.factor(hubble.subs(u, -sp.Rational(1, 2))+sp.Rational(8, 5)/TAU),
            "CD_right_window_H": sp.factor(hubble.subs(u, sp.Rational(1, 2))-sp.Rational(8, 5)/TAU),
            "CD_window_derivative_positive_margin": sp.factor(
                4*(1-w)/(1+w)**2-sp.Rational(48, 25)-(1-4*w)*(52+12*w)/(25*(1+w)**2))}


def control_points():
    eq = bg.equations()
    t = bg.t
    base = {bg.MG2: 1, bg.MF2: 1, bg.M4: 1, bg.a: 1, bg.b: 1, bg.Ng: 1,
            sp.diff(bg.a, t): 0, sp.diff(bg.b, t): 0,
            sp.diff(bg.Ng, t): 0, sp.diff(bg.Nf, t): 0,
            sp.diff(bg.a, t, 2): sp.Rational(1, 10),
            bg.BETAS[2]: 0, bg.BETAS[3]: 0}
    positive_spring = {bg.BETAS[0]: -3, bg.BETAS[1]: 1, bg.BETAS[4]: -1}
    g_bad_nec = {**base, **positive_spring, bg.Nf: sp.Rational(4, 5),
                 sp.diff(bg.b, t, 2): sp.Rational(2, 25),
                 bg.KG: -sp.Rational(2, 5), bg.VG: sp.Rational(1, 5), bg.KF: 0, bg.VF: 0}
    f_bad_nec = {**base, bg.BETAS[0]: 3, bg.BETAS[1]: -1, bg.BETAS[4]: 1,
                 bg.Nf: sp.Rational(4, 5), sp.diff(bg.b, t, 2): sp.Rational(2, 25),
                 bg.KG: 0, bg.VG: 0, bg.KF: -sp.Rational(8, 25), bg.VF: sp.Rational(1, 4)}
    omitted_f = {**base, **positive_spring, bg.Nf: 2, sp.diff(bg.b, t, 2): sp.Rational(1, 5),
                 bg.KG: sp.Rational(4, 5), bg.VG: -sp.Rational(2, 5), bg.KF: 0, bg.VF: 0}
    ratio_two = {**base, bg.BETAS[0]: -6, bg.BETAS[1]: 1, bg.BETAS[4]: -sp.Rational(1, 8),
                 bg.b: 2, bg.Nf: 2, sp.diff(bg.a, t, 2): 0, sp.diff(bg.b, t, 2): 0,
                 bg.KG: 0, bg.VG: 0, bg.KF: 0, bg.VF: 0}
    return {name: {key: sp.factor(eq[key].subs(point, simultaneous=True)) for key in ("EL_Ng", "EL_Nf", "EL_a", "EL_b")}
            for name, point in (("forbidden_g_NEC", g_bad_nec), ("forbidden_f_NEC", f_bad_nec),
                                ("omitted_f_acceleration", omitted_f), ("regular_ratio_two", ratio_two))}


def control_checks():
    data = stationary_equations()
    formal_negative = {bg.MG2: 1, bg.MF2: 1, Y: 1, C: -1, HDOT: 1, NGNULL: 0, NFNULL: 4}
    # Formal continuation only: sqrt|f| on a negative-lapse coordinate chart
    # is not obtained by analytically continuing the positive-root action.
    checks = {"negative_c_formal_weighted_equation": data["dynamic_weighted_null"].subs(formal_negative)}
    for name, values in control_points().items():
        for key, value in values.items():
            expected = sp.Rational(3, 5) if (name, key) == ("omitted_f_acceleration", "EL_b") else 0
            checks[name+"_"+key] = sp.factor(value-expected)
    return checks


def controls():
    data = stationary_equations()
    return {"dropped_f_acceleration_full_equation_error": control_points()["omitted_f_acceleration"]["EL_b"],
            "forbidden_g_null_stress": -sp.Rational(2, 5),
            "forbidden_f_null_stress": -sp.Rational(1, 2),
            "zero_lapse_has_zero_metric_determinant": (-C**2*Y**6).subs(C, 0),
            "formal_negative_c_loses_positive_second_matter_weight": data["second_matter_positive_weight"].subs({Y: 1, C: -1}),
            "formal_negative_c_can_make_weighted_NEC_sum_negative": (NGNULL+C*Y**3*NFNULL).subs({NGNULL: 0, NFNULL: 4, Y: 1, C: -1}),
            "algebraic_branch_dynamic_derivative_inference_error": (HDOT/C**2-HDOT/(C*Y)).subs({HDOT: 1, C: 2, Y: 1}),
            "zero_polynomial_is_a_covered_case_not_a_divisor": data["P"].subs(dict.fromkeys(bg.BETAS[1:4], 0)),
            "regular_stationary_ratio_two_is_not_y_one": sp.Integer(2)-1}
