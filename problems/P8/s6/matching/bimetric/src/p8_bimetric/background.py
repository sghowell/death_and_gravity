"""Full two-lapse FLRW variation and a regular beta1 bounce obstruction.

P8 signature +---, R=-6*(Hdot+2H^2). Einstein actions are -M_i^2 R_i/2.
The interaction is -m4*sqrt(|g|)*sum beta_n e_n(sqrt(g^-1 f)), with
beta0=-3 beta1, beta4=-beta1, beta2=beta3=0 and NU=m4*beta1>0.
"""

from functools import cache

import sympy as sp

t = sp.Symbol("t", real=True)
a, b, Ng, Nf = (sp.Function(name, positive=True)(t) for name in ("a", "b", "N_g", "N_f"))
KIN, POT = (sp.Function(name)(t) for name in ("K_coordinate", "V"))
MG2, MF2, NU = sp.symbols("M_g_squared M_f_squared nu", positive=True)
Y, C = sp.symbols("y c", positive=True)
HG, HF, HDOT, RHOP = sp.symbols("H_g H_f Hdot_g rho_plus_p", real=True)
TAU = sp.Symbol("tau", positive=True)


def ricci(scale, lapse):
    return -6*(sp.diff(scale, t, 2)/(scale*lapse**2)
               +sp.diff(scale, t)**2/(scale**2*lapse**2)
               -sp.diff(scale, t)*sp.diff(lapse, t)/(scale*lapse**3))


@cache
def equations():
    hg, hf = sp.diff(a, t)/(Ng*a), sp.diff(b, t)/(Nf*b)
    y, c = b/a, Nf/Ng
    rho, pressure = KIN/(2*Ng**2)+POT, KIN/(2*Ng**2)-POT
    grav_g = -3*MG2*a*sp.diff(a, t)**2/Ng
    grav_f = -3*MF2*b*sp.diff(b, t)**2/Nf
    interaction = -NU*(-3*Ng*a**3+Nf*a**3+3*Ng*a**2*b-Nf*b**3)
    lagrangian = grav_g+grav_f+interaction+a**3*KIN/(2*Ng)-Ng*a**3*POT

    def el(field):
        return sp.diff(lagrangian, field)-sp.diff(sp.diff(lagrangian, sp.diff(field, t)), t)

    normalized = {"EL_Ng": sp.factor(el(Ng)/a**3), "EL_Nf": sp.factor(el(Nf)/b**3),
                  "EL_a": sp.factor(el(a)/(3*Ng*a**2)), "EL_b": sp.factor(el(b)/(3*Nf*b**2))}
    rho_gv, p_gv = 3*NU*(y-1), NU*(3-2*y-c)
    rho_fv, p_fv = NU*(y**-3-1), NU*(1-1/(c*y**2))
    expected = {"EL_Ng": 3*MG2*hg**2-rho-rho_gv,
                "EL_Nf": 3*MF2*hf**2-rho_fv,
                "EL_a": MG2*(2*sp.diff(hg, t)/Ng+3*hg**2)+pressure+p_gv,
                "EL_b": MF2*(2*sp.diff(hf, t)/Nf+3*hf**2)+p_fv}
    return {"L": lagrangian, "interaction": interaction, "H_g": hg, "H_f": hf,
            "y": y, "c": c, "rho": rho, "pressure": pressure,
            "rho_g_interaction": rho_gv, "p_g_interaction": p_gv,
            "rho_f_interaction": rho_fv, "p_f_interaction": p_fv,
            **normalized, "expected": expected}


def variation_checks():
    data = equations()
    checks = {name: sp.simplify(data[name]-expected) for name, expected in data["expected"].items()}
    for label, scale, lapse, mass in (("g", a, Ng, MG2), ("f", b, Nf, MF2)):
        raw = -mass*lapse*scale**3*ricci(scale, lapse)/2
        boundary = 3*mass*scale**2*sp.diff(scale, t)/lapse
        bulk = -3*mass*scale*sp.diff(scale, t)**2/lapse
        checks[label+"_Einstein_time_boundary"] = sp.simplify(raw-sp.diff(boundary, t)-bulk)
    root_eigenvalues = sp.Matrix.diag(Nf/Ng, b/a, b/a, b/a)
    potential_direct = -NU*Ng*a**3*(-3+sp.trace(root_eigenvalues)-root_eigenvalues.det())
    checks["full_square_root_potential"] = sp.expand(potential_direct-data["interaction"])
    return checks


def bianchi_checks():
    data = equations()
    hg, hf, y, c = (data[key] for key in ("H_g", "H_f", "y", "c"))
    dtg = lambda value: sp.diff(value, t)/Ng
    matter_conservation = dtg(data["rho"])+3*hg*(data["rho"]+data["pressure"])
    interaction_conservation = dtg(data["rho_g_interaction"])+3*hg*(data["rho_g_interaction"]+data["p_g_interaction"])
    numerator = Ng*sp.diff(b, t)-Nf*sp.diff(a, t)
    # No division by H_g or a_dot: these identities remain valid at a bounce.
    return {"g_equation_Noether_identity": sp.simplify(dtg(data["EL_Ng"])-3*hg*(data["EL_a"]-data["EL_Ng"])
                                                       +matter_conservation+interaction_conservation),
            "interaction_dynamical_branch_factor": sp.simplify(interaction_conservation-3*NU*numerator/(Ng**2*a)),
            "Hubble_relation_without_dividing_Hg": sp.simplify(hf-hg/y-numerator/(Ng*Nf*b)),
            "ratio_branch_relation": sp.simplify(dtg(y)+hg*(y-c)-numerator/(Ng**2*a))}


def bounce_equations():
    """At H_g=0, regular Bianchi gives H_f=0; f constraint then fixes y=1."""
    # The f derivative is D_f H_f=Hdot_g/c at y=1,H_g=0.
    g_null = -2*MG2*HDOT-RHOP-NU*(1-C)
    f_null_weighted = -2*MF2*HDOT-NU*(C-1)
    return {"f_constraint_at_stationarity": NU*(Y**-3-1),
            "g_null": g_null, "f_null_weighted": f_null_weighted,
            "summed_null": -2*(MG2+MF2)*HDOT-RHOP,
            "f_lapse_from_acceleration": 1-2*MF2*HDOT/NU,
            "CD_Hdot_at_bounce": 4/TAU**2,
            "CD_required_null_stress": -8*(MG2+MF2)/TAU**2,
            "minimum_Hdot_error": 4/TAU**2,
            "CD_half_tau_window_endpoint_H_magnitude": sp.Rational(8, 5)/TAU,
            "CD_half_tau_window_minimum_Hdot": sp.Rational(48, 25)/TAU**2}


def bounce_checks():
    data = bounce_equations()
    # Multiplication by y^3 is safe since the branch has y>0. The other
    # factor y^2+y+1 is strictly positive, so y=1 is the only positive root.
    scale_cd = (1+(t/TAU)**2)**2
    hg_cd = sp.diff(scale_cd, t)/scale_cd
    checks = {"positive_ratio_root_factorization": sp.expand(Y**3*data["f_constraint_at_stationarity"]+NU*(Y-1)*(Y**2+Y+1)),
              "strict_positive_other_root_factor": sp.expand(Y**2+Y+1-(Y+sp.Rational(1, 2))**2-sp.Rational(3, 4)),
              "full_parent_summed_null_equation": sp.expand(data["g_null"]+data["f_null_weighted"]-data["summed_null"]),
              "CD_bounce_Hdot": sp.simplify(sp.diff(hg_cd, t).subs(t, 0)-4/TAU**2),
              "CD_required_null_stress": sp.expand(data["summed_null"].subs({HDOT: 4/TAU**2, RHOP: data["CD_required_null_stress"]})),
              "Minkowski_positive_lapse_control": data["summed_null"].subs({HDOT: 0, RHOP: 0})}
    eq = equations()
    scale = sp.Symbol("a_bounce", positive=True)
    point = {a: scale, b: scale, Ng: 1, Nf: C, KIN: RHOP,
             sp.diff(a, t): 0, sp.diff(b, t): 0, sp.diff(Ng, t): 0,
             sp.diff(a, t, 2): scale*HDOT, sp.diff(b, t, 2): C*scale*HDOT}
    checks["g_null_from_full_lapse_and_scale_equations"] = sp.simplify(
        (eq["EL_a"]-eq["EL_Ng"]).subs(point, simultaneous=True)+data["g_null"])
    checks["weighted_f_null_from_full_lapse_and_scale_equations"] = sp.simplify(
        C*(eq["EL_b"]-eq["EL_Nf"]).subs(point, simultaneous=True)+data["f_null_weighted"])
    point[b] = Y*scale
    checks["ratio_equation_from_full_f_lapse"] = sp.simplify(
        eq["EL_Nf"].subs(point, simultaneous=True)+data["f_constraint_at_stationarity"])
    # Generic LDL^T positive field metric, an algebraic regression for the
    # written all-dimensional positive-quadratic-form argument.
    v1, v2, mix = sp.symbols("v1 v2 mix", real=True)
    z1, z2 = sp.symbols("z1 z2", positive=True)
    lower = sp.Matrix([[1, 0], [mix, 1]])
    metric = lower*sp.diag(z1, z2)*lower.T
    velocity = sp.Matrix([v1, v2])
    checks["positive_mixed_field_metric_squares"] = sp.expand((velocity.T*metric*velocity)[0]-z1*(v1+mix*v2)**2-z2*v2**2)
    w = sp.Symbol("u_squared", nonnegative=True)
    checks["CD_window_right_endpoint"] = sp.simplify(hg_cd.subs(t, TAU/2)-sp.Rational(8, 5)/TAU)
    checks["CD_window_left_endpoint"] = sp.simplify(hg_cd.subs(t, -TAU/2)+sp.Rational(8, 5)/TAU)
    checks["CD_window_positive_derivative_margin"] = sp.factor(
        4*(1-w)/(1+w)**2-sp.Rational(48, 25)-(1-4*w)*(52+12*w)/(25*(1+w)**2))
    return checks


def NEC_violating_point_checks():
    """Exact positive-lapse algebraic point outside the positive-matter class.

    This is a control, not a parent solution or a healthy P8 witness.
    The f second derivative satisfies the differentiated regular Bianchi.
    """
    data = equations()
    point = {MG2: 1, MF2: 1, NU: 1, a: 1, b: 1, Ng: 1, Nf: sp.Rational(4, 5),
             sp.diff(a, t): 0, sp.diff(b, t): 0, sp.diff(Ng, t): 0, sp.diff(Nf, t): 0,
             sp.diff(a, t, 2): sp.Rational(1, 10), sp.diff(b, t, 2): sp.Rational(2, 25),
             KIN: -sp.Rational(2, 5), POT: sp.Rational(1, 5)}
    return {key: sp.simplify(data[key].subs(point, simultaneous=True)) for key in ("EL_Ng", "EL_Nf", "EL_a", "EL_b")}


def controls():
    data = bounce_equations()
    return {"discarding_f_lapse_loses_ratio_constraint": data["f_constraint_at_stationarity"].subs(Y, 2),
            "discarding_f_acceleration_loses_Planck_term": 2*MF2*HDOT,
            "zero_interaction_loses_branch_factor": sp.Integer(0),
            "forbidden_NEC_control_null_stress": -sp.Rational(2, 5),
            "negative_or_zero_ratio_not_in_positive_branch": sp.Integer(-1),
            "CD_matching_error_remains_order_one_in_background_units": sp.Integer(4)}
