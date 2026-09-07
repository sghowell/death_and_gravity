"""Exact isotropic own-f branch on fixed physical g=eta and theta=T.

This is not an on-shell full-parent/CD background.  The algebraic
desingularization has a regular delta=0 extension, whereas the original
coefficient action does not exist at (u,delta)=(0,0).
"""

from functools import cache

import sympy as sp

M, TAU = sp.symbols("M tau", positive=True)
TIME, U, V = sp.symbols("t u v", real=True)
C = sp.Symbol("c", positive=True)
DELTA = sp.Symbol("delta", positive=True)
ZETA = sp.Symbol("zeta", real=True)
B = sp.Symbol("b_value", positive=True)


@cache
def own_equations():
    """Literal lapse/scale variations and their unfactored consistency.

    E is variation with respect to the common spatial scale b (the sum of
    the three isotropic spatial variations).  A regular solution means
    b>0, N>0, b at least C2 and N at least C1 in physical-g proper time.
    """
    b, lapse = sp.Function("b")(TIME), sp.Function("N")(TIME)
    p, w = sp.Function("beta1")(TIME), sp.Function("beta4")(TIME)
    bd = sp.diff(b, TIME)
    kinetic = -3*M**2*b*bd**2/lapse
    e1, e4 = lapse+3*b, lapse*b**3
    potential = -2*p*e1-2*w*e4
    density = kinetic+potential
    constraint = sp.diff(density, lapse)
    spatial = sp.diff(density, b)-sp.diff(sp.diff(density, bd), TIME)
    consistency = sp.diff(p, TIME)+sp.diff(w, TIME)*b**3-3*p*bd/lapse
    curvature = -6*(sp.diff(b, TIME, 2)/(lapse**2*b)
                    +bd**2/(lapse**2*b**2)-bd*sp.diff(lapse, TIME)/(lapse**3*b))
    boundary = 3*M**2*b**2*bd/lapse
    return {"b": b, "N": lapse, "beta1": p, "beta4": w,
            "e1": e1, "e4": e4, "density": density,
            "kinetic": kinetic, "potential": potential,
            "C": constraint, "E": spatial, "consistency": consistency,
            "unfactored_identity": sp.expand(sp.diff(constraint, TIME)-bd*spatial/lapse+2*consistency),
            "R_B_f": curvature, "covariant_EH": -M**2*lapse*b**3*curvature/2,
            "EH_boundary": boundary,
            "scope": "own-f equations only; g, theta, free chi fixed; beta0 has no own-f variation"}


@cache
def profile():
    """Frozen VARIABLE beta_i=(M²/tau²) B_i(v), v=u², theta=T=tau*u."""
    d = 1+V
    y = 2/d**4
    hp = 4*(1-V)/d**2
    D = C-y
    Q = 32*(1-V)/(C*d**14)
    J = C*d**4*(1-7*V)+12*V
    R = d**8*J/(8*C*(1-V))
    B1 = Q/D
    B4 = 24*V/(C**2*d**2)-B1/y**3
    L = 8/d+sp.diff(J, V)/J+1/(1-V)
    H = D*(sp.diff(J, V)/J-6/d)-sp.diff(D, V)
    return {"v": V, "u": U, "c": C, "delta": DELTA, "d": d,
            "y": y, "h_u": hp, "D": D, "Q": Q, "J": J,
            "R": R, "L": L, "H": H, "B1": B1, "B4": B4,
            "beta1": M**2*B1/TAU**2, "beta4": M**2*B4/TAU**2,
            "algebraic_r_cubed": 1/R,
            "actual_domain": "c=2+delta; 0<delta<=1/100; |u|<=1/100; v=u²",
            "extension_domain": "v in [0,1/10000], delta in [0,1/100]; no action at their common zero"}


@cache
def fixed_c():
    """Un-divided consistency implies an algebraic equation, even for b(u).

    F is defined by b_u/N=u*F(v,b), with no assumption that b or N is
    even.  The nonzero b-Jacobian at u=0 then proves evenness by uniqueness.
    """
    p = profile()
    b1, b4 = p["B1"], p["B4"]
    F = 2*(sp.diff(b1, V)+sp.diff(b4, V)*B**3)/(3*b1)
    equation = 3*B*V*F**2-2*(b1+b4*B**3)

    def at(expression):
        return sp.factor(expression.subs({V: 0, B: 2}))

    jacobian = at(sp.diff(equation, B))
    bv = sp.factor(-at(sp.diff(equation, V))/jacobian)
    bvv = sp.factor(-(at(sp.diff(equation, V, 2))
                     +2*at(sp.diff(equation, V, B))*bv
                     +at(sp.diff(equation, B, 2))*bv**2)/jacobian)
    F0 = at(F)
    N0 = sp.factor(2*bv/F0)
    Nv = sp.factor(2*bvv/F0-2*bv*(at(sp.diff(F, V))+at(sp.diff(F, B))*bv)/F0**2)
    return {"F": F, "equation": equation, "b_jacobian_center": jacobian,
            "F_center": F0, "b_v_center": bv, "b_vv_center": bvv,
            "N_center": N0, "N_v_center": Nv,
            "domain": "fixed c=1 or 2<c<=4; only the positive-delta box has a uniform claim"}


@cache
def joint():
    """Desingularized joint fixed-point map; partials keep zeta independent.

    The positive real cube root is used only on 1-v*D*zeta>0 and R>0.
    In N, b_v is the partial derivative; zeta_v includes the implicit
    derivative T_v/(1-T_zeta).  The separate Fraction certificate bounds
    these expressions on the entire stated box, not a sampled graph.
    """
    p = profile()
    D, Q, R, L, H = (p[name] for name in ("D", "Q", "R", "L", "H"))
    W = V*D*ZETA
    b_cubed = (1-W)/R
    b = b_cubed**sp.Rational(1, 3)
    F = 2*(V*ZETA*H-L)/3
    Fv, Fz = sp.diff(F, V), sp.diff(F, ZETA)
    brv = (-(D+V*sp.diff(D, V))*ZETA/(1-W)-L)/3
    brz = -V*D/(3*(1-W))
    bv, bz = b*brv, b*brz
    T = 3*b*F**2/(2*Q)
    Tv = T*(brv+2*Fv/F+1/(1-V)+14/(1+V))
    Tz = T*(brz+2*Fz/F)
    zv = Tv/(1-Tz)
    lapse = 2*(bv+bz*zv)/F
    return {**{name: p[name] for name in ("D", "Q", "J", "R", "L", "H")},
            "W": W, "b_cubed": b_cubed, "b": b, "F": F,
            "F_v": Fv, "F_zeta": Fz, "b_v_over_b": brv,
            "b_zeta_over_b": brz, "b_v": bv, "b_zeta": bz,
            "T": T, "T_v": Tv, "T_zeta": Tz, "zeta_v": zv, "N": lapse,
            "G": 3*b*F**2-2*Q*ZETA,
            "fixed_point_residual": T-ZETA,
            "box": {V: (0, sp.Rational(1, 10000)), DELTA: (0, sp.Rational(1, 100)),
                    ZETA: (11, 13)},
            "parameter_substitution": {C: 2+DELTA}}


@cache
def center():
    """Center equations with no division by b''; exact profile jets/units."""
    r = sp.Symbol("r0", positive=True)
    # The derivative of r is not assumed positive; r0 alone is positive.
    rdd = sp.Symbol("r_second", real=True)
    p0 = sp.Symbol("beta1_center", real=True, nonzero=True)
    lapse, bdd, ratio = sp.symbols("N0 b_second lapse_ratio", real=True)
    C2 = 3*M**2*r*bdd**2/lapse**2+3*p0*(bdd-rdd)/r
    E0 = 6*M**2*r*bdd/lapse-6*p0+6*p0*lapse/r
    bdd_from_spatial = p0*ratio*(1-ratio)/M**2
    consistency = sp.factor(C2.subs({lapse: r*ratio, bdd: bdd_from_spatial}))
    exact_ratio = 1-M**2*rdd/p0
    f = fixed_c()
    r_uu = -8*(C+2)/C
    exact00 = sp.factor(f["N_center"]**2)
    retained00 = 2*C**2-4
    j = joint()
    point = {V: 0, C: 2, ZETA: 12}
    Gz = sp.simplify(sp.diff(j["G"], ZETA).subs(point))
    z_v = sp.simplify(j["zeta_v"].subs(point))
    z_c = sp.simplify((-sp.diff(j["G"], C)/sp.diff(j["G"], ZETA)).subs(point))
    return {"r0": r, "r_second": rdd, "p0": p0, "N0": lapse,
            "b_second": bdd, "lapse_ratio": ratio,
            "lapse_taylor_order2": C2, "spatial_at_center": E0,
            "b_second_from_spatial": bdd_from_spatial,
            "center_consistency_after_spatial": consistency,
            "exact_lapse_ratio": exact_ratio,
            "profile_r_uu": r_uu, "b_center": sp.Integer(2),
            "b_uu_center": 2*f["b_v_center"], "b_uuuu_center": 12*f["b_vv_center"],
            "N_center": f["N_center"], "N_uu_center": 2*f["N_v_center"],
            "b_TT_center": 2*f["b_v_center"]/TAU**2,
            "N_TT_center": 2*f["N_v_center"]/TAU**2,
            "metric00_center": exact00, "retained_metric00_center": retained00,
            "metric00_center_defect": sp.factor(exact00-retained00),
            "joint_G_zeta_center": Gz, "joint_zeta_v_center": z_v,
            "joint_zeta_c_center": z_c,
            "joint_zeta_center_as_function_c": 3*(C+2)**2/(2*C),
            "joint_point": point}


@cache
def tensors():
    """Literal norm-two homogeneous anisotropies, including physical g EH.

    The common isotropic f is an exact own-f solution; g=eta and the
    prescribed clock need not solve their equations.  Q is a hidden tensor
    amplitude, unrelated to the rational profile numerator named Q.
    """
    e = own_equations()
    b, lapse, beta1 = e["b"], e["N"], e["beta1"]
    q = sp.Function("q_g")(TIME)
    hidden = sp.Function("q_f")(TIME)
    qd, Hd = sp.diff(q, TIME), sp.diff(hidden, TIME)
    K = M**2*b**3/lapse
    e1 = lapse+b*(1+2*sp.cosh((hidden-q)/2))
    e4 = lapse*b**3
    amplitude = sp.Symbol("amplitude", real=True)
    potential = -2*beta1*e1-2*e["beta4"]*e4
    quadratic_potential = sp.diff(potential.subs({q: amplitude*q, hidden: amplitude*hidden}),
                                  amplitude, 2).subs(amplitude, 0)/2
    density = M**2*qd**2/4+K*Hd**2/4+quadratic_potential
    g_equation = -2*(sp.diff(density, q)-sp.diff(sp.diff(density, qd), TIME))
    f_equation = -2*(sp.diff(density, hidden)-sp.diff(sp.diff(density, Hd), TIME))
    nu = 2*beta1*b
    return {"q": q, "Q_tensor": hidden, "K_f": K, "nu": nu,
            "e1": e1, "e4": e4, "literal_potential": potential,
            "quadratic_potential": quadratic_potential, "density": density,
            "g_equation": sp.expand(g_equation), "f_equation": sp.expand(f_equation),
            "g_metric": sp.diag(1, -sp.exp(q), -sp.exp(-q), -1),
            "f_metric": sp.diag(lapse**2, -b**2*sp.exp(hidden), -b**2*sp.exp(-hidden), -b**2),
            "canonical_own_f_pump": sp.diff(sp.sqrt(K), TIME, 2)/sp.sqrt(K),
            "scope": "exact homogeneous restricted action; no full-parent background or boundary/Green choice"}


@cache
def inner():
    """Coefficient limits on each fixed compact x=u/sqrt(delta) interval.

    Joint analyticity makes b,N and their needed v/c derivatives uniformly
    bounded.  ODE solutions converge only for convergent inner-coordinate
    initial data and prescribed inputs; a retarded or other functional
    inverse requires its own state contract.
    """
    x = sp.Symbol("x", real=True)
    p = profile()
    substitution = {V: DELTA*x**2, C: 2+DELTA}
    D_over_delta = sp.factor((p["D"].subs(substitution))/DELTA)
    center_inverse = 8*C/(C-2)
    coupled_inverse = 8*(C**2+16)/(C*(C-2))
    return {"x": x, "D_over_delta": D_over_delta,
            "D_over_delta_limit": 1+8*x**2,
            "K_f_limit": 4*M**2, "N_limit": 2, "b_limit": 2,
            "own_f_inner_coefficient": 16/(1+8*x**2),
            "g_inner_coefficient": 64/(1+8*x**2),
            "coupled_relative_inner_coefficient": 80/(1+8*x**2),
            "own_f_center_mass_squared": center_inverse/TAU**2,
            "coupled_algebraic_center_mass_squared": coupled_inverse/TAU**2,
            "own_f_center_K": 16*M**2/C**2,
            "regular_mass_numerator": 2*p["Q"]*sp.Symbol("N_branch", positive=True)/B**2,
            "canonical_pump_inner_order": "delta*(sqrt(K_f))_uu/sqrt(K_f)=O(delta) on fixed compact x intervals",
            "scope": "fixed compact inner x only; own-f 16 is not coupled-relative 80; no inverse/state or full-window conclusion"}


@cache
def checks():
    e, p, f, j, c, t, i = (own_equations(), profile(), fixed_c(), joint(), center(), tensors(), inner())
    values = {"literal_EH_boundary": e["covariant_EH"]-e["kinetic"]-sp.diff(e["EH_boundary"], TIME),
              "unfactored_consistency": e["unfactored_identity"],
              "profile_beta1": p["B1"]-p["y"]**3*p["h_u"]/(C*p["D"]),
              "profile_root_ratio": p["R"]+p["B4"]/p["B1"],
              "log_root_derivative": p["L"]-sp.diff(p["R"], V)/p["R"],
              "regular_H": p["H"]-p["D"]*(sp.diff(p["Q"], V)/p["Q"]+p["L"])+sp.diff(p["D"], V),
              "fixed_center_equation": f["equation"].subs({V: 0, B: 2}),
              "fixed_center_Jacobian": f["b_jacobian_center"]-96/(C*(C-2)),
              "fixed_center_F": f["F_center"]+4*(C+2)/C,
              "joint_factor_F": f["F"].subs(B**3, j["b_cubed"])-j["F"],
              "joint_F_v": j["F_v"]-sp.diff(j["F"], V),
              "joint_F_zeta": j["F_zeta"]-sp.diff(j["F"], ZETA),
              "joint_b_v_log": sp.diff(j["b_cubed"], V)/(3*j["b_cubed"])-j["b_v_over_b"],
              "joint_b_zeta_log": sp.diff(j["b_cubed"], ZETA)/(3*j["b_cubed"])-j["b_zeta_over_b"],
              "joint_T_v_log": j["T_v"]/j["T"]-j["b_v_over_b"]-2*j["F_v"]/j["F"]+sp.diff(p["Q"], V)/p["Q"],
              "joint_T_zeta_log": j["T_zeta"]/j["T"]-j["b_zeta_over_b"]-2*j["F_zeta"]/j["F"],
              "joint_center_G": j["G"].subs(c["joint_point"]),
              "joint_center_G_zeta": c["joint_G_zeta_center"]+32,
              "joint_center_zeta_v": c["joint_zeta_v_center"]-204,
              "joint_center_zeta_c": c["joint_zeta_c_center"],
              "joint_center_lapse": j["N"].subs(c["joint_point"])-2,
              "center_spatial_solution": c["spatial_at_center"].subs({c["N0"]: c["r0"]*c["lapse_ratio"],
                                             c["b_second"]: c["b_second_from_spatial"]}),
              "center_no_bdd_division": M**2*c["r0"]*c["center_consistency_after_spatial"]/3
                   -c["p0"]*(c["p0"]*(1-c["lapse_ratio"])-M**2*c["r_second"]),
              "profile_center_lapse": c["N_center"]-C**2/2,
              "profile_center_b_uu": c["b_uu_center"]+2*C*(C+2),
              "profile_center_N_uu": c["N_uu_center"]-(27*C**4-12*C**3-236*C**2+592*C-512)/4,
              "profile_exact_metric00": c["metric00_center"]-C**4/4,
              "profile_retained_metric00": c["retained_metric00_center"]-4
                   +2*4*c["profile_r_uu"]/(32/(C*(C-2))),
              "profile_exact_center_defect": c["metric00_center_defect"]-(C**2-4)**2/4,
              "literal_tensor_potential": t["quadratic_potential"]+e["beta1"]*e["b"]*(t["Q_tensor"]-t["q"])**2/2,
              "literal_g_tensor_equation": t["g_equation"]-M**2*sp.diff(t["q"], TIME, 2)-t["nu"]*(t["q"]-t["Q_tensor"]),
              "literal_f_tensor_equation": t["f_equation"]-sp.diff(t["K_f"]*sp.diff(t["Q_tensor"], TIME), TIME)
                   -t["nu"]*(t["Q_tensor"]-t["q"]),
              "inner_D_limit": i["D_over_delta"].subs(DELTA, 0)-i["D_over_delta_limit"],
              "inner_own_vs_relative": i["coupled_relative_inner_coefficient"]-i["own_f_inner_coefficient"]-i["g_inner_coefficient"],
              "inner_own_center_K": i["own_f_center_K"]-M**2*2**3/c["N_center"],
              "inner_own_center_mass": i["own_f_center_mass_squared"]-2*p["beta1"].subs(V, 0)*c["N_center"]/(4*M**2),
              "inner_coupled_center_mass": i["coupled_algebraic_center_mass_squared"]
                   -4*p["beta1"].subs(V, 0)*(1/M**2+1/i["own_f_center_K"])}
    values = {name: sp.factor(sp.cancel(value)) for name, value in values.items()}
    if any(value != 0 for value in values.values()):
        raise ValueError("An exact stationary-branch identity failed")
    return values


def calibration():
    return {"claim": "P8-S6.33.EXACT_STATIONARY",
            "physical_domain": "0<delta<=1/100, |T|<=tau/100; c=2+delta; M,tau>0",
            "normalized_center_mass": "own-f16 and coupled-relative80 are different inner operators",
            "center_lapse": C**2/2, "center_metric00": C**4/4,
            "center_b_second": -2*C*(C+2)/TAU**2,
            "joint_IFT_derivative": -32,
            "regularity": "b C2, N C1 and positive; parity follows from local uniqueness, not an initial assumption",
            "box_proof": "separate exact Fraction interval certificate supplies the uniform Banach/positive-lapse margins",
            "scope": "exact isotropic own-f on fixed g=eta, theta=T; not a full-parent/CD solution, selected Green inverse or controlled matching"}
