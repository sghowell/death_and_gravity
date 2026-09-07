"""Source-preserving, formal own-f derivative expansion of frozen VARIABLE.

No physical metric or matter redefinition is made.  The algebraic inverse
below is not a Green operator or a proof of a controlled derivative expansion.
Curvature is P8(b)'s +--- R_B; both Einstein terms are -M**2 R_B/2.
"""

from functools import cache

import sympy as sp

M, TAU, R = sp.symbols("M tau r", positive=True)
BETA0 = sp.Symbol("beta0", real=True)
BETA1 = sp.Symbol("beta1", real=True, nonzero=True)
C = sp.Symbol("c", positive=True)
V = sp.Symbol("v", real=True)
ETA = sp.diag(1, -1, -1, -1)


def _symmetric(prefix):
    entries = sp.symbols(" ".join(f"{prefix}{i}{j}" for i in range(4) for j in range(i, 4)), real=True)
    matrix = sp.zeros(4)
    index = 0
    for i in range(4):
        for j in range(i, 4):
            matrix[i, j] = matrix[j, i] = entries[index]
            index += 1
    return matrix, entries


def _contract(left, right):
    return sp.trace(ETA*left*ETA*right)


@cache
def potential_hessian():
    """All ten covariant components, including lapse/trace and shifts.

    f=r²(g+h); H=g^-1 h.  beta4*r**4=-r*beta1 is imposed only after
    differentiating the literal square-root/determinant potential series.
    """
    h, entries = _symmetric("h")
    mixed = ETA*h
    t1, t2, t3 = (sp.trace(mixed**j) for j in (1, 2, 3))
    b = R*BETA1
    quadratic = -b*(t2-t1**2)/4
    cubic = b*(t1**3-6*t1*t2+5*t3)/24
    E, source_entries = _symmetric("E")
    linear = M**2*R**2*_contract(E, h)/2
    trace_E = sp.trace(ETA*E)
    solution = M**2*R/BETA1*(E-ETA*trace_E/3)
    substitution = dict(zip(entries, (solution[i, j] for i in range(4) for j in range(i, 4)), strict=True))
    kappa = M**4*R**3/(4*BETA1)
    stationary = sp.expand((linear+quadratic).subs(substitution, simultaneous=True))
    return {"h": h, "entries": entries, "H": mixed, "trace": t1,
            "trace_square": t2, "trace_cube": t3, "b": b,
            "quadratic": quadratic, "cubic": cubic, "linear_EH": linear,
            "E": E, "source_entries": source_entries, "trace_E": trace_E,
            "h_solution": solution, "substitution": substitution,
            "kappa": kappa, "stationary_four_derivative": stationary,
            "hessian": sp.hessian(quadratic, entries)}


@cache
def own_f_correction():
    """Covariant Schouten-type correction, with no g/source equation used."""
    p = potential_hessian()
    ricci, _ = _symmetric("Ric")
    scalar = sp.trace(ETA*ricci)
    E = ricci-ETA*scalar/2
    E_to_ricci = dict(zip(p["source_entries"],
                          (E[i, j] for i in range(4) for j in range(i, 4)), strict=True))
    h = p["h_solution"].subs(E_to_ricci, simultaneous=True).applyfunc(sp.expand)
    # ricci here is R_mu_nu[f0], and scalar=tr_g Ric=r² R[f0].
    expected = M**2*R/BETA1*(ricci-ETA*scalar/6)
    return {"ricci_f0": ricci, "trace_g_ricci_f0": scalar,
            "h": h, "delta_f": R**2*h, "expected_h": expected,
            "coefficient_delta_f": M**2*R**3/BETA1,
            "curvature_order": 2, "first_own_equation_residual_order": 4,
            "retained_action_order": 4, "first_omitted_action_order": 6}


@cache
def conformal_action():
    """Finite tensor polynomial with all log-r clock derivatives retained.

    s=(log r)_phi, t=s_phi.  The pointwise Lorentz frame is only used for
    tensor contraction; Ricci and phi Hessian are unrestricted symmetric
    tensors, and the gradient has all four components.
    """
    ricci, _ = _symmetric("Ric")
    hessian, _ = _symmetric("Phi")
    gradient = sp.Matrix(sp.symbols("p0:4", real=True))
    s, t, kappa, kp = sp.symbols("s t kappa kappa_phi", real=True)
    scalar = sp.trace(ETA*ricci)
    einstein = ricci-ETA*scalar/2
    X = (gradient.T*ETA*gradient)[0]
    box = sp.trace(ETA*hessian)
    sandwich = (gradient.T*ETA*hessian*ETA*gradient)[0]
    ric_pp = (gradient.T*ETA*ricci*ETA*gradient)[0]
    G_pp = (gradient.T*ETA*einstein*ETA*gradient)[0]
    L1, L2 = _contract(hessian, hessian), box**2
    B = einstein-2*s*hessian+2*(s**2-t)*gradient*gradient.T+ETA*(2*s*box+(2*t+s**2)*X)
    trace_B = sp.trace(ETA*B)
    Qg = _contract(ricci, ricci)-scalar**2/3
    Qhat = _contract(B, B)-trace_B**2/3
    expanded = (Qg-4*s*_contract(einstein, hessian)+4*(s**2-t)*ric_pp+2*t*scalar*X
                +4*s**2*(L1-L2)+8*s*(t-s**2)*sandwich
                -4*s*(2*t+s**2)*X*box-12*s**2*t*X**2)
    w = s*gradient
    Hs = s*hessian+t*gradient*gradient.T
    hs, wsq = sp.trace(ETA*Hs), (w.T*ETA*w)[0]
    Vup = ETA*einstein*ETA*w+(ETA*w)*hs-ETA*Hs*ETA*w+wsq*ETA*w
    div_V = (_contract(einstein, Hs)+hs**2-_contract(Hs, Hs)
             -(w.T*ETA*ricci*ETA*w)[0]+2*(w.T*ETA*Hs*ETA*w)[0]+hs*wsq)
    compact_ibp = kappa*Qg+4*kp*(s*G_pp+s**2*(X*box-sandwich)+s**3*X**2)
    rphi = sp.Symbol("r_phi", real=True)
    Y = sp.Symbol("Y", real=True)
    leading = -M**2*(1+R**2)*scalar/2+(sp.Rational(1, 2)-3*M**2*rphi**2)*X+Y/2-2*(BETA0+3*R*BETA1)
    # P[f0]=Ric[f0]-f0 R[f0]/6, not half the Schouten convention.
    P_expected = (ricci-ETA*scalar/6-2*Hs+2*w*w.T-ETA*wsq)
    return {"ricci": ricci, "hessian": hessian, "gradient": gradient,
            "s": s, "t": t, "kappa": kappa, "kappa_phi": kp,
            "R_B": scalar, "X": X, "Y": Y, "box": box,
            "B": B, "trace_B": trace_B, "Qg": Qg, "Qhat": Qhat,
            "expanded_Qhat": expanded, "Vup": Vup, "div_V": div_V,
            "compact_ibp": compact_ibp, "leading_density": leading,
            "P_expected": P_expected, "r_phi": rphi}


@cache
def constant_ratio_calibration():
    """Full two-TT own-f Schur; source stays on g with action j*h_g/4.

    This is a constant-coefficient/off-shell algebra control.  It does not
    assert that the frozen local action has a constant-clock flat vacuum.
    """
    D, hg, hf, j = sp.symbols("D h_g h_f j", real=True)
    G, F, nu = M**2, M**2*R**2, 2*R*BETA1
    action = -(G*D*hg**2+F*D*hf**2+nu*(hg-hf)**2)/8+j*hg/4
    hf_star = nu*hg/(nu+F*D)
    kernel = sp.factor(G*D+nu-nu**2/(nu+F*D))
    retained = (G+F)*D-F**2*D**2/nu
    remainder = F**3*D**3/(nu*(nu+F*D))
    return {"D": D, "h_g": hg, "h_f": hf, "j": j,
            "G": G, "F": F, "nu": nu, "action": action,
            "f_stationary": hf_star, "kernel": kernel,
            "retained_kernel": retained, "kernel_remainder": remainder,
            "c_C": F**2/(4*nu), "kappa": M**4*R**3/(4*BETA1),
            "first_omitted_flat_kernel_coefficient": F**3/nu**2,
            "six_derivative_flat_curvature_coefficient": -F**3/(2*nu**2),
            "own_f_inverse_parameter": F/nu,
            "source_scope": "physical g unchanged; positive canonical T has delta S=+T_mn delta g^mn/2, hence j=+2 Pi_TT for delta g_ij=-h_g e_ij and e:e=1"}


@cache
def profile():
    """Exact even-u profile, v=u².  beta_n=(M²/tau²)*b_n.

    r is the algebraic stationary ratio, not the actual y or lapse c.
    The rational r³ and log-r jets extend to c=2 at v=0 only as limits;
    c=2 is not an admissible member of the frozen action family.
    """
    v, c = V, C
    d = 1+v
    y = 2/d**4
    h2, hp = 16*v/d**2, 4*(1-v)/d**2
    b1 = y**3*hp/(c*(c-y))
    b4 = 3*h2/(2*c**2)-b1/y**3
    theta = 3*h2*(c-y)/(2*c*hp)
    r3 = sp.factor(y**3/(1-theta))
    ell = sp.factor(2*sp.diff(r3, v)/(3*r3))  # (log r)_u=u*ell(v)
    P00 = -2*(hp+ell+2*v*sp.diff(ell, v))-h2+v*ell**2
    kappa_bar = sp.factor(-1/(4*b4))
    f2tt = sp.factor(4*kappa_bar*P00)
    f2space = sp.factor(4*kappa_bar*v*(4/d+ell)**2)
    kbar = 2*(y**3/c-1)*hp-1/(100*d**12)
    return {"v": v, "c": c, "d": d, "y": y, "h_squared": h2,
            "h_u": hp, "b1": b1, "b4": b4, "theta": theta,
            "r_cubed": r3, "log_r_u_over_u": ell,
            "kappa_bar": kappa_bar, "P00_bar": P00,
            "f2_00": f2tt, "f2_ii_over_a_squared": f2space,
            "kbar": kbar}


@cache
def center_jets():
    """Exact physical-clock jets and continuous nonuniformity controls.

    f2 is the additive covariant-metric correction, not h=f2/r².  A spatial
    entry below is divided by a², with f_ii/a²=-y².  u=T/tau and
    varphi=phi/M satisfies varphi_u=sqrt(kbar).
    """
    p, c, v = profile(), C, V

    def jets(expr):
        return tuple(sp.factor(factor*sp.diff(expr, v, order).subs(v, 0))
                     for order, factor in ((0, 1), (1, 2), (2, 12)))

    r3_0, r3_v, r3_vv = (sp.factor(sp.diff(p["r_cubed"], v, n).subs(v, 0)) for n in range(3))
    r2jets = (sp.Integer(4), sp.factor(2*r3_v/3), sp.factor(4*r3_vv-r3_v**2/6))
    f2jets, spacejets, kappajets = (jets(p[name]) for name in
                                  ("f2_00", "f2_ii_over_a_squared", "kappa_bar"))
    rem = (sp.factor(c**2-r2jets[0]-f2jets[0]),
           sp.factor(-r2jets[1]-f2jets[1]), sp.factor(-r2jets[2]-f2jets[2]))
    k0 = sp.factor(p["kbar"].subs(v, 0))
    k2 = sp.factor(2*sp.diff(p["kbar"], v).subs(v, 0))
    clock2 = sp.factor(f2jets[1]/k0)
    clock4 = sp.factor(rem[2]/k0**2-2*rem[1]*k2/k0**3)
    alpha0 = sp.factor(kappajets[0]/2)
    alpha2 = sp.factor(kappajets[1]/2-kappajets[0]*r2jets[1]/8)
    delta = sp.Symbol("delta", positive=True)
    # Numerator certificates hold throughout 2<c<=4; no sampling is used.
    rem4_margin = sp.Poly(sp.cancel(c**2*(rem[2]-10752)).subs(c, 2+delta).expand(), delta)
    first2_margin = sp.Rational(8)-sp.Rational(284, 100)-sp.Rational(184, 10000)-sp.Rational(108, 1000000)
    return {"r_squared_u0_u2_u4": r2jets, "f2_00_u0_u2_u4": f2jets,
            "f2_space_over_a2_u0_u2_u4": spacejets,
            "kappa_bar_u0_u2_u4": kappajets,
            "own_f_inverse_parameter_bar_u0_u2": (alpha0, alpha2),
            "metric00_remainder_u0_u2_u4": rem,
            "r_cubed_center": r3_0, "kbar_center": k0, "kbar_u2_center": k2,
            "f2_00_varphi2": clock2, "remainder00_varphi4": clock4,
            "f2_00_varphi2_limit": sp.factor(clock2.subs(c, 2)),
            "remainder00_varphi4_limit": sp.factor(clock4.subs(c, 2)),
            "remainder_u4_minus_10752_cleared": rem4_margin,
            "f2_u2_above_60_for_delta_le_1over100_margin": first2_margin,
            "clock4_uniform_lower_bound": sp.Rational(107520000, 5755201),
            "lapse_center_exact": c**2, "lapse_center_retained": 4*c-4,
            "inner_metric00_defect_leading": 1+64*sp.Symbol("x", real=True)**2+448*sp.Symbol("x", real=True)**4,
            "domain": "2<c<=4, |u|<=1/10; limits c->2+ are not an action at c=2"}


def domain_margins():
    """Continuous positive-root/clock bounds for the stated positive branch."""
    dmax = sp.Rational(101, 100)
    return {"y_above_19over10": 2/dmax**4-sp.Rational(19, 10),
            "h_u_above_3": 4*sp.Rational(99, 100)/dmax**2-3,
            "theta_upper_bound": sp.Rational(7, 220),
            "one_minus_theta_lower_bound": sp.Rational(213, 220),
            "r_below_21over10_cubed_margin": sp.Rational(21, 10)**3-8*sp.Rational(220, 213),
            "kbar_center_lower_bound": sp.Rational(799, 100),
            "minus_kbar_u2_center_lower_bound": sp.Rational(10794, 25)}


@cache
def checks():
    p, a, f, z, j = (potential_hessian(), conformal_action(), own_f_correction(),
                     constant_ratio_calibration(), center_jets())
    eps = sp.Symbol("epsilon", real=True)
    I = sp.eye(4)
    H = p["H"]
    root = I+eps*H/2-eps**2*H**2/8+eps**3*H**3/16
    volume = (1+eps*p["trace"]/2+eps**2*(p["trace"]**2-2*p["trace_square"])/8
              +eps**3*(p["trace"]**3-6*p["trace"]*p["trace_square"]+8*p["trace_cube"])/48)
    # Exact Newton identities are checked against a literal 4x4 determinant.
    det = (I+eps*H).det(method="domain-ge")
    values = {f"literal_volume_square_order_{n}": sp.expand(volume**2-det).coeff(eps, n) for n in range(4)}
    for n in range(1, 4):
        values[f"literal_square_root_order_{n}"] = sum(
            sp.expand(entry).coeff(eps, n)**2 for entry in root*root-I-eps*H)
    literal = -2*p["b"]*(sp.trace(root)-volume)
    values.update({"literal_linear_potential_vanishes": sp.expand(literal).coeff(eps, 1),
                   "literal_quadratic_potential": sp.expand(literal).coeff(eps, 2)-p["quadratic"],
                   "literal_cubic_potential": sp.expand(literal).coeff(eps, 3)-p["cubic"],
                   "full_algebraic_hessian_determinant": p["hessian"].det()-3*p["b"]**10/16,
                   "stationary_quadratic_action": p["stationary_four_derivative"]-p["kappa"]*(
                       _contract(p["E"], p["E"])-p["trace_E"]**2/3)})
    for n, entry in enumerate(p["entries"]):
        values[f"all_ten_own_f_stationary_equations_{n}"] = sp.diff(p["linear_EH"]+p["quadratic"], entry).subs(p["substitution"], simultaneous=True)
    for n, entry in enumerate(f["h"]-f["expected_h"]):
        values[f"own_f_Schouten_components_{n}"] = entry
    values.update({
        "full_conformal_tensor_square": a["Qhat"]-a["expanded_Qhat"],
        "conformal_Euler_boundary_identity": a["Qhat"]-a["Qg"]+4*a["div_V"],
        "all_clock_derivatives_in_P": sum(sp.expand(entry)**2 for entry in
            a["B"]-ETA*a["trace_B"]/3-a["P_expected"]),
        "variable_kappa_boundary_contraction": a["compact_ibp"]-a["kappa"]*a["Qg"]
            -4*a["kappa_phi"]*(a["gradient"].T*a["Vup"])[0],
        "constant_ratio_own_equation": sp.diff(z["action"], z["h_f"]).subs(z["h_f"], z["f_stationary"]),
        "constant_ratio_stationary_action": z["action"].subs(z["h_f"], z["f_stationary"])
            +z["kernel"]*z["h_g"]**2/8-z["j"]*z["h_g"]/4,
        "constant_ratio_exact_remainder": z["kernel"]-z["retained_kernel"]-z["kernel_remainder"],
        "constant_ratio_all_TT_Weyl_normalization": z["c_C"]-z["kappa"]/2,
        "constant_ratio_first_omitted_six_derivative_symbol": sp.cancel(z["kernel_remainder"]/z["D"]**3).subs(z["D"], 0)
            +2*z["six_derivative_flat_curvature_coefficient"],
        "center_lapse_defect": j["lapse_center_exact"]-j["lapse_center_retained"]-(C-2)**2,
        "center_r2_u2": j["r_squared_u0_u2_u4"][1]+32*(C+2)/C,
        "center_first_correction_u2": j["f2_00_u0_u2_u4"][1]
            +4*(27*C**3-116*C**2+196*C-176)/C,
        "center_remainder_u2": j["metric00_remainder_u0_u2_u4"][1]
            -4*(C-2)*(27*C**2-62*C+80)/C,
        "center_kappa_u2": j["kappa_bar_u0_u2_u4"][1]-(9*C**2-22*C+24)/8,
        "center_own_inverse_parameter_u2": j["own_f_inverse_parameter_bar_u0_u2"][1]
            -(13*C**2-22*C+8)/16,
        "inner_metric_defect_polynomial": j["inner_metric00_defect_leading"]-1
            -sp.cancel(j["metric00_remainder_u0_u2_u4"][1]/(2*(C-2))).subs(C, 2)*sp.Symbol("x", real=True)**2
            -j["metric00_remainder_u0_u2_u4"][2].subs(C, 2)*sp.Symbol("x", real=True)**4/24,
        "center_nonzero_first_clock_derivative_limit": j["f2_00_varphi2_limit"]-sp.Rational(6400, 2399),
        "center_nonzero_remainder_clock_derivative_limit": j["remainder00_varphi4_limit"]-sp.Rational(107520000, 5755201),
    })
    values = {name: sp.factor(sp.expand(value)) for name, value in values.items()}
    if any(value != 0 for value in values.values()):
        raise ValueError("An own-f stationary identity failed")
    return values


def calibration():
    p, j = potential_hessian(), center_jets()
    coefficients = j["remainder_u4_minus_10752_cleared"].all_coeffs()
    if any(value < 0 for value in coefficients) or not any(value > 0 for value in coefficients):
        raise ValueError("The continuous fourth-jet polynomial certificate failed")
    if j["f2_u2_above_60_for_delta_le_1over100_margin"] <= 0:
        raise ValueError("The small-delta second-jet margin failed")
    if any(value <= 0 for value in domain_margins().values()):
        raise ValueError("A continuous positive-root domain margin failed")
    return {"formal_kappa": p["kappa"], "formal_c_C_constant_r": p["kappa"]/2,
            "center_kappa": M**2*TAU**2*j["kappa_bar_u0_u2_u4"][0],
            "center_lapse_defect": (C-2)**2,
            "first_correction_u2_limit": 64, "remainder_u4_limit": 10752,
            "first_correction_varphi2_limit": j["f2_00_varphi2_limit"],
            "remainder_varphi4_lower_bound": j["clock4_uniform_lower_bound"],
            "actual_clock": "varphi=phi/M, varphi_u=sqrt(kbar); d_T=tau^-1 d_u",
            "remainder_scope": "No C2-small first correction or C4-small metric remainder as c->2+; no EFT-wide exclusion",
            "first_omitted_action": "S6=1/2 delta_f2 S_EH,f'' delta_f2 + 1/6 S_pot'''[delta_f2]^3",
            "stationary_branch_scope": "formal local series only; a differential inverse, state, boundary data and a remainder remain unproved",
            "matter_scope": "physical g and literal canonical free chi are unchanged"}
