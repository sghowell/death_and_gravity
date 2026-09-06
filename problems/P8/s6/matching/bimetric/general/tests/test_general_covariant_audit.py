"""Independent arbitrary-beta two-lapse and two-source obstruction audit."""

from functools import cache

import sympy as sp


@cache
def _literal_variation():
    t = sp.Symbol("t", real=True)
    a, b, ng, nf = (sp.Function(name, positive=True)(t) for name in ("a", "b", "Ng", "Nf"))
    kg, kf, vg, vf = (sp.Function(name)(t) for name in ("Kg", "Kf", "Vg", "Vf"))
    mg, mf, interaction_scale = sp.symbols("G F m4", positive=True)
    beta = sp.symbols("beta0:5", real=True)
    z = sp.Symbol("z")
    y, c = b/a, nf/ng
    coefficients = sp.Poly((1+z*c)*(1+z*y)**3, z)
    potential = interaction_scale*ng*a**3*sum(beta[n]*coefficients.nth(n) for n in range(5))
    matter_g = a**3*kg/(2*ng)-ng*a**3*vg
    matter_f = b**3*kf/(2*nf)-nf*b**3*vf
    lagrangian = -3*mg*a*sp.diff(a, t)**2/ng-3*mf*b*sp.diff(b, t)**2/nf-potential+matter_g+matter_f
    euler = lambda field: sp.diff(lagrangian, field)-sp.diff(sp.diff(lagrangian, sp.diff(field, t)), t)
    hg, hf = sp.diff(a, t)/(ng*a), sp.diff(b, t)/(nf*b)
    rho_g, pressure_g = kg/(2*ng**2)+vg, kg/(2*ng**2)-vg
    rho_f, pressure_f = kf/(2*nf**2)+vf, kf/(2*nf**2)-vf
    density_gv = interaction_scale*(beta[0]+3*beta[1]*y+3*beta[2]*y*y+beta[3]*y**3)
    density_fv = interaction_scale*(beta[4]+3*beta[3]/y+3*beta[2]/y**2+beta[1]/y**3)
    pressure_gv = -interaction_scale*(beta[0]+beta[1]*(2*y+c)+beta[2]*(y*y+2*y*c)+beta[3]*c*y*y)
    pressure_fv = -interaction_scale*(beta[4]+beta[3]*(2/y+1/c)+beta[2]*(1/y**2+2/(c*y))+beta[1]/(c*y*y))
    normalized = [euler(ng)/a**3, euler(nf)/b**3,
                  euler(a)/(3*ng*a*a), euler(b)/(3*nf*b*b)]
    expected = [3*mg*hg*hg-rho_g-density_gv, 3*mf*hf*hf-rho_f-density_fv,
                mg*(2*sp.diff(hg, t)/ng+3*hg*hg)+pressure_g+pressure_gv,
                mf*(2*sp.diff(hf, t)/nf+3*hf*hf)+pressure_f+pressure_fv]
    return {"t": t, "a": a, "b": b, "ng": ng, "nf": nf, "beta": beta,
            "mg": mg, "mf": mf, "m4": interaction_scale, "kg": kg, "kf": kf,
            "vg": vg, "vf": vf, "y": y, "c": c, "hg": hg, "hf": hf,
            "matter_g": matter_g, "matter_f": matter_f, "normalized": normalized,
            "expected": expected, "density_gv": density_gv, "density_fv": density_fv,
            "pressure_gv": pressure_gv, "pressure_fv": pressure_fv}


def test_general_four_equations_come_from_generating_potential_before_gauge_fixing():
    data = _literal_variation()
    for actual, expected in zip(data["normalized"], data["expected"], strict=True):
        assert sp.simplify(actual-expected) == 0
    assert sp.diff(data["matter_g"], data["nf"]) == 0
    assert sp.diff(data["matter_f"], data["ng"]) == 0
    # Two separate actions are not one field with simultaneous metric coupling.
    assert sp.diff(data["matter_g"], data["kf"]) == 0
    assert sp.diff(data["matter_f"], data["kg"]) == 0


def test_arbitrary_beta_null_stresses_and_both_interaction_balances():
    data = _literal_variation()
    y, c = data["y"], data["c"]
    beta, m4, t = data["beta"], data["m4"], data["t"]
    polynomial = m4*(beta[1]+2*beta[2]*y+beta[3]*y*y)
    g_null, f_null = (data[f"density_{label}v"]+data[f"pressure_{label}v"] for label in ("g", "f"))
    assert sp.factor(g_null-(y-c)*polynomial) == 0
    assert sp.factor(f_null-(c-y)*polynomial/(c*y**3)) == 0
    difference = y*data["hf"]-data["hg"]
    g_balance = sp.diff(data["density_gv"], t)/data["ng"]+3*data["hg"]*g_null
    f_balance = sp.diff(data["density_fv"], t)/data["nf"]+3*data["hf"]*f_null
    assert sp.factor(g_balance-3*c*polynomial*difference) == 0
    assert sp.factor(f_balance+3*polynomial*difference/(c*y**3)) == 0
    assert sp.factor(g_null+c*y**3*f_null) == 0


def test_general_branch_cancellation_includes_both_independent_NEC_sectors():
    mg, mf, y, c = sp.symbols("G F y c", positive=True)
    hdot, polynomial = sp.symbols("Hdot P", real=True)
    null_g, null_f = sp.symbols("rho_g_plus_p_g rho_f_plus_p_f", nonnegative=True)
    g_equation = -2*mg*hdot-null_g-(y-c)*polynomial
    f_equation = -2*mf*hdot/(c*y)-null_f-(c-y)*polynomial/(c*y**3)
    expected = -2*(mg+mf*y*y)*hdot-null_g-c*y**3*null_f
    assert sp.factor(g_equation+c*y**3*f_equation-expected) == 0
    positive_hdot = sp.Symbol("positive_Hdot", positive=True)
    assert (-expected).subs(hdot, positive_hdot).is_positive is True
    assert sp.factor(g_equation.subs(polynomial, 0)+2*mg*hdot+null_g) == 0
    # This second case includes isolated roots, an identically zero polynomial
    # and branch switching: there is no division by P in the root case.
    assert (-g_equation.subs(polynomial, 0)).subs(hdot, positive_hdot).is_positive is True
    assert sp.factor(expected.subs(null_f, 0)+2*(mg+mf*y*y)*hdot+null_g) == 0


def test_differentiated_regular_branch_has_the_required_lapse_and_ratio_weights():
    t = sp.Symbol("t", real=True)
    y, c, hf = (sp.Function(name)(t) for name in ("y", "c", "Hf"))
    hg = y*hf
    at_stationarity = sp.diff(hg, t).subs(hf, 0)
    # Substitute the value only, not the derivative jet, at the single slice.
    hfdot = sp.Symbol("Hf_coordinate_derivative")
    assert sp.expand(sp.diff(hg, t).subs({sp.diff(hf, t): hfdot, hf: 0}, simultaneous=True)-y*hfdot) == 0
    assert at_stationarity.has(sp.Derivative)
    hdot = sp.Symbol("Hdot")
    assert sp.factor((hfdot/c).subs(hfdot, hdot/y)-hdot/(c*y)) == 0


def test_CD_window_requires_a_forbidden_non_degenerate_zero_even_if_shifted():
    t = sp.Symbol("t", real=True)
    tau = sp.Symbol("tau", positive=True)
    a = (1+t*t/tau**2)**2
    hubble = sp.diff(a, t)/a
    assert sp.simplify(hubble.subs(t, tau/2)-8/(5*tau)) == 0
    assert sp.simplify(hubble.subs(t, -tau/2)+8/(5*tau)) == 0
    w = sp.Symbol("w", nonnegative=True)
    derivative = sp.cancel(tau**2*sp.diff(hubble, t)).subs(t, tau*sp.sqrt(w))
    assert sp.factor(derivative-sp.Rational(48, 25)-(1-4*w)*(52+12*w)/(25*(1+w)**2)) == 0


def test_public_general_action_matches_all_independent_lapse_first_equations():
    from p8_bimetric_general import background as public

    data, actual = _literal_variation(), public.equations()
    replacements = {data["a"]: public.a, data["b"]: public.b,
                    data["ng"]: public.Ng, data["nf"]: public.Nf,
                    data["mg"]: public.MG2, data["mf"]: public.MF2, data["m4"]: public.M4,
                    data["kg"]: public.KG, data["kf"]: public.KF,
                    data["vg"]: public.VG, data["vf"]: public.VF}
    replacements.update(dict(zip(data["beta"], public.BETAS, strict=True)))
    for name, independent in zip(("EL_Ng", "EL_Nf", "EL_a", "EL_b"), data["normalized"], strict=True):
        assert sp.simplify(actual[name]-independent.subs(replacements, simultaneous=True).doit()) == 0


def test_public_root_case_does_not_require_a_stationary_second_metric():
    from p8_bimetric_general import background as public

    t = public.t
    hdot, null_g = sp.symbols("Hdot null_g", real=True)
    y = public.b/public.a
    equation = public.equations()["EL_a"]-public.equations()["EL_Ng"]
    root_equation = equation.subs(public.BETAS[1], -2*public.BETAS[2]*y-public.BETAS[3]*y*y)
    point = {public.Ng: 1, sp.diff(public.Ng, t): 0, sp.diff(public.a, t): 0,
             sp.diff(public.a, t, 2): public.a*hdot, public.KG: null_g}
    # No b_dot, b_ddot, Nf or f-matter condition is inserted in this case.
    assert sp.factor(root_equation.subs(point, simultaneous=True)-2*public.MG2*hdot-null_g) == 0
