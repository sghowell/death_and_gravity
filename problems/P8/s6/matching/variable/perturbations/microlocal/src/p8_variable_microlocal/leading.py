"""Exact leading-action stationary map and candidate principal matrices.

The separate weighted/Cauchy proof is needed to identify these matrices
with a uniform physical principal operator. This module itself proves
only its literal algebra and boundary dictionary.
"""

from functools import cache

import sympy as sp

from . import stueckelberg as st


@cache
def derive():
    d = st.derive()
    u, a, b, c, y, P = (d[key] for key in ("u", "a", "b", "c", "y", "P"))
    L1 = d["boundary_subtracted_L1"]
    U, V = sp.factor(a**2*y*P/2), sp.factor(a**2*P/(2*y))
    first = {st.PSIG: -U*st.PI, st.PSIF: V*st.PI}
    constrained = L1.subs(first, simultaneous=True)
    xf_hessian = sp.factor(sp.diff(constrained, st.XF, 2))
    xf_linear = sp.factor(sp.diff(constrained, st.XF).subs(st.XF, 0))
    xf = sp.factor(-xf_linear/xf_hessian)
    stationary = {**first, st.XF: xf}
    phig = sp.factor(-sp.diff(L1, st.PSIG).subs(
        {**stationary, st.PHIG: 0, st.PHIF: 0}, simultaneous=True)/(2*a))
    phif = sp.factor(-sp.diff(L1, st.PSIF).subs(
        {**stationary, st.PHIG: 0, st.PHIF: 0}, simultaneous=True)/(2*c*b))
    stationary.update({st.PHIG: phig, st.PHIF: phif})
    stationary_L1 = constrained.subs(st.XF, xf)
    # The remaining pi*pi' term is a time boundary with a nonconstant
    # coefficient. Its coefficient can vanish at u=0 while its derivative
    # contributes to the gradient matrix there.
    fpi_coefficient = sp.factor(sp.diff(stationary_L1, st.VPI, st.PI))
    Fpi = fpi_coefficient*st.PI**2/2
    final_L1 = stationary_L1-st.total_time(Fpi)
    gradient = (-sp.hessian(final_L1, st.PHYSICAL)).applyfunc(sp.factor)
    zeroth_substitution = {**stationary,
        st.VPSIG: -U*st.VPI-sp.diff(U, u)*st.PI,
        st.VPSIF: V*st.VPI+sp.diff(V, u)*st.PI,
    }
    stationary_L0 = d["boundary_subtracted_L0"].subs(zeroth_substitution, simultaneous=True)
    kinetic = sp.hessian(stationary_L0, st.PHYSICAL_VELOCITIES).applyfunc(sp.factor)
    weight = sp.diag(1/U, 1, 1)
    normalized_kinetic = (weight.T*kinetic*weight).applyfunc(sp.factor)
    normalized_gradient = (weight.T*gradient*weight).applyfunc(sp.factor)
    return {**d, "U": U, "V": V, "stationary_fields": stationary,
            "xf_hessian": xf_hessian, "stationary_L1": stationary_L1,
            "Fpi_coefficient": fpi_coefficient, "Fpi": Fpi,
            "final_L1": final_L1, "stationary_L0": stationary_L0,
            "kinetic": kinetic, "gradient": gradient,
            "normalized_kinetic": normalized_kinetic,
            "normalized_gradient": normalized_gradient,
            "normalized_field": "Pi=U*pi; time-connection terms are not included in a principal coefficient"}


@cache
def center():
    d = derive()
    u, c = d["u"], d["c"]
    at = {u: 0}
    kinetic = d["kinetic"].subs(at).applyfunc(sp.factor)
    gradient = d["gradient"].subs(at).applyfunc(sp.factor)
    expected_kinetic = sp.diag(30720/(c**2*(c-2)**2), (6400-801*c)/(100*c), 1)
    expected_gradient = sp.diag(256*(9*c**2-22*c+144)/(c**2*(c-2)**2),
                                (6400-801*c)/(100*c), 1)
    return {"c": c, "kinetic": kinetic, "gradient": gradient,
            "expected_kinetic": expected_kinetic, "expected_gradient": expected_gradient,
            "helicity_speed_squared": (9*c**2-22*c+144)/120,
            "Fpi_coefficient_at_center": sp.factor(d["Fpi_coefficient"].subs(at)),
            "Fpi_coefficient_derivative_at_center": sp.factor(sp.diff(d["Fpi_coefficient"], u).subs(at)),
            "weight_second_log_jet": sp.factor((sp.diff(d["U"], u, 2)/d["U"]
                                                  -(sp.diff(d["U"], u)/d["U"])**2).subs(at)),
            "gradient_without_final_boundary": (-sp.hessian(d["stationary_L1"], st.PHYSICAL)).subs(at).applyfunc(sp.factor)}


@cache
def checks():
    d = derive()
    result = {}
    for i, q in enumerate(st.LEADING_AUXILIARIES):
        result[f"stationary_equation_{i}"] = sp.factor(sp.diff(
            d["boundary_subtracted_L1"], q).subs(d["stationary_fields"], simultaneous=True))
    for i, v in enumerate(st.PHYSICAL_VELOCITIES):
        result[f"final_L1_no_velocity_{i}"] = sp.factor(sp.diff(d["final_L1"], v))
    for i, entry in enumerate(d["kinetic"]-d["kinetic"].T):
        result[f"kinetic_symmetric_{i}"] = entry
    c = center()
    for name in ("kinetic", "gradient"):
        for i, entry in enumerate(c[name]-c[f"expected_{name}"]):
            result[f"center_{name}_{i}"] = sp.factor(entry)
    lapse = c["c"]
    result["vanishing_center_boundary_coefficient"] = c["Fpi_coefficient_at_center"]
    result["nonvanishing_center_boundary_derivative"] = sp.factor(
        c["Fpi_coefficient_derivative_at_center"]-256*(9*lapse**2-18*lapse+16)/(lapse**2*(lapse-2)**2))
    result["nonuniform_second_log_weight_jet"] = sp.factor(c["weight_second_log_jet"]+22+8*lapse/(lapse-2))
    result["strict_center_superluminal_excess"] = sp.factor(
        c["helicity_speed_squared"]-sp.Rational(17, 15)-(lapse-2)*(9*lapse-4)/120)
    return result
