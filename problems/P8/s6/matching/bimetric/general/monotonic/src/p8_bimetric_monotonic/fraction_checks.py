"""Independent exact coefficient identities, using the pinned Fraction ring."""

from p8_bimetric_general.polynomial import Laurent


def derive():
    names = ("a", "b", "Ng", "Nf", "ad", "bd", "G", "F")
    a, b, ng, nf, ad, bd, g, f = (Laurent.variable(name) for name in names)
    y, c = b/a, nf/ng
    h, hf = ad/(ng*a), bd/(nf*b)
    inertia = g+f*y*y
    branch = hf-h/y
    dg = lambda value: value.dt()/ng
    df = lambda value: value.dt()/nf
    full = -2*g*dg(h)-2*f*c*y**3*df(hf)
    reduced = -2*inertia*dg(h)+2*f*y*dg(y)*h
    numerator = ng*bd-nf*ad
    checks = {"full_dynamic_weighted_derivative": -2*g*dg(h)-2*f*c*y**3*df(h/y)-reduced,
              "off_branch_derivative_correction": full-reduced+2*f*c*y**3*df(branch),
              "undivided_branch_error": branch-numerator/(ng*nf*b),
              "inertia_chain_rule": dg(inertia)-2*f*y*dg(y),
              "physical_volume_ratio": c*y**3-nf*b**3/(ng*a**3)}
    controls = {"omitting_ratio_derivative": 2*f*y*dg(y)*h,
                "omitting_c_in_weight": -2*g*dg(h)-2*f*y**3*df(h/y)-reduced}
    return {"checks": checks, "controls": controls}


def checks():
    data = derive()
    if any(value.terms for value in data["checks"].values()):
        raise ValueError("An exact Fraction monotonicity identity failed")
    if any(not value.terms for value in data["controls"].values()):
        raise ValueError("A Fraction omission control failed to fire")
    return {"zero_coefficient_residuals": dict.fromkeys(data["checks"], 0),
            "nonzero_omission_term_counts": {name: len(value.terms) for name, value in data["controls"].items()},
            "arithmetic": "coefficientwise Fraction Laurent algebra, reusing only the pinned polynomial primitive"}
