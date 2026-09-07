"""All-sign physical Hubble monotonicity from actual source/EH equations.

No cone, vacuum, infinite-time, auxiliary TT inverse or Hubble-division
premise is used. The index set contains exactly the nonzero constant links;
positive Einstein coefficients and regular positive coframes make K>0.
"""

from functools import cache

import sympy as sp

from .obstruction import exact_rational


@cache
def derive(link_count=2):
    if type(link_count) is not int or link_count not in (1, 2):
        raise ValueError("This trimetric screen needs one or two nonzero links")
    h, hp = sp.symbols("H_u H_u_prime", real=True)
    nh = sp.Symbol("n_h", nonnegative=True)
    links = []
    for name in ("g", "f")[:link_count]:
        g, ratio, cone = sp.symbols(f"G_{name} R_{name} c_{name}", positive=True)
        p = sp.Symbol(f"p_{name}", real=True, nonzero=True)
        endpoint, hi, dhi = sp.symbols(f"b_{name} H_{name} D_{name}_H_{name}", real=True)
        density = 2*p*ratio**3+2*endpoint
        pressure = -2*p*ratio**3/cone-2*endpoint
        di_ratio = ratio*((ratio/cone)*h-hi)
        balance = sp.diff(density, ratio)*di_ratio+3*hi*(density+pressure)
        ratio_prime = ratio*(1-cone)*h
        coefficient = g/ratio**2
        raychaudhuri = coefficient*(hp+(1-cone)*h**2)-p/ratio*(1-cone)
        links.append({"G": g, "R": ratio, "c": cone, "p": p, "b": endpoint,
                      "H_i": hi, "D_i_H_i": dhi,
                      "interaction_density": density, "interaction_pressure": pressure,
                      "D_i_R_kinematic": di_ratio, "interaction_balance": sp.factor(balance),
                      "H_i_on_branch": ratio*h, "R_prime": ratio_prime,
                      "K_i": coefficient, "physical_raychaudhuri_residual": raychaudhuri})
    k = sum(link["K_i"] for link in links)
    kp = sum(sp.diff(link["K_i"], link["R"])*link["R_prime"] for link in links)
    unull = sum(link["p"]/link["R"]*(link["c"]-1) for link in links)-nh/2
    combined = 2*k*hp-h*kp+nh
    zprime = hp/sp.sqrt(k)-h*kp/(2*k**sp.Rational(3, 2))
    return {"link_count": link_count, "links": links, "H_u": h, "H_u_prime": hp, "n_h": nh,
            "K": k, "K_prime": sp.factor(kp), "u_weighted_null_residual": unull,
            "combined_residual": sp.factor(combined),
            "normalized_Hubble": h/sp.sqrt(k), "normalized_Hubble_prime": zprime,
            "normalized_equation_residual": sp.factor(zprime+nh/(2*k**sp.Rational(3, 2)))}


@cache
def checks(link_count=2):
    d = derive(link_count)
    h, hp = d["H_u"], d["H_u_prime"]
    result = {}
    for index, link in enumerate(d["links"]):
        g, r, c, p, hi, dhi = (link[key] for key in ("G", "R", "c", "p", "H_i", "D_i_H_i"))
        lapse = 3*g*hi**2-link["interaction_density"]
        space = g*(2*dhi+3*hi**2)+link["interaction_pressure"]
        lapse_derivative = sp.diff(lapse, hi)*dhi+sp.diff(lapse, r)*link["D_i_R_kinematic"]
        balance_expected = 6*p*r**3*(r*h-hi)/c
        # These equations establish Bianchi conservation as a consequence
        # of the Einstein equations, not as a second independent source axiom.
        result[f"source_balance_factorization_{index}"] = sp.factor(link["interaction_balance"]-balance_expected)
        result[f"full_Einstein_Bianchi_residual_{index}"] = sp.factor(lapse_derivative+3*hi*(lapse-space)+link["interaction_balance"])
        result[f"ordered_ratio_derivative_{index}"] = sp.cancel((c/r)*link["D_i_R_kinematic"].subs(hi, r*h)-link["R_prime"])
        di_hi = (r/c)*(link["R_prime"]*h+r*hp)
        old_null = di_hi-p*r**3*(1/c-1)/g
        result[f"proper_to_physical_Raychaudhuri_{index}"] = sp.cancel(g*c*old_null/r**4-link["physical_raychaudhuri_residual"])
    raychaudhuri_sum = sum(link["physical_raychaudhuri_residual"] for link in d["links"])
    result["full_combined_equation_residual"] = sp.factor(2*raychaudhuri_sum-2*d["u_weighted_null_residual"]-d["combined_residual"])
    result["normalized_monotonic_equation"] = sp.simplify(d["normalized_equation_residual"]-d["combined_residual"]/(2*d["K"]**sp.Rational(3, 2)))
    result["regular_zero_Hubble_residual"] = sp.simplify(d["combined_residual"].subs(h, 0)-2*d["K"]*hp-d["n_h"])
    return result


def bounce_residual_budget(k_lower, hubble_derivative_lower, null_lower=0,
                          hubble_abs=0, kprime_abs=0):
    """Necessary actual-equation threshold, not an error estimate.

    K>=k>0,H_u'>=a>0,n_h>=nu>=0,|H_u|<=eta,|K'|<=D imply that
    2K H_u'-H_u K'+n_h is at least 2ka+nu-eta D. Units are physical:
    K mass²,H' mass²,n_h mass4,H mass,K' mass3. No CD dictionary assumed.
    """
    k, acceleration, null, ha, kd = map(exact_rational,
        (k_lower, hubble_derivative_lower, null_lower, hubble_abs, kprime_abs))
    if not (k > 0 and acceleration > 0 and null >= 0 and ha >= 0 and kd >= 0):
        raise ValueError("Require positive K/acceleration and nonnegative remaining bounds")
    margin = 2*k*acceleration+null-ha*kd
    return {"combined_residual_lower_bound": margin,
            "status": "POSITIVE_ACTUAL_EQUATION_RESIDUAL" if margin > 0 else "INCONCLUSIVE",
            "scope": "Necessary cancellation magnitude, not a supplied omitted-operator or loop bound"}


def cd_endpoint_test(length, normalized_hubble_error):
    """Conditional ACTUAL physical H comparison at T=+-L*tau.

    The map must identify this parent's actual h=u^T eta u and its proper
    time T with the prescribed CD metric/time. An error bound e means
    tau*abs(H_parent-H_CD)<=e at BOTH endpoints. No derivative, field-content
    or omitted-operator dictionary is supplied by this endpoint test.
    """
    length, error = map(exact_rational, (length, normalized_hubble_error))
    if not (length > 0 and error >= 0):
        raise ValueError("Require L>0 and normalized Hubble error>=0")
    threshold = 4*length/(1+length**2)
    return {"L": length, "normalized_hubble_error": error,
            "necessary_normalized_endpoint_error": threshold,
            "strict_sign_margin": threshold-error,
            "status": "EXCLUDED_BY_ACTUAL_PARENT_HUBBLE_MONOTONICITY" if error < threshold else "INCONCLUSIVE",
            "dictionary": "Actual h metric and its proper T at both identified endpoints; not an approximate reduced equation alone"}


def controls():
    d = derive()
    return {"positive_bounce_H0_residual": sp.factor(d["combined_residual"].subs({d["H_u"]: 0, d["H_u_prime"]: 1})),
            "exact_positive_margin": bounce_residual_budget(1, 1, 1)["combined_residual_lower_bound"],
            "near_zero_Hubble_positive_margin": bounce_residual_budget(1, 1, 1, sp.Rational(1, 10), 5)["combined_residual_lower_bound"],
            "uncontrolled_connection_is_inconclusive": bounce_residual_budget(1, 1, 1, 1, 5)["combined_residual_lower_bound"],
            "degenerate_crossing_required_negative_null": -6*sp.Symbol("T", real=True)**2,
            "CD_half_duration_endpoint_threshold": cd_endpoint_test(sp.Rational(1, 2), 0)["necessary_normalized_endpoint_error"],
            "CD_unit_duration_endpoint_threshold": cd_endpoint_test(1, 0)["necessary_normalized_endpoint_error"]}
