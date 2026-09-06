"""Exact limiting forward-time flow and zero-charge homogeneous comparison."""

from functools import cache

import sympy as sp

from . import scaled


@cache
def general_jets():
    v, a = sp.symbols("v a", positive=True)
    vp, ap = 1+a-v**2, -a*v

    def derivative(value):
        return sp.diff(value, v)*vp+sp.diff(value, a)*ap

    tr = (3*v**2-1-a)/(2*v)
    nr = sp.factor(derivative(tr)+tr**2)
    mass = v**2-v
    return {"v": v, "a": a, "vprime": vp, "aprime": ap,
            "vsecond": derivative(vp), "theta_relative": tr, "N_relative": nr,
            "mass_squared": mass, "mass_squared_prime": derivative(mass),
            "mass_squared_second": sp.factor(derivative(derivative(mass))),
            "heavy_growth_coefficient": sp.factor(nr-mass),
            "general_hierarchy_margin": sp.factor(nr-sp.Rational(3, 4)*mass)}


@cache
def solution():
    u = sp.Symbol("u", nonnegative=True)
    z = 12*sp.cosh(u)+2*sp.sinh(u)-11
    zp = sp.diff(z, u)
    v, a, j = zp/z, 11/z, 11*u/z
    rr = v-j
    return {"u": u, "z": z, "zprime": zp, "v": v, "a": a, "j": j,
            "R": rr, "eta": 3*(rr**2-1),
            "positive_R_numerator": zp-11*u}


def checks():
    d = general_jets()
    v, a = (d[key] for key in ("v", "a"))
    nr = (3*(a+1)**2+3*v**4-2*v**2)/(4*v**2)
    data = solution()
    u, z, zp, vv, aa, jj, rr = (data[key] for key in ("u", "z", "zprime", "v", "a", "j", "R"))
    zz = sp.Symbol("z", positive=True)
    energy_v2 = 1+22/zz-19/zz**2
    return {
        "general_relative_normalization_curvature": sp.factor(d["N_relative"]-nr),
        "general_mass_led_positive_margin": sp.factor(d["general_hierarchy_margin"]
             -(3*(a+1)**2+v**2*(3*v-2))/(4*v**2)),
        "full_time_heavy_coefficient": sp.factor(d["heavy_growth_coefficient"]
             -3*(a+1)**2/(4*v**2)-(-v**2+4*v-2)/4),
        "tuned_second_derivative_survives": sp.factor(nr.subs(a, 3*v**2-1)-(15*v**2-1)/2),
        "integrable_scale_equation": sp.simplify(sp.diff(z, u, 2)-z-11),
        "integrable_scale_initial": z.subs(u, 0)-1,
        "integrable_scale_velocity_initial": zp.subs(u, 0)-2,
        "exact_scale_first_integral": sp.trigsimp(zp**2-z**2-22*z+19),
        "full_time_v_flow": sp.simplify(sp.diff(vv, u)-1-aa+vv**2),
        "full_time_a_flow": sp.simplify(sp.diff(aa, u)+aa*vv),
        "full_time_j_flow": sp.simplify(sp.diff(jj, u)-aa+jj*vv),
        "full_time_R_flow": sp.simplify(sp.diff(rr, u)-1+rr*vv),
        "positive_R_numerator_derivative": sp.simplify(sp.diff(data["positive_R_numerator"], u)-z),
        "density_sign_crossing_identity": sp.simplify((rr-1-(11*(1-u)-10*sp.exp(-u))/z).rewrite(sp.exp)),
        "uniform_v_squared_upper_square": sp.factor(sp.Rational(140, 19)-energy_v2-(11*zz-19)**2/(19*zz**2)),
        "uniform_v_squared_lower": sp.factor(energy_v2-1-(22*zz-19)/zz**2),
        "uniform_heavy_coefficient_quarter_margin": sp.expand((-v**2+4*v-2)-1-(v-1)*(3-v)),
        "comparison_cosine_hyperbolic_solution": sp.diff(sp.cosh(u/2), u, 2)-sp.cosh(u/2)/4,
    }


def canonical_limit_checks():
    d = scaled.canonical_jets()
    g = general_jets()
    v, a, vp, vpp = (g[key] for key in ("v", "a", "vprime", "vsecond"))
    e = d["e"]
    base = {d["kappa"]: 1/v, d["kappa_prime"]: -vp/v**2,
            d["kappa_second"]: 2*vp**2/v**3-vpp/v**2,
            d["lambda"]: -v, d["lambda_prime"]: -vp, d["h"]: 0, d["hprime"]: 0}
    targets = {"theta_sum": 0, "N_sum": 0, "omega": 0,
               "theta_relative": g["theta_relative"], "N_relative": g["N_relative"],
               "mass_squared": g["mass_squared"], "mass_squared_prime": g["mass_squared_prime"],
               "mass_squared_second": g["mass_squared_second"],
               "physical_source_relative": -sp.sqrt(v)}
    checks = {name: sp.simplify(d[name].subs(base, simultaneous=True).subs(e, 0)-target)
              for name, target in targets.items()}
    # First nonzero light/mixing orders at the centre: h=0, hprime=a*e.
    centre = dict(base)
    centre[d["hprime"]] = a*e
    checks["theta_sum_over_e"] = sp.simplify(
        sp.limit(d["theta_sum"].subs(centre, simultaneous=True)/e, e, 0)-(1+a-3*v)/2)
    checks["N_sum_over_e"] = sp.simplify(
        sp.limit(d["N_sum"].subs(centre, simultaneous=True)/e, e, 0)-(v-3)/2)
    checks["omega_over_e"] = sp.simplify(
        sp.limit(d["omega"].subs(centre, simultaneous=True)/e, e, 0)+g["theta_relative"]/sp.sqrt(v))
    checks["omega_prime_over_e"] = sp.simplify(
        sp.limit(d["omega_prime"].subs(centre, simultaneous=True)/e, e, 0)+g["N_relative"]/sp.sqrt(v))
    return checks


def physical_projection_checks():
    d = general_jets()
    v, a = d["v"], d["a"]
    at = {v: 2, a: 11}
    vp, vpp = (d[key].subs(at) for key in ("vprime", "vsecond"))
    growth = d["heavy_growth_coefficient"].subs(at)
    # gamma_eff=-(2/M)*sqrt(v)*H in the limiting zero-common-field gauge.
    # Report derivatives of M*gamma_eff/2 to avoid an arbitrary Planck scale.
    first = -vp/(2*sp.sqrt(2))
    second = -vpp/(2*sp.sqrt(2))+vp**2/(4*2**sp.Rational(3, 2))-sp.sqrt(2)*growth
    return {"initial_heavy_growth_coefficient": growth-sp.Rational(55, 2),
            "initial_physical_tensor_time_dependence": sp.simplify(first+2*sp.sqrt(2)),
            "initial_physical_tensor_tidal_derivative": sp.simplify(second+10*sp.sqrt(2))}
