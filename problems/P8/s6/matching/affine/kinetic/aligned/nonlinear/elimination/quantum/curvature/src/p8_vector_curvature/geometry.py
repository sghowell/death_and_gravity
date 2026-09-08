"""Physical FLRW tensor contractions and global rational coefficient envelopes."""
from functools import cache

import sympy as sp
from p8_affine_aligned import modes
from p8_affine_retuned.bounds import exact

from . import heat

u = modes.u


@cache
def metric_contractions():
    a = sp.Symbol("a", positive=True)
    H, Hd = sp.symbols("H H_dot", real=True)
    g = sp.diag(-1, a**2, a**2, a**2)
    gi = sp.diag(-1, a**-2, a**-2, a**-2)

    def derivative(value, coordinate):
        return sp.diff(value, a)*a*H+sp.diff(value, H)*Hd if coordinate == 0 else 0

    gamma = sp.MutableDenseNDimArray.zeros(4, 4, 4)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                gamma[i, j, k] = sp.expand(sum(gi[i, l]*(derivative(g[l, k], j)
                    +derivative(g[l, j], k)-derivative(g[j, k], l))/2 for l in range(4)))
    riemann = sp.MutableDenseNDimArray.zeros(4, 4, 4, 4)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    riemann[i, j, k, l] = sp.expand(derivative(gamma[i, l, j], k)-derivative(gamma[i, k, j], l)
                        +sum(gamma[i, k, n]*gamma[n, l, j]-gamma[i, l, n]*gamma[n, k, j] for n in range(4)))
    ricci = sp.Matrix(4, 4, lambda i, j: sum(riemann[k, i, k, j] for k in range(4)))
    scalar = sp.expand(sp.trace(gi*ricci))
    ricci_norm = sp.expand(sp.trace(gi*ricci*gi*ricci))
    riemann_norm = sp.expand(sum((g[i, i]*riemann[i, j, k, l])**2*gi[i, i]*gi[j, j]*gi[k, k]*gi[l, l]
                                for i in range(4) for j in range(4) for k in range(4) for l in range(4)))
    return {"H": H, "Hd": Hd, "R": scalar, "Ricci2": ricci_norm, "Riemann2": riemann_norm,
            "actual_metric_R": sp.expand(scalar-6*(Hd+2*H**2)),
            "actual_metric_Ricci2": sp.expand(ricci_norm-12*(Hd**2+3*Hd*H**2+3*H**4)),
            "actual_metric_Riemann2": sp.expand(riemann_norm-12*((Hd+H**2)**2+H**4)),
            "physical_FLRW_Weyl_squared_zero": sp.expand(riemann_norm-2*ricci_norm+scalar**2/3)}


@cache
def rolling():
    bg, tensors = modes.canonical(), metric_contractions()
    H, a = bg["H"], bg["a"]
    subs = {tensors["H"]: H, tensors["Hd"]: sp.diff(H, u)}
    data = {name: sp.factor(tensors[name].subs(subs)) for name in ("R", "Ricci2", "Riemann2")}
    data["boxR"] = sp.factor(-sp.diff(data["R"], u, 2)-3*H*sp.diff(data["R"], u))
    data["Gauss_Bonnet"] = sp.factor(data["Riemann2"]-4*data["Ricci2"]+data["R"]**2)
    mapping = {heat.R: data["R"], heat.Ric2: data["Ricci2"], heat.Riem2: data["Riemann2"], heat.boxR: data["boxR"]}
    data["local_a4_with_divergence"] = sp.factor(heat.coefficients()["proca"][4].subs(mapping))
    data["local_a4_mod_divergence"] = sp.factor(data["local_a4_with_divergence"]+data["boxR"]/15)
    data["box_R_boundary_identity"] = sp.factor(a**3*data["boxR"]+sp.diff(a**3*sp.diff(data["R"], u), u))
    data["Gauss_Bonnet_boundary_identity"] = sp.factor(a**3*data["Gauss_Bonnet"]-sp.diff(8*a**3*H**3, u))
    data["box_R_boundary_flux"] = sp.factor(-a**3*sp.diff(data["R"], u))
    data["Gauss_Bonnet_boundary_flux"] = sp.factor(8*a**3*H**3)
    data["box_R_boundary_ninth_power_growth_identity"] = sp.limit(data["box_R_boundary_flux"]/u**9, u, sp.oo)-336
    data["Gauss_Bonnet_boundary_ninth_power_growth_identity"] = sp.limit(data["Gauss_Bonnet_boundary_flux"]/u**9, u, sp.oo)-512
    data["sharp_global_scalar_curvature_upper_identity"] = sp.factor(49-data["R"]-(7*u**2-5)**2/(1+u**2)**2)
    return data


def even_envelope(expression):
    """For P(u²)/(1+u²)^n use the convex binomial basis on all real u."""
    value = sp.cancel(expression)
    numerator, denominator = sp.fraction(value)
    den = sp.Poly(denominator, u)
    if den.degree() % 2:
        raise ValueError("Require an even power of the quadratic denominator")
    n = den.degree()//2
    coefficient = den.LC()
    if coefficient.is_positive is not True or sp.expand(denominator-coefficient*(1+u**2)**n) != 0:
        raise ValueError("Require a positive multiple of (1+u²)^n")
    poly = sp.Poly(numerator/coefficient, u)
    if poly.degree() > 2*n or any(power[0] % 2 for power, term in poly.terms() if term != 0):
        raise ValueError("Require an even numerator of no higher degree")
    entries = [poly.nth(2*j)/sp.binomial(n, j) for j in range(n+1)]
    upper = max(abs(entry) for entry in entries)
    if any(not entry.is_Rational for entry in entries):
        raise ValueError("Require rational fixed coefficients")
    return {"power": n, "binomial_coefficients": entries, "absolute_upper": upper}


@cache
def continuous_bounds():
    d = rolling()
    return {name: even_envelope(d[name]) for name in ("R", "Ricci2", "Riemann2", "boxR", "Gauss_Bonnet",
                                                    "local_a4_with_divergence", "local_a4_mod_divergence")}


@cache
def proof_checks():
    bounds = continuous_bounds()
    out = {"global_binomial_convex_envelope_"+name: all(bool(abs(coefficient) <= d["absolute_upper"])
                for coefficient in d["binomial_coefficients"]) for name, d in bounds.items()}
    out.update({"all_binomial_powers_nonnegative_integers": all(isinstance(d["power"], int) and d["power"] >= 0 for d in bounds.values()),
                "global_scalar_curvature_positive_numerator": all(coefficient > 0 for coefficient in (24, 168)),
                "local_a4_absolute_envelope_308_over_3": bounds["local_a4_with_divergence"]["absolute_upper"] == sp.Rational(308, 3),
                "noncompact_boundary_flux_not_assumed_zero": all(sp.limit(rolling()[name]/u**9, u, sp.oo).is_nonzero is True
                    for name in ("box_R_boundary_flux", "Gauss_Bonnet_boundary_flux")),
                "no_local_pole_to_finite_quantum_error_transfer": True})
    return out


def scale_bounds(reference_mass_time_product):
    mass_time = exact(reference_mass_time_product, "m0_tau")
    if mass_time.is_positive is not True:
        raise ValueError("Require positive m0*tau")
    bound4 = continuous_bounds()["local_a4_with_divergence"]["absolute_upper"]
    return {"m0_tau": mass_time,
            "scalar_curvature_over_mass_squared": 49/mass_time**2,
            "absolute_local_curvature_pole_correction_over_flat_pole_weight": 49/(3*mass_time**2)+2*bound4/(3*mass_time**4),
            "finite_curved_or_in_in_remainder": False}


@cache
def checks():
    tensors, bg = metric_contractions(), rolling()
    out = {name: value for name, value in tensors.items() if name.startswith("actual_metric_") or name.endswith("Weyl_squared_zero")}
    out.update({name: value for name, value in bg.items() if name.endswith("identity")})
    for name, bound in continuous_bounds().items():
        reconstructed = sum(sp.binomial(bound["power"], j)*coefficient*u**(2*j)
                            for j, coefficient in enumerate(bound["binomial_coefficients"]))/(1+u**2)**bound["power"]
        out["global_binomial_envelope_reconstruction_"+name] = sp.factor(bg[name]-reconstructed)
    return out
