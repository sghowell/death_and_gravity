"""Candidate covariant continuation and finite on-clock first variations."""
from functools import cache

import sympy as sp
from p8_vector_dimensional import local, metric


@cache
def scalar_curvature_tensor():
    D, H, Hd = local.dimension, local.H, local.Hd
    N, Nd = sp.symbols("N N_dot", real=True)
    acc, expansion2 = (Hd+H**2)/N**2-H*Nd/N**3, H**2/N**2
    R = 2*D*acc+D*(D-1)*expansion2
    Ricci2 = D**2*acc**2+D*(acc+(D-1)*expansion2)**2
    Riemann2 = 4*D*acc**2+2*D*(D-1)*expansion2**2
    # The four-dimensional scalar heat coefficients are continued as
    # covariant invariants in D+1 dimensions, not re-ranked with D.
    scalar_a4 = R**2/72-Ricci2/180+Riemann2/180
    density = 2*N*scalar_a4
    first = sp.diff(density, N).subs({N: 1, Nd: 0})
    rate = sp.diff(density, Nd).subs({N: 1, Nd: 0})
    energy = sp.factor(-first+local.derivative(rate)+D*H*rate)
    pressure = sp.factor(-energy-local.derivative(energy)/(D*H))
    clockR = sp.factor(R.subs({N: 1, Nd: 0}))
    boxR = sp.factor(-local.derivative(local.derivative(clockR))-D*H*local.derivative(clockR))
    return {"R_D": clockR, "box_R_D": boxR, "energy": energy, "pressure": pressure,
            "gradient_normal": sp.factor(energy-boxR/6),
            "gradient_spatial_average": sp.factor(pressure+boxR/6)}


@cache
def mass_pole_extension():
    D, H, Hd, alpha, beta = local.dimension, local.H, local.Hd, local.alpha, local.beta
    d = scalar_curvature_tensor()
    ricci_normal, ricci_spatial = -D*(Hd+H**2), Hd+D*H**2
    C2_normal = -d["R_D"]/12-sp.Rational(5, 6)*ricci_normal
    C2_spatial = d["R_D"]/12-sp.Rational(5, 6)*ricci_spatial
    poles = {0: -sp.Rational(3, 2)*alpha-sp.Rational(9, 2)*beta,
             1: sp.factor(2*(-alpha*C2_normal+3*beta*C2_spatial)),
             2: sp.factor((-alpha*d["gradient_normal"]+3*beta*d["gradient_spatial_average"])/2)}
    return {"mass_poles_in_D": poles,
            "mass_counterterm_evanescent_energy": {n: sp.factor(-2*sp.diff(value, D).subs(D, 3)) for n, value in poles.items()}}


@cache
def finite_coefficients():
    ordinary = metric.variations()
    mass = mass_pole_extension()
    out = {}
    for n, raw in local.integrated().items():
        total_evanescent = ordinary["energy_evanescent_counterterm"][n]+mass["mass_counterterm_evanescent_energy"][n]
        out[n] = {"energy": sp.factor(raw["radial_finite_part"]-total_evanescent),
                  "pressure": metric.ordinary_matching()["ordinary_covariantly_subtracted_local_terms"][n]["pressure"],
                  "total_energy_counterterm_evanescent_coefficient": total_evanescent}
    return out


@cache
def checks():
    D, H, Hd, alpha, beta = local.dimension, local.H, local.Hd, local.alpha, local.beta
    mass = mass_pole_extension()
    finite = finite_coefficients()
    ordinary = {alpha: 0, beta: 0}
    out = {}
    for n, value in mass["mass_poles_in_D"].items():
        actual = local.integrated()[n]["pole"]
        out[str(n)+"_actual_clock_pole_is_retained"] = sp.factor(value.subs(D, 3)-actual+actual.subs(ordinary))
        out[str(n)+"_ordinary_finite_matching_limit"] = sp.factor(finite[n]["energy"].subs(ordinary)
                     -metric.ordinary_matching()["ordinary_covariantly_subtracted_local_terms"][n]["energy"])
    out["frozen_flat_finite_lapse_term_not_shifted"] = sp.factor(finite[0]["energy"]-local.integrated()[0]["radial_finite_part"])
    out["second_order_mass_evanescent_coefficient"] = sp.factor(mass["mass_counterterm_evanescent_energy"][1]
                                    -alpha*(H**2+sp.Rational(8, 3)*Hd)-beta*(3*H**2-2*Hd))
    out["second_order_finite_clock_energy_at_mu_mass"] = sp.factor(finite[1]["energy"].subs(local.ell, 0)
                                    +(10+alpha-3*beta)*H**2+sp.Rational(2, 3)*beta*Hd)
    return out
