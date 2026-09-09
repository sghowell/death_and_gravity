"""Actual physical matter density/pressure jets in the natural exponential chart."""

from functools import cache

import sympy as sp

from . import model, source


@cache
def data():
    u = model.u
    a = (1 + u * u) ** 2
    bg = model.parent.coefficients()["background"]
    ell, h = bg["ell"], bg["h"]
    v, chi, _Pv, Pc = model.FIELDS
    c = Pc / a**3
    n1 = sp.factor(
        -model.on_clock()["linear_inhomogeneous_lapse_force"]
        / model.on_clock()["independent_mean_lapse_pivot"]
    )
    dc1 = c - 3 * ell * v
    dc2 = sp.Rational(9, 2) * ell * v * v - 3 * v * c
    common = (
        ell * dc2
        + dc1**2 / 2
        - 3 * ell * dc1 * n1 / h
        + 3 * ell**2 * (1 + 3 * h) * n1 * n1 / (4 * h * h)
    )
    gradient = model.k**2 * chi * chi / a**2
    intrinsic_rho = sp.factor(common + gradient / 2)
    intrinsic_p = sp.factor(common - gradient / 6)
    mean_n, xi = sp.symbols(
        "induced_homogeneous_lapse induced_homogeneous_log_hat_scale", real=True
    )
    mean = -3 * ell * ell * (xi + mean_n / (2 * h))
    rhoM = sp.hessian(intrinsic_rho, model.FIELDS).applyfunc(sp.factor)
    pM = sp.hessian(intrinsic_p, model.FIELDS).applyfunc(sp.factor)
    return {
        "linear_inhomogeneous_lapse": n1,
        "intrinsic_quadratic_matter_density": intrinsic_rho,
        "intrinsic_quadratic_matter_pressure": intrinsic_p,
        "intrinsic_density_Hessian": rhoM,
        "intrinsic_pressure_Hessian": pM,
        "additional_mean_density_and_pressure": mean,
        "mean_lapse_not_counted_again_inside_intrinsic_quadratic_observable": True,
        "checks": {
            "physical_density_pressure_difference_is_actual_spatial_gradient": sp.factor(
                intrinsic_rho - intrinsic_p - sp.Rational(2, 3) * gradient
            ),
        },
    }


@cache
def center():
    d = data()
    u = model.u
    induced = (
        -sp.Rational(3, 200)
        * model.on_clock()["actual_mean_lapse_source"].subs(u, 0)
        * sp.Rational(80, 243)
    )
    total = sp.factor(d["intrinsic_quadratic_matter_density"].subs(u, 0) + induced)
    old = source.density_center()
    v, chi, Pv, Pc = model.FIELDS
    symbols = {
        str(v): v for v in old["natural_exponential_physical_density_jet"].free_symbols
    }
    mapping = {
        symbols["v"]: v,
        symbols["chi"]: chi,
        symbols["p"]: Pv,
        symbols["pc"]: Pc,
        symbols["positive_wave_scale"]: model.k,
    }
    expected = old["natural_exponential_physical_density_jet"].subs(
        mapping, simultaneous=True
    )
    return {
        "anchor_full_physical_density_jet": total,
        "anchor_positive_I4_physical_density_trace": sp.factor(
            sp.trace(sp.hessian(total, model.FIELDS)) / 2
        ),
        "checks": {
            "all_center_physical_density_pairs_match_independent_nonlinear_observable": sp.factor(
                total - expected
            ),
        },
    }
