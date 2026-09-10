"""Global positive Euclidean multiplier of the actual on-shell self-energy."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_light_pole import kernel as parent_kernel
from p8_vacuum_spin2 import triangles


@cache
def data():
    old = parent_kernel.data()
    p = model.data()
    x = old["Feynman_parameter"]
    b = old["normalized_parameter_weight"]
    pref = old["self_energy_prefactor"]
    y = sp.Symbol("nonnegative_Euclidean_momentum_squared", nonnegative=True)
    u, t = sp.symbols("positive_increment unit_auxiliary_parameter", positive=True)
    f = 1 - sp.log(1 + u) / u
    primitive = t - sp.log(1 + u * t) / u
    integrand = u * t / (1 + u * t)
    moment = b * (1 - sp.log(1 + b * (y + 1)) / (b * (y + 1)))
    ratio = pref * sp.Integral(moment, (x, 0, 1))
    alpha = pref * sp.Integral(b, (x, 0, 1))
    parent_euclidean = old["local_analytic_on_shell_kernel"].subs(
        old["Minkowski_invariant"], -y
    )
    tri = triangles.data()
    tri0 = (tri["light_stress_triangle"] + tri["heavy_stress_triangle"]).subs(
        tri["t"], 0
    )
    tri0 = sp.integrate(tri0, (tri["z"], 0, 1))
    tri0 = (tri["g"] * tri0 / (16 * sp.pi**2)).subs(
        {
            tri["M"]: p["heavy_mass_squared"],
            tri["g"]: p["cubic_coupling_squared"],
            tri["x"]: x,
        }
    )
    a = sp.Symbol("relative_uniform_self_energy_bound", positive=True)
    q = sp.Symbol("actual_relative_self_energy", nonnegative=True)
    return {
        "x": x,
        "y": y,
        "b": b,
        "prefactor": pref,
        "universal_positive_function": f,
        "positive_integral_representation": sp.Integral(integrand, (t, 0, 1)),
        "relative_Euclidean_self_energy": ratio,
        "sharp_uniform_supremum": alpha,
        "same_light_and_stress_normalization": {
            "parent_Pi_prime_at_one": pref * sp.Integral(b, (x, 0, 1)),
            "parent_unsubtracted_stress_charge_integrand": tri0,
            "once_fixed_kinetic_counterterm": -alpha,
            "renormalized_stress_charge_not_alpha": sp.Integer(1),
        },
        "displayed_Euclidean_inverse": (y + 1) * (1 - q),
        "scope": "Real Euclidean momenta at the fixed one-loop on-shell subtraction. Positive quadratic-form bounds and repetitions of this computed kernel are not timelike analytic-contour, reflection-positivity, full higher-loop or UV-completion theorems.",
        "checks": {
            "actual_parent_analytic_kernel_Euclidean_continuation": sp.simplify(
                parent_euclidean - (b * (y + 1) - sp.log(1 + b * (y + 1)))
            ),
            "relative_kernel_factorization": sp.simplify(
                moment - parent_euclidean / (y + 1)
            ),
            "positive_integral_primitive_derivative": sp.factor(
                sp.diff(primitive, t) - integrand
            ),
            "positive_integral_primitive_anchor": primitive.subs(t, 0),
            "positive_integral_endpoint_equals_universal_function": sp.simplify(
                primitive.subs(t, 1) - f
            ),
            "universal_continuous_zero_limit": sp.limit(f, u, 0, dir="+"),
            "universal_large_increment_limit": sp.limit(f, u, sp.oo) - 1,
            "strict_monotonicity_log_gap_derivative": sp.factor(
                sp.diff(sp.log(1 + u) - u / (1 + u), u) - u / (1 + u) ** 2
            ),
            "strict_monotonicity_log_gap_anchor": (sp.log(1 + u) - u / (1 + u)).subs(
                u, 0
            ),
            "universal_function_derivative": sp.factor(
                sp.diff(f, u) - (sp.log(1 + u) - u / (1 + u)) / u**2
            ),
            "positive_lower_integrand_gap": sp.factor(
                integrand
                - u * t / (1 + u)
                - u * u * t * (1 - t) / ((1 + u) * (1 + u * t))
            ),
            "linear_upper_integrand_gap": sp.factor(
                u * t - integrand - u * u * t * t / (1 + u * t)
            ),
            "unit_upper_gap": sp.simplify(1 - f - sp.log(1 + u) / u),
            "actual_light_and_stress_charge_same_integrand": sp.factor(tri0 - pref * b),
            "same_parent_sharp_alpha_majorant": sp.factor(
                sp.integrate((1 - x) / p["heavy_mass_squared"], (x, 0, 1))
                - 1 / (2 * p["heavy_mass_squared"])
            ),
            "Euclidean_inverse_lower_gap": sp.expand(
                (y + 1) * (1 - q) - (1 - a) * (y + 1) - (a - q) * (y + 1)
            ),
            "Euclidean_covariance_upper_gap": sp.factor(
                1 / (1 - a) - 1 / (1 - q) - (a - q) / ((1 - a) * (1 - q))
            ),
        },
    }
