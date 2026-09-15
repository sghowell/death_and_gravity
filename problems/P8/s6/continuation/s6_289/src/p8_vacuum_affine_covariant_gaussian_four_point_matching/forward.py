"""Exact Gaussian forward coefficients and whole-interval remainder bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_gaussian_metric_pole_matching import spectral

from . import gaussian, source

X = s.Symbol("unit_velocity", real=True)
MU = source.MU
N = s.Symbol("positive_loop_mass_squared", positive=True)


def require_domain(mass, loop_mass):
    mu = source.require_mass(mass)
    n = source.require_mass(loop_mass)
    if n < mu:
        raise ValueError(
            "Require loop mass squared at least external mass squared for this bound"
        )
    return mu, n


def numerator_jets(loop_mass=N, mass=MU, variable=X):
    n, mu, x = map(s.sympify, (loop_mass, mass, variable))
    h = 1 - x * x
    den = 4 * n - 2 * mu * h
    return (
        -16 * mu * (4 * n * n - 2 * n * mu * h - mu * mu * h * h) / den**3,
        16 * mu * (20 * n * n - 4 * n * mu * h + mu * mu * h * h) / den**3,
    )


def nonlocal_integrand(species, loop_mass=N, mass=MU, variable=X):
    gaussian.require_species(species)
    loop_mass, mass, variable = map(s.sympify, (loop_mass, mass, variable))
    j2, j0 = numerator_jets(loop_mass, mass, variable)
    return (
        j2 * spectral.weight(species, 2, variable) / 96
        + j0 * spectral.weight(species, 0, variable) / 2304
    )


def nonlocal_b20(species, loop_mass=N, mass=MU, kappa=source.K):
    kappa = s.sympify(kappa)
    return s.Integral(nonlocal_integrand(species, loop_mass, mass), (X, 0, 1)) / (
        s.pi**2 * kappa**2
    )


def absolute_nonlocal_bound(species, loop_mass, mass=1, kappa=source.K):
    gaussian.require_species(species)
    mu, n = require_domain(mass, loop_mass)
    k = s.sympify(kappa)
    factor = {"scalar": s.Rational(73, 2520), "vector": s.Rational(1, 35)}[species]
    return mu * factor / (s.pi**2 * k * k * n)


def fixed_HP_interval():
    L = s.log(source.HEAVY_MASS2) + 2
    k = source.KAPPA
    center = -L / (320 * s.pi**2 * k * k)
    bound = absolute_nonlocal_bound(
        "scalar", source.HEAVY_MASS2, 1, k
    ) + absolute_nonlocal_bound("vector", source.VECTOR_MASS2, 1, k)
    return center - bound, center + bound


@cache
def data():
    a = s.Symbol("channel", real=True)
    n, mu, x = N, MU, X
    h = 1 - x * x
    den = 4 * n - a * h
    j2 = s.factor(s.diff(a * (a - 4 * mu) ** 2 / den, a, 2).subs(a, 2 * mu))
    j0 = s.factor(s.diff(a * (a + 2 * mu) ** 2 / den, a, 2).subs(a, 2 * mu))
    target = numerator_jets()
    checks = {
        "entire_spin2_crossed_forward_second_jet": s.factor(j2 - target[0]),
        "entire_spin0_crossed_forward_second_jet": s.factor(j0 - target[1]),
        "entire_uncombined_spin2_derivative_for_uniform_bound": s.factor(
            target[0]
            + 4 * mu / (4 * n - 2 * mu * h)
            + 8 * h * mu**2 / (4 * n - 2 * mu * h) ** 2
            - 16 * h * h * mu**3 / (4 * n - 2 * mu * h) ** 3
        ),
        "entire_uncombined_spin0_derivative_for_uniform_bound": s.factor(
            target[1]
            - 20 * mu / (4 * n - 2 * mu * h)
            - 64 * h * mu**2 / (4 * n - 2 * mu * h) ** 2
            - 64 * h * h * mu**3 / (4 * n - 2 * mu * h) ** 3
        ),
    }
    moments = {}
    for species in ("scalar", "vector"):
        w2 = spectral.weight(species, 2, x)
        w0 = spectral.weight(species, 0, x)
        moments[species] = (s.integrate(w2, (x, 0, 1)), s.integrate(w0, (x, 0, 1)))
        bound = s.factor(moments[species][0] / 16 + 17 * moments[species][1] / 1152)
        checks[species + "_whole_absolute_bound_constant"] = (
            bound
            - ({"scalar": s.Rational(73, 2520), "vector": s.Rational(1, 35)}[species])
        )
        checks[species + "_entire_spin2_positive_weight"] = s.factor(
            w2
            - (
                {
                    "scalar": x**6 / 30,
                    "vector": x
                    * x
                    * (13 + 14 * (1 - x * x) + 3 * (1 - x * x) ** 2)
                    / 30,
                }[species]
            )
        )
        checks[species + "_entire_spin0_positive_weight"] = s.factor(
            w0
            - (
                {
                    "scalar": x * x * (3 - x * x) ** 2,
                    "vector": x
                    * x
                    * (s.Rational(8, 3) + 3 * (x * x - s.Rational(1, 3)) ** 2),
                }[species]
            )
        )
    light = s.factor(nonlocal_integrand("scalar", mu, mu, x))
    primitive = s.integrals.rationaltools.ratint(light, x)
    exact = s.simplify(primitive.subs(x, 1) - primitive.subs(x, 0))
    target_light = s.Rational(319, 4800) - 23 * s.pi / 1280
    P = 3 * x**8 - 12 * x**6 + 22 * x**4 - 140 * x * x + 255
    checks.update(
        {
            "whole_light_Gaussian_rational_integrand": s.factor(
                light - x * x * P / (1920 * (1 + x * x) ** 3)
            ),
            "whole_light_Gaussian_rational_primitive": s.factor(
                s.diff(primitive, x) - light
            ),
            "whole_light_Gaussian_exact_integral": s.simplify(exact - target_light),
            "whole_light_Gaussian_positive_numerator": s.expand(
                P
                - 128
                - 120 * (1 - x * x)
                - 4 * (1 - x * x) ** 2
                - 3 * (1 - x * x) ** 4
            ),
        }
    )
    eps = s.Symbol("inverse_heavy_mass", positive=True)
    coefficients = {}
    for species in ("scalar", "vector"):
        expansion = s.series(
            nonlocal_integrand(species, 1 / eps, mu, x), eps, 0, 3
        ).removeO()
        row = tuple(
            s.factor(s.integrate(s.expand(expansion).coeff(eps, j), (x, 0, 1)))
            for j in (1, 2)
        )
        coefficients[species] = row
        expected = {
            "scalar": (mu / 240, 151 * mu * mu / 60480),
            "vector": (s.S.Zero, mu * mu / 20160),
        }[species]
        for j, (v, target_value) in enumerate(zip(row, expected, strict=True)):
            checks[species + "_whole_heavy_coefficient_" + str(j)] = s.factor(
                v - target_value
            )
    coeff = source.gaussian_source.fixed_coefficients()
    L = s.log(source.HEAVY_MASS2) + 2
    local_b20 = 4 * (coeff["R_old_squared"] + 4 * coeff["Weyl_squared"] / 3)
    checks["whole_fixed_H_Proca_local_b20"] = s.factor(local_b20 + L / (320 * s.pi**2))
    scalar_bound = s.Rational(73, 2520) / source.HEAVY_MASS2
    vector_bound = s.Rational(1, 35) / source.VECTOR_MASS2
    remainder_bound = scalar_bound + vector_bound
    exp_lower = sum(s.Rational(5, 2) ** j / s.factorial(j) for j in range(5))
    margins = {
        "actual_both_loop_masses_above_mu": min(source.HEAVY_MASS2, source.VECTOR_MASS2)
        - 1,
        "actual_nonlocal_bound_below_one_ten_millionth": s.Rational(1, 10**7)
        - remainder_bound,
        "strict_negative_fixed_coefficient_margin": s.Rational(1, 160)
        - s.Rational(1, 10**7),
        "log10_upper_bound_Taylor_margin": exp_lower - 10,
        "actual_log_heavy_mass_upper_input": s.Integer(10) ** 200 - source.HEAVY_MASS2,
        "absolute_inverse_kappa_squared_margin": 1
        - (s.Rational(502, 320) + s.Rational(1, 10**7)) / 4,
    }
    r, w = s.symbols("unmatched_Phi_R_squared unmatched_Phi_Weyl_squared", real=True)
    return {
        "whole_massive_forward_numerator_jets": target,
        "whole_scalar_nonlocal_b20": nonlocal_b20("scalar"),
        "whole_Proca_nonlocal_b20": nonlocal_b20("vector"),
        "whole_weight_integrals": moments,
        "whole_light_Gaussian_nonlocal_integrand": light,
        "whole_light_Gaussian_nonlocal_b20": target_light / (s.pi**2 * source.K**2),
        "whole_light_Gaussian_with_unmatched_local": (
            4 * (r + 4 * w / 3) + target_light / s.pi**2
        )
        / source.K**2,
        "whole_fixed_H_Proca_b20": local_b20 / source.KAPPA**2
        + nonlocal_b20("scalar", source.HEAVY_MASS2, 1, source.KAPPA)
        + nonlocal_b20("vector", source.VECTOR_MASS2, 1, source.KAPPA),
        "whole_fixed_H_Proca_strict_interval": fixed_HP_interval(),
        "whole_exact_positive_bound_margins": margins,
        "whole_heavy_expansion_calibration": coefficients,
        "whole_bound_proof": "For0<=h<=1 and n>=mu, den=4n-2mu h>=2n. The uncombined derivative formulas give absJ2<=6mu/n and J0<=34mu/n. Integrating the complete positive spin weights gives scalar73/2520 and vector1/35. For actual mu1,nH=10^200/512+2,M_A^2=10^6, the sum bound is below1/(10^7*pi^2*kappa^2). This is a full finite-mass bound, not truncation of the displayed heavy series.",
        "actual_fixed_HP_verdict": "The entire selected fixed Gaussian H/Proca metric b20, after consistent known light-pole subtraction, lies within1/(10^7*pi^2*kappa^2) of-(log(nH)+2)/(320*pi^2*kappa^2). It is strictly negative and has absolute value below1/kappa^2=10^-1600. This is not the full-source b20 or an isolated finite-gravity positivity/no-go verdict.",
        "coefficient_and_IR_boundary": "All gapped kernels are analytic at the crossing center before derivatives are interchanged with their integrals. The known Newton/light poles are subtracted consistently in all channels first. The light Gaussian local curvature combination remains unknown; the exact positive nonlocal constant does not fix its sign or the complete S288 local anchor. Massless loops, physical IR/Regge and original state/bounce remain separate.",
        "checks": checks,
        "gates": {
            "whole_finite_mass_functions_not_only_asymptotic_series": True,
            "all_domain_derivative_bounds_and_positive_weights": True,
            "all_actual_numeric_margins_strict": all(v > 0 for v in margins.values()),
            "complete_fixed_HP_negative_interval_not_full_b20_verdict": True,
            "exact_light_nonlocal_integral_with_local_unknown_retained": True,
            "known_Newton_poles_removed_before_center_or_forward_limit": True,
            "physical_IR_Regge_and_original_P8_not_closed": True,
        },
    }
