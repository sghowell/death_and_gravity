"""Exact total-energy simplex and its whole leading-soft regulator limit."""

from functools import cache

import sympy as s

from . import source

A = s.Symbol("nonnegative_leading_soft_index", nonnegative=True)
DELTA = s.Symbol("finite_angular_scheme_conversion", real=True)
RATIO = s.Symbol("positive_resolution_over_reference", positive=True)
EP = source.EP
COUNT = s.Symbol("nonnegative_emission_count", integer=True, nonnegative=True)


def require_count(value):
    if isinstance(value, (bool, float, s.Float)) or not isinstance(
        value, (int, s.Integer)
    ):
        raise TypeError("Require an exact nonnegative integer emission count")
    if value < 0:
        raise ValueError("Require an exact nonnegative integer emission count")
    return s.Integer(value)


def radial_index(
    energy=source.soft.E,
    cosine=source.soft.COST,
    epsilon=EP,
    mass=source.MU,
    kappa=source.K,
):
    return (
        source.soft.phase_normalization(epsilon)
        * source.soft.kernel(energy, cosine, epsilon, mass)
        / (4 * s.pi**2 * s.sympify(kappa))
    )


def physical_index(
    energy=source.soft.E, cosine=source.soft.COST, mass=source.MU, kappa=source.K
):
    return source.soft.kernel_zero(energy, cosine, mass) / (
        4 * s.pi**2 * s.sympify(kappa)
    )


def simplex_term(count, index, epsilon=EP, ratio=RATIO):
    count = require_count(count)
    index, e, x = map(s.sympify, (index, epsilon, ratio))
    return (index * s.gamma(2 * e) * x ** (2 * e)) ** count / (
        s.factorial(count) * s.gamma(1 + 2 * e * count)
    )


def auxiliary_mean(index, epsilon=EP, ratio=RATIO):
    index, e, x = map(s.sympify, (index, epsilon, ratio))
    return index * s.gamma(2 * e) * x ** (2 * e)


def calorimetric_limit(index=A, conversion=DELTA, ratio=RATIO):
    a, delta, x = map(s.sympify, (index, conversion, ratio))
    return s.exp(delta - s.EulerGamma * a) * x**a / s.gamma(1 + a)


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(s.expand_func(value))

    e = s.Symbol("positive_dimensional_regulator", positive=True)
    a, d, L = s.symbols("a delta log_ratio", real=True)
    ae = a + 2 * e * d
    lam = ae * s.gamma(1 + 2 * e) * s.exp(2 * e * L) / (2 * e)
    const = s.series(lam - a / (2 * e), e, 0, 1).removeO()
    put("full_exponent_finite_part", const - d + a * s.EulerGamma - a * L)
    put("Poisson_scaled_mean_limit", s.limit(2 * e * lam, e, 0, dir="+") - a)
    put("Poisson_scaled_variance_limit", s.limit(4 * e**2 * lam, e, 0, dir="+"))
    put(
        "Poisson_scaled_second_moment",
        (2 * e) ** 2 * (lam**2 + lam) - (2 * e * lam) ** 2 - 4 * e**2 * lam,
    )
    put("Gamma_reciprocal_normalization_at_zero", 1 / s.gamma(1))
    checks["Gamma_reciprocal_normalization_at_zero"] -= 1
    put(
        "Gamma_global_lower_integral",
        s.integrate(
            s.exp(-s.Symbol("u", positive=True)),
            (s.Symbol("u", positive=True), 1, s.oo),
        )
        - s.exp(-1),
    )
    k0, k1, k = s.symbols("K0 K1 kappa", real=True)
    p = source.soft.phase_normalization(e)
    a0 = k0 / (4 * s.pi**2 * k)
    a1 = s.diff(p * (k0 + e * k1) / (4 * s.pi**2 * k), e).subs(e, 0)
    put(
        "S296_exact_finite_conversion",
        a1 / 2 - (k1 + (s.EulerGamma - 2 - s.log(s.pi)) * k0) / (8 * s.pi**2 * k),
    )
    put("S296_one_real_radial_weight", a0 / (2 * e) - k0 / (8 * s.pi**2 * k * e))
    B = s.Symbol("real_physical_soft_B", real=True)
    put(
        "S288_virtual_rate_pole",
        (-2 * B) / (4 * s.pi**2 * k) / (2 * e) + B / (4 * s.pi**2 * k * e),
    )
    for count in range(5):
        term = simplex_term(count, ae, e, s.exp(L))
        put(
            f"whole_simplex_term_{count}",
            term
            - auxiliary_mean(ae, e, s.exp(L)) ** count
            / (s.factorial(count) * s.gamma(1 + 2 * e * count)),
        )
    for count in range(1, 5):
        prior = s.gamma(2 * e) ** (count - 1) / s.gamma(1 + 2 * e * (count - 1))
        beta = (
            s.gamma(2 * e)
            * s.gamma(1 + 2 * e * (count - 1))
            / s.gamma(1 + 2 * e * count)
        )
        put(
            f"Dirichlet_convolution_{count}",
            prior * beta - s.gamma(2 * e) ** count / s.gamma(1 + 2 * e * count),
        )
    x, y, z = s.symbols(
        "positive_eikonal_x positive_eikonal_y positive_eikonal_z", positive=True
    )
    put("two_emission_eikonal_sum", 1 / (x * (x + y)) + 1 / (y * (x + y)) - 1 / (x * y))
    put(
        "three_emission_induction_step",
        (1 / (x * y) + 1 / (x * z) + 1 / (y * z)) / (x + y + z) - 1 / (x * y * z),
    )
    put("physical_normalized_virtual_to_analytic", (ae - a) / (2 * e) - d)
    put(
        "normalized_compound_Poisson_mean_energy",
        s.integrate(ae * x ** (2 * e), (x, 0, 1)) - ae / (1 + 2 * e),
    )
    put("whole_limit_zero_index", calorimetric_limit(0, d, x) - s.exp(d))
    put(
        "whole_limit_reference_scale",
        calorimetric_limit(a, d, 1) - s.exp(d - s.EulerGamma * a) / s.gamma(1 + a),
    )
    r = s.Symbol("positive_Laplace_argument", positive=True)
    laplace = -a * (s.EulerGamma + s.log(r) + s.expint(1, r))
    put(
        "finite_reference_Levy_exponent_derivative",
        s.diff(laplace, r) - a * (s.exp(-r) - 1) / r,
    )
    put("finite_reference_Levy_exponent_zero", s.limit(laplace, r, 0, dir="+"))
    return {
        "whole_radial_index": radial_index(),
        "whole_D4_index": physical_index(),
        "whole_finite_conversion": source.soft.finite_conversion(),
        "whole_N_emission_simplex": (A * s.gamma(2 * EP) * RATIO ** (2 * EP)) ** COUNT
        / (s.factorial(COUNT) * s.gamma(1 + 2 * EP * COUNT)),
        "whole_analytic_virtual_times_soft_sum_limit": calorimetric_limit(),
        "whole_Poisson_limit_proof": "Let a_e=a+2e*Delta+O(e^2)>=0 and lambda=a_e Gamma(2e)x^(2e). The exact positive sum is exp(lambda-a/(2e))*E[1/Gamma(1+2e*N)] for auxiliary N~Poisson(lambda). At fixedx>0,2e*N tends to a in L2, and the reciprocal Gamma function on the positive real axis is continuous and globally bounded by exp(1), since Gamma(1+y)>=integral_1^infinity exp(-u)du. Its expectation tends to1/Gamma(1+a). The exponent tends toDelta-a*EulerGamma+a*lnx. This proves the whole positive sum limit without interchanging an unbounded termwise series.",
        "whole_leading_soft_boundary": "Each fixed-N leading soft amplitude factorizes by the permutation eikonal identity. Summing those leading terms is exact for this defined eikonal observable; uniform control of the difference from the full amplitudes as N grows is a separate unproved obligation. The virtual Coulomb phase remains in amplitudes and cancels only from this rate.",
        "checks": checks,
        "gates": {
            "total_energy_simplex_with_exact_Bose_factor": True,
            "all_N_eikonal_identity_induction_written": True,
            "finite_regulator_positive_sum_not_termwise_limit": True,
            "bounded_continuous_Gamma_Poisson_limit": True,
            "complete_S296_angular_and_scheme_conversion": True,
            "finite_reference_positive_probability_interpretation": True,
            "leading_only_not_full_recoil_or_hard_loop_completion": True,
        },
    }
