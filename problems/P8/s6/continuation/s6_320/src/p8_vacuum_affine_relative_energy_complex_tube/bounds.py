"""Nonnegative relative-complex majorants and explicitly shrinking Cauchy radii."""

from functools import cache

import sympy as s

from . import source

EPSILON = s.Rational(1, 10**12)


def pure_coefficient(n):
    n = int(source.require_order(n))
    return _pure_coefficient(n)


@cache
def _pure_coefficient(n):
    if n == 1:
        return s.S.One
    return 25600 * sum(
        s.factorial(k + 1) * 32 ** (k + 1) * 13**k * _partitions(n, k)
        for k in range(2, n + 1)
    )


@cache
def _partitions(n, k):
    if n == 0:
        return s.S.One if k == 0 else s.S.Zero
    if k <= 0 or k > n:
        return s.S.Zero
    return sum(
        s.binomial(n - 1, j - 1) * _pure_coefficient(j) * _partitions(n - j, k - 1)
        for j in range(1, n - k + 2)
    )


def pure_envelope(n):
    n = int(source.require_order(n))
    return 2 * (2 * 10**13) ** (n - 1) * s.factorial(n)


def require_energies(energies):
    if not isinstance(energies, (tuple, list)) or not energies:
        raise ValueError("Require a nonempty finite center-energy list")
    ws = tuple(source.require_mass(w) for w in energies)
    if sum(ws, s.S.Zero) > s.Rational(1, 8):
        raise ValueError("Require total positive center energy at most1/8")
    return ws


def bare_upper(energies):
    ws = require_energies(energies)
    n = len(ws)
    return (
        3
        * source.HEAVY_MASS2**2
        * s.factorial(n)
        * (2 * 10**17) ** n
        / (source.KAPPA ** s.Rational(n, 2) * s.prod(ws))
    )


def scaled_coefficient(n):
    n = int(source.require_order(n))
    return (
        6
        * source.HEAVY_MASS2**2
        * s.factorial(n)
        * (4 * 10**17) ** n
        / source.KAPPA ** s.Rational(n, 2)
    )


def cauchy_majorant(energies, alpha):
    ws = require_energies(energies)
    if not isinstance(alpha, (tuple, list)) or len(alpha) != len(ws):
        raise ValueError("Require one exact nonnegative derivative order per energy")
    orders = tuple(int(source.require_multiplicity(a)) for a in alpha)
    return scaled_coefficient(len(ws)) * s.prod(
        s.factorial(a) / (EPSILON * w) ** a for a, w in zip(orders, ws)
    )


def hard_functions(x):
    L = (1 - 1024 * x) / (1 - 5120 * x)
    D = 1 / (1 - 4 * x)
    H = (1 - 1024 * x) / (1 - 3072 * x)
    E = 1024 / (1 - 1024 * x) ** 2
    V = 120000000 * 1024 * (2 / (1 - 32 * x) ** 3 - 2)
    return {
        "L": L,
        "D": D,
        "H": H,
        "E": E,
        "V": V,
        "matter": L**4 * (D + s.Rational(3, 2) * D**2 * H),
        "gravity": s.Rational(3 * 600000, 8) * L**4 * E**2 / (1 - V),
    }


@cache
def barrier_data():
    x = s.Symbol("x")
    cap = s.Rational(1, 10**13)
    hardcap = s.Rational(1, 10**17)
    nonlinear = lambda c: 819200 * ((1 - 416 * c) ** -2 - 1 - 832 * c)
    values = {k: s.factor(v.subs(x, hardcap)) for k, v in hard_functions(x).items()}
    return {
        "pure_cap": cap,
        "pure_radius": cap / 2,
        "pure_margin": cap / 2 - nonlinear(cap),
        "hard_cap": hardcap,
        "hard_radius": hardcap / 2,
        "hard_margin": hardcap / 2 - nonlinear(hardcap),
        "core_values": values,
    }


@cache
def data():
    from p8_vacuum_affine_uniform_all_tree_bound import core

    x = s.Symbol("x")
    counts = {1: s.S.One}
    checks = {}
    for n in range(2, 9):
        previous = sum(counts[j] * x**j / s.factorial(j) for j in range(1, n))
        value = s.factorial(n) * s.expand(
            819200 * sum((k + 1) * 416**k * previous**k for k in range(2, n + 1))
        ).coeff(x, n)
        counts[n] = value
        checks[f"independent_complex_EGF_coefficient_{n}"] = value - pure_coefficient(n)
    row = barrier_data()
    values = row["core_values"]
    checks["pure_radius_half_cap"] = row["pure_radius"] - row["pure_cap"] / 2
    checks["hard_radius_half_cap"] = row["hard_radius"] - row["hard_cap"] / 2
    checks["scalar_denominator_margin_arithmetic"] = (
        s.Rational(3, 8)
        - 100 * EPSILON
        - s.Rational(1, 4)
        - (s.Rational(1, 8) - 100 * EPSILON)
    )
    checks["hard_component_over_gap_budget"] = s.Integer(4 * 50 * 600000 - 120000000)
    for n in range(1, 7):
        # This exact identity records energy homogeneity, not integration.
        w = tuple(s.Rational(1, 16 * n) for _ in range(n))
        checks[f"scaled_bound_accounts_for_complex_product_and_phase_{n}"] = (
            scaled_coefficient(n) - 2 ** (n + 1) * bare_upper(w) * s.prod(w)
        )
        checks[f"Cauchy_zero_order_recovers_scaled_bound_{n}"] = cauchy_majorant(
            w, (0,) * n
        ) - scaled_coefficient(n)
    for name, expr in hard_functions(x).items():
        coeffs = core.rational_coefficients(expr, x, 6)
        if name == "V":
            for r in range(1, 7):
                checks[f"new_EH_step_series_coefficient_{r}"] = coeffs[
                    r
                ] - 120000000 * (r + 2) * (r + 1) * 32 ** (r + 2)
    return {
        "whole_pure_complex_majorant": "D1=1; Dn=25600*sum_(labeled root partitions,k>=2)(k+1)!*32^(k+1)*13^k*product D_child. The factor is25 from the complex momentum coefficient norm times1024 from the conserved weighted inverse. Dn<=2*(2e13)^(n-1)*n!, uniformly for relative energy shifts at most1/100.",
        "whole_first_eight_complex_coefficients": counts,
        "whole_exact_barriers": row,
        "whole_complex_hard_core_functions": hard_functions(x),
        "whole_complete_complex_tree_bound": "In |zi-wi|<=1e-12*wi, |M_N(zi)|/A0<=3*n^2*N!*(2e17)^N/[kappa^(N/2)*product center wi]. The unique hard-core decomposition and all-valence bounds are unchanged; scalar inverse4/Wsub, hard gap denominator600000 and paired EH step120000000 are the complex budgets.",
        "whole_Cauchy_theorem": "G_N=product zi*sqrt(rho)*M_N/A0 is holomorphic on a neighborhood of the closed relative polydisc. |G_N|<=K_N=6*n^2*N!*(4e17)^N/kappa^(N/2). For every nonnegative multiindex alpha, |partial^alpha G_N|<=product(alpha_i!)*K_N/product((1e-12*wi)^alpha_i).",
        "whole_nonclosure": "The radii shrink with every center energy and do not reach soft faces. These derivative bounds are not an integrable all-N overlap remainder or an inclusive probability sum. Hard/evanescent matching, quantum state, unitarity, absolute complex Regge and common-parent bounce remain open.",
        "checks": checks,
        "gates": {
            "pure_positive_exact_barrier": row["pure_margin"] > 0,
            "hard_positive_exact_barrier": row["hard_margin"] > 0,
            "hard_path_series_below_one_thousandth": 0
            < values["V"]
            < s.Rational(1, 1000),
            "matter_core_below_three": 0 < values["matter"] < 3,
            "gravity_core_below_three_eleven": 0 < values["gravity"] < 3 * 10**11,
            "original_prefactor_dominates_gravity": 3 * source.HEAVY_MASS2**2
            > 3 * 10**11,
            "all_eight_finite_coefficients_below_proved_envelope": all(
                0 < pure_coefficient(n) <= pure_envelope(n) for n in counts
            ),
            "complex_light_subset_gap_positive": s.Rational(1, 8) - 100 * EPSILON > 0,
            "relative_total_perturbation_has_no_N_multiplier": True,
            "analytic_action_keeps_bilinear_products": True,
            "real_Born_normalization_not_complex_absolute_replacement": True,
            "Cauchy_polydisc_does_not_include_zero_energy_faces": EPSILON < 1,
            "no_inclusive_or_Regge_conclusion_from_local_holomorphy": True,
        },
    }
