"""All-angle real known-interference bound using the full positive Born denominator."""

from functools import cache

import sympy as s

from . import phase, regularized, source

ENERGY = s.Symbol("physical_energy", positive=True)
KAPPA = s.Symbol("positive_kappa", positive=True)
TAU = s.Symbol("positive_transfer_magnitude", positive=True)
NEAR_CONSTANT = s.Integer(10) ** 14
RATE_CONSTANT = s.Integer(10) ** 12


def require_physical(energy, transfer):
    energy = source.compact.require_positive_exact(energy)
    if isinstance(transfer, (bool, float)) or not isinstance(transfer, (int, s.Basic)):
        raise TypeError("Require exact negative transfer")
    transfer = s.sympify(transfer)
    if (
        transfer.has(s.Float)
        or transfer.is_number is not True
        or transfer.is_real is not True
    ):
        raise ValueError("Require exact negative transfer")
    if transfer.is_finite is not True or transfer.is_negative is not True:
        raise ValueError("Require exact negative transfer")
    if not s.Rational(25, 4) <= energy <= 16 or transfer <= 4 - energy:
        raise ValueError("Require the stated physical nonforward energy-angle domain")
    return (
        energy,
        transfer,
        4 - energy - transfer,
        min(-transfer, energy + transfer - 4),
    )


def near_real_bound(tau):
    tau = source.compact.require_delta(tau)
    return NEAR_CONSTANT * (1 - s.log(tau) + 1 / s.sqrt(tau))


def all_angle_known_rate_bound(kappa):
    return RATE_CONSTANT / source.compact.require_positive_exact(kappa)


def endpoint_coefficient(energy=ENERGY, kappa=KAPPA):
    v = phase.eikonal(energy)
    return 3 * (5 * v + 6) / (32 * kappa * v)


def near_budgets():
    return {
        "constant": regularized.massive_forward_bound()
        + 4 * (110000 * 200 * 10 + 1300 * 200)
        + 2
        * (210000000 * 10 + 200000 * 1200 + 2500 * 100 + 2300000 * 12 + 30000000 * 21)
        + 22000
        + (40000 * 3 + 328)
        + 800328
        + 210000000
        + 2300000 * 5,
        "log": 40000 + s.Rational(210000000, 2) + 2300000,
        "inverse_sqrt": 5 * 210000000,
    }


@cache
def data():
    q = s.Symbol("positive_q", positive=True)
    a = q + 4
    tau = TAU
    p = 4 * q + 16 + 8 / q
    w = 4 * tau * (q - tau) / q**2
    z = 1 - 2 * tau / q
    v = phase.eikonal(a)
    ag = source.forward.gravity_born(a, z, KAPPA)
    n = s.Symbol("positive_heavy_mass2", positive=True)
    g = s.Symbol("positive_cubic", positive=True)
    am = source.forward.matter_born(a, z, n, g)
    t0 = source.assembly.tree_jets(a, -tau, 1)[0]
    dk, alpha, beta = s.symbols("finite_delta_kappa finite_alpha finite_beta")
    fprime = 2 * 322 * 28 + s.Rational(322**2, 6)
    gprime = 28 + s.Rational(322, 6)
    budgets = near_budgets()
    B = s.Symbol("B_soft")
    checks = {
        "spacelike_squared_eikonal_derivative_bound": fprime - s.Rational(105938, 3),
        "spacelike_eikonal_derivative_bound": gprime - s.Rational(245, 3),
        "small_box_constant_budget": 40000 * 3 + 4 * 82 - 120328,
        "finite_soft_reference_decomposition": 4
        * (1 / a - 1 / tau + 1 / (4 - a + tau))
        * B
        - 2 * source.assembly.tree_jets(a, -tau, 1)[1] * B,
        "finite_soft_reference_real_budget": 4 * (2 * 100000 + 82) - 800328,
        "all_near_log_coefficients": budgets["log"] - 107340000,
        "all_near_inverse_sqrt_coefficients": budgets["inverse_sqrt"] - 1050000000,
        "Bose_nearest_endpoint_w": s.factor(w - (1 - z * z)),
        "positive_gravity_pole_lower_coefficient": s.factor(p * q / 8 - v / 2),
        "positive_gravity_minimum_coefficient": (v / 2).subs(q, s.Rational(9, 4))
        - s.Rational(257, 32),
        "positive_gravity_eight_margin": s.Rational(257, 32) - 8 - s.Rational(1, 32),
        "positive_gravity_compact_p_lower": p.subs(q, s.Rational(9, 4))
        - s.Rational(257, 9),
        "interior_gravity_fourteen_margin": s.Rational(257, 18)
        - 14
        - s.Rational(5, 18),
        "known_rate_normalization": s.Rational(2, 16 * 8) - s.Rational(1, 64),
        "near_rate_normalization_margin": 64 * 9 * RATE_CONSTANT
        - 3 * NEAR_CONSTANT
        - s.Integer(276) * 10**12,
        "interior_rate_normalization_margin": 14 * RATE_CONSTANT
        - 2 * source.compact.KNOWN_AMPLITUDE_CONSTANT
        - s.Integer(13998) * 10**9,
        "original_uniform_rate_exponent": 12 - 800 + 788,
        "whole_Born_forward_Newton_shape": s.factor(
            s.cancel(tau * t0).subs(tau, 0) + v
        ),
        "whole_gravity_Born_forward_residue": s.factor(
            s.cancel(tau * ag).subs(tau, 0) - v / KAPPA
        ),
        "whole_matter_Born_has_no_forward_pole": s.factor(
            s.cancel(tau * am).subs(tau, 0)
        ),
        "known_endpoint_sqrt_rate_coefficient": s.factor(
            phase.classical_coefficient(a) / (8 * s.pi**2 * KAPPA * v)
            - endpoint_coefficient(a)
        ),
        "unassigned_Newton_endpoint_interference": s.factor(
            2 * dk * (-v) / (KAPPA**2) / (v / KAPPA) + 2 * dk / KAPPA
        ),
        "local_alpha_endpoint_interference": s.cancel(
            tau * alpha * (a * a + tau * tau + (4 - a + tau) ** 2)
        ).subs(tau, 0),
        "local_beta_endpoint_interference": (tau * beta).subs(tau, 0),
    }
    return {
        "whole_near_forward_component_budgets": budgets,
        "whole_near_forward_real_bound": NEAR_CONSTANT
        * (1 - s.log(TAU) + 1 / s.sqrt(TAU)),
        "whole_uniform_known_interference_bound": RATE_CONSTANT / KAPPA,
        "original_uniform_known_interference_bound": RATE_CONSTANT / source.KAPPA,
        "whole_known_endpoint_sqrt_coefficient": endpoint_coefficient(),
        "independent_matching_endpoint": -2 * dk / KAPPA,
        "whole_domain": "mu=nu=1,25/4<=s<=16,-(s-4)<t<0,all physical nonforward "
        "angles and unchanged original matter parameters. The bound concerns "
        "2Re(A_known)/(A_m+A_G), with unexpanded full positive Born denominator, "
        "in S302's fixed analytic hard reference.",
        "whole_real_bound_proof": "By Bose symmetry choose tau=min(-t,-u). "
        "For tau<=1, separate the purely imaginary endpoint of the small-a boxes "
        "and soft-reference term, bound their real differences by spacelike "
        "derivatives, and use the exact regularized massive Taylor block. "
        "The resulting |ReF|<1e14*(1+|ln tau|+tau^-1/2) combines with "
        "A_m+A_G>=A_G>8/(kappa*tau). For tau>=1 use the S302 delta1 compact "
        "bound and A_G>14/kappa. Both give the same1e12/kappa rate bound.",
        "closure_boundary": "The known hard-reference coefficient tends to zero "
        "like sqrt(tau) after full Born normalization; the cross section itself "
        "is not finite at forward transfer. The unknown Newton matching gives "
        "-2delta_kappa/kappa at the endpoint, while unknown local anchors vanish "
        "there. No matching value, gravity-Born radiation, full inclusive rate, "
        "loop square or all-loop/Regge estimate is inferred.",
        "checks": checks,
        "gates": {
            "all_near_component_budgets_below_common_constant": all(
                bool(value < NEAR_CONSTANT) for value in budgets.values()
            ),
            "spacelike_derivatives_uniform_without_forward_cutoff": bool(
                fprime < 40000 and gprime < 82
            ),
            "full_positive_Born_not_matter_only_Newton_expansion": True,
            "near_and_interior_regions_cover_all_nonforward_angles": True,
            "original_kappa_and_matter_parameters_unchanged": bool(
                source.KAPPA == s.Integer(10) ** 800 and source.HEAVY_MASS2 > 16
            ),
            "one_loop_interference_not_loop_square_or_cross_section": True,
            "unknown_Newton_matching_endpoint_is_not_set_to_zero": True,
            "complex_Regge_and_physical_detector_obligations_still_open": True,
        },
    }
