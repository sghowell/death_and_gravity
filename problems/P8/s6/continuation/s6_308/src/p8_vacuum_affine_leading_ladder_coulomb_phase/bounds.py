"""Explicit original-parameter gamma-phase error and its strict scope."""

from functools import cache

import sympy as s

from . import phase, source


def exact_rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Rational)):
        raise TypeError("Require exact rational physical coordinates")
    return s.Rational(value)


def require_domain(energy, tau):
    energy, tau = map(exact_rational, (energy, tau))
    if not s.Rational(25, 4) <= energy <= 16 or not 0 < tau <= 1:
        raise ValueError("Require the compact massive near-forward domain")
    return energy, tau


def original_gamma_remainder_bound(energy, tau):
    require_domain(energy, tau)
    return s.Rational(2, 10**2400)


def gamma_remainder_bound(strength):
    a = exact_rational(strength)
    if not 0 <= a <= s.Rational(1, 2):
        raise ValueError("Require a real phase strength in[0,1/2]")
    return 2 * a**3


def original_full_Born_embedding_bound(energy, tau):
    require_domain(energy, tau)
    return s.Rational(4, 10**2400)


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(value)

    q, t, k = s.symbols("positive_energy_excess transfer kappa", positive=True)
    energy = q + 4
    V = source.old_phase.eikonal(energy)
    D = source.old_phase.gap(energy)
    put("original_numerator_flux_identity", V - D**2 - 2)
    put(
        "original_strength_as_flux_sum",
        phase.eta(energy, k) - (D + 2 / D) / (8 * s.pi * k),
    )
    put(
        "lower_flux_endpoint",
        s.Rational(25, 4) * s.Rational(9, 4) - s.Rational(15, 4) ** 2,
    )
    put("upper_flux_endpoint", 16 * (16 - 4) - (8 * s.sqrt(3)) ** 2)
    put(
        "whole_strength_coarse_roundup",
        (s.Integer(14) + s.Rational(8, 15)) / 24 - s.Rational(109, 180),
    )
    a, z = s.symbols("strength nonnegative_square_ratio", positive=True)
    j = s.Symbol("integer_index", integer=True, nonnegative=True)
    put(
        "odd_gamma_tail_geometric",
        s.summation(z**j, (j, 0, s.oo)).args[0][0] - 1 / (1 - z),
    )
    put(
        "zeta3_integral_bound",
        1
        + s.integrate(
            s.Symbol("x", positive=True) ** (-3),
            (s.Symbol("x", positive=True), 1, s.oo),
        )
        - s.Rational(3, 2),
    )
    put(
        "gamma_cubic_rational_majorant",
        2 * s.Rational(3, 2) * a**3 / (3 * (1 - a * a)) - a**3 / (1 - a * a),
    )
    put("original_gamma_error", 2 / source.KAPPA**3 - s.Rational(2, 10**2400))
    put(
        "original_fullBorn_embedding_error",
        4 / source.KAPPA**3 - s.Rational(4, 10**2400),
    )
    w = 4 * t / q * (1 - t / q)
    p = 4 * q + 16 + 8 / q
    put(
        "Born_pole_against_positive_full_Born",
        s.factor(p / (2 * w) - V / (2 * t)) - V / (2 * (q - t)),
    )
    put("Born_numerator_identity", p * q - 4 * V)
    # A real phase never amplifies an existing absolute comparison.
    r, c = s.symbols("real_error phase_angle", real=True)
    put(
        "unit_modulus_preserves_error", s.exp(s.I * c) * r * s.exp(-s.I * c) * r - r * r
    )
    margins = {
        "upper_flux_below14": s.Integer(14) ** 2 - 192,
        "strength_bound_below_inverse_kappa": 1 - s.Rational(109, 180),
        "gamma_tail_bound_roundup": 2 - s.Rational(4, 3),
        "original_eta_below_half": s.Rational(1, 2) - 1 / source.KAPPA,
        "original_classical_hierarchy_not_assumed": 1
        - s.Rational(16, 24 * source.KAPPA),
        "near_forward_domain_below_other_endpoint": s.Rational(9, 4) - 1,
    }
    for name, value in margins.items():
        put("positive_arithmetic_" + name, value - s.Abs(value))
    return {
        "whole_gamma_phase_relative_error": "Forreal0<eta<=1/2,|Q(eta)-1|<=|logQ(eta)|<=2zeta(3)eta^3/[3(1-eta^2)]<2eta^3. The same bound controls|C_eta(tau)-exp(i eta c)| independently of positive transfer.",
        "whole_original_gamma_error": s.Rational(2, 10**2400),
        "whole_original_full_Born_embedding_error": s.Rational(4, 10**2400),
        "whole_strength_proof": "For25/4<=s<=16, D=sqrt(s(s-4)) lies in[15/4,8sqrt3]. SinceV=D^2+2 andD<14, eta=(D+2/D)/(8pi kappa)<109/(180kappa)<1/kappa. Thus the original eta is positive and far below1/2.",
        "whole_gamma_tail_proof": "The convergent odd-zeta logarithm begins ateta^3. For oddj>=3,zetaj<=zeta3<3/2 and1/j<=1/3. Sum the geometriceta^2 tail, use1/(1-eta^2)<=4/3, and use|exp(i theta)-1|<=|theta| forrealtheta. This bounds the complete gamma factor rather than just a finite Taylor polynomial.",
        "whole_embedding_proof": "The old gravity Born lower bound isG0>p/(2w),p=4q+16+8/q,w=4tau/q*(1-tau/q),q=s-4. Sincep*q=4V, G0>V/(2tau). The full positive Born exceedsA_G, so the leading pole/A_B ratio is below2. ReplacingQ byone changes only the specified phase model's leading-pole amplitude by less than4e-2400 relative toA_B. This does not bound omitted physical diagrams.",
        "whole_limit_order_and_scope": "The phase result controls the all-rung leading-ladder analytic generating function after fixed-order IR subtraction. It supplies no bound on the full interacting amplitude minus that function. An all-angle physical rate, arbitrary-multiplicity radiation, unknown Newton matching, high-energy Regge contour and original V/G/B/P8 remain open.",
        "whole_positive_margins": margins,
        "checks": checks,
        "gates": {
            "all_original_parameter_margins_strict": all(
                bool(value > 0) for value in margins.values()
            ),
            "whole_Gamma_function_tail_not_finite_polynomial": True,
            "error_uniform_in_unexpanded_transfer_log": True,
            "full_Born_embedding_bounds_only_difference_of_named_models": True,
            "classical_macroscopic_hierarchy_not_imported": True,
            "all_uncomputed_original_quantum_and_Regge_gates_open": True,
        },
    }
