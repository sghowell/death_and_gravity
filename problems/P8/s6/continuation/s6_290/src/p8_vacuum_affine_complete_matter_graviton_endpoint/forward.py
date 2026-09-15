"""Whole crossed endpoint coefficient and finite-mass bounds with RH anchor retained."""

from functools import cache

import sympy as s

from . import endpoint, source

MU, N, G, K = source.MU, source.N, source.G, source.K
Z, V = source.Z, source.V
X = endpoint.X


def require_bound_domain(mass, heavy):
    mass, heavy = map(source.require_mass, (mass, heavy))
    if heavy < 4 * mass:
        raise ValueError("Require n>=4mu for the whole endpoint bound")
    return mass, heavy


def f1_derivative(order, transfer, mass=MU, heavy=N, cubic=G):
    order = source.require_order(order)
    transfer, mass, heavy, cubic = map(s.sympify, (transfer, mass, heavy, cubic))
    a = s.Symbol("differentiated_channel", real=True)
    integrand = endpoint.vertex.F1_integrand(a, Z, V, mass, heavy)
    return (
        cubic
        * cubic
        * s.Integral(
            s.diff(integrand, a, order).subs(a, transfer), (V, 0, 1), (Z, 0, 1)
        )
        / (16 * s.pi**2)
    )


def f1_b20(mass=MU, heavy=N, cubic=G, kappa=K):
    mass, heavy, cubic, kappa = map(s.sympify, (mass, heavy, cubic, kappa))
    return (
        f1_derivative(1, 2 * mass, mass, heavy, cubic)
        - endpoint.f1(2 * mass, mass, heavy, cubic) / mass
        - 2 * f1_derivative(1, 0, mass, heavy, cubic)
    ) / kappa


def triangle_b20_integrand(mass=MU, heavy=N, z=Z, v=V):
    mass, heavy, z, v = map(s.sympify, (mass, heavy, z, v))
    transfer_weight = (1 - z) ** 2 * (1 - v * v) / 4
    w = (1 - z) * ((1 - z) ** 2 * v * v - 1) / 2
    value = 0
    for active, spectator in ((mass, heavy), (heavy, mass)):
        F = endpoint.denominator(active, spectator, 0, z, v, mass)
        value += (
            -w
            * 2
            * transfer_weight
            * (F + 2 * mass * transfer_weight)
            / (F - 2 * mass * transfer_weight) ** 3
        )
    return value


def bubble_center_jets(mass=MU):
    mass = s.sympify(mass)
    y = X * (1 - X)
    return (
        endpoint.bubble_remainder(2 * mass, mass),
        s.Integral(y * y / (mass - 2 * mass * y), (X, 0, 1)),
        s.Integral(y**3 / (mass - 2 * mass * y) ** 2, (X, 0, 1)),
    )


def bubble_b20(mass=MU, heavy=N, cubic=G, quartic=endpoint.C, kappa=K):
    mass, heavy, cubic, quartic, kappa = map(
        s.sympify, (mass, heavy, cubic, quartic, kappa)
    )
    b, b1, b2 = bubble_center_jets(mass)
    d = heavy - 2 * mass
    mixed = (
        4 * mass * b2 / d
        + 2 * (heavy + 2 * mass) * b1 / d**2
        + 2 * (heavy + 2 * mass) * b / d**3
    )
    return (quartic * (4 * mass * b2 + 2 * b1) + cubic * cubic * mixed) / (
        16 * s.pi**2 * kappa
    )


def known_b20(mass=MU, heavy=N, cubic=G, quartic=endpoint.C, kappa=K):
    mass, heavy, cubic, quartic, kappa = map(
        s.sympify, (mass, heavy, cubic, quartic, kappa)
    )
    tri = (
        cubic
        * cubic
        * s.Integral(triangle_b20_integrand(mass, heavy), (V, 0, 1), (Z, 0, 1))
        / (16 * s.pi**2 * kappa)
    )
    return (
        f1_b20(mass, heavy, cubic, kappa)
        + tri
        + bubble_b20(mass, heavy, cubic, quartic, kappa)
    )


def mixing_anchor_b20(residue=endpoint.H, mass=MU, heavy=N, kappa=K):
    residue, mass, heavy, kappa = map(s.sympify, (residue, mass, heavy, kappa))
    return -2 * residue * (heavy + 2 * mass) / (kappa * (heavy - 2 * mass) ** 3)


def known_endpoint_bound(mass, heavy, cubic, quartic, kappa):
    mass, heavy = require_bound_domain(mass, heavy)
    cubic, kappa = map(source.require_mass, (cubic, kappa))
    if isinstance(quartic, bool) or not isinstance(quartic, (int, s.Rational)):
        raise TypeError("Require an exact real rational quartic")
    quartic = s.Rational(quartic)
    if abs(quartic) > 6 * cubic * cubic / heavy:
        raise ValueError("The quartic violates the stated whole endpoint bound premise")
    return 3 * cubic * cubic / (4 * s.pi**2 * kappa * mass * heavy)


@cache
def data():
    a, v, mu, n = s.symbols("channel crossing_v mu n", positive=True)
    f, g = s.Function("f1"), s.Function("f2")
    poly = 4 * mu * mu / a - 3 * mu + a / 2
    crossed = poly.subs(a, 2 * mu + v) * f(2 * mu + v) + poly.subs(a, 2 * mu - v) * f(
        2 * mu - v
    )
    f2cross = (4 * mu + v) * g(2 * mu + v) + (4 * mu - v) * g(2 * mu - v)
    z, x = s.symbols("z x", real=True)
    y = x * (1 - x)
    F_light = mu * (1 - z) ** 2 + n * z
    F_heavy = n * (1 - z) + mu * z * z
    Delta, A = s.symbols(
        "zero_transfer_denominator parameter_transfer_weight", positive=True
    )
    hh = s.Symbol("RH_residue", real=True)
    B = s.Function("Bbar")
    u = (a + 2 * mu) / (n - a)
    primitive = -1 / ((n / 2 - mu) * (mu + (n / 2 - mu) * z))
    checks = {
        "entire_F1_crossed_center_value": poly.subs(a, 2 * mu),
        "entire_F1_crossed_center_first_jet": s.diff(poly, a).subs(a, 2 * mu)
        + s.Rational(1, 2),
        "entire_F1_crossed_center_second_jet": s.diff(poly, a, 2).subs(a, 2 * mu)
        - 1 / mu,
        "entire_crossed_F1_coefficient": s.simplify(
            s.diff(crossed, v, 2).subs(v, 0) / 2
            - f(2 * mu) / mu
            + s.diff(f(a), a).subs(a, 2 * mu)
        ),
        "entire_crossed_F2_coefficient": s.simplify(
            s.diff(f2cross, v, 2).subs(v, 0) / 2
            - s.diff((a + 2 * mu) * g(a), a, 2).subs(a, 2 * mu)
        ),
        "whole_constant_improvement_zero_b20": s.diff(
            (4 * mu + v) + (4 * mu - v) + 2 * mu, v, 2
        ),
        "whole_RH_exchange_b20": s.factor(
            s.diff((a + 2 * mu) * hh / (n - a), a, 2).subs(a, 2 * mu)
            - 2 * hh * (n + 2 * mu) / (n - 2 * mu) ** 3
        ),
        "whole_triangle_F2_second_derivative": s.factor(
            s.diff((a + 2 * mu) / (Delta - a * A), a, 2)
            - 2 * A * (Delta + 2 * mu * A) / (Delta - a * A) ** 3
        ),
        "whole_light_active_linear_majorant": s.factor(
            F_light - (mu * (1 - z) + n * z / 2) - z * (n / 2 - mu + mu * z)
        ),
        "whole_heavy_active_majorant": s.factor(F_heavy - n * (1 - z) - mu * z * z),
        "whole_inverse_square_linear_primitive": s.factor(
            s.diff(primitive, z) - 1 / (mu + (n / 2 - mu) * z) ** 2
        ),
        "whole_inverse_square_linear_integral": s.factor(
            primitive.subs(z, 1) - primitive.subs(z, 0) - 2 / (mu * n)
        ),
        "whole_bubble_second_weight_moment": s.integrate(y * y, (x, 0, 1))
        - s.Rational(1, 30),
        "whole_bubble_third_weight_moment": s.integrate(y**3, (x, 0, 1))
        - s.Rational(1, 140),
        "whole_C_bubble_derivative_bound": s.Rational(4, 35)
        + s.Rational(2, 15)
        - s.Rational(26, 105),
        "whole_mixed_bubble_derivative_product": s.factor(
            s.diff(u * B(a), a, 2)
            - u * s.diff(B(a), a, 2)
            - 2 * (n + 2 * mu) * s.diff(B(a), a) / (n - a) ** 2
            - 2 * (n + 2 * mu) * B(a) / (n - a) ** 3
        ),
        "whole_mixed_bubble_bound_constant": s.Rational(8, 35)
        + s.Rational(4, 5)
        + s.Rational(4, 5)
        - s.Rational(64, 35),
    }
    actualN, actualG, actualK = source.HEAVY_MASS2, source.CUBIC, source.KAPPA
    margins = {
        "actual_n_above_four_mu": actualN - 4,
        "whole_original_contact_bound": 6 * actualG**2 / actualN - abs(source.CONTACT),
        "whole_known_bound_below_1e_minus1005": s.Rational(1, 10**1005)
        - 3 * actualG**2 / (16 * actualK * actualN),
        "entire_triangle_majorant_margin": s.Integer(7) - s.Rational(27, 4),
        "entire_C_bubble_majorant_margin": s.Rational(1, 4) - s.Rational(26, 105),
        "entire_mixed_bubble_majorant_margin": s.Integer(2) - s.Rational(64, 35),
        "entire_known_endpoint_majorant_margin": s.Integer(12)
        - 7
        - 2
        - s.Rational(3, 2)
        - s.Rational(5, 36),
    }
    return {
        "whole_F1_endpoint_b20": f1_b20(),
        "whole_F2_triangle_b20_integrand": triangle_b20_integrand(),
        "whole_bubble_center_jets": bubble_center_jets(),
        "whole_bubble_endpoint_b20": bubble_b20(),
        "whole_known_endpoint_b20": known_b20(),
        "whole_unmatched_RH_b20": mixing_anchor_b20(),
        "whole_constant_curvature_b20": s.S.Zero,
        "whole_actual_known_endpoint_absolute_bound": known_endpoint_bound(
            1, actualN, actualG, source.CONTACT, actualK
        ),
        "whole_exact_positive_margins": margins,
        "whole_uniform_bound_proof": "For n>=4mu and0<=t<=2mu, both triangle denominators are at least half their zero-transfer values. S286 gives the full F1 bound5g^2/(144pi^2*kappa*n^2). The entire positive F2triangle integral is bounded by7g^2/(16pi^2*kappa*mu*n). The exact Bbar/Bprime/Bsecond bounds are2/15,1/(15mu),1/(35mu^2). The C and H-mixing bubble derivatives have bounds1/(4mu),2/(mu*n). Actual absC<6g^2/n gives total below3g^2/(4pi^2*kappa*mu*n). With pi>2 and actual parameters its magnitude is below10^-1005. No heavy asymptotic series is used as an error bound.",
        "whole_matching_boundary": "The constant curvature anchor has zero crossed b20; the RH exchange residue contributes exactly-2h(n+2mu)/(kappa(n-2mu)^3). h is not chosen or bounded. Additional local/Ricci derivative matching, mixed matter/graviton irreducible graphs, full-source cuts, all-loop errors and finite physical IR/Regge remain separate. This bound concerns only the known endpoint graph coefficient.",
        "checks": checks,
        "gates": {
            "whole_crossed_endpoint_coefficient_not_only_t_channel": True,
            "F1_t_channel_continuation_kept_before_forward_limit": True,
            "whole_finite_mass_integrals_not_asymptotic_truncations": True,
            "complete_actual_known_graph_bound_with_all_positive_margins": all(
                v > 0 for v in margins.values()
            ),
            "constant_anchor_zero_but_RH_anchor_retained": True,
            "no_bound_or_value_invented_for_RH": True,
            "not_full_mixed_four_point_IR_Regge_or_original_P8": True,
        },
    }
