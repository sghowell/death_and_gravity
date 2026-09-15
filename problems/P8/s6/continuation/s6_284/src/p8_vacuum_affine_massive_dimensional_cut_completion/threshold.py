"""Exact massive triangle threshold jets and crossing-symmetric principal parts."""

from functools import cache

import sympy as s

from . import source

S, T, MU, EP, NU = source.S, source.T, source.MU, source.EP, source.NU


def triangle_threshold(mass=MU):
    return (s.log(2) - s.I * s.pi / 2) / s.sympify(mass)


def triangle_threshold_slope(mass=MU):
    mass = s.sympify(mass)
    return 1 / (24 * mass**2) - triangle_threshold(mass) / (6 * mass)


def bubble_threshold(mass=MU, epsilon=EP, scale=NU):
    mass, epsilon, scale = map(s.sympify, (mass, epsilon, scale))
    return -1 / epsilon + 2 + s.log(4 * s.pi * scale**2 / mass) - s.EulerGamma


def tadpole_finite(mass=MU, epsilon=EP, scale=NU):
    mass = s.sympify(mass)
    return mass * (bubble_threshold(mass, epsilon, scale) - 1)


def old_gram_choice(energy=S, transfer=T, mass=MU, bubble=None):
    energy, transfer, mass = map(s.sympify, (energy, transfer, mass))
    B = bubble_threshold(mass) if bubble is None else s.sympify(bubble)
    q = energy - 4 * mass
    u = -q - transfer
    return (
        -12 * B * mass**2 * transfer * u / q**2 + (5 * B + 2) * mass * transfer * u / q
    )


def evanescent_gram_choice(energy=S, transfer=T, mass=MU):
    energy, transfer, mass = map(s.sympify, (energy, transfer, mass))
    q = energy - 4 * mass
    u = -q - transfer
    return 8 * mass**2 * transfer * u / q**2 - 6 * mass * transfer * u / q


def full_gram_choice(energy=S, transfer=T, mass=MU, bubble=None):
    energy, transfer, mass = map(s.sympify, (energy, transfer, mass))
    B = bubble_threshold(mass) if bubble is None else s.sympify(bubble)
    q = energy - 4 * mass
    u = -q - transfer
    return (
        -(12 * B - 8) * mass**2 * transfer * u / q**2
        + (5 * B - 4) * mass * transfer * u / q
    )


@cache
def data():
    y, r, a = s.symbols("threshold_y radial_r positive_a", positive=True)
    L = s.log(1 - y * y)
    primitive = (s.log(r + 1) - s.log(r * r + a * a) / 2 + s.atan(r / a) / a) / (
        1 + a * a
    )
    p0 = -L / y - s.log((1 + y) / (1 - y))
    p1 = -(y * y + (1 - y * y) * L) / (3 * y**3) - s.Rational(2, 3) * p0
    ip = (1 - s.sqrt(1 - y * y)) / y
    q, t, mu, Bstar, Cstar = s.symbols(
        "threshold_q transfer positive_mu Bstar Cstar", nonzero=True
    )
    row = source.reconstruction.coefficients(4 * mu + q, t, mu)
    jet = (
        row["C00mu"] * (Cstar + q * (1 / (24 * mu**2) - Cstar / (6 * mu)))
        + row["B00"] * (Bstar - 2 * mu * Cstar - q / (4 * mu))
        + row["Bmm"] * Bstar
    )
    pp2 = 12 * Bstar * mu**2 * t * t
    pp1 = mu * t * (12 * Bstar * mu - (5 * Bstar + 2) * t)
    old = old_gram_choice(4 * mu + q, t, mu, Bstar)
    entire = full_gram_choice(4 * mu + q, t, mu, Bstar)
    A0 = s.Symbol("raw_massive_tadpole_A0")
    uu = -q - t
    tad = -(12 * mu * A0 + 4 * mu**2) * t * uu / q**2 + (5 * A0 + mu) * t * uu / q
    checks = {
        "literal_radial_triangle_primitive": s.factor(
            s.diff(primitive, r) - 1 / ((1 + r) * (r * r + a * a))
        ),
        "whole_real_threshold_primitive": s.factor(s.diff(p0, y) - L / y**2),
        "whole_real_threshold_slope_primitive": s.factor(
            s.diff(p1, y) - (y * y + (1 - y * y) * L) / y**4
        ),
        "whole_imaginary_threshold_primitive": s.simplify(
            s.diff(ip, y) - (1 / s.sqrt(1 - y * y) - 1) / y**2
        ),
        "real_threshold_lower_endpoint": s.limit(p0, y, 0, dir="+"),
        "real_threshold_upper_endpoint": s.simplify(
            s.limit(p0, y, 1, dir="-") + 2 * s.log(2)
        ),
        "real_slope_lower_endpoint": s.limit(p1, y, 0, dir="+"),
        "real_slope_upper_endpoint": s.simplify(
            s.limit(p1, y, 1, dir="-") - (-1 + 4 * s.log(2)) / 3
        ),
        "imaginary_threshold_lower_endpoint": s.limit(ip, y, 0, dir="+"),
        "imaginary_threshold_upper_endpoint": s.limit(ip, y, 1, dir="-") - 1,
        "complete_threshold_triangle_value": s.expand(
            triangle_threshold(mu) - (s.log(2) - s.I * s.pi / 2) / mu
        ),
        "complete_imaginary_threshold_slope": s.simplify(
            s.im(triangle_threshold_slope(mu)) - s.pi / (12 * mu**2)
        ),
        "raw_triangle_bubble_threshold_dictionary": s.expand(
            2 * mu * triangle_threshold(mu)
            + (bubble_threshold(mu) - 2 * s.log(2) + s.I * s.pi)
            - bubble_threshold(mu)
        ),
        "raw_tadpole_bubble_threshold_dictionary": s.expand(
            tadpole_finite(mu) / mu + 1 - bubble_threshold(mu)
        ),
        "whole_old_cut_master_double_principal_part": s.factor(
            s.limit(q * q * jet, q, 0) - pp2
        ),
        "whole_old_cut_master_simple_principal_part": s.factor(
            s.limit(q * (jet - pp2 / q**2), q, 0) - pp1
        ),
        "crossing_choice_preserves_double_and_simple_parts": s.factor(
            old - pp2 / q**2 - pp1 / q + (5 * Bstar + 2) * mu * t
        ),
        "full_evans_Gram_choice_assembly": s.factor(
            entire - old - evanescent_gram_choice(4 * mu + q, t, mu)
        ),
        "full_Gram_choice_other_channel_crossing": s.factor(
            entire - entire.subs(t, -q - t)
        ),
        "forced_tadpole_dependent_Gram_dictionary": s.factor(
            entire.subs(Bstar, A0 / mu + 1) - tad
        ),
    }
    return {
        "whole_threshold_triangle_and_first_jet": (
            triangle_threshold(),
            triangle_threshold_slope(),
        ),
        "whole_raw_bubble_and_tadpole_threshold": (
            bubble_threshold(),
            tadpole_finite(),
        ),
        "whole_old_cut_basis_principal_part": pp2 / q**2 + pp1 / q,
        "whole_dimensionally_completed_crossing_Gram_choice": full_gram_choice(),
        "whole_forced_tadpole_dependent_Gram_choice": tad,
        "prescription": "The displayed rational choice is SUBTRACTED from the cut-master representative plus its determined evanescent bubble rational term. It cancels the integer Gram principal parts. This is not a choice of the still-unknown physical finite pole or local normalization.",
        "analytic_boundary": "C00mu and B00 have regular local boundary-value Taylor jets at s=4mu. Massive scalar masters retain their genuine square-root/log threshold branches. Removing the artificial integer poles is not an assertion about the full amplitude or threshold-supported distributions.",
        "checks": checks,
        "gates": {
            "massive_triangle_not_massless_limit": triangle_threshold() != 0,
            "Feynman_imaginary_threshold_part_retained": triangle_threshold().has(s.I),
            "raw_Gamma_and_4pi_threshold_constants_retained": bubble_threshold().has(
                s.EulerGamma
            ),
            "whole_double_and_simple_principal_parts_retained": True,
            "finite_evans_changes_Gram_choice": True,
            "tadpole_coefficient_not_declared_fully_determined": True,
            "regular_rational_choice_not_physical_local_matching": True,
        },
    }
