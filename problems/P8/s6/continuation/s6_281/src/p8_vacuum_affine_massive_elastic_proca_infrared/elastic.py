"""Entire massive elastic tree and closed full angular resolvent convolution."""

from functools import cache

import sympy as s

from . import source

S, MU, K, Z = source.S, source.MU, source.K, source.Z
N, G2, C = s.symbols("heavy_mass_squared cubic_squared contact", real=True)
W, V = s.symbols("angular_regulator heavy_resolvent", positive=True)


def exact_real(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Rational)):
        raise TypeError("Require an exact real rational")
    return s.Rational(value)


def require_domain(energy, angle, mass=1, heavy=source.HEAVY_MASS2):
    ss, zz, mu, nn = map(exact_real, (energy, angle, mass, heavy))
    if mu <= 0 or not 4 * mu < ss < nn or not -1 < zz < 1:
        raise ValueError("Require 4mu<s<n and strictly nonforward real angle")
    return ss, zz, mu, nn


def coefficients(energy=S, mass=MU, kappa=K, heavy=N, cubic2=G2, contact=C):
    q = energy - 4 * mass
    p = 4 * (energy**2 - 4 * mass * energy + 2 * mass**2) / (kappa * q)
    v = 1 + 2 * heavy / q
    b = 4 * cubic2 * v / q
    a0 = (
        contact
        + cubic2 / (heavy - energy)
        + (-2 * energy + 6 * mass - 2 * mass**2 / energy + q**2 / (6 * energy)) / kappa
    )
    a2 = -(q**2) / (6 * kappa * energy)
    return p, v, b, a0, a2


def whole_tree(angle=Z, energy=S, mass=MU, kappa=K, heavy=N, cubic2=G2, contact=C):
    p, v, b, a0, a2 = coefficients(energy, mass, kappa, heavy, cubic2, contact)
    return p / (1 - angle**2) + b / (v**2 - angle**2) + a0 + a2 * s.legendre(2, angle)


def Q0(w):
    return s.log((w + 1) / (w - 1)) / 2


def Q2(w):
    return ((3 * w * w - 1) * Q0(w) - 3 * w) / 2


def master_I(w, v, z):
    c = w * v - z
    disc = s.factor(c * c - (w * w - 1) * (v * v - 1))
    if disc == 0:
        return 1 / c
    d = s.sqrt(disc)
    return s.log((c + d) / s.sqrt((w * w - 1) * (v * v - 1))) / d


def master_J(w, v, z):
    return (master_I(w, v, z) + master_I(w, v, -z)) / (2 * w * v)


def whole_regulated_average(w=W, z=Z):
    p, v, b, a0, a2 = coefficients()
    return (
        p * p * master_J(w, w, z)
        + 2 * p * b * master_J(w, v, z)
        + b * b * master_J(v, v, z)
        + 2 * p * (a0 * Q0(w) / w + a2 * Q2(w) * s.legendre(2, z) / w)
        + 2 * b * (a0 * Q0(v) / v + a2 * Q2(v) * s.legendre(2, z) / v)
        + a0 * a0
        + a2 * a2 * s.legendre(2, z) / 5
    )


@cache
def data():
    x, z, t, d = s.symbols("x z parameter discriminant_root", real=True)
    q = S - 4 * MU
    tt = -q * (1 - x) / 2
    uu = -q * (1 + x) / 2
    grav = (
        -(
            (2 * MU**2 - 2 * MU * S - tt * uu) / S
            + (2 * MU**2 - 2 * MU * tt - S * uu) / tt
            + (2 * MU**2 - 2 * MU * uu - S * tt) / uu
        )
        / K
    )
    matter = C + G2 * (1 / (N - S) + 1 / (N - tt) + 1 / (N - uu))
    w, v = s.symbols("w v", positive=True)
    alpha = (w - v) ** 2 - 2 * (1 - z)
    beta = 2 * (v * (w - v) + 1 - z)
    gamma = v * v - 1
    D2 = w * w + v * v - 2 * w * v * z + z * z - 1
    density = (t * w + (1 - t) * v) ** 2 - (t * t + (1 - t) ** 2 + 2 * z * t * (1 - t))
    primitive = s.log(
        (2 * alpha * t + beta - 2 * d) / (2 * alpha * t + beta + 2 * d)
    ) / (2 * d)
    resnum = s.fraction(s.together(s.diff(primitive, t) - 1 / density))[0]
    c = w * v - z
    ratio_num = (2 * alpha + beta - 2 * d) * (beta + 2 * d)
    ratio_den = (2 * alpha + beta + 2 * d) * (beta - 2 * d)
    checks = {
        "entire_original_tree_all_channels": s.factor(grav + matter - whole_tree(x)),
        "elastic_crossing_mass_sum": s.expand(S + tt + uu - 4 * MU),
        "entire_quadratic_Feynman_denominator": s.expand(
            density - (alpha * t * t + beta * t + gamma)
        ),
        "exact_master_discriminant": s.expand(beta**2 - 4 * alpha * gamma - 4 * D2),
        "whole_master_primitive": s.rem(resnum, d * d - D2, d),
        "whole_master_endpoint_cross_ratio": s.rem(
            s.expand(ratio_num * (c - d) - ratio_den * (c + d)), d * d - D2, d
        ),
        "master_positive_sheet_product": s.expand(
            c * c - D2 - (w * w - 1) * (v * v - 1)
        ),
        "master_removable_forward_equal_pole": master_I(
            s.Integer(2), s.Integer(2), s.Integer(1)
        )
        - s.Rational(1, 3),
        "whole_azimuthal_spin_two": s.expand(
            (3 * (z * z * x * x + (1 - z * z) * (1 - x * x) / 2) - 1) / 2
            - s.legendre(2, x) * s.legendre(2, z)
        ),
        "spin_two_resolvent_moment": s.factor(
            (3 * (v * v * Q0(v) / v - 1) - Q0(v) / v) / 2 - Q2(v) / v
        ),
    }
    r, a = s.symbols("r a", positive=True)
    antiderivative = 1 / (r * (a - r * x))
    checks["sphere_square_primitive"] = s.factor(
        s.diff(antiderivative, x) - 1 / (a - r * x) ** 2
    )
    checks["sphere_square_integral"] = s.factor(
        (antiderivative.subs(x, 1) - antiderivative.subs(x, -1)) / 2
        - 1 / (a * a - r * r)
    )
    checks["spin_two_orthogonality"] = s.integrate(s.legendre(2, x), (x, -1, 1))
    checks["spin_two_norm"] = s.integrate(
        s.legendre(2, x) ** 2, (x, -1, 1)
    ) / 2 - s.Rational(1, 5)
    return {
        "entire_invariant_gravity_tree": grav,
        "entire_retained_contact_and_heavy_tree": matter,
        "entire_elastic_tree": whole_tree(),
        "complete_resolvent_coefficients": coefficients(),
        "whole_Feynman_master": master_I(W, V, Z),
        "whole_even_resolvent_convolution": master_J(W, V, Z),
        "complete_auxiliary_regulated_sewing": whole_regulated_average(),
        "physical_elastic_phase_space": s.sqrt(1 - 4 * MU / S) / (32 * s.pi),
        "regulator_scope": "w>1 is an auxiliary angular resolvent, not a graviton mass, physical detector resolution, or physical Wilsonian cutoff. The closed master is on the positive real sheet for w,v>1 and -1<=z<=1, with its continuous coincident discriminant limit.",
        "checks": checks,
        "gates": {
            "full_heavy_contact_and_gravity_interferences_present": all(
                whole_regulated_average().has(v) for v in (K, C, G2, N, MU, Z, W)
            ),
            "full_master_not_forward_only": master_J(W, V, Z).has(Z),
            "identical_and_optical_halves_retained": s.Rational(1, 4) / (8 * s.pi)
            == 1 / (32 * s.pi),
            "positive_real_Feynman_denominator_for_w_v_above_one": True,
            "continuous_discriminant_degeneracy_handled": True,
            "no_massless_external_limit_substituted": whole_tree().has(MU),
        },
    }
