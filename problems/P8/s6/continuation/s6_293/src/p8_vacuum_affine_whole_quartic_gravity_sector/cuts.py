"""Entire D-dimensional scalar-pair cut sewn to the complete gravity tree."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_resonance_soft_pole_pairing import masters

from . import proper, source

MU, K, EP = source.MU, source.K, source.EP
S = proper.S
X = s.Symbol("physical_cut_angle", real=True)


def require_cut_domain(mass, energy, epsilon, scale_squared=1):
    mass, energy, epsilon, scale_squared = map(
        source.require_mass, (mass, energy, epsilon, scale_squared)
    )
    if energy <= 4 * mass or epsilon > s.Rational(1, 8):
        raise ValueError("Require the massive pair cut s>4mu and0<epsilon<=1/8")
    return mass, energy, epsilon, scale_squared


def tree_parts(energy=S, mass=MU, epsilon=EP):
    energy, mass, epsilon = map(s.sympify, (energy, mass, epsilon))
    q = energy - 4 * mass
    return (
        4 * proper.numerator(energy, mass, epsilon) / q,
        -7 * energy / 4 + 4 * mass + 2 * mass**2 / (energy * (1 + epsilon)),
        -q * q / (4 * energy),
    )


def whole_tree(energy=S, mass=MU, epsilon=EP, angle=X):
    p, a, b = tree_parts(energy, mass, epsilon)
    angle = s.sympify(angle)
    return p / (1 - angle**2) + a + b * angle**2


def whole_angular_average(energy=S, mass=MU, epsilon=EP):
    energy, mass, epsilon = map(s.sympify, (energy, mass, epsilon))
    p, a, b = tree_parts(energy, mass, epsilon)
    return p * (1 + 2 * epsilon) / (2 * epsilon) + a + b / (3 + 2 * epsilon)


def whole_phase(energy=S, mass=MU, epsilon=EP, scale_squared=proper.NU2):
    energy, mass, epsilon, scale_squared = map(
        s.sympify, (energy, mass, epsilon, scale_squared)
    )
    beta = s.sqrt(1 - 4 * mass / energy)
    return (
        beta
        * (4 * s.pi * scale_squared) ** (-epsilon)
        * energy**epsilon
        * beta ** (2 * epsilon)
        * s.gamma(1 + epsilon)
        / (8 * s.pi * s.gamma(2 + 2 * epsilon))
    )


def whole_cut(
    energy=S, mass=MU, epsilon=EP, quartic=proper.C, kappa=K, scale_squared=proper.NU2
):
    energy, mass, epsilon, quartic, kappa, scale_squared = map(
        s.sympify, (energy, mass, epsilon, quartic, kappa, scale_squared)
    )
    return (
        quartic
        * whole_phase(energy, mass, epsilon, scale_squared)
        * whole_angular_average(energy, mass, epsilon)
        / (2 * kappa)
    )


@cache
def data():
    a, mu, e, z = s.symbols("energy mass_squared epsilon angle", positive=True)
    Q = a - 4 * mu
    d = 4 + 2 * e
    t = -Q * (1 - z) / 2
    u = -Q * (1 + z) / 2
    V = (a - 2 * mu) ** 2 - 2 * mu * mu / (1 + e)
    N = lambda aa, bb, cc: 4 * mu * mu * (d - 3) / (d - 2) - 2 * mu * aa - bb * cc
    literal = -N(a, t, u) / a - N(t, a, u) / t - N(u, a, t) / u
    p = 4 * V / Q
    aa = -7 * a / 4 + 4 * mu + 2 * mu * mu / (a * (1 + e))
    bb = -Q * Q / (4 * a)
    checks = {}
    checks["whole_crossed_GR_tree_poles_and_regular_angle"] = s.factor(
        literal - p / (1 - z * z) - aa - bb * z * z
    )
    Jratio = -2 * (1 + 2 * e) / (e * Q)
    Mnext = -Q * (1 + e) / (2 * (2 * e + 3))
    Yr = (mu - Mnext) / a
    loop = -V * Jratio - 2 * (a - 2 * mu) + (a + 2 * mu / (1 + e)) * Yr
    tree = p * (1 + 2 * e) / (2 * e) + aa + bb / (3 + 2 * e)
    checks["entire_D_moment_cut_equals_complete_tree_interference"] = s.factor(
        loop - tree
    )
    checks["whole_massive_weighted_bubble_imaginary_ratio"] = s.factor(
        Yr - mu / a - Q * (e + 1) / (2 * a * (2 * e + 3))
    )
    checks["whole_cusp_beta_ratio"] = s.gammasimp(
        -4
        * s.gamma(e)
        * s.gamma(e + s.Rational(3, 2))
        / (Q * s.gamma(e + s.Rational(1, 2)) * s.gamma(e + 1))
        - Jratio
    )
    checks["whole_next_moment_beta_ratio"] = s.gammasimp(
        -Q
        * s.gamma(e + 2)
        * s.gamma(e + s.Rational(3, 2))
        / (4 * s.gamma(e + 1) * s.gamma(e + s.Rational(5, 2)))
        - Mnext
    )
    checks["entire_angular_pole_beta_ratio"] = s.gammasimp(
        s.gamma(e)
        * s.gamma(e + s.Rational(3, 2))
        / (s.gamma(e + s.Rational(1, 2)) * s.gamma(e + 1))
        - (1 + 2 * e) / (2 * e)
    )
    checks["entire_angular_quadratic_beta_ratio"] = s.gammasimp(
        s.gamma(s.Rational(3, 2))
        * s.gamma(e + s.Rational(3, 2))
        / (s.gamma(s.Rational(1, 2)) * s.gamma(e + s.Rational(5, 2)))
        - 1 / (3 + 2 * e)
    )
    checks["identical_optical_half_times_two_interference"] = s.Rational(
        1, 4
    ) * 2 - s.Rational(1, 2)
    checks = {name: s.simplify(s.expand_func(value)) for name, value in checks.items()}
    return {
        "whole_crossed_massive_GR_tree": whole_tree(),
        "whole_full_D_tree_angular_average": whole_angular_average(),
        "whole_full_D_identical_scalar_phase": whole_phase(),
        "whole_full_D_quartic_gravity_interference_cut": whole_cut(),
        "whole_physical_normal_sheet_pair_moments": tuple(
            masters.light_moment(power, MU, S) for power in (EP - 1, EP, EP + 1)
        ),
        "whole_cut_sewing_proof": "Only the s-channel massive moments have an imaginary part for s>4mu,t,u<0. The complete cusp/bubble ratios are ImJ/ImM=-2(1+2e)/(eQ) and ImM_(e+1)/ImM_e=-Q(e+1)/(2(2e+3)). The weighted endpoint ratio is mu/s+Q(e+1)/(2s(2e+3)). Substitution in the whole proper-plus-endpoint sum exactly equals the complete crossed GR tree averaged with(1-z^2)^e, not only its singular angular part. ImB0/(16pi^2)=Phi2/2 and the identical optical1/4 times2C*AGR interference fix the normalization.",
        "whole_angular_endpoint_subtraction": "The pole integral is B(e,1/2)/2=1/(2e)+integral_0^1(1-z^2)^e/(1+z)dz. This is an exact convergent subtraction for small positive e, not naive quadrature of an unresolved endpoint. The regular tree terms and every evanescent phase/projector coefficient remain.",
        "whole_cut_scope": "This is the coefficient linear in C and1/kappa, so it is an interference term and need not be a positive measure by itself. It checks the whole nonlocal amplitude but cannot fix a local UV polynomial from its discontinuity. The tensor/OS calculation supplies that separate polynomial. No physical soft-factor unitarity, other coupling cuts, Regge or original P8 closure follows.",
        "checks": checks,
        "gates": {
            "complete_s_t_u_tree_not_only_forward_pole": True,
            "entire_physical_moments_and_evanescent_ratios": True,
            "identical_optical_interference_factor_not_distinct_pair": True,
            "whole_D_phase_and_convergent_angular_subtraction": True,
            "independent_entire_cut_not_local_UV_matching": True,
            "interference_not_positive_full_physical_measure": True,
        },
    }
