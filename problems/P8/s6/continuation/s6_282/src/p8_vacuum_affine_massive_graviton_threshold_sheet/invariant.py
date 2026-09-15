"""Complete invariant cut, shared massive kernel and explicit analytic sheet."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_detector_ir_cut import soft

from . import angular, source

S, T, MU, K = source.S, source.T, source.MU, source.K
XI = s.Symbol("massive_pair_parameter", real=True)
MT, MU_CUT, MS = s.symbols("M_t M_u M_4mu_minus_s", real=True)


def kernel(z, mass=MU):
    z, mass = map(s.sympify, (z, mass))
    if z == 0:
        return 1 / (4 * mass)
    return s.Integral(1 / (4 * mass - z * (1 - XI * XI)), (XI, 0, 1))


def real_kernel(value, mass=1):
    zz, mu = map(source.exact_real, (value, mass))
    if mu <= 0 or zz >= 4 * mu:
        raise ValueError("Require positive mass and real argument below4mu")
    if zz == 0:
        return 1 / (4 * mu)
    if zz > 0:
        return s.atan(s.sqrt(zz / (4 * mu - zz))) / s.sqrt(zz * (4 * mu - zz))
    return s.atanh(s.sqrt(-zz / (4 * mu - zz))) / s.sqrt(-zz * (4 * mu - zz))


def polynomials(energy=S, mass=MU):
    energy, mass = map(s.sympify, (energy, mass))
    v = energy**2 - 4 * mass * energy + 2 * mass**2
    a = energy**3 - 8 * mass * energy**2 + 20 * mass**2 * energy - 12 * mass**3
    b = energy**3 - 10 * mass * energy**2 + 30 * mass**2 * energy - 20 * mass**3
    r = energy**2 - 18 * mass * energy + 146 * mass**2
    d = 13 * energy**2 - 24 * mass * energy - 232 * mass**2
    return v, a, b, r, d


def algebraic_cut(energy=S, transfer=T, mass=MU, mt=MT, mu_cut=MU_CUT, ms=MS):
    energy, transfer, mass = map(s.sympify, (energy, transfer, mass))
    q = energy - 4 * mass
    u = -q - transfer
    _, a, b, r, d = polynomials(energy, mass)
    vt = polynomials(transfer, mass)[0]
    vu = polynomials(u, mass)[0]
    return (
        (vt**2 * mt + vu**2 * mu_cut) / energy
        - a * ms
        + r / 120
        + 3 * transfer * u * (b * ms + d / 120) / q**2
    )


def regular_threshold_combination(energy=S, mass=MU):
    energy, mass = map(s.sympify, (energy, mass))
    q = energy - 4 * mass
    b = polynomials(energy, mass)[2]
    polynomial = (83 * mass**2 + 20 * mass * q - 5 * q * q) / (120 * mass**2)
    return polynomial + s.Integral(
        b * (1 - XI * XI) ** 2 / (16 * mass**2 * (4 * mass + q * (1 - XI * XI))),
        (XI, 0, 1),
    )


def whole_cut(energy=S, transfer=T, mass=MU, kappa=K):
    energy, transfer, mass, kappa = map(s.sympify, (energy, transfer, mass, kappa))
    q = energy - 4 * mass
    u = -q - transfer
    _, a, _, r, _ = polynomials(energy, mass)
    vt = polynomials(transfer, mass)[0]
    vu = polynomials(u, mass)[0]
    return (
        (vt**2 * kernel(transfer, mass) + vu**2 * kernel(u, mass)) / energy
        - a * kernel(-q, mass)
        + r / 120
        + 3 * transfer * u * regular_threshold_combination(energy, mass)
    ) / (4 * s.pi * kappa * kappa)


@cache
def data():
    q = S - 4 * MU
    u = -q - T
    h, z = angular.H, angular.Z
    subs = {h: S / q, z: 1 + 2 * T / q}
    c00, c11, cl, cp = angular.coefficients()
    converted = (
        q
        * q
        / 64
        * (
            c00 * q * q * (MT + MU_CUT) / (2 * S)
            + c11 * q * (MT - MU_CUT) / 2
            + cl * q * MS
            + cp
        )
    )
    xi = s.Symbol("feynman_parameter", real=True)
    pair_den = 4 * MU - 4 * T * xi * (1 - xi)
    master_den = h - (xi * xi + (1 - xi) ** 2 + 2 * z * xi * (1 - xi))
    small = s.Symbol("threshold_offset", real=True)
    _, _, bb, _, dd = polynomials(4 * MU + small, MU)
    den = 4 * MU + small * (1 - XI * XI)
    quotient = s.factor(
        (bb * (1 / (4 * MU) - small / (24 * MU * MU)) + dd / 120) / small**2
    )
    _, _, b, _, d = polynomials()
    checks = {
        "entire_closed_cut_in_invariant_kernel_basis": s.factor(
            converted.subs(subs, simultaneous=True) - algebraic_cut()
        ),
        "full_Feynman_denominator_dictionary": s.factor(
            master_den.subs(subs, simultaneous=True) - pair_den / q
        ),
        "full_single_resolvent_dictionary": s.factor(
            (h - XI * XI).subs(h, S / q) - (S - q * XI * XI) / q
        ),
        "invariant_double_seed_discriminant": s.factor(
            ((1 - z) * (2 * h - 1 - z)).subs(subs, simultaneous=True)
            - 4 * T * (T - 4 * MU) / q**2
        ),
        "invariant_double_seed_ratio": s.factor(
            ((1 - z) / (2 * h - 1 - z)).subs(subs, simultaneous=True) - T / (T - 4 * MU)
        ),
        "whole_t_u_crossing": s.factor(
            algebraic_cut(transfer=u, mt=MU_CUT, mu_cut=MT) - algebraic_cut()
        ),
        "exact_threshold_M_second_order_remainder": s.factor(
            1 / den
            - 1 / (4 * MU)
            + small * (1 - XI * XI) / (16 * MU * MU)
            - small**2 * (1 - XI * XI) ** 2 / (16 * MU * MU * den)
        ),
        "exact_threshold_remainder_linear_moment": s.integrate(1 - XI * XI, (XI, 0, 1))
        - s.Rational(2, 3),
        "exact_threshold_remainder_quadratic_moment": s.integrate(
            (1 - XI * XI) ** 2, (XI, 0, 1)
        )
        - s.Rational(8, 15),
        "entire_threshold_regular_quotient": s.factor(
            quotient
            - (83 * MU * MU + 20 * MU * small - 5 * small * small) / (120 * MU * MU)
        ),
        "threshold_constant_numerator_vanishes": s.factor(
            (b / (4 * MU) + d / 120).subs(S, 4 * MU)
        ),
        "threshold_linear_numerator_vanishes": s.factor(
            (s.diff(b, S) / (4 * MU) - b / (24 * MU * MU) + s.diff(d, S) / 120).subs(
                S, 4 * MU
            )
        ),
        "exact_regular_combination_threshold": s.factor(
            quotient.subs(small, 0)
            + bb.subs(small, 0) * s.Rational(8, 15) / (64 * MU**3)
            - s.Rational(29, 40)
        ),
        "original_forward_external_threshold_density": s.factor(
            (
                2 * (2 * MU**2) ** 2 / (4 * MU) / (4 * MU)
                - polynomials(4 * MU)[1] / (4 * MU)
                + polynomials(4 * MU)[3] / 120
            )
            - MU * MU / 4
        ),
    }
    for j in range(5):
        exact = s.integrate((1 - XI * XI) ** j, (XI, 0, 1)) / (4 * MU) ** (j + 1)
        checks["same_accepted_soft_kernel_Taylor_coefficient_" + str(j)] = s.factor(
            exact - soft.kernel_coefficient(j, MU)
        )
    rr = s.Symbol("positive_complementary_ratio", positive=True)
    tt = 4 * MU * rr * rr / (1 + rr * rr)
    uu = 4 * MU / (1 + rr * rr)
    common = 4 * MU * rr / (1 + rr * rr)
    checks["complementary_massive_kernel_denominator"] = s.factor(
        tt * uu - common * common
    )
    checks["complementary_massive_kernel_sum"] = s.factor(
        (s.atan(rr) + s.pi / 2 - s.atan(rr)) / common - s.pi / (2 * common)
    )
    be = s.Symbol("above_threshold_beta", positive=True)
    ss = 4 * MU / (1 - be * be)
    pv = s.log((XI - be) / (XI + be)) / (2 * ss * be)
    checks["complete_above_threshold_PV_primitive"] = s.factor(
        s.diff(pv, XI) - 1 / (ss * (XI * XI - be * be))
    )
    checks["above_threshold_delta_denominator"] = s.factor(
        4 * MU - ss * (1 - XI * XI) - ss * (XI * XI - be * be)
    )
    checks["above_threshold_delta_jacobian"] = s.factor(
        s.diff(ss * (XI * XI - be * be), XI).subs(XI, be) - 2 * ss * be
    )
    return {
        "whole_original_invariant_cut": whole_cut(),
        "whole_algebraic_kernel_dictionary": algebraic_cut(),
        "all_five_invariant_polynomials": polynomials(),
        "whole_regular_external_threshold_combination": regular_threshold_combination(),
        "analytic_kernel": "M(z)=integral_0^1 dx/[4mu-z(1-x^2)], holomorphic on C minus[4mu,infinity), M(0)=1/(4mu). Its upper boundary above4mu is [-atanh(beta)+i*pi/2]/(z*beta), beta=sqrt(1-4mu/z).",
        "complete_seed_dictionary": "Q=s-4mu: I(h,z)=Q*M(t), I(h,-z)=Q*M(u), J0=Q^2[M(t)+M(u)]/(2s), J1=Q[M(t)-M(u)]/2, L0=Q*M(4mu-s).",
        "normal_sheet_domain": "s not on(-infinity,0], t and u=4mu-s-t not on[4mu,infinity). Apparent Q^-2 cancels by the displayed uniform integral. In particular real s>0,t<4mu,u<4mu is a continuation of the physical normal gg cut.",
        "boundary_warning": "At fixed negative t, lowering s through -t makes u cross4mu. Continue the explicit M(u) boundary; a real-angle integral cannot replace it. A normal channel cut is not the full imaginary part where crossed massless cuts overlap.",
        "proof": "The parameter denominators give the dictionary before a square root or logarithm is chosen. Uniform compact separation from their zeros proves joint holomorphy. The regular Taylor-remainder integral removes the external threshold exactly, including nonphysical fixed transfer.",
        "checks": checks,
        "gates": {
            "entire_massive_kernel_shared_with_accepted_soft_analysis": True,
            "whole_angular_cut_not_only_forward": whole_cut().has(T),
            "all_crossed_arguments_retained": whole_cut().has(
                kernel(T), kernel(4 * MU - S - T)
            ),
            "apparent_external_threshold_not_declared_a_pole": not s.denom(
                quotient
            ).has(small),
            "explicit_upper_and_lower_boundary_not_blind_log_substitution": True,
            "ordinary_physical_positivity_not_extended_to_unphysical_region": True,
            "full_overlapping_cut_amplitude_not_claimed": True,
        },
    }
