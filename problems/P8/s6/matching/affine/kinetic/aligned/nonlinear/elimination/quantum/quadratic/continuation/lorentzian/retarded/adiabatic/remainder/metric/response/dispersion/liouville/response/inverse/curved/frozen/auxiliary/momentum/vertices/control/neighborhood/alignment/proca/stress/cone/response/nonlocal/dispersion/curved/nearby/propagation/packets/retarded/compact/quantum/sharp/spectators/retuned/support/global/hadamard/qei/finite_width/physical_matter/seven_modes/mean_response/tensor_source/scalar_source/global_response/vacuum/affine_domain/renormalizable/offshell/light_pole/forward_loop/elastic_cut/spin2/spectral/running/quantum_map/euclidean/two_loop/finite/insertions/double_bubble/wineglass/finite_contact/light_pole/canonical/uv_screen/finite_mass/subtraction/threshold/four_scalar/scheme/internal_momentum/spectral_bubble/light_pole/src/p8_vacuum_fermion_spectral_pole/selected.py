"""Explicit selected quadratic functional, not a complete two-loop propagator."""

from functools import cache

import sympy as sp


@cache
def data():
    s, h, r, fp, r2 = sp.symbols("s h scalar_slope fermion_slope spectral_family_slope")
    Fs, Ps, P2s = sp.symbols(
        "fermion_OS_remainder scalar_OS_remainder spectral_OS_remainder"
    )
    f1, p1, p21 = sp.symbols(
        "whole_fermion_mass_anchor whole_scalar_mass_anchor whole_spectral_family_mass_anchor"
    )
    kappa = 1 + h * (r - fp) + h * h * r2
    mass = 1 + h * (p1 - f1) + h * h * p21
    f = f1 + (s - 1) * fp + Fs
    p = p1 + (s - 1) * r + Ps
    p2 = p21 + (s - 1) * r2 + P2s
    raw = mass - s + h * (f - p) - h * h * p2
    normalized = raw / kappa
    target = 1 - s + (h * (Fs - Ps) - h * h * P2s) / kappa
    bF, bS, b2 = sp.symbols("bF bS b2")
    local = target.subs(
        {Fs: bF * (s - 1) ** 2, Ps: bS * (s - 1) ** 2, P2s: b2 * (s - 1) ** 2}
    )
    t, dF, dS, d2 = sp.symbols(
        "t nonnegative_fermion_defect nonnegative_scalar_defect nonnegative_spectral_defect"
    )
    spacelike = target.subs({s: -t, Fs: -dF, Ps: dS, P2s: d2, h: 1}, simultaneous=True)
    kactual = 1 + r - fp + r2
    checks = {
        "explicit_mass_and_kinetic_reference_identity": sp.factor(normalized - target),
        "selected_mass_one_anchor": sp.simplify(local.subs(s, 1)),
        "selected_unit_residue_anchor": sp.simplify(sp.diff(local, s).subs(s, 1) + 1),
        "selected_spacelike_all_sector_signs": sp.factor(
            (1 + t - spacelike) * kactual - dF - dS - d2
        ),
        "one_loop_limit_preserved": sp.factor(sp.diff(target, h).subs(h, 0) - Fs + Ps),
        "named_second_order_field_normalization_retained": sp.factor(
            sp.diff(target, h, 2).subs(h, 0) / 2 + P2s + (r - fp) * (Fs - Ps)
        ),
        "higher_orders_not_declared_absent": sp.factor(sp.diff(kappa, h, 2) / 2 - r2),
    }
    return {
        "formal_loop_parameter": h,
        "selected_local_mass_reference": mass,
        "selected_kinetic_normalization": kappa,
        "selected_normalized_inverse": target,
        "selected_spacelike_inverse": spacelike,
        "scope": "Complete one-loop quadratic functional plus only the first scalar covariance-insertion family, with its explicit local references and exact constant normalization. It is not the complete new-model two-loop inverse or a later-loop error estimate; the formal h expansion retains the displayed field-normalization products.",
        "checks": checks,
    }
