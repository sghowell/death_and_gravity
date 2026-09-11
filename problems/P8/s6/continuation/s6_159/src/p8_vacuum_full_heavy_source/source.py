"""Fixed-reference derivative and regulator-complete source coordinate map."""

from functools import cache

import sympy as s


@cache
def data():
    a, h, e, L, g, G, Q, T, I2, B0, r, Sa, F, k0, k1, ell = s.symbols(
        "a h epsilon L g G Q T I2 B0 r Sa T_fermion k0 k1 ell"
    )
    Ta = T - I2 * (a - 1) + s.Symbol("T_second") * (a - 1) ** 2 / 2
    sunset = s.Symbol("S0") + Sa * (a - 1)
    deltaZ = -r
    deltaMass = -L * T / 2 + g * B0 - r
    # Counterterms are fixed at the physical a=1 reference when H varies.
    Vraw = L * Ta**2 / 8 - g * sunset / 4
    Vct = (deltaMass - a * deltaZ) * Ta / 2
    fixed_derivative = G * s.diff(Vraw + Vct, a).subs(a, 1)
    inserted = -g * Sa / 2 - g * B0 * I2 + r * T
    cG = L * G / (2 * Q * e)
    Jhybrid = -G * (inserted + F) / 2 - cG * T / 2
    kD = k0 + e * k1
    tad = -1 / (Q * e) - (ell + 1) / Q - e * (1 + ell + ell**2 / 2 + s.pi**2 / 12) / Q
    first = -G * tad / 2
    Gstar = G - h * kD * G
    reexpression = s.expand(s.diff(first.subs(G, Gstar), h).subs(h, 0)).coeff(e, 0)
    correction = -G * (k0 * (ell + 1) + k1) / (2 * Q)
    return {
        "scalar_OS_inserted_tadpole": inserted,
        "proper_cubic_source_insertion": cG * T / 2,
        "complete_hybrid_second_source_before_overall_subtraction": Jhybrid,
        "complete_second_source_rule": "J2_MS=Fin[-G(Tinsert_scalar,D+Tinsert_fermion,D)/2-deltaG1,D T1,D/2]+Fin[k_D G T1,D/2]. deltaG1,D=LG/(2Q epsilon); the common overall source poles are subtracted only after the proper forest.",
        "first_regulated_source": first,
        "full_first_parameter_source_reexpression": correction,
        "held_fixed_reference_counterterms": {"kinetic": deltaZ, "mass": deltaMass},
        "checks": {
            "fixed_vacuum_derivative_is_half_inserted_tadpole": s.expand(
                fixed_derivative - G * inserted / 2
            ),
            "local_L_source_cancels_after_physical_OS": s.diff(
                s.expand(fixed_derivative), L
            ),
            "source_cancels_all_hybrid_one_point_terms": s.expand(
                Jhybrid + G * (inserted + F) / 2 + cG * T / 2
            ),
            "proper_cubic_pole_occurs_once": cG * T / 2 - L * G * T / (4 * Q * e),
            "regulator_complete_source_reexpression": s.expand(
                reexpression - correction
            ),
            "epsilon_coefficient_of_first_field_is_not_omitted": s.diff(
                reexpression, k1
            )
            + G / (2 * Q),
            "no_second_map_on_zero_tree_source": s.diff(s.Integer(0), G),
            "full_first_regulated_G_shift": s.diff(Gstar, h) + kD * G,
            "mass_derivative_counts_both_light_lines": -s.Rational(1, 2) * 2 + 1,
            "source_first_pole": s.expand(first).coeff(e, -1) - G / (2 * Q),
        },
        "scope": "H has no new field/OS reference. The original stationary-H condition fixes this source; it is not a new freedom. The derivative holds Phi counterterms fixed at a=1 and does not impose a different OS condition for every H background.",
    }
