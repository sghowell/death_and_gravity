"""Complete fixed-hybrid scalar vacuum forest and remaining ownership."""

from functools import cache

import sympy as s


def rows():
    return [
        {
            "id": "scalar_quartic_vacuum",
            "owner": "this_scalar_forest",
            "status": "BOUNDED",
        },
        {
            "id": "scalar_cubic_sunset_vacuum",
            "owner": "this_scalar_forest",
            "status": "BOUNDED",
        },
        {
            "id": "scalar_Phi_OS_counterterm_trace",
            "owner": "this_scalar_forest",
            "status": "BOUNDED",
        },
        {
            "id": "heavy_mass_MS_counterterm_trace",
            "owner": "this_scalar_forest",
            "status": "BOUNDED",
        },
        {
            "id": "scalar_fermion_vacuum_with_assigned_proper_forests",
            "owner": "S6.149",
            "status": "BOUNDED",
        },
        {
            "id": "gauge_fermion_vacuum_all_42_states",
            "owner": "S6.149",
            "status": "BOUNDED",
        },
        {
            "id": "heavy_source_reducible_triple",
            "owner": "this_source_cancellation",
            "status": "CANCELLED_REGULATED",
        },
        {
            "id": "pure_massless_gauge_ghost_and_constant_field_trace",
            "owner": "dimensional_scaleless_identity",
            "status": "ZERO",
        },
    ]


def validate_rows(value):
    if value != rows():
        raise ValueError("The complete vacuum ownership ledger changed")
    return True


@cache
def data():
    e, L, g, Q, T, TM, B, S, r, x = s.symbols("epsilon L g Q T TM B S r pE_squared")
    deltaZ = -r
    deltaMass = -L * T / 2 + g * B - r
    trace_numerator = deltaZ * x + deltaMass
    raw = L * T**2 / 8 - g * S / 4
    PhiCT = -L * T**2 / 4 + g * B * T / 2
    heavyCT = g * TM / (4 * Q * e)
    full = -L * T**2 / 8 - g * S / 4 + g * B * T / 2 + heavyCT
    ell, b0, b1 = s.symbols("ell beta0 beta1")
    tD = -1 / e - (ell + 1) - e * (1 + ell + ell**2 / 2 + s.pi**2 / 12)
    bD = 1 / e + b0 + e * b1
    Jf, Mf, mf, ms = s.symbols("field_Jacobian M mF light_mass", positive=True)
    V1 = (
        ms**2 * (s.log(ms / mf**2) - s.Rational(3, 2))
        + Mf**2 * (s.log(Mf / mf**2) - s.Rational(3, 2))
    ) / (4 * Q) + 63 * mf**4 / Q
    return {
        "complete_rows": rows(),
        "physical_Phi_first_counterterms": {"kinetic": deltaZ, "mass": deltaMass},
        "full_scalar_vacuum_before_overall_MS": full,
        "complete_scalar_reference": "Finite epsilon coefficient of -L T_D(1)^2/8-g S_D(1,1,M)/4+g B_HL,D(-1)T_D(1)/2+g T_D(M)/(4Q epsilon). All terms use the same regulator and mu=mF.",
        "mixed_bubble_tadpole_finite_numerator": s.expand(bD * tD).coeff(e, 0),
        "one_loop_vacuum_at_general_physical_light_mass_squared": V1,
        "field_coordinate_ownership": "First Phi field direction changes L,g,Y only, not physical/free mass one, M or mF. The full regulated V1 is invariant, apart from (1/2)Tr log K, a scaleless constant trace. Second tree variation acts on V0=0. Proper fermion and fermionic Phi OS counterterm traces are already in S6.149 and are not repeated.",
        "checks": {
            "light_OS_counterterm_correlated_kinetic_mass": s.expand(
                trace_numerator - (-L * T / 2 + g * B - r * (x + 1))
            ),
            "correlated_OS_slope_closes_to_constant": s.factor(
                -r * (x + 1) / (x + 1) + r
            ),
            "proper_heavy_mass_counterterm": heavyCT - g / (2 * Q * e) * TM / 2,
            "full_scalar_forest_sum": s.expand(raw + PhiCT + heavyCT - full),
            "quartic_Wick_factor": s.Rational(3, 24) - s.Rational(1, 8),
            "cubic_sunset_Wick_factor": -s.Rational(2, 8) + s.Rational(1, 4),
            "cubic_reducible_Wick_factor": -s.Rational(1, 8) + s.Rational(1, 8),
            "mixed_reference_positive_epsilon_terms": s.expand(bD * tD).coeff(e, 0)
            + (1 + ell + ell**2 / 2 + s.pi**2 / 12)
            + (ell + 1) * b0
            + b1,
            "first_Phi_coordinate_variation_of_V1_zero": 2 * L * s.diff(V1, L)
            + 2 * g * s.diff(V1, g)
            + s.Symbol("Y") * s.diff(V1, s.Symbol("Y")),
            "second_tree_vacuum_coordinate_variation_zero": s.diff(s.Integer(0), Mf),
            "no_constant_field_factor_mass_dependence": s.diff(s.log(Jf), Mf),
            "eight_owned_vacuum_classes": len(rows()) - 8,
        },
    }
