"""Independent identical-pair unitarity and a positive Q2 expansion."""

from functools import cache

import sympy as sp

from . import cut


@cache
def data():
    d = cut.data()
    Z = d["Z"]
    z = sp.Symbol("scattering_cosine", real=True)
    n = sp.Symbol("positive_integer_index", integer=True, positive=True)
    P2 = (3 * z * z - 1) / 2
    P2Z = (3 * Z * Z - 1) / 2
    primitive = -P2Z * sp.log(Z - z) - 3 * z * z / 4 - 3 * Z * z / 2
    anchored = (primitive.subs(z, 1) - primitive.subs(z, -1)) / 2
    logQ = (3 * Z * Z - 1) * (sp.log(Z + 1) - sp.log(Z - 1)) / 4 - 3 * Z / 2
    coefficient = 2 * n / ((2 * n + 1) * (2 * n + 3))
    pe, pa, T, g = sp.symbols(
        "positive_external_pair_momentum positive_internal_pair_momentum positive_transfer positive_cubic_squared",
        positive=True,
    )
    beta = 2 * pa / sp.sqrt(T)
    angularK = 2 * pe * pa
    pref = beta / (64 * sp.pi) * (pa * pa / (pe * pe)) * 4 * g / angularK
    A, C = sp.symbols(
        "positive_external_threshold_gap positive_internal_threshold_gap", positive=True
    )
    generic_pref = g * C / (8 * sp.pi * sp.sqrt(T) * A ** sp.Rational(3, 2))
    contact, B, k = sp.symbols(
        "angle_independent_contact positive_B positive_k", real=True
    )
    full_vertex = contact + g / (B - k * z) + g / (B + k * z)
    return {
        "Legendre_P2": P2,
        "full_even_crossed_tree_vertex": full_vertex,
        "anchored_Q2_angular_primitive": primitive,
        "Q2_positive_series_coefficient": coefficient,
        "Q2_positive_series": "sum n>=1 2n/((2n+1)(2n+3))*Z^(-2n-1), Z>1",
        "generic_identical_pair_spin_two_cut_prefactor": generic_pref,
        "normalization": "2 Im and the identical-intermediate-particle 1/2! are separate. dPhi2/dcos=beta/(16pi); hence Im F=beta/(64pi)*(pa^2/pe^2)*integral A P2 dcos.",
        "scope": "Partial waves here are ordinary integer-spin two-particle projections of the specified tree vertices. No complex-spin continuation or Regge tower is asserted.",
        "checks": {
            "exact_Legendre_angular_division": sp.factor(
                P2 / (Z - z) - (P2Z / (Z - z) - 3 * (z + Z) / 2)
            ),
            "literal_angular_primitive_derivative": sp.factor(
                sp.diff(primitive, z) - P2 / (Z - z)
            ),
            "anchored_Q2_angular_integral": sp.expand(anchored - logQ),
            "positive_real_log_ratio_identity_derivative": sp.factor(
                sp.diff(sp.log(Z + 1) - sp.log(Z - 1) - sp.log((Z + 1) / (Z - 1)), Z)
            ),
            "positive_real_log_ratio_identity_anchor": sp.simplify(
                (logQ - d["universal_positive_Q2"]).subs(Z, 2)
            ),
            "local_contact_has_no_spin_two_projection": sp.integrate(P2, (z, -1, 1)),
            "even_crossed_tree_vertex": sp.factor(
                full_vertex.subs(z, -z) - full_vertex
            ),
            "even_moment_positive_series_coefficient": sp.factor(
                (3 / (2 * n + 3) - 1 / (2 * n + 1)) / 2 - coefficient
            ),
            "positive_series_coefficient_upper_gap": sp.factor(
                sp.Rational(2, 15)
                - coefficient
                - 2 * (4 * n - 3) * (n - 1) / (15 * (2 * n + 1) * (2 * n + 3))
            ),
            "geometric_majorant_to_threshold_cube_gap": sp.expand(
                Z * (Z * Z - 1) - (Z - 1) ** 3 - (Z - 1) * (3 * Z - 1)
            ),
            "two_Im_and_identical_phase_factors": sp.factor(
                beta / (16 * sp.pi) / 2 / 2 - beta / (64 * sp.pi)
            ),
            "internal_external_spin_two_momentum_ratio": sp.factor(
                pref - g * pa * pa / (16 * sp.pi * sp.sqrt(T) * pe**3)
            ),
            "generic_partial_wave_matches_triangle_cut_prefactor": sp.simplify(
                pref.subs({pe: sp.sqrt(A) / 2, pa: sp.sqrt(C) / 2}) - generic_pref
            ),
        },
    }
