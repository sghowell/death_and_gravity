"""Positive-pole coefficient diagnostic, explicitly not an exact unitary completion."""

from functools import cache

import sympy as s

from . import coefficients as c
from . import moment

D = 2 * c.LAMBDA / c.GAMMA
MASS2 = D + 2
COUPLING2 = c.GAMMA * D**4


@cache
def data():
    t, v = c.transfer, c.v
    ss, uu = 2 - t / 2 + v, 2 - t / 2 - v
    amplitude = COUPLING2 * (1 / (MASS2 - ss) + 1 / (MASS2 - t) + 1 / (MASS2 - uu))
    b2 = 2 * COUPLING2 / (D + t / 2) ** 3
    b4 = 2 * COUPLING2 / D**5
    J4 = 2 * COUPLING2 / D**4
    return {
        "diagnostic_rational_crossing_amplitude": amplitude,
        "diagnostic_positive_s_channel_atom": "rho_atom(S,0)=pi*g^2 delta(S-M_H^2); the new heavy pole is retained in the positive spectral measure, not subtracted as a known original light pole. Its t-channel term is independent of v and does not alter b2.",
        "diagnostic_mass_squared": MASS2,
        "diagnostic_coupling_squared": COUPLING2,
        "diagnostic_b2_transfer_function": b2,
        "diagnostic_nonzero_v4_coefficient": b4,
        "diagnostic_full_atom_low_moment": J4,
        "scope": "This separately named rational positive-spectral diagnostic matches only b20 and b21 and has a nonzero v4 coefficient. It is not the original parent, an exact unitary S matrix, a physical narrow-width estimate, a matched full tree/function family, a renormalizable UV completion, a prescribed state or a common-parent bounce construction. The real tree rational amplitude alone lacks the light elastic cut required by exact unitarity; no parent action is modified.",
        "purpose": "It demonstrates that the required extra absorptive weight is compatible with these two low coefficients at the level of the positive dispersion data. The mandatory enhancement is not a no-go or a proof that the original first elastic coefficient dominates the full measure.",
        "checks": {
            "diagnostic_half_second_derivative": s.cancel(
                s.diff(amplitude, v, 2).subs(v, 0) / 2 - b2
            ),
            "diagnostic_exact_original_b20": s.cancel(b2.subs(t, 0) - 4 * c.LAMBDA),
            "diagnostic_exact_original_b21": s.cancel(
                s.diff(b2, t).subs(t, 0) + 3 * c.GAMMA
            ),
            "diagnostic_nonzero_v4_coefficient": s.cancel(
                s.diff(amplitude, v, 4).subs({v: 0, t: 0}) / s.factorial(4) - b4
            ),
            "diagnostic_v4_value": s.cancel(b4 - c.GAMMA**2 / c.LAMBDA),
            "diagnostic_positive_atom_low_moment": J4 - 2 * c.GAMMA,
            "actual_rational_trilinear_coupling_squared": COUPLING2
            - s.Rational(1, 2**26),
            "positive_atom_matches_full_transfer_dispersion": -s.Rational(3, 2) * J4
            + 3 * c.GAMMA,
            "positive_atom_saturates_moment_Gram_bound_only": 4 * c.LAMBDA * b4 - J4**2,
        },
        "gates": {
            "diagnostic_heavy_pole_above_original_light_threshold": MASS2 > 4,
            "diagnostic_heavy_pole_below_named_spectral_split": MASS2 < moment.SPLIT,
            "diagnostic_coupling_and_higher_coefficient_positive": COUPLING2 > 0
            and b4 > 0,
            "positive_atom_exceeds_required_full_low_moment": J4
            > moment.NAMED_REQUIRED,
            "rational_diagnostic_not_claimed_as_exact_unitary_completion": True,
            "moment_bound_attainment_not_exact_S_matrix_optimality": True,
        },
    }
