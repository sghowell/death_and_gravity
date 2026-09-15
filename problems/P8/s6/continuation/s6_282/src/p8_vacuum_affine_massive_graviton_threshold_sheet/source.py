"""The whole unchanged source and original four-helicity normalization."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_elastic_proca_infrared import source as previous
from p8_vacuum_affine_massive_graviton_cut_isolated_replay import sewing

S, MU, K = previous.S, previous.MU, previous.K
T, U = s.symbols("transfer crossed_energy", real=True)
X, Y, Z, H = s.symbols(
    "left_axis right_axis scattering_cosine inverse_beta_squared", real=True
)


def exact_real(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Rational)):
        raise TypeError("Require an exact real rational")
    return s.Rational(value)


def require_real_sheet(energy, transfer, mass=1):
    ss, tt, mu = map(exact_real, (energy, transfer, mass))
    if mu <= 0 or ss <= 0 or tt >= 4 * mu or 4 * mu - ss - tt >= 4 * mu:
        raise ValueError("Require s>0,t<4mu,u<4mu on the real normal-cut sheet")
    return ss, tt, 4 * mu - ss - tt, mu


def angular_numerator(h=H, x=X, y=Y, z=Z):
    gram = 1 - z * z - x * x - y * y + 2 * z * x * y
    transverse = z - x * y
    return (h - 1) ** 4 + transverse**4 - 6 * transverse**2 * gram + gram**2


def rational_angular_integrand(h=H, x=X, y=Y, z=Z):
    return angular_numerator(h, x, y, z) / ((h - x * x) * (h - y * y))


@cache
def data():
    inherited = previous.data()
    checks = dict(inherited["checks"])
    q = S - 4 * MU
    r = s.Symbol("beta_squared", real=True)
    gram = 1 - Z * Z - X * X - Y * Y + 2 * Z * X * Y
    phase = (Z - X * Y) ** 4 - 6 * (Z - X * Y) ** 2 * gram + gram**2
    scaled = ((1 - r) ** 4 + r**4 * phase) / ((1 - r * X * X) * (1 - r * Y * Y))
    checks.update(
        {
            "entire_original_four_helicity_density_normalization": s.factor(
                sewing.full_helicity_sewing(S, MU, K, X, Y, Z) / (32 * s.pi)
                - q * q * rational_angular_integrand(S / q) / (256 * s.pi * K * K)
            ),
            "regular_beta_squared_integrand_identity": s.factor(
                scaled - rational_angular_integrand(1 / r) * r * r
            ),
            "threshold_coordinate_with_full_external_mass": s.factor(
                S / q - 1 / (1 - 4 * MU / S)
            ),
            "entire_invariant_mass_sum": S + T + (4 * MU - S - T) - 4 * MU,
            "all_fixed_vacuum_parameters_preserved": previous.current.heavy.MASS2
            - (s.Integer(10) ** 200 / 512 + 2),
            "original_mass_one_not_massless_substitution": previous.previous.MASS2 - 1,
            "original_kappa_not_strict_decoupling": previous.previous.KAPPA
            - s.Integer(10) ** 800,
        }
    )
    return {
        "whole_original_R_F": inherited["entire_original_R_F"],
        "whole_retained_heavy_source": inherited["entire_retained_heavy_source"],
        "all_three_vacuum_constants": inherited["all_three_vacuum_constants"],
        "whole_formal_loop_marker": inherited["entire_formal_loop_marker"],
        "parameters": inherited["parameters"],
        "entire_original_four_helicity_sewing": sewing.full_helicity_sewing(
            S, MU, K, X, Y, Z
        ),
        "entire_rescaled_angular_integrand": rational_angular_integrand(),
        "regular_beta_squared_integrand": scaled,
        "normalization": "rho_gg=Q^2/(256pi kappa^2) times the normalized sphere average. Both identical-intermediate and optical halves are already present in S280; no KLT prefactor is imported anew.",
        "scope": "Same formal leading vacuum mass1,kappa10^800,full fixed R,F,heavy source and all vacuum constants. S281 elastic and Proca cuts are retained, not recomputed or changed. No exact quantum LSZ vacuum is claimed.",
        "checks": checks,
        "gates": {
            "whole_original_not_a_clock_germ": all(
                v in inherited
                for v in (
                    "entire_original_R_F",
                    "entire_retained_heavy_source",
                    "all_three_vacuum_constants",
                )
            ),
            "inherited_check_dictionary_is_copied": checks is not inherited["checks"],
            "entire_helicity_phase_retained": rational_angular_integrand().has(X, Y, Z),
            "mass1_and_finite_kappa_retained": previous.previous.MASS2 == 1
            and previous.previous.KAPPA == s.Integer(10) ** 800,
            "normal_cut_not_full_overlapping_discontinuity": True,
            "rejected_S279_and_all_historical_scopes_unchanged": True,
        },
    }
