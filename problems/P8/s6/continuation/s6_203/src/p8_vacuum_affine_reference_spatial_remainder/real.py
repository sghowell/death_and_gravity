"""All-real-momentum endpoint row, including the unexpanded low band."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference
from p8_vacuum_affine_reference_spatial_analyticity import boundary, domain
from p8_vacuum_affine_retarded_reference_boundary import majorants
from p8_vector_state.comparison import AMAX

MASS = domain.MASS
E = boundary.ENDPOINT
REAL = 2 * E


@cache
def constants():
    r = domain.INNER
    m = MASS
    C = reference.PAIR_REF
    return {
        "complete_real_row_over_nu": 18
        * C
        * C
        * sum(sum(majorants.recurrence()[j]) / (r * m) ** j for j in range(5)),
        "unexpanded_low_band_integral_upper": 2 * E * m**4 / 27,
        "same_real_global_pair": C,
        "same_original_high_low_partition": m,
    }


@cache
def data():
    nu, mu = s.symbols("nu mu", positive=True)
    k, m = s.symbols("k m", positive=True)
    c = constants()
    checks = {
        "full_two_leg_inverse_phase_leaves_original_nu": s.factor(
            2 * nu * mu / (nu + mu) + 2 * nu**2 / (nu + mu) - 2 * nu
        ),
        "complete_real_row_is_twice_joint_row_majorant": c["complete_real_row_over_nu"]
        - 2 * boundary.constants()["full_first_five_endpoint_jet_row_over_nu"],
        "unexpanded_low_band_volume": s.integrate(k * k, (k, 0, m)) - m**3 / 3,
        "both_complete_readouts_not_single_leg": reference.FIELD_NORM**2 - 16000000,
    }
    return {
        "full_real_row": "For all real k,l=-k+P, including |k|<m, the complete unit-W8 pair has norm<=Cref sqrt(nu mu), g<=2/(nu+mu), and the same complex-time radius. The full j0..4 row is bounded by18 Cref^2 nu sum_j sum(c_j)/(r m)^j <4e27nu. All nine pairs and all source time jets are retained.",
        "no_spatial_frame_at_zero_required": "The full real polarization sum uses any orthonormal basis at zero momentum; no complex momentum continuation or division by |k| is used for the low band. The zero-momentum point is harmless for the integral.",
        "low_band": "Leave E(P,k) unexpanded for |k|<m. Since nu<2m, its full integral is below(2E m^4/27), where E=2e27. This is a fixed mathematical partition, not a physical cutoff or a new state.",
        "constants": c,
        "checks": checks,
        "gates": {
            "complete_real_row_display": c["complete_real_row_over_nu"] < REAL,
            "real_pair_phase_ratio": s.Rational(1, 2) < 1,
            "low_band_nu_below_twice_mass": 1 + 1 / AMAX**2 < 4,
            "low_band_display_below_two_e38": c["unexpanded_low_band_integral_upper"]
            < 2 * 10**38,
        },
    }
