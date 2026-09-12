"""Joint analytic endpoint coefficients and a controlled far-momentum tail."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_retarded_reference_boundary import majorants
from p8_vector_state.comparison import AMAX

from . import domain

ENDPOINT = s.Integer(2) * 10**27
TAYLOR = s.Integer(4) * 10**57
FAR = s.Integer(10) ** 54
TAIL = s.Integer(10) ** 58


@cache
def constants():
    rows = majorants.recurrence()
    C = reference.PAIR_REF
    r = domain.INNER
    m = domain.MASS
    norm = 9 * C * C * sum(sum(rows[j]) / (r * m) ** j for j in range(5))
    return {
        "endpoint_time_rows": tuple(rows[j] for j in range(5)),
        "full_first_five_endpoint_jet_row_over_nu": norm,
        "fourth_spatial_Taylor_remainder_numerator": 2 * ENDPOINT / domain.DELTA**5,
        "far_full_integral_upper": TAYLOR * AMAX**3 / (24 * m),
        "far_removed_both_leg_tail_numerator": TAYLOR * AMAX**4 / 9,
    }


@cache
def data():
    x = s.Symbol("x")
    z = s.Symbol("z", real=True)
    p, nu, K = s.symbols("p nu K", positive=True)
    c = constants()
    checks = {
        "all_five_endpoint_coefficient_rows": sum(
            sum(row) for row in c["endpoint_time_rows"]
        )
        - (1 + 3 + 14 + 91 + 765),
        "spatial_Cauchy_geometric_remainder": s.factor(
            1 / (1 - x) - sum(x**j for j in range(5)) - x**5 / (1 - x)
        ),
        "actual_spatial_Taylor_numerator": c[
            "fourth_spatial_Taylor_remainder_numerator"
        ]
        - TAYLOR,
        "fifth_spatial_power_produces_integrable_internal_fourth_power": s.factor(
            2 * ENDPOINT * nu * (p / (domain.DELTA * nu)) ** 5 - TAYLOR * p**5 / nu**4
        ),
        "far_half_band_radial_tail": s.integrate(z**-2, (z, K / 2, s.oo)) - 2 / K,
        "canonical_far_finite_piece": 4 * FAR / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 746,
        "canonical_far_regulator_tail": 4 * TAIL / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 742,
    }
    return {
        "pre_current_row": "Apply the S198 product-rule recursion to the COMPLETE pre-current endpoint coefficients, before the final imaginary part and before multiplying test Fourier transforms. For j0..4, the coefficient of Gamma^(r) is bounded by9 Cref^2 c_jr INNER^(-(j-r)) nu^(1-j). The complete row from source time jets0..4 has norm below2e27nu jointly in complex p.",
        "holomorphy": "The conjugate detector amplitude is replaced by its analytic Schwarz continuation. The W8 amplitudes and inverse phase are jointly holomorphic; time differentiation preserves that property. No assertion that Im itself is holomorphic, and no spatial derivative of the tests or Borel state, is used.",
        "Cauchy": "For any real unit direction v, the pre-current row E(lambda v) is holomorphic for |lambda|<=delta nu and bounded by2e27nu. Its degree-four Taylor error at |P|<=delta nu/2 is<=4e57 |P|^5 nu^-4, in the operator norm from the five source time jets to the detector tensor.",
        "far_region": "Define FAR by real |k|>=m and |P|<=delta nu(k)/2. Integrate ONLY the degree-four spatial Taylor remainder over this region, retaining the original both-created-mode regulator. This is not the full endpoint kernel or its local Taylor polynomial.",
        "far_bound": "Since integral nu^-4<=Amax^3/(24m), the far remainder is below1e54||D||L2 X46[Gamma], where X46^2=sum_(j=0)^4||(1-Delta)^3 partial_t^j Gamma||L2^2. The factor |P|^5 is controlled by the six-spatial-derivative weight.",
        "far_tail": "On FAR, |P|<=delta |k| because nu<=2|k|. Hence |l|<=(1+delta)|k|, and a removed both-leg pair implies |k|>K/(1+delta)>K/2. The far remainder tail is below1e58||D||L2 X46[Gamma]/K for K>=1000. Both canonical displays are4e-746 and4e-742/K.",
        "remaining": "The complementary near-momentum region, the full finite/divergent Taylor coefficients, complete reference contact and original fixed covariant spatial matching remain unresolved. This far Taylor piece is not a full response inverse, interacting background, physical cutoff or original V/G/B closure.",
        "constants": c,
        "checks": checks,
        "gates": {
            "complete_first_five_endpoint_bound": c[
                "full_first_five_endpoint_jet_row_over_nu"
            ]
            < ENDPOINT,
            "Cauchy_half_disc_factor": s.Rational(1, 2) < 1,
            "full_far_integral_display": c["far_full_integral_upper"] < FAR,
            "removed_far_tail_display": c["far_removed_both_leg_tail_numerator"] < TAIL,
            "same_two_leg_far_tail_geometry": 1 + domain.DELTA < 2,
            "far_internal_power_is_absolutely_integrable": s.integrate(
                z**-2, (z, 1, s.oo)
            )
            == 1,
        },
    }
