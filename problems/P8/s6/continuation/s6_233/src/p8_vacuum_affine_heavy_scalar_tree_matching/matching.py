"""Exact complete original on-shell target plus a uniform massive rational remainder."""

from functools import cache

import sympy as s
from p8_vacuum_canonical_affine_decoupling import amplitude as source

from . import model as m

S, T = s.symbols("s t", real=True)
U = 4 - S - T
v = s.Symbol("crossing_v", real=True)
D, g = m.D, m.g
GAM = g * g / D**4
LAM = g * g / (2 * D**3)
TREE = m.C + g * g * (1 / (D + 2 - S) + 1 / (D + 2 - T) + 1 / (D + 2 - U))
TARGET = (
    2 * LAM * ((S - 2) ** 2 + (T - 2) ** 2 + (U - 2) ** 2)
    + 3 * GAM * S * T * U
    - 8 * GAM
)
REMAINDER = GAM * sum((z - 2) ** 4 / (D - (z - 2)) for z in (S, T, U))
ENERGY = s.Integer(10) ** 98
RATIO_UPPER = s.Rational(32, 625)
RELATIVE_BOUND = s.Rational(1, 60)


@cache
def data():
    xi = s.Symbol("shifted_channel")
    finite = sum(xi**n / D ** (n + 1) for n in range(4))
    original = source.data()["complete_on_shell_tree_amplitude"]
    symbols = {str(z): z for z in original.free_symbols}
    literal = original.subs(
        {
            symbols["s"]: S,
            symbols["t"]: T,
            symbols["u"]: U,
            symbols["lambda"]: LAM,
            symbols["gamma"]: GAM,
        },
        simultaneous=True,
    )
    forward = REMAINDER.subs({S: 2 + v, T: 0})
    r = s.Symbol("positive_relative_s_gap", positive=True)
    x = s.Symbol("angle_cosine", real=True)
    tp, up = -(S - 4) * (1 - x) / 2, -(S - 4) * (1 + x) / 2
    lower = 2 * LAM * (S - 2) ** 2 - 8 * GAM
    lambda_part = 2 * LAM * ((S - 2) ** 2 + (tp - 2) ** 2 + (up - 2) ** 2)
    r_actual = (ENERGY**2 - 2) / m.GAP
    return {
        "complete_named_model_tree": TREE,
        "complete_original_target_tree": TARGET,
        "exact_full_angular_remainder": REMAINDER,
        "exact_matching_identity": "For s+t+u=4, A_V2S_T1=A_original+gamma*sum(channel-2)^4/[D-(channel-2)]. The exact finite geometric identity retains the original mass and potential. It is not an identification of full-field actions or a derivative expansion with an omitted remainder.",
        "actual_matching_energy": ENERGY,
        "uniform_tree_relative_error": "On4<=s<=10^196 and the complete physical angular interval, 0<A_V2S_T1-A_original<(1/60)A_original. For r=(s-2)/D<1, the relative error is at most6r^2/(1-r). The actual r<32/625 gives6144/370625<1/60.",
        "uniform_bound_inputs": "All |channel-2|<=s-2, every heavy denominator is positive below its s-channel pole, and A_original>=lambda(s-2)^2 because the Galileon term is nonnegative physically and4lambda>8gamma. No massless approximation or forward-only estimate is used.",
        "forward_remainder": 16 * GAM / (D + 2) + 2 * GAM * D * v**4 / (D**2 - v**2),
        "forward_coefficients": "The named model has exactly b20_tree=4lambda and b21_tree=-3gamma, b40_tree=gamma^2/lambda, and its forward constant differs from the original tree by16gamma/(D+2). Its higher coefficients do not vanish.",
        "boundary": "The uniform error compares two complete TREE amplitudes. It is not the full angular quantum error premise of S231, a physical cutoff, exact unitarity or common-parent bounce matching. The heavy pole and all its higher terms remain.",
        "checks": {
            "literal_original_source_tree": s.expand(literal - TARGET),
            "massive_shifted_channel_sum": sum(z - 2 for z in (S, T, U)) + 2,
            "massive_shifted_cubic_identity": s.expand(
                sum((z - 2) ** 3 for z in (S, T, U)) - 3 * S * T * U + 8
            ),
            "finite_geometric_identity_with_remainder": s.cancel(
                1 / (D - xi) - finite - xi**4 / (D**4 * (D - xi))
            ),
            "complete_tree_matching_identity": s.cancel(TREE - TARGET - REMAINDER),
            "actual_lambda_parameter_match": m.G2 / (2 * m.GAP**3) - m.LAMBDA,
            "actual_gamma_parameter_match": m.G2 / m.GAP**4 - m.GAMMA,
            "exact_forward_remainder": s.cancel(
                forward - 16 * GAM / (D + 2) - 2 * GAM * D * v**4 / (D**2 - v**2)
            ),
            "forward_second_coefficient_remainder_zero": s.cancel(
                s.diff(forward, v, 2).subs(v, 0)
            ),
            "forward_transfer_second_coefficient_remainder_zero": s.cancel(
                s.diff(REMAINDER.subs(S, 2 - T / 2 + v), v, 2, T).subs({v: 0, T: 0})
            ),
            "forward_fourth_coefficient_remainder": s.cancel(
                s.diff(forward, v, 4).subs(v, 0) / 24 - GAM**2 / LAM
            ),
            "complete_lambda_lower_polynomial": s.expand(
                lambda_part
                - lower
                - 8 * GAM
                - 2 * LAM * ((tp - 2) ** 2 + (up - 2) ** 2)
            ),
            "relative_bound_from_gap_relation": s.cancel(
                3 * GAM * (r * D) ** 2 / (LAM * (D - r * D)) - 6 * r * r / (1 - r)
            ),
            "relative_error_monotonicity": s.cancel(
                s.diff(6 * r * r / (1 - r), r) - 6 * r * (2 - r) / (1 - r) ** 2
            ),
            "named_rational_error_upper": 6 * RATIO_UPPER**2 / (1 - RATIO_UPPER)
            - s.Rational(6144, 370625),
        },
        "gates": {
            "original_massive_threshold_lower_margin_positive": 4 * m.LAMBDA
            - 8 * m.GAMMA
            > 0,
            "named_energy_above_scalar_threshold": ENERGY**2 > 4,
            "named_energy_below_original_heavy_pole": ENERGY**2 < m.MASS2,
            "actual_ratio_below_exact_rational_majorant": 0
            < r_actual
            < RATIO_UPPER
            < 1,
            "uniform_relative_tree_error_below_one_sixtieth": s.Rational(6144, 370625)
            < RELATIVE_BOUND,
            "actual_forward_constant_remainder_positive": 16 * m.GAMMA / (m.GAP + 2)
            > 0,
        },
    }
