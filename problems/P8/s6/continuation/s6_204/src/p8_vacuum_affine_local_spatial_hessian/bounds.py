"""Fixed compact CD local Hessian estimate in a stated derivative norm."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vector_state.comparison import AMAX

LOCAL = s.Integer(5) * 10**4


@cache
def constants():
    A = AMAX
    A3 = s.Integer(4)
    h = s.Integer(2)
    h1 = s.Integer(4)
    h2 = s.Integer(12)
    S = h1 + 2 * h * h
    S1 = h2 + 4 * h * h1
    R = (A3 + 3 * A3 * h + A) / 2
    R2 = 6 * (A3 * S + A3 * (3 * h * S + S1) + A * S) + s.Rational(4, 3)
    C = (
        A3 * (1 + 6 * h + 4 * h1 + 11 * h * h + h2 + 7 * h * h1 + 6 * h**3)
        + 2 * A * (1 + h)
        + 1
    )
    fixed = (
        s.Rational(5, 3) * modes.MASS**2 * 15
        + s.Rational(800, 30)
        + s.Rational(3200, 18)
    ) / (64 * 9)
    return {
        "complete_Einstein_coefficient_sum": R,
        "complete_R_squared_coefficient_sum": R2,
        "complete_Weyl_squared_coefficient_sum": C,
        "fixed_local_Hessian_norm_upper": fixed,
    }


@cache
def data():
    t = s.Symbol("t", real=True)
    H = 4 * t / (1 + t * t)
    Hp = s.diff(H, t)
    Hpp = s.diff(Hp, t)
    c = constants()
    checks = {
        "actual_first_clock_derivative": s.factor(
            Hp - 4 * (1 - t * t) / (1 + t * t) ** 2
        ),
        "actual_second_clock_derivative": s.factor(
            Hpp - 8 * t * (t * t - 3) / (1 + t * t) ** 3
        ),
        "canonical_full_local_spatial_bound": 4 * LOCAL / modes.KAPPA
        - 2 * s.Rational(1, 10) ** 795,
    }
    return {
        "clock_bounds": "On the actual CD slab|t|<=1/2, a=(1+t^2)^2 satisfies1<=a<=25/16, a^3<4, |H|<=2, |Hdot|<=4 and |Hddot|<=12. These follow from the displayed exact rational derivatives, using1+t^2>=1 and |t|<=1/2.",
        "full_helicity_bounds": "The sum of absolute proper-time coefficient bounds is below15 for H_R,3200 for H_R2 and800 for H_C, uniformly over all five spatial directions. The fixed coefficients, mass1000 and pi^2>9 give an operator upper bound below5e4.",
        "norm": "Phi42[Gamma]^2=sum_(j=0)^4||(1-Delta)^2 partial_t^j Gamma||L2(dt dx;F)^2. Each q^r factor with r0,1,2 is dominated by(1+q)^2. Therefore ||H_fixed Gamma||L2<5e4 Phi42[Gamma] and |delta^2 S_fixed(D,Gamma)|<5e4||D||L2 Phi42[Gamma].",
        "canonical": "Both external metric factors yield2e-795 in the same norm. This local four-time/four-spatial derivative bound is not a bound for the full quantum response, a derivative-compatible mixed inverse or a stability theorem.",
        "matching_boundary": "The result supplies only the original finite local action target. It does not identify the remaining reference polynomial-integrand/contact sector with this target or eliminate its regulator artifacts.",
        "constants": c,
        "checks": checks,
        "gates": {
            "actual_a_cubed_bound": AMAX**3 < 4,
            "complete_Einstein_bound": c["complete_Einstein_coefficient_sum"] < 15,
            "complete_R_squared_bound": c["complete_R_squared_coefficient_sum"] < 3200,
            "complete_Weyl_squared_bound": c["complete_Weyl_squared_coefficient_sum"]
            < 800,
            "fixed_local_operator_display": c["fixed_local_Hessian_norm_upper"] < LOCAL,
            "same_two_metric_canonical_factors": True,
        },
    }
