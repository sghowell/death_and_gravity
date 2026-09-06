"""Fraction-only reconstruction of the mode and inverse norm constants.

No symbolic engine, primary bound helper, quadrature or floating-point
coefficient is used here. The proofs of the integral representations and
infinite-series domination remain the accompanying written lemmas.
"""

from fractions import Fraction as Q
from math import factorial


def dyson(order, duration, split):
    # Sum all allocations between the two actual complex-conjugate modes.
    allocations = sum((Q(1, factorial(left)*factorial(order-left))
                       for left in range(order+1)), Q(0))
    infrared = allocations*duration**(2*order)*split**2/2
    ultraviolet = allocations*duration**order*split**(2-order)/(order-2)
    return {"infrared": infrared, "ultraviolet": ultraviolet,
            "total": infrared+ultraviolet}


def logarithm_upper(ratio, count):
    v = (ratio-1)/(ratio+1)
    power, subtotal = v, Q(0)
    for j in range(count):
        subtotal += 2*power/(2*j+1)
        power *= v*v
    return subtotal+2*power/((2*count+1)*(1-v*v))


def mode_calibration(a7):
    first_jet = Q(a7["derived_finite_geometry"]["U_derivative_bounds_through_5_per_delta"][1])
    cap = 2*Q(1, 10**14)*first_jet
    duration, future = Q(3), Q(1, 4)
    amplitude = cap*duration
    z = amplitude*duration**2
    exp_upper = 1/(1-2*z)
    # 3/2 bounds the radial IR+UV weights; 3/2 bounds n^2/n!
    # relative to 1/(n-3)!; all three base factors contribute 2^3.
    tail_coefficient = Q(3, 2)*Q(3, 2)*2**3
    full_derivative = 2*z+tail_coefficient*z*z*exp_upper
    log_upper = logarithm_upper(duration/future, 8)
    shared_quadratic = amplitude*future**2*(Q(5, 4)+log_upper/2)
    shared_higher = tail_coefficient*z*z*exp_upper*future/duration
    return {"prepared_U1": first_jet, "derivative_cap": cap, "strength": z,
            "full_derivative": full_derivative,
            "full_value": duration*full_derivative,
            "log_ratio_upper": log_upper,
            "shared_quadratic": shared_quadratic,
            "shared_higher": shared_higher,
            "shared_derivative": shared_quadratic+shared_higher,
            "shared_value": future*(shared_quadratic+shared_higher),
            "full_rounding_gap": Q(1, 10**8)-full_derivative,
            "shared_rounding_gap": Q(1, 10**10)-shared_quadratic-shared_higher}


def inverse_calibration():
    length, split = Q(1, 2**1024), Q(2**512)
    pole_argument = 3*length
    pole = 80*pole_argument/(1-pole_argument)
    lower_cut = 8/split
    upper_cut = Q(80, 255)
    norm = pole+lower_cut+upper_cut
    return {"length": length, "split": split, "pole_norm_upper": pole,
            "low_cut_norm_upper": lower_cut, "high_cut_norm_upper": upper_cut,
            "norm_upper": norm, "rounded_gap": Q(1, 3)-norm,
            "abstract_q": Q(1, 3), "abstract_distance": Q(1, 10)/(1-Q(1, 3)),
            "abstract_ball_margin": 1-Q(1, 10)-Q(1, 3)}


def exact_polynomial_checks():
    # Divide x^4(1-x)^4 by 1+x^2 in an independent coefficient array.
    numerator = [Q(0)]*4+[Q(1), Q(-4), Q(6), Q(-4), Q(1)]
    remainder = numerator[:]
    quotient = [Q(0)]*7
    for degree in range(8, 1, -1):
        coefficient = remainder[degree]
        quotient[degree-2] = coefficient
        remainder[degree] -= coefficient
        remainder[degree-2] -= coefficient
    if any(remainder[degree] for degree in range(2, 9)) or remainder[:2] != [Q(-4), Q(0)]:
        raise ValueError("The independent positive pi-integral division failed")
    integral_polynomial = sum((coefficient/Q(degree+1)
                               for degree, coefficient in enumerate(quotient)), Q(0))
    if integral_polynomial != Q(22, 7):
        raise ValueError("The independent positive pi-integral constant failed")
    # 3/2 - n/((n-1)(n-2)) has positive numerator
    # (n-3)(3n-2)/2. Form the left difference and right product separately.
    left_numerator = [Q(3, 2)*2, Q(3, 2)*(-3)-1, Q(3, 2)]
    right_numerator = [Q(0)]*3
    for j, first in enumerate((-3, 1)):
        for k, second in enumerate((-2, 3)):
            right_numerator[j+k] += Q(first*second, 2)
    if left_numerator != right_numerator:
        raise ValueError("The factorial majorant coefficient identity failed")
    return {"pi_integral_polynomial_part": str(integral_polynomial),
            "pi_integral_remainder": "-4/(1+x^2)",
            "pi_squared_upper_gap": str(10-Q(22, 7)**2),
            "factorial_gap_factored": "(n-3)*(3*n-2)/(2*(n-1)*(n-2))"}


def serialize(value):
    if isinstance(value, dict):
        return {name: serialize(item) for name, item in value.items()}
    return str(value)


def replay(a7):
    modes, inverse = mode_calibration(a7), inverse_calibration()
    if min(modes["full_rounding_gap"], modes["shared_rounding_gap"], inverse["rounded_gap"]) <= 0:
        raise ValueError("An independent response or inverse rounded margin failed")
    return {"mode_calibration": serialize(modes),
            "inverse_calibration": serialize(inverse),
            "polynomial_checks": exact_polynomial_checks(),
            "dyson_allocation_checks": {str(n): serialize(dyson(n, Q(3), Q(1, 3)))
                                        for n in range(3, 10)}}
