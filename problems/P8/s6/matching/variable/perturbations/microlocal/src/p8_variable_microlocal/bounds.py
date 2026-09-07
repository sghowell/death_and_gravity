"""Constructive uniform high-K Cauchy denominator bound at fixed delta.

The threshold is a sufficient mathematical symbol threshold, not a
physical cutoff or an admissible EFT window. It is intentionally not
uniform in c->2 without a fixed positive delta_min.
"""

from fractions import Fraction
from functools import cache

from . import rational as ra
from . import symbol


def threshold(delta_min):
    if isinstance(delta_min, bool) or not isinstance(delta_min, (int, Fraction)):
        raise TypeError("delta_min must be an exact rational number")
    delta_min = Fraction(delta_min)
    if not 0 < delta_min <= 2:
        raise ValueError("Require 0<delta_min<=2 and 2+delta_min<=c<=4")
    return Fraction(2_000_000)/delta_min**2


@cache
def derive():
    d = symbol.derive()
    u, c, _K = map(ra.Rational, ra.CONTEXT.gens())
    numerator, denominator = d["det_C"].numerator, d["det_C"].denominator
    if {int(power[2]) for power in denominator.to_dict()} != {3}:
        raise ValueError("The generic Cauchy denominator is not K³ times a background polynomial")
    if max(int(power[2]) for power in numerator.to_dict()) != 3:
        raise ValueError("The generic Cauchy numerator degree is not three")
    leading = ra.Rational(numerator).leading()
    expected = c**4*(1+u*u)**18*(c*(1+u*u)**4-2)**2*(u*u-1)/344064
    if leading != expected:
        raise ValueError("The literal leading numerator factorization failed")
    norms = [Fraction(0), Fraction(0), Fraction(0)]
    for raw_power, coefficient in numerator.to_dict().items():
        power = tuple(map(int, raw_power))
        if power[2] < 3:
            norms[power[2]] += abs(Fraction(str(coefficient)))*Fraction(1, 10)**power[0]*4**power[1]
    rounded = (16, 14, 1)
    if any(value >= bound for value, bound in zip(norms, rounded, strict=True)):
        raise ValueError("The continuous numerator coefficient majorant failed")
    # c>=2, 1-u²>=99/100, d>=1, c*d^4-2>=delta_min.
    leading_lower_factor = Fraction(99, 2_150_400)
    relative_error = Fraction(sum(rounded), 1)/(2_000_000*leading_lower_factor)
    if not relative_error < Fraction(1, 2):
        raise ValueError("The strict high-K denominator margin failed")
    return {"leading_numerator": leading, "expected_leading_numerator": expected,
            "coefficient_norms": tuple(norms), "rounded_coefficient_bounds": rounded,
            "leading_absolute_lower_factor_times_delta_squared": leading_lower_factor,
            "relative_determinant_error_upper": relative_error,
            "threshold_numerator": 2_000_000,
            "threshold_formula": "K >= 2000000/delta_min², 0<delta_min<=2, 2+delta_min<=c<=4",
            "interpretation": "Sufficient exact Cauchy-symbol threshold; not a physical cutoff"}
