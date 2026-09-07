"""Separately constructed Fraction norm and polynomial checks, no SymPy."""

from fractions import Fraction
from math import comb

F = Fraction
R_MAX = F(1, 1000)
DEFICIT_MAX = F(1, 100)


def calibration(r=R_MAX, deficit=DEFICIT_MAX):
    for value in (r, deficit):
        if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
            raise TypeError("Only integers and Fractions are exact inputs")
    r, deficit = F(r), F(deficit)
    if not 0 < r <= F(1, 1000) or not 0 <= deficit <= F(1, 100):
        raise ValueError("Outside the declared pulse/radius domain")
    moment = F(5, 28)-deficit-4*r*r
    observable = 1-9*r*r
    source_times_map_over_mu = F(8, 5)/F(13, 4)
    amplitude = source_times_map_over_mu*observable*moment
    extra = 2*F(101, 75)*(10+60*r*r)*5
    return {"r": r, "deficit_over_r": deficit,
            "moment_over_r_three_halves": moment,
            "amplitude_over_r_squared_lower": amplitude,
            "additional_feedback_error_over_r_fourth": extra,
            "simplified_response_error_over_r_fourth": 1000+extra,
            "rounded_amplitude_over_r_squared_lower": F(2, 25),
            "rounded_response_error_over_r_fourth": F(1200),
            "phase_separation_over_r_squared_lower": 2*F(2, 25)-2*1200*r*r,
            "claimed_separation_over_r_squared": F(3, 20),
            "smooth_transition_width_over_r": F(1, 800),
            "smooth_pulse_deficit_over_r_upper": 4*F(1, 800)}


def checks():
    # Coefficient proof of (1+t)^-9>=1-9t for every t>=0:
    # multiply by (1+t)^9; all remaining nonconstant coefficients of
    # 1-(1-9t)(1+t)^9 are nonnegative.
    polynomial = [F(0)]
    for j in range(1, 11):
        value = 9*comb(9, j-1)-(comb(9, j) if j <= 9 else 0)
        assert value >= 0
        polynomial.append(F(value))
    c = calibration()
    assert c["amplitude_over_r_squared_lower"] > F(2, 25)
    assert c["simplified_response_error_over_r_fourth"] < 1200
    assert c["phase_separation_over_r_squared_lower"] > F(3, 20)
    return {"inverse_ninth_nonnegative_coefficients": polynomial,
            "rectangular_moment_lower": (1-F(3, 8))/F(7, 2),
            "smooth_deficit_fraction": 4*F(1, 800),
            "all_continuous_margins_at_max_radius": True}
