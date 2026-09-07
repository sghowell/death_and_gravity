"""Independent rational norm and matrix-positivity replay, no Arb or SymPy."""
from fractions import Fraction as Q


def norm_constants():
    ell, radius, root, mulo, muhi = Q(1, 100), Q(1, 50), Q(1, 7), Q(3), Q(13, 4)
    residual = sum((ell*ell*28, ell*ell*root*35, ell*ell*root*14*muhi,
                    ell*root*160*radius*radius/mulo,
                    radius*root*14/mulo, ell*radius*200/mulo), Q(0))
    derivative = (ell*ell*(Q(101, 100)+root/100)
                  +ell/mulo*(root/100+radius*Q(101, 100)))
    kappa = Q(1, 2)+Q(1, 20)
    difference = Q(2)*Q(1, 2000)*Q(3, 2)/(5000*Q(1, 10))
    return {
        "ell": ell, "radius_upper": radius, "sqrt_radius_upper": root,
        "energy_remainder_entry_sum_upper": residual,
        "logarithmic_energy_exponent": kappa,
        "duration_one_exponential_upper": 1+kappa+kappa*kappa*3/(6-2*kappa),
        "momentum_generator_derivative_entry_sum_upper": derivative,
        "input_squared_upper": ell**4+ell**2*radius/mulo**2,
        "input_norm_upper": Q(1, 2000),
        "momentum_generator_derivative_upper": Q(1, 5000),
        "midpoint_to_any_K_distance": Q(3, 2),
        "source_map_K_difference_upper": difference,
        "midpoint_scaled_outer_singular_lower": Q(6, 10**6),
        "uniform_scaled_outer_gramian_lower": Q(9, 10**12),
        "uniform_Frobenius_gramian_lower": Q(9, 4*10**12),
        "uniform_Frobenius_gramian_upper": Q(16, 25),
        "generic_minimum_source_L2_lower": Q(5, 4),
        "generic_minimum_source_L2_upper": Q(2*10**6, 3),
    }


class Interval:
    """Unrounded exact rational intervals for the 4x4 Schur replay."""

    def __init__(self, low, high=None):
        if isinstance(low, bool) or not isinstance(low, (Q, int)):
            raise TypeError("Exact interval endpoints only")
        high = low if high is None else high
        if isinstance(high, bool) or not isinstance(high, (Q, int)):
            raise TypeError("Exact interval endpoints only")
        self.lo, self.hi = Q(low), Q(high)
        if self.lo > self.hi:
            raise ValueError("Reversed exact interval")

    def __sub__(self, other):
        return Interval(self.lo-other.hi, self.hi-other.lo)

    def __mul__(self, other):
        products = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(products), max(products))

    def __truediv__(self, other):
        if other.lo <= 0 <= other.hi:
            raise ValueError("Zero-containing Schur denominator")
        return self*Interval(1/other.hi, 1/other.lo)


def positive_gramian(matrix):
    if (not isinstance(matrix, (tuple, list)) or len(matrix) != 4
            or any(not isinstance(row, (tuple, list)) or len(row) != 4 for row in matrix)):
        raise ValueError("Require four by four exact interval matrix")
    values = [[Interval(*pair) for pair in row] for row in matrix]
    # Exact symmetry is known from the matrix ODE; require the independently
    # enclosed opposite entries to have a nonempty common interval.
    for i in range(4):
        for j in range(i+1, 4):
            lo, hi = max(values[i][j].lo, values[j][i].lo), min(values[i][j].hi, values[j][i].hi)
            values[i][j] = values[j][i] = Interval(lo, hi)
        values[i][i] = values[i][i]-Interval(Q(36, 10**12))
    pivots = []
    for j in range(4):
        pivot = values[j][j]
        if pivot.lo <= 0:
            raise ValueError("A strictly positive rational pivot was not proved")
        pivots.append((pivot.lo, pivot.hi))
        for i in range(j+1, 4):
            for k in range(j+1, 4):
                values[i][k] = values[i][k]-values[i][j]*values[j][k]/pivot
    return {"strict_positive_pivots": len(pivots),
            "pivot_lower_bounds": tuple(pair[0] for pair in pivots),
            "shift": Q(36, 10**12), "arithmetic": "exact unrounded Fraction intervals"}


def checks():
    c = norm_constants()
    if not (c["energy_remainder_entry_sum_upper"] < Q(1, 20)
            and c["duration_one_exponential_upper"] < 2
            and c["momentum_generator_derivative_entry_sum_upper"] < Q(1, 5000)
            and c["input_squared_upper"] < Q(1, 2000)**2):
        raise ValueError("An independent continuous norm inequality failed")
    if (c["midpoint_scaled_outer_singular_lower"]-c["source_map_K_difference_upper"])**2 != c["uniform_scaled_outer_gramian_lower"]:
        raise ValueError("An independent singular-value margin failed")
    if c["uniform_Frobenius_gramian_lower"]*c["generic_minimum_source_L2_upper"]**2 != 1:
        raise ValueError("An independent source cost inverse failed")
    return c
