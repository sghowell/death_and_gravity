"""Independent Fraction whole-box proof for the desingularized own-f branch.

No binary floats, SymPy, Arb, or frozen implementation modules enter this
engine.  The delta=0 face belongs only to the regularized equations: it is
not a member of the literal VARIABLE action family.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction


def rational(value):
    """Admit exact integers/Fractions only; in particular reject bool/float."""
    if type(value) is int or isinstance(value, Fraction):
        return Fraction(value)
    raise TypeError("an exact int or Fraction is required")


@dataclass(frozen=True, slots=True)
class Interval:
    """Closed rational interval with inclusion-preserving arithmetic."""

    lo: Fraction
    hi: Fraction | None = None

    def __post_init__(self):
        low = rational(self.lo)
        high = low if self.hi is None else rational(self.hi)
        if low > high:
            raise ValueError("interval endpoints are reversed")
        object.__setattr__(self, "lo", low)
        object.__setattr__(self, "hi", high)

    def __add__(self, other):
        other = as_interval(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -as_interval(other)

    def __rsub__(self, other):
        return as_interval(other) + -self

    def __mul__(self, other):
        other = as_interval(other)
        products = (self.lo * other.lo, self.lo * other.hi,
                    self.hi * other.lo, self.hi * other.hi)
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def inverse(self):
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError("a reciprocal interval may not contain zero")
        return Interval(1 / self.hi, 1 / self.lo)

    def __truediv__(self, other):
        return self * as_interval(other).inverse()

    def __rtruediv__(self, other):
        return as_interval(other) * self.inverse()

    def __pow__(self, exponent):
        if type(exponent) is not int:
            raise TypeError("interval powers require an exact integer exponent")
        if exponent < 0:
            return self.inverse() ** (-exponent)
        if exponent == 0:
            return Interval(1)
        endpoints = (self.lo ** exponent, self.hi ** exponent)
        low = 0 if exponent % 2 == 0 and self.lo <= 0 <= self.hi else min(endpoints)
        return Interval(low, max(endpoints))

    def contains(self, value):
        value = rational(value)
        return self.lo <= value <= self.hi

    def subset_of(self, other):
        other = as_interval(other)
        return other.lo <= self.lo and self.hi <= other.hi

    @property
    def max_abs(self):
        return max(abs(self.lo), abs(self.hi))


def as_interval(value):
    return value if isinstance(value, Interval) else Interval(value)


def outward(value, denominator=1_000_000):
    """Enclose by multiples of 1/denominator, including negative endpoints."""
    if type(denominator) is not int:
        raise TypeError("the display denominator must be an integer")
    if denominator <= 0:
        raise ValueError("the display denominator must be positive")
    value = as_interval(value)
    low = (value.lo.numerator * denominator) // value.lo.denominator
    high = -((-value.hi.numerator * denominator) // value.hi.denominator)
    return Interval(Fraction(low, denominator), Fraction(high, denominator))


def fixed_box():
    return {"v": Interval(0, Fraction(1, 10_000)),
            "delta": Interval(0, Fraction(1, 100)),
            "zeta": Interval(11, 13)}


def _box(box):
    domain = fixed_box()
    if box is None:
        return domain
    if not isinstance(box, Mapping):
        raise TypeError("the box must be a mapping of named exact intervals")
    if set(box) != set(domain):
        raise ValueError("the box needs exactly v, delta, and zeta")
    for name, value in box.items():
        if not isinstance(value, Interval):
            raise TypeError("box entries must be Interval instances")
        if not value.subset_of(domain[name]):
            raise ValueError(f"{name} lies outside the proved whole box")
    return dict(box)


def enclose(box=None):
    """Continuous bounds, including all partials needed for implicit lapse.

    v and zeta derivatives below hold c=2+delta fixed.  The cube-root
    bracket is validated against b_cubed rather than assumed or rounded.
    The same bracket is sufficient for any admitted sub-box.
    """
    box = _box(box)
    v, delta, zeta = (box[name] for name in ("v", "delta", "zeta"))
    c, d = 2 + delta, 1 + v
    J = c * d**4 * (1 - 7*v) + 12*v
    J_v = c * d**3 * (-3 - 35*v) + 12
    J_vv = c * d**2 * (-44 - 140*v)
    Q = 32 * (1 - v) / (c * d**14)
    Q_v_over_Q = -1/(1 - v) - 14/d
    R = d**8 * J / (8*c*(1 - v))
    L = 8/d + J_v/J + 1/(1 - v)
    L_v = -8/d**2 + J_vv/J - (J_v/J)**2 + 1/(1 - v)**2
    D = c - 2/d**4
    D_v, D_vv = 8/d**5, -40/d**6
    H = D*(J_v/J - 6/d) - D_v
    H_v = (D_v*(J_v/J - 6/d)
           + D*(J_vv/J - (J_v/J)**2 + 6/d**2) - D_vv)
    F = 2*(v*zeta*H - L)/3
    F_v = 2*(zeta*(H + v*H_v) - L_v)/3
    F_zeta = 2*v*H/3
    A = 1 - v*D*zeta
    b_cubed = A/R
    b = Interval(Fraction(199, 100), Fraction(201, 100))
    if not b.lo**3 < b_cubed.lo or not b_cubed.hi < b.hi**3:
        raise ArithmeticError("the rational cube-root bracket was not proved")
    b_v = b*(-zeta*(D + v*D_v)/A - L)/3
    b_zeta = -b*v*D/(3*A)
    T = 3*b*F**2/(2*Q)
    T_v = 3*(b_v*F**2 + 2*b*F*F_v - b*F**2*Q_v_over_Q)/(2*Q)
    T_zeta = 3*(b_zeta*F**2 + 2*b*F*F_zeta)/(2*Q)
    zeta_v = T_v/(1 - T_zeta)
    b_total_v = b_v + b_zeta*zeta_v
    N = 2*b_total_v/F
    lapse_implicit_correction = 2*b_zeta*zeta_v/F
    return {"v": v, "delta": delta, "zeta": zeta, "c": c, "d": d,
            "J": J, "J_v": J_v, "J_vv": J_vv, "Q": Q,
            "Q_v_over_Q": Q_v_over_Q, "R": R, "L": L, "L_v": L_v,
            "D": D, "D_v": D_v, "D_vv": D_vv, "H": H, "H_v": H_v,
            "F": F, "F_v": F_v, "F_zeta": F_zeta, "A": A,
            "b_cubed": b_cubed, "b": b, "b_v": b_v, "b_zeta": b_zeta,
            "T": T, "T_v": T_v, "T_zeta": T_zeta, "zeta_v": zeta_v,
            "b_total_v": b_total_v, "N": N,
            "lapse_implicit_correction": lapse_implicit_correction}


def calibration():
    """Concise exact enclosures and strictly positive proof margins."""
    values = enclose()
    names = ("J", "Q", "R", "A", "b_cubed", "F", "b_v", "b_total_v",
             "T", "T_v", "T_zeta", "zeta_v", "N")
    bounds = {name: outward(values[name]) for name in names}
    b = values["b"]
    margins = {
        "cube_root_lower": bounds["b_cubed"].lo - b.lo**3,
        "cube_root_upper": b.hi**3 - bounds["b_cubed"].hi,
        "self_map_lower": bounds["T"].lo - 11,
        "self_map_upper": 13 - bounds["T"].hi,
        "contraction_below_one_hundredth": Fraction(1, 100) - bounds["T_zeta"].max_abs,
        "implicit_derivative_denominator": 1 - bounds["T_zeta"].hi,
        "negative_F": -bounds["F"].hi,
        "negative_total_scale_derivative": -bounds["b_total_v"].hi,
        "scale_above_19over10": b.lo - Fraction(19, 10),
        "scale_below_21over10": Fraction(21, 10) - b.hi,
        "lapse_above_19over10": bounds["N"].lo - Fraction(19, 10),
        "lapse_below_21over10": Fraction(21, 10) - bounds["N"].hi,
    }
    margins.update({f"positive_{name}": bounds[name].lo for name in ("J", "Q", "R", "A")})
    return {"domain": fixed_box(), "b": b, "outward_enclosures": bounds,
            "strict_margins": margins, "u_half_width": Fraction(1, 100),
            "contraction_upper": Fraction(1, 100),
            "literal_delta_positive": True,
            "delta_zero_is_only_desingularized_limit": True,
            "full_parent_solution": False, "tensor_Green_inverse_selected": False}


def checks():
    """Return Boolean proof checks; every entry must be True."""
    data, values = calibration(), enclose()
    results = {name: value > 0 for name, value in data["strict_margins"].items()}
    results.update({f"outward_contains_{name}": values[name].subset_of(bound)
                    for name, bound in data["outward_enclosures"].items()})
    results["D_nonnegative_on_closed_box"] = values["D"].lo >= 0
    results["all_checks_are_boolean"] = all(type(value) is bool for value in results.values())
    return results


def _encode(value):
    if isinstance(value, Interval):
        return [str(value.lo), str(value.hi)]
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {key: _encode(item) for key, item in value.items()}
    return value


def report():
    """Read-only JSON-compatible exact data; no report file is written."""
    data = calibration()
    data["checks"] = checks()
    return _encode(data)
