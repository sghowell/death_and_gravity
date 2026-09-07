"""Exact finite-delta own-f coefficient bounds inherited from frozen S6.33.

The background enclosure is explicitly imported, not claimed as a second
independent proof.  The domain inclusion, polynomial gap identity, and
monotone coefficient inequalities below use separate Fraction arithmetic.
No binary floats or numerical cube roots enter the proof.
"""

from fractions import Fraction

from p8_exact_stationary import intervals as parent

DELTA_MAX = Fraction(1, 625)
WIDTH = Fraction(1, 4)
V_MAX = Fraction(1, 10_000)
K_MIN, K_MAX = Fraction(19, 5), Fraction(21, 5)
S_MIN, S_MAX = Fraction(42), Fraction(65)


def rational(value):
    if type(value) is int or isinstance(value, Fraction):
        return Fraction(value)
    raise TypeError("an exact int or Fraction is required")


def _point(delta, x, extension):
    delta, x = rational(delta), rational(x)
    if type(extension) is not bool:
        raise TypeError("extension must be an explicit Boolean")
    if not -WIDTH <= x <= 0:
        raise ValueError("x lies outside [-1/4,0]")
    if not 0 <= delta <= DELTA_MAX or delta == 0 and not extension:
        raise ValueError("literal delta must be in (0,1/625]; zero is extension-only")
    return delta, x, delta*x*x


def polynomial(coefficients, value):
    """Ascending exact coefficients, evaluated by a separate Horner rule."""
    value = rational(value)
    coefficients = tuple(rational(item) for item in coefficients)
    if not coefficients:
        raise ValueError("a polynomial needs at least one coefficient")
    result = Fraction(0)
    for coefficient in reversed(coefficients):
        result = result*value + coefficient
    return result


def gap_polynomials():
    """Exact positive-coefficient proof valid for every v>=0.

    P(v)=4+6v+4v²+v³; d⁴-1=vP; 4d⁴-P has no constant term and
    strictly positive other coefficients.  Thus 0<P/d⁴<=4 for v>=0.
    """
    d4 = (1, 4, 6, 4, 1)
    P = (4, 6, 4, 1)
    shifted = (0, *P)
    identity = tuple(d4[i] - (1 if i == 0 else 0) - shifted[i] for i in range(5))
    padded = (*P, 0)
    upper = tuple(4*d4[i] - padded[i] for i in range(5))
    return {"d_fourth": d4, "P": P, "d4_minus_one_minus_vP": identity,
            "four_d4_minus_P": upper}


def denominator_ratio(delta, x, *, extension=False):
    """Exact desingularized D/delta; delta=0 needs explicit opt-in.

    Unlike evaluating D and dividing afterward, this formula remains
    defined at the excluded literal action endpoint.
    """
    _, x, v = _point(delta, x, extension)
    return 1 + 2*x*x*polynomial(gap_polynomials()["P"], v)/(1 + v)**4


def denominator_identity(delta, x, *, extension=False):
    """Return D-delta*(desingularized ratio), exactly zero on the domain."""
    delta, x, v = _point(delta, x, extension)
    D = 2 + delta - 2/(1 + v)**4
    return D - delta*denominator_ratio(delta, x, extension=extension)


def inherited_enclosures():
    """Use only the explicit sub-box justified by v=delta*x².

    These are frozen S6.33 background enclosures, including the actual
    implicit lapse derivative.  No replacement by a partial lapse is made.
    """
    return parent.enclose({"v": parent.Interval(0, V_MAX),
                           "delta": parent.Interval(0, DELTA_MAX),
                           "zeta": parent.Interval(11, 13)})


def _outward(low, high, denominator=1_000_000):
    """Separate exact directed rounding of a pair of rational endpoints."""
    low, high = rational(low), rational(high)
    if low > high:
        raise ValueError("coefficient endpoints are reversed")
    if type(denominator) is not int:
        raise TypeError("rounding denominator must be an integer")
    if denominator <= 0:
        raise ValueError("rounding denominator must be positive")
    lower = (low.numerator*denominator)//low.denominator
    upper = -((-high.numerator*denominator)//high.denominator)
    return Fraction(lower, denominator), Fraction(upper, denominator)


def _raw_coefficients():
    data = inherited_enclosures()
    # k=b³/N is rational in the regularized background: use the tighter
    # b_cubed enclosure directly, not the cube of its coarse root bracket.
    k_lower = data["b_cubed"].lo/data["N"].hi
    k_upper = data["b_cubed"].hi/data["N"].lo
    # 1<=D/delta<=3/2, independently proved by gap_polynomials().
    s_lower = 2*data["Q"].lo*data["b"].lo/Fraction(3, 2)
    s_upper = 2*data["Q"].hi*data["b"].hi
    return {"k_lower": k_lower, "k_upper": k_upper,
            "s_lower": s_lower, "s_upper": s_upper}


def coefficient_box():
    """Concise rational coefficient enclosures and the safe theorem bounds."""
    raw = _raw_coefficients()
    kl, ku = _outward(raw["k_lower"], raw["k_upper"])
    sl, su = _outward(raw["s_lower"], raw["s_upper"])
    return {"k_lower": kl, "k_upper": ku, "s_lower": sl, "s_upper": su,
            "k_min": K_MIN, "k_max": K_MAX, "s_min": S_MIN, "s_max": S_MAX,
            "delta_max": DELTA_MAX, "width": WIDTH, "v_max": V_MAX,
            "x_min": -WIDTH, "x_max": Fraction(0),
            "gap_min": Fraction(1), "gap_max": Fraction(3, 2)}


def calibration():
    data = coefficient_box()
    source = inherited_enclosures()
    inherited = {name: _outward(source[name].lo, source[name].hi)
                 for name in ("Q", "b", "b_cubed", "N")}
    margins = {
        "k_above_19over5": data["k_lower"] - K_MIN,
        "k_below_21over5": K_MAX - data["k_upper"],
        "s_above_42": data["s_lower"] - S_MIN,
        "s_below_65": S_MAX - data["s_upper"],
        "delta_below_parent_cap": Fraction(1, 100) - DELTA_MAX,
        "positive_parent_Q": inherited["Q"][0],
        "positive_parent_b_cubed": inherited["b_cubed"][0],
        "positive_parent_N": inherited["N"][0],
    }
    return {**data, "strict_margins": margins,
            "inherited_rational_enclosures": inherited,
            "gap_polynomials": gap_polynomials(),
            "domain_v_identity": DELTA_MAX*WIDTH**2 - V_MAX,
            "spatial_fourier_mode": "homogeneous tensor probe",
            "physical_time": "T=tau*sqrt(delta)*x",
            "flux_normalization": "K_f=M²*k; s=delta*tau²*nu/M²",
            "equation": "(k*Q_x)_x+s*(Q-q)=0",
            "literal_delta_zero_allowed": False,
            "x_zero_is_an_allowed_literal_observation": True,
            "full_parent_background_claim": False,
            "fixed_physical_low_frequency_band_claim": False}


def checks():
    """All values are Boolean; strict inequalities use exact rationals."""
    data, raw = calibration(), _raw_coefficients()
    polynomials = data["gap_polynomials"]
    results = {name: value > 0 for name, value in data["strict_margins"].items()}
    results.update({
        "domain_v_identity": data["domain_v_identity"] == 0,
        "gap_numerator_positive": all(value > 0 for value in polynomials["P"]),
        "gap_identity": all(value == 0 for value in polynomials["d4_minus_one_minus_vP"]),
        "gap_upper_constant_zero": polynomials["four_d4_minus_P"][0] == 0,
        "gap_upper_other_coefficients_positive": all(value > 0 for value in polynomials["four_d4_minus_P"][1:]),
        "gap_max_from_width": 1 + 8*WIDTH**2 == Fraction(3, 2),
        "outward_k_lower": data["k_lower"] <= raw["k_lower"],
        "outward_k_upper": data["k_upper"] >= raw["k_upper"],
        "outward_s_lower": data["s_lower"] <= raw["s_lower"],
        "outward_s_upper": data["s_upper"] >= raw["s_upper"],
    })
    return results
