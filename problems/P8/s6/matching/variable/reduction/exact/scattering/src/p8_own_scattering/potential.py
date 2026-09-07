"""Finite uniform canonical own-f pole bound, with exact interval jet proof.

The delta=0 corner is only the desingularized background limit.  The
bounded pole remainder is not jointly continuous there.  No derivative
bound on that combined remainder or fixed-window transfer is asserted.
"""

from fractions import Fraction

from p8_exact_stationary import intervals as parent

from .jets import constant, variable

Interval = parent.Interval
RADIUS = Fraction(1, 100)
DELTA_MAX = Fraction(1, 100)
AV_CAP = Fraction(129)
PUMP_CAP = Fraction(22)
D4_CAP = Fraction(1001, 1000)
POTENTIAL_CAP = Fraction(44)


def partial_jets(box=None, *, root_bracket=None):
    """All partials are at independent (v,zeta), holding c fixed.

    The frozen parent validates the entire input box and its positive-root
    chart.  A narrower root bracket must itself pass a cubed-endpoint test.
    """
    inherited = parent.enclose(box)
    v = variable(0, inherited["v"])
    zeta = variable(1, inherited["zeta"])
    c, d = constant(2 + inherited["delta"]), 1 + v
    J = c*d**4*(1 - 7*v) + 12*v
    J_v = c*d**3*(-3 - 35*v) + 12
    Q = 32*(1 - v)/(c*d**14)
    R = d**8*J/(8*c*(1 - v))
    L = 8/d + J_v/J + 1/(1 - v)
    D = c - 2/d**4
    H = D*(J_v/J - 6/d) - 8/d**5
    F = 2*(v*zeta*H - L)/3
    b_cubed = (1 - v*D*zeta)/R
    bracket = inherited["b"] if root_bracket is None else root_bracket
    b = b_cubed.cube_root(bracket)
    T = 3*b*F**2/(2*Q)
    return {"v": v, "delta": constant(inherited["delta"]), "zeta": zeta,
            "c": c, "d": d, "J": J, "J_v": J_v, "Q": Q, "R": R,
            "L": L, "D": D, "H": H, "F": F, "b_cubed": b_cubed,
            "b": b, "T": T}


def implicit_enclosures(box=None, *, root_bracket=None):
    """Differentiate zeta=T(v,zeta,c) through order three exactly."""
    jets = partial_jets(box, root_bracket=root_bracket)
    T, b, F = (jets[name] for name in ("T", "b", "F"))
    denominator = 1 - T.partial(0, 1)
    z1 = T.partial(1, 0)/denominator
    z2 = (T.partial(2, 0) + 2*T.partial(1, 1)*z1 + T.partial(0, 2)*z1**2)/denominator
    z3 = (T.partial(3, 0) + 3*T.partial(2, 1)*z1 + 3*T.partial(1, 2)*z1**2
          + T.partial(0, 3)*z1**3
          + 3*(T.partial(1, 1) + T.partial(0, 2)*z1)*z2)/denominator
    b1 = b.partial(1, 0) + b.partial(0, 1)*z1
    b2 = (b.partial(2, 0) + 2*b.partial(1, 1)*z1 + b.partial(0, 2)*z1**2
          + b.partial(0, 1)*z2)
    b3 = (b.partial(3, 0) + 3*b.partial(2, 1)*z1 + 3*b.partial(1, 2)*z1**2
          + b.partial(0, 3)*z1**3
          + 3*(b.partial(1, 1) + b.partial(0, 2)*z1)*z2 + b.partial(0, 1)*z3)
    f1 = F.partial(1, 0) + F.partial(0, 1)*z1
    f2 = (F.partial(2, 0) + 2*F.partial(1, 1)*z1 + F.partial(0, 2)*z1**2
          + F.partial(0, 1)*z2)
    return {"v": jets["v"].partial(0, 0), "delta": jets["delta"].partial(0, 0),
            "Q": jets["Q"].partial(0, 0), "Q_v": jets["Q"].partial(1, 0),
            "T_zeta": T.partial(0, 1), "implicit_denominator": denominator,
            "zeta_v": z1, "zeta_vv": z2, "zeta_vvv": z3,
            "b": b.partial(0, 0), "b_v": b1, "b_vv": b2, "b_vvv": b3,
            "F": F.partial(0, 0), "F_v": f1, "F_vv": f2}


def potential_enclosures(box=None, *, root_bracket=None):
    """A_v and the exact canonical pump; no derivative of V-16/Z is taken."""
    data = implicit_enclosures(box, root_bracket=root_bracket)
    b, b1, b2, b3 = (data[name] for name in ("b", "b_v", "b_vv", "b_vvv"))
    F, f1, f2 = (data[name] for name in ("F", "F_v", "F_vv"))
    A = 4*data["Q"]*b1/(F*b**2)
    Av = A*(data["Q_v"]/data["Q"] + b2/b1 - f1/F - 2*b1/b)
    ell = 3*b1/b - b2/b1 + f1/F
    ell_v = 3*(b2/b - (b1/b)**2) - (b3/b1 - (b2/b1)**2) + (f2/F - (f1/F)**2)
    pump = ell + 2*data["v"]*ell_v + data["v"]*ell**2
    return {**data, "A": A, "A_v": Av, "ell": ell, "ell_v": ell_v, "pump": pump}


def center_data(delta, *, extension=False):
    """Exact rational jets at v=0 on the selected fixed point.

    The known b=2 root is verified by the same cube-root guard.  No interval
    width or approximate limit evaluation is used for these center data.
    """
    delta = parent.rational(delta)
    if type(extension) is not bool:
        raise TypeError("extension must be explicitly Boolean")
    if not 0 <= delta <= DELTA_MAX or delta == 0 and not extension:
        raise ValueError("literal delta is in (0,1/100]; zero requires extension=True")
    c = 2 + delta
    zeta = 3*(c + 2)**2/(2*c)
    box = {"v": Interval(0), "delta": Interval(delta), "zeta": Interval(zeta)}
    values = potential_enclosures(box, root_bracket=Interval(2))
    if any(value.lo != value.hi for value in values.values()):
        raise ArithmeticError("a purported exact center jet has nonzero interval width")
    return {name: value.lo for name, value in values.items()}


def pole_polynomials():
    d4 = (1, 4, 6, 4, 1)
    E = (20, 40, 30, 8)
    P = (4, 6, 4, 1)
    E_upper = tuple(20*d4[i] - (E[i] if i < len(E) else 0) for i in range(5))
    z_minus_d = (0, 0, *E)
    difference = tuple((8*d4[i-1] if i >= 1 else 0)
                       - (2*P[i-1] if 1 <= i <= len(P) else 0)
                       - z_minus_d[i] for i in range(6))
    return {"E_numerator": E, "E_upper": E_upper,
            "Z_minus_D_identity": difference,
            "D_lower_delta_coefficients": (0, 4, 6, 4, 1),
            "D_lower_other_coefficients": (0, 0, 12, 8, 2)}


def directional_limit(ratio):
    """Corner limit for finite delta/v -> ratio>=0; not a literal c=2 action."""
    ratio = parent.rational(ratio)
    if ratio < 0:
        raise ValueError("the admitted directional ratio must be nonnegative")
    w = ratio + 8
    return 24 - 112/w + 320/w**2


def calibration():
    values = potential_enclosures()
    names = ("T_zeta", "zeta_v", "zeta_vv", "zeta_vvv", "b_v", "b_vv", "b_vvv",
             "F", "F_v", "F_vv", "A", "A_v", "ell", "ell_v", "pump")
    outward = {name: parent.outward(values[name]) for name in names}
    bound = D4_CAP*(max(Fraction(8), AV_CAP/8) + 5) + PUMP_CAP
    center = center_data(0, extension=True)
    margins = {
        "A_v_below_129": AV_CAP - outward["A_v"].max_abs,
        "pump_below_22": PUMP_CAP - outward["pump"].max_abs,
        "d4_upper": D4_CAP - (1 + RADIUS**2)**4,
        "potential_below_44": POTENTIAL_CAP - bound,
        "implicit_denominator_positive": values["implicit_denominator"].lo,
        "negative_b_v": -outward["b_v"].hi,
        "negative_F": -outward["F"].hi,
    }
    return {"outward_enclosures": outward, "strict_margins": margins,
            "radius": RADIUS, "delta_max": DELTA_MAX, "A_v_cap": AV_CAP,
            "pump_cap": PUMP_CAP, "d4_cap": D4_CAP,
            "potential_cap": POTENTIAL_CAP, "derived_potential_bound": bound,
            "pole_polynomials": pole_polynomials(), "corner_data": center,
            "vertical_corner_limit": 8 - center["ell"],
            "horizontal_corner_limit": center["A_v"]/8 + 5 - center["ell"],
            "combined_remainder_jointly_continuous": False,
            "combined_remainder_derivative_bound": False,
            "literal_delta_zero_allowed": False,
            "fixed_physical_window_transfer_proved": False,
            "low_frequency_EFT_or_UV_verdict": False}


def checks():
    """Boolean enclosure checks and exact corner/polynomial controls."""
    data = calibration()
    polynomials = data["pole_polynomials"]
    center = data["corner_data"]
    results = {name: value > 0 for name, value in data["strict_margins"].items()}
    results.update({
        "pole_difference_identity": all(value == 0 for value in polynomials["Z_minus_D_identity"]),
        "E_positive_coefficients": all(value > 0 for value in polynomials["E_numerator"]),
        "E_upper_constant_zero": polynomials["E_upper"][0] == 0,
        "E_upper_nonconstant_positive": all(value > 0 for value in polynomials["E_upper"][1:]),
        "D_lower_coefficients_nonnegative": all(value >= 0 for name in
             ("D_lower_delta_coefficients", "D_lower_other_coefficients") for value in polynomials[name]),
        "corner_A": center["A"] == 16,
        "corner_A_v": center["A_v"] == -48,
        "corner_log_k_v": center["ell"] == -16,
        "vertical_corner": data["vertical_corner_limit"] == 24,
        "horizontal_corner": data["horizontal_corner_limit"] == 15,
        "directional_formula": directional_limit(0) == data["horizontal_corner_limit"],
    })
    return results
