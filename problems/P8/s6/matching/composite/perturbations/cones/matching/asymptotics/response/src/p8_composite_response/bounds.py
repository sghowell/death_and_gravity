"""Explicit finite-epsilon continuation and source-response error majorants.

Large constants are stored as exponents: ``2**(-N)`` is an exact rational
number, but materializing its denominator is neither needed nor useful.
There are no floating-point estimates in this certificate interface.
"""

from fractions import Fraction
from functools import cache

from . import model
from .intervals import Interval, Jet, exact

BOX = {"epsilon": (Fraction(0), Fraction(1, 10000)),
       "u": (Fraction(0), Fraction(1)),
       "z": (Fraction(1, 2), Fraction(12)),
       "R": (Fraction(1, 10), Fraction(3)),
       "q": (Fraction(1, 256), Fraction(1, 64))}


def ceiling(value):
    value = exact(value)
    return -((-value.numerator)//value.denominator)


def log2_majorant(value):
    """Nonnegative integer n with exact nonnegative value <=2**n."""
    value = exact(value)
    if value < 0:
        raise ValueError("A majorant must be nonnegative")
    n = max(1, ceiling(value))
    return (n-1).bit_length()


def integer_enclosure(value):
    return {"lower": str(value.lo), "upper": str(value.hi),
            "absolute_integer_majorant": ceiling(value.magnitude())}


@cache
def enclosure():
    eps = Jet.variable(Interval(*BOX["epsilon"]), 0)
    z = Jet.variable(Interval(*BOX["z"]), 1)
    root = Jet.variable(Interval(*BOX["R"]), 2)
    u, q = (Jet.constant(Interval(*BOX[key])) for key in ("u", "q"))
    d = model.expressions(eps, u, z, root, q)
    d["initial_root"] = model.initial_root(eps)
    return d


@cache
def build():
    d = enclosure()
    positive = ("X_bar", "CD_denominator", "lapse_numerator", "lapse_denominator",
                "root_denominator", "kappa", "A_g", "A_f", "C_lock", "null_over_e")
    if any(d[key].value.lo <= 0 for key in positive):
        raise ValueError("The outward coefficient box lost a required positive denominator")
    if not d["e"].value.hi < Fraction(1, 100):
        raise ValueError("The exact e<1/100 constraint failed")
    initial_slope = d["initial_root"].gradient[0].magnitude()
    if not initial_slope < 3:
        raise ValueError("The finite-family initial root offset exceeds 3epsilon")

    flow = [d["zprime"], d["Rprime"]]
    kb = ceiling(max(value.derivative_norm((1, 2)) for value in flow))
    bb = ceiling(max(value.gradient[0].magnitude() for value in flow))
    # ||b_eps-b_0|| <=(3+bb) exp(kb) eps <=2**background_exponent eps.
    bgexp = log2_majorant(3+bb)+2*kb

    matrix, forcing = d["matrix"], d["forcing"]
    kp = ceiling(max(sum(entry.value.magnitude() for entry in row) for row in matrix))
    bp = ceiling(max(entry.value.magnitude() for entry in forcing))
    da = ceiling(max(sum(entry.derivative_norm() for entry in row) for row in matrix))
    db = ceiling(max(entry.derivative_norm() for entry in forcing))
    # p<=1; same bound for the limiting and finite zero-data phase columns.
    phase_exp = log2_majorant(bp)+2*kp
    dxexp = log2_majorant(da+db)+bgexp+phase_exp+2*kp
    c = ceiling(max(1, d["output_g"].value.magnitude()+d["output_f"].value.magnitude()))
    dc = ceiling(d["output_g"].derivative_norm()+d["output_f"].derivative_norm())
    out_exp = 1+max(log2_majorant(c)+dxexp, log2_majorant(dc)+bgexp+phase_exp)

    # Prepared source-free data have h_g=h_f and their u derivatives equal,
    # hence H=H'=0 in the exact relative canonical coordinate initially.
    # |L0|,|Lprime0|<=1. Include the nonconstant initial coefficients.
    initial_difference = ceiling(max(1, d["A_g"].value.magnitude(),
                                    3*d["A_f"].derivative_norm(),
                                    3*d["C_lock"].derivative_norm()))
    initial_size = ceiling(max(1, *(d[key].value.magnitude()
                                   for key in ("A_g", "A_f", "C_lock"))))
    prep_phase_exp = log2_majorant(initial_size)+2*kp
    prep_dx_exp = 1+max(log2_majorant(initial_difference)+2*kp,
                        log2_majorant(da)+bgexp+4*kp)
    prep_output_exp = 1+max(log2_majorant(c)+prep_dx_exp,
                           log2_majorant(dc)+bgexp+2*kp)
    # One more factor2 bounds full-minus-locked, both compared with their
    # identical limit. epsilon<=2**(-N) implies all advertised error gates.
    eps_exp = max(14, bgexp+7, out_exp+10, prep_output_exp+11)
    amplitude_exp = max(phase_exp, prep_phase_exp)+7
    exact_margins = {
        "parameter_box": Fraction(1, 10000)-Fraction(1, 1 << 14),
        "z_lower_tube": Fraction(1)-Fraction(1, 128)-BOX["z"][0],
        "z_upper_tube": BOX["z"][1]-Fraction(11)-Fraction(1, 128),
        "R_lower_tube": Fraction(2, 11)-Fraction(1, 128)-BOX["R"][0],
        "R_upper_tube": BOX["R"][1]-Fraction(2)-Fraction(1, 128),
        "error_gate": Fraction(1, 600)-Fraction(1, 1024),
        "small_metric": Fraction(1, 100)-Fraction(1, 128),
    }
    if not all(value > 0 for value in exact_margins.values()):
        raise ValueError("A continuation, error or amplitude margin failed")
    if not (eps_exp >= bgexp+7 and eps_exp >= out_exp+10
            and eps_exp >= prep_output_exp+11 and amplitude_exp >= phase_exp+7):
        raise ValueError("An exponent gate is inconsistent")

    records = {}
    for key, value in d.items():
        if isinstance(value, Jet):
            records[key] = {"value": integer_enclosure(value.value),
                            "partials_epsilon_z_R": [integer_enclosure(x) for x in value.gradient]}
    return {
        "box": {key: [str(x) for x in value] for key, value in BOX.items()},
        "rounding": "Exact Fraction operations rounded OUTWARD after each operation to 32 binary places; integer-isqrt radical enclosure",
        "positive_box_coefficients": list(positive),
        "records": records,
        "background_state_Lipschitz": kb, "background_parameter_derivative": bb,
        "initial_R_offset_majorant": 3, "background_error_power_of_two": bgexp,
        "phase_matrix_norm": kp, "phase_source_norm": bp,
        "phase_matrix_parameter_state_derivative": da,
        "phase_source_parameter_state_derivative": db,
        "phase_bound_power_of_two": phase_exp,
        "phase_error_power_of_two": dxexp,
        "output_norm": c, "output_parameter_state_derivative": dc,
        "output_error_power_of_two": out_exp,
        "prepared_initial_difference": initial_difference,
        "prepared_initial_size": initial_size,
        "prepared_phase_bound_power_of_two": prep_phase_exp,
        "prepared_phase_error_power_of_two": prep_dx_exp,
        "prepared_output_error_power_of_two": prep_output_exp,
        "epsilon_range": {"strict_lower": "0", "inclusive_upper": f"2^(-{eps_exp})",
                          "negative_binary_exponent": eps_exp},
        "sigma_amplitude_range": {"inclusive_lower": "0", "inclusive_upper": f"2^(-{amplitude_exp})",
                                  "negative_binary_exponent": amplitude_exp},
        "rational_margins": {key: str(value) for key, value in exact_margins.items()},
    }
