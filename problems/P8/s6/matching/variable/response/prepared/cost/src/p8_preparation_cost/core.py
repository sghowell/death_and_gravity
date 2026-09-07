"""Literal action/chart identities and continuous source-map energy bounds."""
from functools import cache

import sympy as sp
from p8_variable_response import operator, source

from . import construction

Q = sp.Rational


@cache
def calibration():
    ell, r, root_r, mu_lo, mu_hi = Q(1, 100), Q(1, 50), Q(1, 7), Q(3), Q(13, 4)
    energy_remainder = (ell**2*28+ell**2*root_r*(28+7)
                        +ell**2*root_r*14*mu_hi+ell*root_r*160*r**2/mu_lo
                        +r*root_r*14/mu_lo+ell*r*200/mu_lo)
    momentum_derivative = (ell**2*Q(101, 100)+ell**2*root_r/100
                           +ell*root_r/(100*mu_lo)+ell*r*Q(101, 100)/mu_lo)
    input_squared = ell**4+ell**2*r/mu_lo**2
    exponent = Q(11, 20)
    exponential = 1+exponent+exponent**2/(2*(1-exponent/3))
    difference = 2*Q(1, 2000)*Q(1, 5000)*Q(3, 2)/sp.sqrt(ell)
    midpoint_sigma = Q(6, 10**6)
    return {
        "ell": ell, "radius_upper": r, "sqrt_radius_upper": root_r,
        "energy_remainder_entry_sum_upper": energy_remainder,
        "logarithmic_energy_exponent": exponent,
        "duration_one_exponential_upper": exponential,
        "momentum_generator_derivative_entry_sum_upper": momentum_derivative,
        "input_squared_upper": input_squared,
        "input_norm_upper": Q(1, 2000),
        "momentum_generator_derivative_upper": Q(1, 5000),
        "midpoint_to_any_K_distance": Q(3, 2),
        "source_map_K_difference_upper": difference,
        "midpoint_scaled_outer_singular_lower": midpoint_sigma,
        "uniform_scaled_outer_gramian_lower": (midpoint_sigma-difference)**2,
        "uniform_Frobenius_gramian_lower": (midpoint_sigma-difference)**2/4,
        "uniform_Frobenius_gramian_upper": Q(16, 25),
        "generic_minimum_source_L2_lower": Q(5, 4),
        "generic_minimum_source_L2_upper": Q(2*10**6, 3),
    }


@cache
def checks():
    d, old, c = construction.derive(), operator.derive(), calibration()
    substitutions = {old["u"]: d["u"], old["c"]: 2, old["kbar"]**2: d["K"]}
    mapping = operator.physical_map().subs(substitutions)
    canonical_generator = operator.first_order().subs(substitutions)
    physical_residual = (sp.diff(mapping, d["u"])+mapping*canonical_generator
                         -d["phase"]*mapping)
    actual_source = source.coefficients()
    forcing = sp.Matrix([0, actual_source["j_light"], 0,
                         actual_source["j_heavy"]]).subs(substitutions)
    source_residual = mapping*forcing-d["physical_input"]
    residuals = {
        "literal_punctured_spring": sp.factor(
            (old["mass"]*old["K_relative"]).subs(substitutions)-d["spring"]),
        "physical_g_friction": sp.factor(d["phase"][1, 1]+12*d["u"]/d["d"]),
        "physical_f_friction": sp.factor(d["phase"][3, 3]-12*d["u"]/d["d"]),
        "indicial_frequency": operator.MU**2+Q(1, 4)-10,
    }
    for index, entry in enumerate(physical_residual):
        residuals[f"full_physical_generator_{index}"] = sp.factor(entry)
    for index, entry in enumerate(source_residual):
        residuals[f"actual_physical_g_input_{index}"] = sp.factor(entry)
    r = sp.Symbol("r", positive=True)
    aa, cc, dc, ee, fc, bb, jl, jh = sp.symbols("A C d_c E f_c b j_l j_h", real=True)
    ell, mu = Q(1, 100), operator.MU
    w = operator.outer_map(r, -1)*sp.diag(1, 1/ell, 1, 1)
    generic = sp.Matrix([[0, 1, 0, 0], [-aa, 0, -cc, r*dc],
                         [0, 0, 0, 1], [-ee, r*fc, -10/r**2-bb, 0]])
    direct = w.inv()*(ell*generic*w+ell*sp.diff(w, r))
    expected = sp.Matrix([
        [0, 1, 0, 0],
        [-ell**2*aa, 0, -ell**2*sp.sqrt(r)*(cc+dc/2), -ell**2*sp.sqrt(r)*dc*mu],
        [0, 0, 0, -ell*mu/r],
        [ell*sp.sqrt(r)*ee/mu, -r**Q(3, 2)*fc/mu, ell*mu/r+ell*r*bb/mu, 0]])
    expected_input = sp.Matrix([0, ell**2*jl, 0, -ell*sp.sqrt(r)*jh/mu])
    for index, entry in enumerate(direct-expected):
        residuals[f"scaled_outer_clock_chain_{index}"] = sp.simplify(entry)
    for index, entry in enumerate(w.inv()*ell*sp.Matrix([0, jl, 0, jh])-expected_input):
        residuals[f"scaled_outer_actual_source_{index}"] = sp.simplify(entry)
    kappa = c["logarithmic_energy_exponent"]
    margins = {
        "sqrt_radius": Q(1, 49)-Q(1, 50),
        "d_fourth_continuous": Q(101, 100)-(1+Q(1, 50)**2)**4,
        "energy_remainder": Q(1, 20)-c["energy_remainder_entry_sum_upper"],
        "energy_exponent_equals_half_plus_margin": kappa-Q(1, 2),
        "propagator_norm": 2-c["duration_one_exponential_upper"],
        "momentum_derivative": Q(1, 5000)-c["momentum_generator_derivative_entry_sum_upper"],
        "source_norm_squared": Q(1, 2000)**2-c["input_squared_upper"],
        "singular_value_perturbation": c["midpoint_scaled_outer_singular_lower"]
        -c["source_map_K_difference_upper"],
    }
    if any(value != 0 for value in residuals.values()):
        raise ValueError("A literal full-action or scaled-chart identity failed")
    if any(value.is_positive is not True for value in margins.values()):
        raise ValueError("A continuous source-map energy margin failed")
    return {"residuals": residuals, "strict_margins": margins}

@cache
def coefficient_bridges():
    """Exact rational coefficients versus Arb jet enclosures, not trajectory samples."""
    from fractions import Fraction

    from flint import arb, ctx

    from . import gramian

    d = construction.derive()
    matrices = [d["scaled_phase"].applyfunc(
        lambda entry, n=n: sp.factor(sp.diff(entry, d["u"], n)*d["ell"]**n/sp.factorial(n)))
        for n in range(4)]
    old_precision, old_cap = ctx.prec, ctx.cap
    rows = []
    try:
        ctx.prec, ctx.cap = 256, 5
        for x in (Q(0), Q(1, 3), Q(2, 3), Q(1)):
            actual = gramian._generator(arb(int(x.p))/int(x.q), 3)
            count = 0
            for n, matrix in enumerate(matrices):
                reference = matrix.subs({d["u"]: d["left"]+d["ell"]*x, d["K"]: Q(5, 2)})
                for i in range(4):
                    for j in range(4):
                        value = reference[i, j]
                        if not isinstance(value, sp.Rational):
                            raise TypeError("An exact midpoint generator coefficient is not rational")
                        fraction = Fraction(int(value.p), int(value.q))
                        lower = gramian._fraction(actual[n][i, j].lower())
                        upper = gramian._fraction(actual[n][i, j].upper())
                        if not lower <= fraction <= upper:
                            raise ValueError("A literal rational-to-Arb coefficient bridge failed")
                        count += 1
            rows.append({"x": x, "K": Q(5, 2), "Taylor_orders": 4, "comparisons": count})
    finally:
        ctx.prec, ctx.cap = old_precision, old_cap
    return rows
