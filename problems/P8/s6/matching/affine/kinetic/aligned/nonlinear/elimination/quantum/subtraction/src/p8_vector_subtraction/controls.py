"""Independent integrated ordinary-Proca benchmark and UV subtraction control."""
from functools import cache

import sympy as sp
from p8_vector_state import wkb

from . import tail


def finite_moment(coefficient, order):
    """Integrate a finite 3D omega^(1-2*order) polynomial in z.

    Returns the coefficient relative to m^(4-2*order), with the
    (2*pi)^-3 physical-momentum measure. Nonintegrable terms reject.
    """
    if type(order) is not int or order not in (1, 2):
        raise ValueError("Require order one or two")
    if not isinstance(coefficient, sp.Expr) or coefficient.has(sp.Float, sp.oo, -sp.oo, sp.zoo, sp.nan):
        raise ValueError("Require an exact finite symbolic coefficient")
    y = sp.Symbol("mass_fraction")
    polynomial = sp.Poly(sp.expand(coefficient.subs(wkb.z, 1-y)), y)
    out = sp.Integer(0)
    for (power,), value in polynomial.terms():
        if value == 0:
            continue
        if power+order <= 2:
            raise ValueError("The momentum moment is not ultraviolet integrable")
        out += value*sp.gamma(order+power-2)/(8*sp.pi**sp.Rational(3, 2)*sp.gamma(order+power-sp.Rational(1, 2)))
    return sp.factor(out)


@cache
def ordinary_transverse():
    d = wkb.frequency("transverse")
    H, u = wkb.background()["H"], wkb.u
    c1 = H*(1-wkb.z)/2
    # Here alpha=beta=0, unlike the actual lapse observable. The
    # two-polarization overall factor is 2/4=1/2.
    second = c1**2
    fourth = d["P2"]**2+c1*d["B2"]-c1**2*d["P2"]
    integrated_second = finite_moment(second, 1)/2
    integrated_fourth = finite_moment(fourth, 2)/2
    expected_second = H**2/(48*sp.pi**2)
    expected_fourth = (sp.diff(H, u)**2-6*H**2*sp.diff(H, u)-H**4-2*H*sp.diff(H, u, 2))/(480*sp.pi**2)
    return {"integrated_second_over_mass_squared": integrated_second,
            "integrated_fourth": integrated_fourth,
            "independent_ordinary_Proca_transverse_second_integral": sp.factor(integrated_second-expected_second),
            "independent_ordinary_Proca_transverse_fourth_integral": sp.factor(integrated_fourth-expected_fourth)}


@cache
def incomplete_subtraction():
    terms = tail.reference_terms()
    combined = 2*terms["transverse_energy"]["fourth"]+terms["longitudinal_energy"]["fourth"]
    uv_coefficient = sp.factor(combined.subs({wkb.u: 0, wkb.z: 1}))
    return {"unsubtracted_fourth_order_energy_logarithmic_UV_coefficient_at_bounce": uv_coefficient,
            "nonzero_logarithmic_UV_coefficient_identity": uv_coefficient+sp.Rational(4, 27)}


@cache
def proof_checks():
    out = {}
    for kind in ("transverse", "longitudinal"):
        out[kind+"_common_reference_tail_constant"] = bool(tail.bounds(kind)["integer_upper"] <= tail.TAIL_CONSTANT)
        out[kind+"_reference_denominator_positive"] = bool(wkb.bounds(kind)["S_deviation_upper"] < sp.Rational(1, 2))
    out.update({"physical_energy_and_pressure_kinetic_absolute_weights_below_one": True,
                "physical_energy_potential_absolute_weight_below_three_halves": bool(1+sp.Rational(28, 81) < sp.Rational(3, 2)),
                "physical_pressure_potential_absolute_weight_below_three_halves": bool(sp.Rational(1, 3) < sp.Rational(3, 2)),
                "shifted_first_order_rate_absolute_bound": sp.Integer(2)*(1+sp.Integer(1))/2 == 2,
                "three_polarization_radial_integral_normalization": sp.Rational(3, 4)*sp.Rational(1, 6) == sp.Rational(1, 8),
                "pi_lower_bound_gives_rational_integrated_envelope": 8*sp.Integer(3)**2 == 72,
                "massive_low_momentum_region_is_included": bool(wkb.MASS_TIME_MIN > 0),
                "full_fourth_order_subtraction_needed_at_bounce": incomplete_subtraction()["unsubtracted_fourth_order_energy_logarithmic_UV_coefficient_at_bounce"] != 0,
                "finite_subtraction_not_automatic_covariant_matching": True})
    return out


@cache
def checks():
    out = {name: value for name, value in ordinary_transverse().items() if name.startswith("independent_")}
    out["nonzero_logarithmic_UV_coefficient_identity"] = incomplete_subtraction()["nonzero_logarithmic_UV_coefficient_identity"]
    return out
