"""Continuous exact-mode/WKB energy comparison, with all momenta integrated."""
from functools import cache

import sympy as sp
from p8_affine_retuned.bounds import exact

from . import wkb

RESIDUAL_CONSTANT = sp.Integer(2_000_000)
AMAX = sp.Rational(25, 16)


@cache
def oscillator_identities():
    Wd, rate, phase = sp.symbols("W_dot physical_rate phase", real=True)
    positiveW = sp.Symbol("positive_W", positive=True)
    f = sp.exp(-sp.I*phase)/sp.sqrt(2*positiveW)
    fd = (-Wd/(2*positiveW)-sp.I*positiveW)*f
    wronskian = sp.expand(f*sp.conjugate(fd)-sp.conjugate(f)*fd)
    rho = sp.Symbol("rho", real=True)
    B = sp.Matrix([[-sp.I*rho*sp.conjugate(f)*f, -sp.I*rho*sp.conjugate(f)**2],
                   [sp.I*rho*f**2, sp.I*rho*f*sp.conjugate(f)]])
    metric = sp.diag(1, -1)
    pseudo_unitarity = (B.conjugate().T*metric+metric*B).applyfunc(sp.simplify)
    return {"unit_canonical_Wronskian": sp.simplify(wronskian-sp.I),
            "variation_of_constants_preserves_CCR": pseudo_unitarity,
            "variation_of_constants_readout_constraint": (sp.Matrix([[f, sp.conjugate(f)]])*B).applyfunc(sp.simplify),
            "variation_of_constants_exact_forced_equation": (sp.Matrix([[fd, sp.conjugate(fd)]])*B
                                        +rho*sp.Matrix([[f, sp.conjugate(f)]])).applyfunc(sp.simplify),
            "reference_physical_derivative_norm": sp.simplify((fd-rate*f)*sp.conjugate(fd-rate*f)
                          -(positiveW**2+(rate+Wd/(2*positiveW))**2)/(2*positiveW))}


@cache
def momentum_integral():
    y = sp.Symbol("dimensionless_momentum", nonnegative=True)
    # k=Amax*m*y. After angular integration, 1/(2*pi²) times
    # integral k²/(m²+k²/Amax²)² dk = Amax³/(8*pi*m).
    actual = sp.integrate(y**2/(1+y**2)**2, (y, 0, sp.oo))
    return {"radial_resolvent_square_integral": sp.simplify(actual-sp.pi/4)}


@cache
def proof_checks():
    data = {kind: wkb.bounds(kind) for kind in ("transverse", "longitudinal")}
    out = {}
    for kind, d in data.items():
        out[kind+"_positive_WKB_reference"] = bool(d["S_deviation_upper"] < sp.Rational(1, 2))
        out[kind+"_reference_log_rate_below_four"] = bool(d["reference_log_rate_absolute_upper"] < 4)
        out[kind+"_common_residual_constant"] = bool(d["residual_over_inverse_frequency_fourth_upper"] <= RESIDUAL_CONSTANT)
    out.update({"initial_to_final_time_length_one": sp.Rational(1, 2)-sp.Rational(-1, 2) == 1,
                "positive_mass_floor": bool(wkb.MASS_TIME_MIN**2 > 15),
                "maximum_variation_of_constants_exponent_below_quarter": bool(2*RESIDUAL_CONSTANT/wkb.MASS_TIME_MIN**5 < sp.Rational(1, 4)),
                "actual_lapse_kinetic_weight_lower": 1-sp.Rational(4, 9) == sp.Rational(5, 9),
                "actual_lapse_potential_weight_upper_below_three_halves": bool(1+sp.Rational(28, 81) < sp.Rational(3, 2)),
                "reference_rate_plus_half_log_rate_below_five": 3+sp.Integer(4)/2 == 5,
                "reference_energy_frequency_numerator_upper_four": bool(sp.Rational(9, 4)+sp.Rational(3, 2)+25/wkb.MASS_TIME_MIN**2 < 4),
                "global_scale_factor_on_compact_interval": AMAX == (1+sp.Rational(1, 4))**2,
                "three_polarization_integrated_energy_constant": bool(2*AMAX**4 < 12),
                "pressure_absolute_weight_ratio": 1/sp.Rational(5, 9) == sp.Rational(9, 5),
                "finite_evolution_difference_not_full_renormalized_energy": True})
    return out


def physical_bounds(planck_time_product, reference_mass_time_product):
    L = exact(planck_time_product, "M_tau")
    R = exact(reference_mass_time_product, "m0_tau")
    if L.is_positive is not True or not bool(R >= wkb.MASS_TIME_MIN):
        raise ValueError("Require positive M*tau and m0*tau>=1000")
    return {"M_tau": L, "m0_tau": R,
            "physical_lapse_energy_evolution_difference_over_reference_density": 12*RESIDUAL_CONSTANT/(R*L**2),
            "pressure_evolution_difference_over_reference_density": sp.Rational(108, 5)*RESIDUAL_CONSTANT/(R*L**2),
            "full_renormalized_energy_or_all_order_Hadamard_claim": False}


@cache
def checks():
    return {**oscillator_identities(), **momentum_integral()}
