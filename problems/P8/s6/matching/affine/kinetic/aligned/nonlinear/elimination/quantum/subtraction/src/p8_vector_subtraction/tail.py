"""Generic exact reference subtraction identity and continuous tail bound."""
from functools import cache

import sympy as sp
from p8_vector_state import comparison, energy, wkb

p, r, b2, b4, c = sp.symbols("P2 P4 B2 B4 c1", real=True)
wa, wb = sp.symbols("kinetic_weight potential_weight", real=True)
t = sp.Symbol("inverse_frequency_squared", positive=True)
TAIL_CONSTANT = sp.Integer(400_000)


@cache
def algebra():
    S = 1+p*t+r*t**2
    shifted_rate = c+(b2*t+b4*t**2)/(2*S)
    actual = wa*S+wb/S+wa*t*shifted_rate**2/S
    zero = wa+wb
    second = (wa-wb)*p+wa*c**2
    fourth = (wa-wb)*r+wb*p**2+wa*(c*b2-c**2*p)
    inverse_tail = (2*p*r-p**3+r*(r-p**2)*t)/S
    inverse_second = (p**2-r+p*r*t)/S
    squared_inverse_first = (b4-b2*(2*p+(p**2+2*r)*t+2*p*r*t**2+r**2*t**3))/S**2
    remainder = (wb*inverse_tail+wa*c**2*inverse_second+wa*c*squared_inverse_first
                 +wa*(b2+b4*t)**2/(4*S**3))
    epsilon = sp.Symbol("adiabatic_order", real=True)
    scaling = {p: epsilon**2*p, r: epsilon**4*r, b2: epsilon**3*b2,
               b4: epsilon**5*b4, c: epsilon*c}
    formal = sp.series(actual.subs(scaling, simultaneous=True), epsilon, 0, 5).removeO().expand()
    expected = zero+epsilon**2*t*second+epsilon**4*t**2*fourth
    return {"S": S, "shifted_rate": shifted_rate, "actual": actual,
            "zero": zero, "second": second, "fourth": fourth,
            "remainder": remainder,
            "exact_reference_subtraction_tail": sp.factor(actual-zero-t*second-t**2*fourth-t**3*remainder),
            "formal_derivative_order_not_naive_mass_order": sp.expand(formal-expected),
            "inverse_through_second_order_tail": sp.factor(1/S-1+p*t-(p**2-r)*t**2-t**3*inverse_tail)}


@cache
def bounds(kind):
    b = wkb.bounds(kind)["coefficient_envelopes"]
    P, R, B2, B4 = [b[name] for name in ("P2", "P4", "B2", "B4")]
    tmax = 1/wkb.MASS_TIME_MIN**2
    # Both energy and pressure in the common omega/(4*a_s^3)
    # normalization have |wa|<=1, |wb|<=3/2, |c1|<=2.
    inverse = 2*(2*P*R+P**3+R*(R+P**2)*tmax)
    inverse_second = 2*(P**2+R+P*R*tmax)
    squared_inverse = 4*(B4+B2*(2*P+(P**2+2*R)*tmax+2*P*R*tmax**2+R**2*tmax**3))
    constant = sp.Rational(3, 2)*inverse+4*inverse_second+2*squared_inverse+2*(B2+B4*tmax)**2
    return {"inverse_tail_bound": inverse,
            "inverse_second_order_tail_bound": inverse_second,
            "squared_inverse_first_order_tail_bound": squared_inverse,
            "reference_subtracted_bracket_over_inverse_frequency_sixth_upper": constant,
            "integer_upper": sp.ceiling(constant)}


@cache
def physical_weights():
    z, H = wkb.z, wkb.background()["H"]
    h = (1+wkb.u**2)**3
    alpha, beta = 4/(9*h), 28/(81*h)
    return {"transverse_energy": {"mode": "transverse", "A": sp.Integer(1), "B": 1+beta*(1-z), "c1": H*(1-z)/2},
            "longitudinal_energy": {"mode": "longitudinal", "A": 1-alpha*z, "B": 1+beta, "c1": H*(1+z)/2},
            "transverse_pressure": {"mode": "transverse", "A": sp.Rational(1, 3), "B": (2*z-1)/3, "c1": H*(1-z)/2},
            "longitudinal_pressure": {"mode": "longitudinal", "A": (1+2*z)/3, "B": -sp.Rational(1, 3), "c1": H*(1+z)/2}}


@cache
def physical_weight_checks():
    canonical = energy.canonical()
    clock = energy.actual_coefficients()
    h = (1+wkb.u**2)**3
    mapping = {energy.alpha: clock["a_N"].subs(clock["h"], h),
               energy.beta: clock["b_N"].subs(clock["h"], h),
               energy.H: wkb.background()["H"],
               energy.q: wkb.z*energy.m2/(1-wkb.z)}
    out = {}
    for name, key in (("transverse_energy", "rhoT"), ("longitudinal_energy", "rhoL"),
                      ("transverse_pressure", "pT"), ("longitudinal_pressure", "pL")):
        weights = physical_weights()[name]
        rate = weights["c1"]-wkb.background()["lambda"]/2
        expected = (weights["A"]*(energy.dv-rate*energy.v)**2
                    +weights["B"]*energy.m2/(1-wkb.z)*energy.v**2)/(2*energy.a**3)
        out[name+"_original_physical_quadratic_form"] = sp.factor(canonical[key].subs(mapping)-expected)
    return out


@cache
def reference_terms():
    generic = algebra()
    out = {}
    for name, weights in physical_weights().items():
        data = wkb.frequency(weights["mode"])
        mapping = {p: data["P2"], r: data["P4"], b2: data["B2"], b4: data["B4"],
                   c: weights["c1"], wa: weights["A"], wb: weights["B"]}
        out[name] = {key: sp.factor(generic[key].subs(mapping)) for key in ("zero", "second", "fourth")}
    return out


@cache
def momentum_integral():
    y = sp.Symbol("dimensionless_momentum", nonnegative=True)
    actual = sp.integrate(y**2/(1+y**2)**sp.Rational(5, 2), (y, 0, sp.oo))
    return {"radial_inverse_frequency_fifth_integral": sp.simplify(actual-sp.Rational(1, 3))}


def physical_bounds(planck_time_product, reference_mass_time_product):
    prior = comparison.physical_bounds(planck_time_product, reference_mass_time_product)
    L, R = prior["M_tau"], prior["m0_tau"]
    reference = TAIL_CONSTANT/(72*R**2*L**2)
    return {"M_tau": L, "m0_tau": R,
            "reference_energy_or_pressure_subtracted_tail_over_reference_density": reference,
            "full_subtracted_energy_integral_over_reference_density":
                prior["physical_lapse_energy_evolution_difference_over_reference_density"]+reference,
            "full_subtracted_pressure_integral_over_reference_density":
                prior["pressure_evolution_difference_over_reference_density"]+reference,
            "full_covariant_finite_matching_or_all_order_Hadamard_claim": False}


@cache
def checks():
    data = algebra()
    return {**{name: data[name] for name in
               ("exact_reference_subtraction_tail", "formal_derivative_order_not_naive_mass_order", "inverse_through_second_order_tail")},
            **momentum_integral(), **physical_weight_checks()}
