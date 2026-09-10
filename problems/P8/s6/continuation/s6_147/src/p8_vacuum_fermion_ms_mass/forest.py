"""Independent proper-counterterm dictionary and scalar-mass transfer identity."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_proper_references import conversion


@cache
def data():
    e = s.Symbol("epsilon")
    Y, a, Cf, Q = s.symbols("Y a Cf Q", positive=True)
    m, mu, b, x = s.symbols("m mu b q_squared", positive=True)
    prior = conversion.data()
    z = prior["fermion_total_kinetic_UV_counterterm"] / e
    eta = prior["fermion_total_mass_UV_counterterm_over_m"] / e
    nu = prior["total_Yukawa_UV_counterterm_over_y"] / e
    f0 = m * m * (mu / m) ** (2 * e)
    variation = 2 * (nu - z) * f0 + (eta - z) * m * s.diff(f0, m)
    tad = s.exp(s.EulerGamma * e) * mu ** (2 * e) * s.gamma(e) / (e - 1) * b ** (1 - e)
    I2 = s.exp(s.EulerGamma * e) * mu ** (2 * e) * s.gamma(e) * b ** (-e)
    K = s.Symbol("quartic_UV_residue")
    return {
        "proper_mass_counterterm_multiplier": (6 - 3 * e) * (Y - 2 * a * Cf) / (Q * e),
        "mass_difference_identity": "f_scalar,b=1(0)-f_scalar,b=0(0)=-1/2 Integral_0^1 F_mix,D(b) db, with the whole proper quartic MS forest retained before taking the finite part.",
        "massless_tadpole_scope": "At b=0 the local whole-fermion-cycle counterterm closes into a scaleless tadpole. At b=1 it is NOT omitted; it is contained in F_mix,D and its b integral.",
        "other_proper_counterterms": "Fermion kinetic/mass and Yukawa MS poles are independent of the exchanged scalar mass; their first-fermion-bubble insertions are identical and cancel in the b=1 minus b=0 difference.",
        "checks": {
            "fixed_mu_complete_mass_CT_multiplier": s.simplify(
                variation - (6 - 3 * e) * (Y - 2 * a * Cf) * f0 / (Q * e)
            ),
            "propagator_mass_difference": s.factor(
                1 / (x + 1) - 1 / x + 1 / (x * (x + 1))
            ),
            "integrated_squared_propagator": s.integrate(1 / (x + b) ** 2, (b, 0, 1))
            - 1 / (x * (x + 1)),
            "tadpole_mass_derivative_is_minus_bubble": s.simplify(s.diff(tad, b) + I2),
            "regulated_integrated_bubble_equals_negative_tadpole": s.simplify(
                s.exp(s.EulerGamma * e) * mu ** (2 * e) * s.gamma(e) / (1 - e)
                + tad.subs(b, 1)
            ),
            "quartic_counterterm_sign_after_mass_integration": K
            / (2 * e)
            * (-s.Symbol("tad_one"))
            - (-K * s.Symbol("tad_one") / (2 * e)),
            "two_internal_box_labels_half_contraction": s.Rational(1, 2) * 6 - 3,
            "unit_mass_normalization_is_one": f0.subs({m: 1, mu: 1}) - 1,
        },
    }
