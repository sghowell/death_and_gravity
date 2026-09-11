"""Absolute continuum bound for the complete time-ordered memory remainder."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference
from p8_vector_state.comparison import AMAX

CREF = reference.PAIR_REF
CERR = reference.PAIR_ERROR
DISPLAY = s.Integer(10) ** 23


@cache
def constants():
    A, m = AMAX, reference.MASS
    J4 = A**3 / (24 * m)
    J10 = 5 * A**3 / (1536 * m**7)
    return {
        "J4_pi_lower_bound": J4,
        "J10_pi_lower_bound": J10,
        "complete_memory_H1_coefficient": 72 * (CREF * CERR * J4 + CERR**2 * J10),
    }


@cache
def data():
    nu, mu = s.symbols("nu mu", positive=True)
    x, y, z, w = s.symbols(
        "reference_D reference_Gamma error_D error_Gamma", complex=True
    )
    c = constants()
    h = nu**-6 + mu**-6
    L = s.Symbol("L", positive=True)
    checks = {
        "complete_pair_product_difference": s.expand(
            s.conjugate(x + z) * (y + w)
            - s.conjugate(x) * y
            - s.conjugate(x) * w
            - s.conjugate(z) * y
            - s.conjugate(z) * w
        ),
        "linear_internal_weight_symmetry": s.expand(
            nu * mu * h - (mu / nu) * nu**-4 - (nu / mu) * mu**-4
        ),
        "quadratic_weight_split_gap": s.factor(
            2 * (nu**-12 + mu**-12) - h * h - (nu**-6 - mu**-6) ** 2
        ),
        "both_mixed_pair_terms_retained": 2 * CREF * CERR - 2 * 10**24,
        "full_polarization_and_current_factor": (s.Rational(1, 4) * 2 * 2) * 3**2 - 9,
        "linear_and_quadratic_integrals": s.expand(
            9
            * (
                2 * CREF * CERR * 2 * L * s.Symbol("J4")
                + CERR**2 * 4 * L * s.Symbol("J10")
            )
            - 36 * L * (CREF * CERR * s.Symbol("J4") + CERR**2 * s.Symbol("J10"))
        ),
    }
    return {
        "comparison": "At each finite spherical regulator subtract only the same mode formula with u_actual replaced by alpha_initial u_W8. Keep the complete actual preparation and its nonzero error. The reference expression is a proof device, not a new physical state or renormalization prescription.",
        "domain": "Real smooth compact spatial tracefree D,Gamma inside the open unit CD slab; Gamma is zero near the common preparation. Their times may overlap. The retarded theta(t-s) remains in both memory terms.",
        "pointwise_pair_inputs": "|a_ref|<=1e9 sqrt(nu mu); |a_actual-a_ref|<=1e15 sqrt(nu mu)(nu^-6+mu^-6). No time derivative is applied to the oscillatory actual mixing.",
        "complete_product_bound": "nu mu [2 Cref Cerr (nu^-6+mu^-6)+Cerr^2(nu^-6+mu^-6)^2]. Both mixed terms and the error-error term are retained.",
        "internal_bound": "For P=k+l, mu<=L nu and nu<=L mu, L=1+|P|/(Amax m). The linear integral is <=2L J4 and the quadratic integral<=4L J10. All nine polarization pairs give36L(Cref Cerr J4+Cerr^2 J10).",
        "time_and_space": "Use |theta|<=1 without changing theta. Unit-time Cauchy and L<2(1+|P|^2), followed by weighted Fourier Cauchy, bound the memory remainder by1e23 M[D]M[Gamma], where M[f]^2=||f||L2(dt dx;F)^2+||grad_x f||L2(dt dx;F)^2.",
        "constants": c,
        "checks": checks,
        "gates": {
            "full_J4_below_inverse_mass": c["J4_pi_lower_bound"] < 1 / reference.MASS,
            "full_J10_below_inverse_mass_seventh": c["J10_pi_lower_bound"]
            < 1 / reference.MASS**7,
            "complete_memory_strict_display": c["complete_memory_H1_coefficient"]
            < DISPLAY,
            "same_actual_reference_pair_coefficients": CREF == 10**9 and CERR == 10**15,
            "mass_gap_kept": reference.MASS == 1000,
        },
    }
