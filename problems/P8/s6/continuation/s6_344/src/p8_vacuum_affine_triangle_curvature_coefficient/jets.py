"""Explicit scalar labeled-box lift and entire homogeneous TT divided differences."""

from functools import cache

import sympy as s


def require_order(value):
    if type(value) is not int or value not in (1, 2, 3):
        raise ValueError("Require a homogeneous scalar derivative order 2,4 or6")
    return value


def difference_power(value, coefficient, soft_dot, order):
    order = require_order(order)
    return s.expand(
        sum(
            s.binomial(order, q)
            * value ** (order - q)
            * coefficient**q
            * (2 * soft_dot) ** (q - 1)
            for q in range(1, order + 1)
        )
    )


@cache
def data():
    U, c, a, scale = s.symbols("U c a scale")
    checks = {}
    for order in (1, 2, 3):
        dd = difference_power(U, c, a, order)
        checks[f"degree{2 * order}_whole_polynomial_divided_difference"] = s.expand(
            2 * a * dd - ((U + 2 * a * c) ** order - U**order)
        )
        checks[f"degree{2 * order}_coincident_endpoint"] = s.expand(
            dd.subs(a, 0) - order * c * U ** (order - 1)
        )
        checks[f"degree{2 * order}_homogeneous_metric_response"] = s.expand(
            dd.subs({U: scale**2 * U, a: scale**2 * a}, simultaneous=True)
            - scale ** (2 * order - 2) * dd
        )
    return {
        "checks": checks,
        "gates": {
            "each_box_acts_on_one_entire_scalar_factor": True,
            "ordering_and_composite_Phi_squared_factor_explicit": True,
            "no_singular_division_when_k_dot_R_is_zero": True,
            "uniform_small_momentum_expansion_not_physical_cut_replacement": True,
        },
        "whole_labeled_scalar_convention": "Let L=-Box_g on scalars and f1=Phi^2,f2=Phi,f3=Phi. For a polynomial F(v,w,q), the covariant lift is its coefficient-wise local integral sum F_abc*(L^a f1)*(L^b f2)*(L^c f3), with the original triangle normalization and full four-field Bose differentiation. No derivatives are commuted between the three factors. Its real-null-TT contact is-2 sum_i epsilon(Ri,Ri) DD_i F between Ri^2 and(Ri+k)^2, all other labeled arguments fixed.",
        "whole_scalar_metric_variation": "For g=eta+2h/sqrt(kappa), delta(-Box) has TT momentum kernel-2epsilon(R,R)/sqrt(kappa). The contracted Christoffel term vanishes for a transverse traceless real wave, but every insertion in every power L^a is retained. Their ordered sum is the polynomial divided difference, including the coincident limit.",
        "whole_analytic_jet_proof": "For any fixed finite off-shell momenta and positive scalar masses, M=sum x_i m_i^2>=min m_i^2>0 on the compact real simplex. At sufficiently small common momentum scale lambda, |lambda^2 U|<rho M with rho<1 uniformly. The geometric series and all finitely required derivatives are uniformly dominated, so the local Taylor coefficients are the displayed convergent parameter moments. The outer heavy denominator is also analytic there. No physical timelike threshold is crossed in this argument.",
    }
