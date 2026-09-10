"""Co-scaled canonical light-action limit with its DHOST interaction retained."""

from functools import cache

import sympy as sp

from . import analytic, family


@cache
def data():
    eps = sp.Symbol("inverse_sqrt_kappa", positive=True)
    phi, Y = sp.symbols(
        "fixed_canonical_field fixed_canonical_gradient_square", real=True
    )
    gamma, lam = sp.symbols("fixed_positive_gamma fixed_positive_lambda", positive=True)
    u, X = family.u, family.X
    # q^(gamma/eps^2)=exp[(gamma/eps^2) log(1-eps^4 Y^2)].
    # This exact expansion through eps^4 suffices for every displayed limit.
    S = 1 - gamma * eps**2 * Y**2 + gamma**2 * eps**4 * Y**4 / 2
    q = 1 - eps**4 * Y**2
    h = (1 + eps**2 * phi**2) ** 3
    R = 1 + (1 - S) * (eps**2 * Y - 1) / h
    RX = ((1 - S) + 2 * gamma * Y * (eps**2 * Y - 1) * S / q) / h
    Ftree = family.data()["original_retuned_tree_scalar"]
    Fv = analytic.vacuum_lower_function(eps**-2, gamma / eps**2, lam)
    Ftree = Ftree.subs({u: eps * phi, X: eps**2 * Y}, simultaneous=True)
    Fv = Fv.subs({u: eps * phi, X: eps**2 * Y}, simultaneous=True)
    L = 1 - eps**4 * phi**4
    full_lower = (Fv + (1 - L * S) * (Ftree - Fv)) / eps**2
    limiting_lower = sp.limit(full_lower, eps, 0, dir="+")
    coefficients = {
        "canonical_lower": limiting_lower,
        "canonical_L3": sp.limit(RX / Y, eps, 0, dir="+"),
        "canonical_L4": sp.limit(
            -RX / Y - sp.Rational(7, 4) * eps**2 * RX**2 / R, eps, 0, dir="+"
        ),
        "canonical_L5": sp.limit(eps**2 * RX**2 / (R * Y), eps, 0, dir="+"),
        "tensor_deviation_over_eps_squared": sp.limit(
            (R - 1) / eps**2, eps, 0, dir="+"
        ),
    }
    # Polynomial integration-by-parts identity in arbitrary Lorentzian jets:
    # div[X(v^mu box(phi)-H^{mu nu}v_nu)] =
    # 2(L3-L4)+X[(box(phi))^2-H_mn H^mn].
    eta = sp.diag(1, -1, -1, -1)
    v = sp.Matrix(sp.symbols("v0:4", real=True))
    entries = {
        (i, j): sp.Symbol("H" + str(i) + str(j), real=True)
        for i in range(4)
        for j in range(i, 4)
    }
    HH = sp.Matrix(4, 4, lambda i, j: entries[min(i, j), max(i, j)])
    xx = (v.T * eta * v)[0]
    box = sp.trace(eta * HH)
    hsq = sp.trace(eta * HH * eta * HH)
    L3 = (v.T * eta * HH * eta * v)[0] * box
    L4 = (v.T * eta * HH * eta * HH * eta * v)[0]
    first_gradient = 2 * (HH * eta * v)
    bracket = eta * v * box - eta * HH * eta * v
    div = (first_gradient.T * bracket)[0] + xx * (box**2 - hsq)
    return {
        "scaling": "even n -> infinity, kappa_n=n/gamma, M_n=sqrt(kappa_n)/tau; tau and canonical masses, lambda and gamma fixed",
        "canonical_mass": sp.Integer(1),
        "limits": coefficients,
        "constant_canonical_domain": "fixed compact canonical field/jet sets with Y>-1/(4 gamma)",
        "canonical_tensor_light_mixing": "R-1=O(eps^2); multiplication by a canonical graviton curvature O(eps) leaves O(eps) mixing",
        "retained_flat_quartic_L3_minus_L4": -2 * gamma,
        "local_scalar_potential": "Phi^2/2+gamma Phi^4/3",
        "no_finite_M_gravity_diagrams_removed": True,
        "checks": {
            "coscaled_lower_has_fixed_mass_and_both_interactions": sp.factor(
                limiting_lower - (Y / 2 - phi**2 / 2 + lam * Y**2 - gamma * phi**4 / 3)
            ),
            "coscaled_L3_retained": sp.factor(coefficients["canonical_L3"] + 2 * gamma),
            "coscaled_L4_retained": sp.factor(coefficients["canonical_L4"] - 2 * gamma),
            "coscaled_L5_vanishes": coefficients["canonical_L5"],
            "tensor_deviation_is_suppressed_before_graviton_scaling": sp.factor(
                coefficients["tensor_deviation_over_eps_squared"] + gamma * Y**2
            ),
            "flat_Galileon_integration_by_parts_polynomial": sp.expand(
                div - 2 * (L3 - L4) - xx * (box**2 - hsq)
            ),
            "switch_log_generator_leading_term": sp.limit(
                gamma * sp.log(q) / eps**4, eps, 0, dir="+"
            )
            + gamma * Y**2,
            "selected_order_is_even": sp.Integer(analytic.ORDER) % 2,
            "selected_action_normalization_on_same_scaling_curve": sp.Integer(
                analytic.ORDER
            )
            / analytic.FIXED_GAMMA
            - analytic.KAPPA,
        },
    }
