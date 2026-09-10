"""Independent local action jets and regular affine coefficient limits."""

from functools import cache

import sympy as sp
from p8_exceptional_vacuum import analytic as previous

from . import family


@cache
def data():
    u, X = family.u, family.X
    d = family.data()
    kap = sp.Symbol("positive_canonical_kappa", positive=True)
    lam = sp.Symbol("fixed_derivative_quartic", positive=True)
    phi, Y = sp.symbols("canonical_field canonical_gradient_square", real=True)
    eps = sp.Symbol("independent_field_amplitude", real=True)
    tree = family.original.data()["original_retuned_tree_scalar"]
    fv = previous.vacuum_lower_function(kap, family.N, lam)
    full = tree + sp.exp(-(u**4)) * d["S"] * (fv - tree)
    transformed = kap * full.subs(
        {u: eps * phi / sp.sqrt(kap), X: eps**2 * Y / kap}, simultaneous=True
    )
    expansion = sp.series(transformed, eps, 0, 5).removeO().expand()
    h = (1 + u * u) ** 3
    RXX = sp.factor(sp.diff(d["R"], X, 2).subs(X, 0))
    RuXX = sp.factor(sp.diff(d["R"], X, 2, u).subs(X, 0))
    return {
        "actual_quadratic": expansion.coeff(eps, 2),
        "actual_quartic": expansion.coeff(eps, 4),
        "source_c_over_X_at_vacuum": RXX / 4,
        "source_F4_at_vacuum": -RXX / 4,
        "source_J3_over_X_at_vacuum": RuXX / 8,
        "regular_auxiliary_q_order": sp.Integer(4),
        "checks": {
            "actual_new_analytic_vacuum_quadratic": sp.expand(
                expansion.coeff(eps, 2) - (Y - phi**2) / 2
            ),
            "actual_new_analytic_vacuum_quartic": sp.expand(
                expansion.coeff(eps, 4) - lam * Y**2 + family.N * phi**4 / (3 * kap)
            ),
            "no_actual_scalar_cubic": expansion.coeff(eps, 3),
            "source_c_removable_coefficient": sp.factor(RXX / 4 + family.N / (2 * h)),
            "source_F4_removable_value": sp.factor(-RXX / 4 - family.N / (2 * h)),
            "source_J3_removable_coefficient": sp.factor(
                RuXX / 8 - family.N * sp.diff(h, u) / (4 * h * h)
            ),
            "new_scalar_full_clock_value": sp.factor(full.subs(X, 1) - tree.subs(X, 1)),
            "same_finite_vacuum_tensor_and_DHOST_jets": sp.factor(
                RXX + 2 * family.N / h
            ),
        },
    }
