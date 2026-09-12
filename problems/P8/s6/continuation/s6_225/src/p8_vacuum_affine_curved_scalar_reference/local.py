"""Complete original finite scalar Hessian and its retained lower-order terms."""

from functools import cache

import sympy as s
from p8_vacuum_affine_isolated_shear_resolvent import normalization as shear
from p8_vacuum_affine_local_tensor_response import prescription

from . import geometry as g


@cache
def data():
    raw = g.literal()
    SD, WD = g.channels(g.wd, g.cd)
    SG, WG = g.channels(g.wg, g.cg)
    TD, TG = 3 * g.wd - g.cd, 3 * g.wg - g.cg
    d = lambda f: s.diff(f, g.eta)
    rDG = (
        24 * d(g.wd) * d(g.wg)
        - 8 * (d(g.wd) * d(g.cg) + d(g.cd) * d(g.wg))
        + 4 * d(g.cd) * d(g.cg)
        - 20 * g.q * g.wd * g.wg
        + 4 * g.q * (g.cd * g.wg + g.cg * g.wd)
    )
    mixed = lambda value: g.jet(value).coeff(g.ed, 1).coeff(g.eg, 1)
    volume = mixed(raw["volume"])
    Einstein = mixed(raw["volume"] * raw["Rold"])
    R2 = mixed(raw["volume"] * raw["Rold"] ** 2)
    fixed = prescription.data()["literal_fixed_finite_density_before_64_pi_squared"]
    names = {str(v): v for v in fixed.free_symbols}
    fixed = s.expand(fixed.subs(names["mass"], g.bridge.MASS))
    Cw = fixed.coeff(names["Weyl_squared"])
    Cr = fixed.coeff(names["R_old"], 2)
    Ce = fixed.coeff(names["R_old"], 1)
    Cv = fixed.subs(
        {names["Weyl_squared"]: 0, names["Euler_density"]: 0, names["R_old"]: 0}
    )
    Weyl = s.Rational(8, 3) * Cw * WD * WG
    complete = Cv * volume + Ce * Einstein + Cr * R2 + Weyl
    principal = 72 * Cr * SD * SG + Weyl
    remainder = Cr * (
        12 * g.U * rDG + 72 * g.U * (TD * SG + TG * SD) + 36 * g.U**2 * TD * TG
    )
    remainder += (
        Ce * g.a**2 * (rDG + 6 * (TD * SG + TG * SD) + 6 * g.U * TD * TG)
        + Cv * g.a**4 * TD * TG
    )
    variables = []
    weights = []
    mapping = {}
    for field in (g.wd, g.cd, g.wg, g.cg):
        for order in range(3):
            symbol = s.Symbol(str(field.func) + "_jet" + str(order))
            variables.append(symbol)
            weights.append(order)
            mapping[s.diff(field, g.eta, order)] = symbol
    polynomial = s.Poly(s.expand(remainder.xreplace(mapping)), *variables, g.k)
    orders = tuple(
        sorted(
            {
                sum(e * w for e, w in zip(mon, weights + [1]))
                for mon, coef in polynomial.terms()
                if coef != 0
            }
        )
    )
    TF = s.simplify(
        complete.subs({g.wd: g.cd / 3, g.wg: g.cg / 3, g.k: 0}).doit()
        / s.Rational(8, 3)
    )
    kinetic = s.expand(TF).coeff(s.diff(g.cd, g.eta)).coeff(s.diff(g.cg, g.eta))
    fourth = s.expand(TF).coeff(s.diff(g.cd, g.eta, 2)).coeff(s.diff(g.cg, g.eta, 2))
    bg = g.profile()
    aa = bg["a"]
    actual_kinetic = s.factor(
        kinetic.subs(
            {
                s.diff(g.a, g.eta, 2): aa * bg["U"],
                s.diff(g.a, g.eta): aa * bg["h"],
                g.a: aa,
            },
            simultaneous=True,
        )
    )
    checks = {
        "literal_mixed_curvature_second_variation": s.simplify(
            mixed(raw["Rold"]) * g.a**2 - rDG
        ),
        "complete_R2_density_Hessian": s.simplify(
            R2
            - 72 * SD * SG
            - 12 * g.U * rDG
            - 72 * g.U * (TD * SG + TG * SD)
            - 36 * g.U**2 * TD * TG
        ),
        "complete_Einstein_density_Hessian": s.simplify(
            Einstein - g.a**2 * (rDG + 6 * (TD * SG + TG * SD) + 6 * g.U * TD * TG)
        ),
        "complete_volume_density_Hessian": s.simplify(volume - g.a**4 * TD * TG),
        "unchanged_Weyl_finite_coefficient": Cw
        - shear.data()["actual_fixed_fourth_coefficient"],
        "unchanged_trace_curvature_square_coefficient": 72 * Cr + 4,
        "complete_original_local_Hessian_not_only_curvature_square": s.simplify(
            complete - principal - remainder
        ),
        "zero_transfer_unit_Frobenius_TF_fourth_bridge": s.factor(
            fourth + s.Rational(1, 30)
        ),
        "zero_transfer_unit_Frobenius_TF_kinetic_bridge": s.factor(
            actual_kinetic
            - 2 * aa * aa * g.tensor_input.A.subs(g.tensor_input.t, bg["t"])
        ),
        "zero_transfer_TF_full_local_reconstruction": s.simplify(
            TF
            - kinetic * s.diff(g.cd, g.eta) * s.diff(g.cg, g.eta)
            - fourth * s.diff(g.cd, g.eta, 2) * s.diff(g.cg, g.eta, 2)
        ),
    }
    return {
        "mixed_rescaled_scalar_curvature": rDG,
        "trace_volume_jets": (TD, TG),
        "curvature_square_factor": principal,
        "complete_lower_order_local_remainder": s.expand(remainder),
        "actual_fixed_coefficients": {
            "Weyl_squared": Cw,
            "R_old_squared": Cr,
            "R_old": Ce,
            "constant": Cv,
        },
        "weighted_derivative_orders_in_remainder": orders,
        "zero_transfer_TF_normalization": s.Rational(8, 3),
        "boundary": "This is the complete original finite local Hessian on the exponential synchronous spatial scalar quotient, before64pi^2 and in conformal Lebesgue measure. Compact Euler and box-R variations vanish as variational boundaries. The background-curvature, second-metric-variation, volume and Einstein terms are retained.",
        "order_not_norm": "The remainder has total differential order at most2 in conformal time and comoving space. This does not make it a bounded or small same-space perturbation of the complete quantum-force system.",
        "checks": checks,
        "gates": {
            "all_lower_order_terms_retained": s.expand(remainder) != 0,
            "remainder_has_no_fourth_order_monomial": max(orders) <= 2,
            "remainder_has_actual_second_order_terms": max(orders) == 2,
            "original_mass_and_volume_coefficients_nonzero": Ce > 0 and Cv > 0,
            "complete_trace_curvature_coefficient_negative": Cr < 0,
            "complete_Weyl_curvature_coefficient_negative": Cw < 0,
        },
    }
