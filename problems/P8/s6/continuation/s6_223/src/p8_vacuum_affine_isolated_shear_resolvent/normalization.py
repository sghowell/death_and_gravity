"""Physical cut normalization and the unchanged flat finite shear Hessian."""

from functools import cache

import sympy as s
from p8_proca_rank_one_inverse import spectral as old_trace
from p8_vacuum_affine_flat_tensor_cut import projectors
from p8_vacuum_affine_local_tensor_response import prescription

y, m, p, sigma = s.symbols("y mass p spectral_mass_squared", positive=True)
W2 = y * y * (30 - 20 * y * y + 3 * y**4) / 30
W0 = y * y * (3 - 2 * y * y + 3 * y**4)
A0 = s.Rational(1, 30)


@cache
def data():
    cut2 = projectors.above_threshold_density(2).subs(
        {projectors.S: sigma, projectors.MASS: m}
    )
    cut0 = projectors.above_threshold_density(0).subs(
        {projectors.S: sigma, projectors.MASS: m}
    )
    change = 4 * m * m / (1 - y * y)
    jacobian = s.diff(change, y)
    radial = {}
    for spin, rho, factor in ((2, cut2, 64), (0, cut0, 768)):
        value = (factor * s.pi**2 * p * rho / (sigma**3 * (sigma + p))).subs(
            sigma, change
        ) * jacobian
        radial[spin] = s.cancel(s.refine(value, s.Q.positive(1 - y * y)))
    fixed = prescription.data()["literal_fixed_finite_density_before_64_pi_squared"]
    names = {str(v): v for v in fixed.free_symbols}
    weyl = s.expand(fixed).coeff(names["Weyl_squared"])
    t = s.Symbol("time", real=True)
    shear = s.Function("unit_Frobenius_shear")(t)
    # Flat homogeneous TF curvature: C^2_(2)=tr(gamma''^2)/2.
    # Euler has zero compact variation and R_(1)=0 in this channel.
    quadratic = weyl * s.diff(shear, t, 2) ** 2 / 2
    euler = s.diff(s.diff(quadratic, s.diff(shear, t, 2)), t, 2)
    fourth = s.expand(euler).coeff(s.diff(shear, t, 4))
    checks = {
        "actual_old_trace_weight_bridge": s.expand(
            old_trace.data()["positive_radial_weight"].subs(old_trace.y, y) - W0
        ),
        "actual_old_trace_finite_coefficient_bridge": old_trace.data()["H_at_zero"] - 4,
        "original_Weyl_coefficient": weyl + s.Rational(1, 30),
        "unit_TF_quadratic_density": quadratic + s.diff(shear, t, 2) ** 2 / 60,
        "Hessian_not_quadratic_action_coefficient": fourth + s.Rational(1, 30),
        "physical_spin2_radial_conversion": s.factor(
            radial[2] - p * W2 / (4 * m * m + p * (1 - y * y))
        ),
        "physical_trace_normalization12_radial_conversion": s.factor(
            radial[0] - p * W0 / (4 * m * m + p * (1 - y * y))
        ),
        "A2_threshold": A0 - s.integrate(W2 / y**2, (y, 0, 1)) + s.Rational(172, 225),
        "A2_dimensionless_slope": s.integrate(W2, (y, 0, 1)) / 4 - s.Rational(3, 56),
        "old_trace_threshold_gap": 4
        - s.integrate(W0 / y**2, (y, 0, 1))
        - s.Rational(16, 15),
        "positive_weight_decomposition": s.expand(
            30 - 20 * y * y + 3 * y**4 - 13 - (1 - y * y) * (17 - 3 * y * y)
        ),
    }
    return {
        "cut_densities": {2: cut2, 0: cut0},
        "normalization": "For unit-Frobenius homogeneous TF metric Q, F2(p)=-1/30+64pi^2 D2(-p)/p^2, with the lower-order polynomial kept outside this isolated fourth-order factor. The full source Q=2vI has squared norm12 and factor768pi^2 in the trace comparison.",
        "actual_fixed_fourth_coefficient": fourth,
        "radial_weights": {2: W2, 0: W0},
        "A2_definition": "A2=-F2=1/30+integral_0^1 p W2(y)/[4m^2+p(1-y^2)]dy.",
        "scope": "This is the isolated flat-vacuum reference fourth-order factor with the original finite Hessian. It is not the full tree-plus-loop symbol, a curved state kernel, or a changed subtraction prescription.",
        "checks": checks,
        "gates": {
            "finite_fourth_has_opposite_sign_to_positive_cut": fourth < 0,
            "zero_p_value_positive": A0 > 0,
            "threshold_value_negative": -s.Rational(172, 225) < 0,
            "dimensionless_zero_p_slope_positive": s.Rational(3, 56) > 0,
        },
    }
