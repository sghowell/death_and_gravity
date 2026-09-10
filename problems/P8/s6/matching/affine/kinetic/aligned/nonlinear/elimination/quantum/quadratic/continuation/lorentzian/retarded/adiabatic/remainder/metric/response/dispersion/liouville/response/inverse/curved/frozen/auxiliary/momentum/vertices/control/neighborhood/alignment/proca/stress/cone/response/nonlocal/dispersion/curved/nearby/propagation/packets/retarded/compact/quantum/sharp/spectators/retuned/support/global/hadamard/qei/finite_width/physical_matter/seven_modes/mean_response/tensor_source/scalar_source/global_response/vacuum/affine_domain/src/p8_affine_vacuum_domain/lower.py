"""Regular zero-gradient lower dictionary and the inherited boundary obstruction."""

from functools import cache

import sympy as sp
from p8_affine import lower as original

from . import family


@cache
def data():
    X = family.X
    R, RX, Ru, q, qu, F = sp.symbols("positive_R R_X R_u q q_u target_F", real=True)
    p = sp.sqrt(R) / 2
    px = -RX / (4 * sp.sqrt(R))
    pu = Ru / (4 * sp.sqrt(R))
    fphi = Ru / 2
    forcing = 3 * RX * Ru / (4 * R)
    coefficient = 1 / (2 * X) - 3 * RX / (4 * R)
    qX = forcing - coefficient * q
    J3 = (q + fphi) / (4 * p * X)
    J2 = F - X * qu + 3 * X * (q + 2 * fphi) ** 2 / (16 * p * p)
    reduced = original.source_lower(p, px, pu, fphi, J3, -X, 1, J2)
    r = sp.Symbol("positive_R_variable", positive=True)
    W = (
        sp.Rational(4, 3) * r ** -sp.Rational(3, 4)
        + 4 * r ** sp.Rational(1, 4)
        - sp.Rational(16, 3)
    )
    n, hu, h = sp.symbols("positive_order h_u positive_h", positive=True)
    # R=1-nX^2/h+O(X^3); R_u=n h_u X^2/h^2+O(X^3).
    leading_forcing = -sp.Rational(3, 2) * n * n * hu / h**3
    leading_q = -n * n * hu / (3 * h**3)
    y = sp.Symbol("positive_fourth_root_R", positive=True)
    xp = sp.Symbol("positive_integral_X", positive=True)
    zz = sp.Symbol("positive_dummy_X", positive=True)
    RR = sp.Function("positive_smooth_R")(xp)
    gg = sp.Function("smooth_forcing")(xp)
    II = sp.sqrt(xp) * RR ** -sp.Rational(3, 4)
    q_integral = sp.Integral((II * gg).subs(xp, zz), (zz, 0, xp)) / II
    return {
        "global_regular_q_absolute_bound": "|q(u,X)|<32 n |h_u|/h^3 <=96 n/h^2 on the whole stated domain",
        "global_R_X_bound": "|R_X|<7n/h",
        "global_R_u_bound": "|R_u|<5|h_u|/(4h^2)",
        "physical_X_ODE_coefficient": coefficient,
        "physical_X_ODE_forcing": forcing,
        "affine_cubic": J3,
        "affine_lower_scalar": J2,
        "regular_integral": "q(u,X)=X R(u,X)^(3/4) int_0^1 sqrt(s) R(u,sX)^(-3/4) forcing(u,sX) ds",
        "regular_q_leading_X_four_coefficient": leading_q,
        "boundary_obstruction_nonnegative_function": W,
        "old_q_vacuum_divergence_coefficient": "C(u)=-(3h_u/(8h)) int_0^1 W(R(u,z))/sqrt(z) dz",
        "regular_clock_boundary": "q_reg(u,1)=(3h_u/(8h)) int_0^1 W(R(u,z))/sqrt(z) dz",
        "checks": {
            "literal_regular_integral_solves_ODE": sp.factor(
                sp.diff(q_integral, xp) + sp.diff(II, xp) * q_integral / II - gg
            ),
            "forcing_total_derivative_W_sign": sp.factor(
                sp.sqrt(xp)
                * r ** -sp.Rational(3, 4)
                * 3
                * RX
                * (hu / h)
                * (1 - r)
                / (4 * r)
                + 3 * hu * sp.sqrt(xp) * sp.diff(W, r) * RX / (4 * h)
            ),
            "global_integral_norm_constant": sp.Rational(6, 5)
            * 3
            * sp.Rational(2, 3)
            * sp.Rational(105, 8)
            - sp.Rational(63, 2),
            "literal_lower_Q1": sp.factor(reduced["Q1"] - q),
            "literal_lower_Q2": sp.factor(reduced["Q2"] + 2 * qX),
            "literal_scalar_after_boundary": sp.factor(reduced["P"] + X * qu - F),
            "W_derivative": sp.factor(
                sp.diff(W, r) - (r - 1) * r ** -sp.Rational(7, 4)
            ),
            "W_at_one": W.subs(r, 1),
            "W_positive_factorization": sp.factor(
                W.subs(r, y**4)
                - 4 * (y - 1) ** 2 * (3 * y * y + 2 * y + 1) / (3 * y**3)
            ),
            "regular_q_indicial_coefficient": sp.factor(
                sp.Rational(9, 2) * leading_q - leading_forcing
            ),
            "original_inverse_factor_derivative": sp.factor(
                sp.diff(sp.sqrt(X) * sp.Function("R")(X) ** -sp.Rational(3, 4), X)
                / (sp.sqrt(X) * sp.Function("R")(X) ** -sp.Rational(3, 4))
                - (
                    1 / (2 * X)
                    - 3 * sp.diff(sp.Function("R")(X), X) / (4 * sp.Function("R")(X))
                )
            ),
        },
    }
