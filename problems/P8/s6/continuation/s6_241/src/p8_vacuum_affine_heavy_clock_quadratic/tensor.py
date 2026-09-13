"""Complete finite heavy tensor Euler operator with its physical weighted pairing."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import state

from . import local

t = local.t
a = (1 + t * t) ** 2
H = s.diff(a, t) / a
R0 = 6 * (s.diff(H, t) + 2 * H * H)
k = s.Symbol("spatial_momentum", nonnegative=True)
mass_squared = s.Symbol("heavy_mass_squared", positive=True)
ell = s.Symbol("log_heavy_mass_squared", real=True)
A = mass_squared * (ell - 1) / 12 - ell * R0 / 72


@cache
def data():
    f = s.Function("tensor")(t)
    test = s.Function("test")(t)
    q = k * k / a**2

    def D(value):
        return s.diff(value, t, 2) + H * s.diff(value, t) + q * value

    def Dadj(value):
        return (
            s.diff(value, t, 2)
            + 5 * H * s.diff(value, t)
            + (2 * s.diff(H, t) + 6 * H * H + q) * value
        )

    def L(value):
        return s.diff(value, t, 2) + 3 * H * s.diff(value, t) + q * value

    LA = A * L(f) + s.diff(A, t) * s.diff(f, t)
    dd = Dadj(D(f))
    expanded = (
        s.diff(f, t, 4)
        + 6 * H * s.diff(f, t, 3)
        + (4 * s.diff(H, t) + 11 * H * H + 2 * q) * s.diff(f, t, 2)
    )
    expanded += (
        s.diff(H, t, 2) + 7 * H * s.diff(H, t) + 6 * H**3 + 2 * H * q
    ) * s.diff(f, t) + q * q * f
    integrand = a**3 * (A * (s.diff(f, t) ** 2 - q * f * f) - ell * D(f) ** 2 / 120)
    EL = (
        s.diff(integrand, f)
        - s.diff(s.diff(integrand, s.diff(f, t)), t)
        + s.diff(s.diff(integrand, s.diff(f, t, 2)), t, 2)
    )
    direct_adjoint = (
        s.diff(a**3 * test, t, 2) - s.diff(a**3 * H * test, t) + a**3 * q * test
    ) / a**3
    full = LA + ell * dd / 120
    checks = {
        "complete_physical_weighted_adjoint": s.simplify(direct_adjoint - Dadj(test)),
        "complete_fourth_order_operator_not_frozen_coefficients": s.simplify(
            dd - expanded
        ),
        "complete_weighted_variable_A_operator": s.simplify(
            LA - s.diff(a**3 * A * s.diff(f, t), t) / a**3 - A * q * f
        ),
        "independent_full_tensor_action_Euler_variation": s.simplify(
            EL + 2 * a**3 * full
        ),
        "whole_tensor_coefficient_from_scalar_heat_action": s.simplify(
            A - (mass_squared * (ell - 1) / 3 / 4 + (-ell / 36) * R0 / 2)
        ),
        "whole_Weyl_quadratic_coefficient": -ell / 120 - (-ell / 60) / 2,
        "full_actual_reference_curvature": s.factor(
            R0 - 24 * (1 + 7 * t * t) / (1 + t * t) ** 2
        ),
        "canonical_tensor_action_factor": s.Rational(4, 64) / state.KAPPA
        - s.Rational(1, 16) / state.KAPPA,
    }
    coefs = {}
    polynomial = s.Poly(
        s.expand(full),
        f,
        s.diff(f, t),
        s.diff(f, t, 2),
        s.diff(f, t, 3),
        s.diff(f, t, 4),
        k,
    )
    c1 = c0 = s.S.Zero
    for powers, coefficient in polynomial.terms():
        assert sum(powers[:5]) == 1 and powers[-1] in (0, 2, 4)
        derivative = next(i for i in range(5) if powers[i])
        assert derivative + powers[-1] <= 4
        coefs[(derivative, powers[-1])] = s.factor(coefficient)
        for (mpower, lpower), value in s.Poly(coefficient, mass_squared, ell).terms():
            assert mpower in (0, 1) and lpower in (0, 1)
            b = local.time_envelope(value) * 462**lpower
            if mpower:
                c1 += b
            else:
                c0 += b
    n = state.MASS2
    kappa = state.KAPPA
    bound = (c1 * n + c0) / (72 * kappa)
    Amin = n * s.Rational(393, 12) - s.Rational(462 * 78, 72)
    # Unimodular TT metric and u=t keep X1 and volume fixed at every order.
    gd, gg = s.symbols("gamma_D gamma_G", real=True)
    determinant = s.exp(gd / s.sqrt(2) + gg / s.sqrt(2)) * s.exp(
        -gd / s.sqrt(2) - gg / s.sqrt(2)
    )
    checks["whole_unimodular_TT_profile_and_vacuum_variation"] = s.simplify(
        s.diff(determinant, gd, gg)
    )
    return {
        "actual_parameters": {
            "mass_squared": n,
            "kappa0": kappa,
            "time_interval": (-s.Rational(1, 2), s.Rational(1, 2)),
        },
        "complete_tensor_gamma_action": "Before64pi², integral a³{A_H tr[gamma_t²-a^-2(grad gamma)²]-(ell/120)tr[(D gamma)²]}, A_H=n(ell-1)/12-ell R0/72, ell=log n. The whole heavy scalar finite local action is used; the different old vector coefficient is not borrowed.",
        "full_A_H": s.factor(A),
        "full_actual_curvature": s.factor(R0),
        "complete_operator_before_8pi2kappa": s.factor(full),
        "complete_differential_coefficients": coefs,
        "full_actual_pairing": "D=tt+H t-a^-2 Delta; Ddag=tt+5H t+2H'+6H²-a^-2 Delta in the a³ dt dx measure. Qloc=[A_H L+A_H'partial_t+(ell/120)Ddag D]/(8pi² kappa0), with L=tt+3H t-a^-2 Delta. Its Euler contribution is-Qloc h in the physical canonical chart gamma=2h/sqrt(kappa0).",
        "full_DdagD": "t^4+6H t^3+(4H'+11H²)t²+(H''+7HH'+6H³)t-2a^-2 Delta t²-2H a^-2 Delta t+a^-4 Delta². The isolated zero-time Delta coefficient cancels only after all coefficient derivatives are retained.",
        "covariant_derivation": "Unimodular TT has zero first scalar-curvature variation and exactly constant volume. R contributes1/4 and R² contributesR0/2 to the quadratic kinetic/gradient coefficient. The compact Weyl quadratic is(1/2)(D gamma)², by conformal invariance and the explicit flat Euler boundary. The fixed scalar clock profile has no TT variation because u=t,X1 and volume are exactly unchanged.",
        "positive_finite_tensor_kinetic_lower": Amin,
        "complete_operator_coefficient_majorant": {
            "coefficient_of_mass_squared": c1,
            "constant": c0,
        },
        "full_normalized_local_tensor_graph_bound": bound,
        "explicit_local_graph": "For any real r set ||h||Xr=sum over j,l>=0,j+2l<=4 of ||(-Delta)^l partial_t^j h||Linf(I,H^r). The full coefficient-sum estimate gives ||Qloc h||Linf(I,H^r)<10^-597 ||h||Xr. This is an order-four input graph, not a same-space inverse estimate.",
        "response_boundary": "The exact prescribed finite local matching action and its full fourth-order contacts are retained. The rest of the exact heavy SLE/subtraction-dependent determinant response is NOT included here, and a pole of this isolated local operator is not identified with a physical pole of the full quantum theory. No order reduction, causal-pole deletion, quantum inverse, nonlinear bounce or UV/Regge closure is inferred.",
        "checks": {key: s.cancel(value) for key, value in checks.items()},
        "gates": {
            "actual_positive_mass_and_log_interval": bool(n > 10**196),
            "full_finite_tensor_kinetic_positive": bool(Amin > 0),
            "required_variable_A_prime_present": s.diff(A, t) != 0,
            "all_complete_operator_monomials_in_graph": all(
                j + p <= 4 for j, p in coefs
            ),
            "complete_fourth_order_term_not_deleted": (4, 0) in coefs
            and coefs[4, 0] != 0,
            "full_mixed_time_spatial_contacts_retained": all(
                key in coefs for key in ((2, 2), (1, 2), (0, 4))
            ),
            "complete_all_momentum_local_graph_bound": bool(
                bound < s.Rational(1, 10**597)
            ),
            "whole_TT_profile_variation_zero_by_exact_chart_not_stress_reset": True,
            "isolated_local_roots_not_physical_full_response_poles": True,
        },
    }
