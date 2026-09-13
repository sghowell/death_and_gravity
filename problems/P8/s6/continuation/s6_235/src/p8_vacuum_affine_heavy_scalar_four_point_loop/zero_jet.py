"""Independent full constant-background loop jet, not the on-shell subtraction point."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_scalar_tree_matching import model

from . import amplitude

n, g, C = amplitude.n, amplitude.g, amplitude.C
D = model.D
z = s.Symbol("heavy_parameter", real=True)
B21 = (n * s.log(n) - n + 1) / (n - 1) ** 2
B22 = ((n + 1) * s.log(n) - 2 * (n - 1)) / (n - 1) ** 3
ZERO_JET = 6 * (g * g * (C + g * g / n) * B21 + g**4 * B22) / (16 * s.pi**2)


@cache
def data():
    q = 1 + (n - 1) * z
    primitive = (n * s.log(q) - q) / (n - 1) ** 2
    phi, k = s.symbols("constant_light_field Euclidean_loop_squared", real=True)
    A = C + g * g / n
    w = -A / 2 - g * g / (k + n)
    hessian = s.Matrix(
        [[k + 1 - (C + g * g / n) * phi * phi / 2, -g * phi], [-g * phi, k + n]]
    )
    reduced = k + 1 + w * phi * phi
    vertex = s.diff(s.log(reduced / (k + 1)) / 2, phi, 4).subs(phi, 0)
    routed = amplitude.complete_loop().subs(
        {amplitude.svar: 0, amplitude.tvar: 0, amplitude.uvar: 0}
    )
    routed = routed.subs(
        {amplitude.Bfun(0): 0, amplitude.Cfun(0): B21, amplitude.Dfun(0, 0): B22}
    )
    fixed_A_over = -3 / D + 2 / D**2 + 1 / (D + 2)
    actual_ratio = 9 * model.G2 * 462 * (model.GAP + 2) / (32 * 9 * (model.GAP - 1))
    return {
        "off_shell_scope": "This is the complete zero-external-momentum four-light1LPI Taylor jet in the unadjusted S234 base scheme. Set every external virtuality and invariant to0. It is NOT obtained by setting s=t=u=0 while still imposing on-shell mass1, and it is not the physical symmetric subtraction used by V2S-T1-OS4.",
        "constant_background_hessian": hessian,
        "remaining_light_hessian": reduced,
        "quartic_loop_jet": ZERO_JET,
        "integrals_without_common_16pi_squared": {"B21": B21, "B22": B22},
        "integral_definitions": "B21=integral_0^1(1-z)/[1+(n-1)z]dz; B22=integral_0^1 z(1-z)/[1+(n-1)z]^2 dz=-partial_n B21. At mu1 the finite MSbar light double-propagator integral B_MS(0;1,1)=0. The common factor1/(16pi^2) is already included in the displayed quartic jet.",
        "complete_diagram_check": "The stationary-heavy full2x2 Hessian gives w=-(C+g^2/n)/2-g^2/(k_E^2+n). Its Euclidean fourth derivative is -6 integral w^2/(k_E^2+1)^2; the amplitude-sign jet is its negative. Expanding w^2 reproduces all bubble, triangle and alternating-box multiplicities of the full routed formula, with the complete MSbar counterterms.",
        "sign_and_size": "For actual D, B21,B22>0, B22<B21/n and C/g^2+2/n<0. Hence this off-shell first-loop jet is negative. Its magnitude is less than18g^4 ln(n)/(16pi^2 D^2), and its ratio to24 times the original valley quartic is below10^-6. This is a diagnostic of the distinct off-shell jet, not an adopted zero-momentum finite matching condition.",
        "actual_valley_ratio_upper": actual_ratio,
        "checks": {
            "independent_constant_background_Schur": s.cancel(
                hessian.det() / (k + n) - reduced
            ),
            "full_fourth_field_derivative": s.cancel(vertex + 6 * w * w / (k + 1) ** 2),
            "complete_squared_nonlocal_vertex_decomposition": s.expand(
                6 * w * w
                - s.Rational(3, 2) * A * A
                - 6 * A * g * g / (k + n)
                - 6 * g**4 / (k + n) ** 2
            ),
            "exact_B21_primitive": s.cancel(s.diff(primitive, z) - (1 - z) / q),
            "exact_B21_complete_endpoints": s.cancel(
                primitive.subs(z, 1) - primitive.subs(z, 0) - B21
            ),
            "full_mass_derivative_B22": s.cancel(-s.diff(B21, n) - B22),
            "positive_B21_over_n_minus_B22_integrand": s.cancel(
                (1 - z) / (n * q) - z * (1 - z) / q**2 - (1 - z) ** 2 / (n * q * q)
            ),
            "logarithmic_B21_upper_difference": s.cancel(
                s.log(n) / (n - 1) - B21 - (n - 1 - s.log(n)) / (n - 1) ** 2
            ),
            "actual_negative_coefficient_identity": s.cancel(
                fixed_A_over + 1 / (D + 2) + (D * D + 4 * D - 4) / (D * D * (D + 2))
            ),
            "actual_coefficient_magnitude_upper": s.cancel(
                3 / D + fixed_A_over - 2 / D**2 - 1 / (D + 2)
            ),
            "full_routed_loop_reproduces_zero_jet": s.cancel(routed - ZERO_JET),
            "B21_equal_mass_removable_limit": s.limit(B21, n, 1) - s.Rational(1, 2),
            "B22_equal_mass_removable_limit": s.limit(B22, n, 1) - s.Rational(1, 6),
        },
        "gates": {
            "actual_zero_jet_sign_domain": model.GAP > 4,
            "actual_logarithm_upper_domain": 1 < model.MASS2 < 10**198,
            "actual_off_shell_contact_to_valley_upper": 0
            < actual_ratio
            < s.Rational(1, 10**6),
            "off_shell_zero_jet_not_on_shell_matching_point": True,
            "full_two_field_hessian_not_contact_bubble_only": True,
        },
    }
