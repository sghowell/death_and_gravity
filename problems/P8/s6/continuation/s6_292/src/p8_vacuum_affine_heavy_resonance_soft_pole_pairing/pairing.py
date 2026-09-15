"""Complete known formal real/virtual threshold pairing, not physical IR closure."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_graviton_production_threshold import (
    production as real_production,
)
from p8_vacuum_affine_heavy_graviton_production_threshold import (
    threshold as real_threshold,
)

from . import masters, proper, source

MU, N, G, K, EP = source.MU, source.N, source.G, source.K, source.EP
Y = s.Symbol("dimensionless_threshold_coordinate", positive=True)
L = s.Symbol("compact_dimensionless_window", positive=True)
W = s.Function("fixed_C1_test_weight")


def require_window(mass, heavy, length, scale_squared=1):
    return real_threshold.require_window(mass, heavy, length, scale_squared)


def angular_log(mass=MU, heavy=N):
    mass, heavy = map(s.sympify, (mass, heavy))
    b2 = 1 - 4 * mass / heavy
    return s.Integral(
        real_production.angular_kernel(b2, masters.X) * s.log(1 - masters.X**2),
        (masters.X, 0, 1),
    )


def finite_local_coefficient(mass=MU, heavy=N, cubic=G, kappa=K):
    mass, heavy, cubic, kappa = map(s.sympify, (mass, heavy, cubic, kappa))
    R3 = proper.radiation_coefficient(mass, heavy)
    F0 = real_production.angular_closed(1 - 4 * mass / heavy)
    K0 = real_threshold.leading_coefficient(heavy, mass, cubic, kappa)
    return (
        K0
        * (
            1
            + s.log(heavy / 4)
            + angular_log(mass, heavy) / F0
            + s.re(proper.B1(mass, heavy)) / R3
        )
        / 2
    )


def paired_finite_functional(length=L, mass=MU, heavy=N, cubic=G, kappa=K):
    length, mass, heavy, cubic, kappa = map(
        s.sympify, (length, mass, heavy, cubic, kappa)
    )
    coeff = lambda y: (
        real_threshold.leading_coefficient(heavy * (1 + y), mass, cubic, kappa)
        * W(heavy * (1 + y))
    )
    return (
        s.Integral((coeff(Y) - coeff(0)) / Y, (Y, 0, length))
        + coeff(0) * s.log(length)
        + finite_local_coefficient(mass, heavy, cubic, kappa) * W(heavy)
    )


def paired_remainder(
    epsilon, length, M2, M10, M11, Mvirt2, weight_at_threshold, cubic=G, kappa=K
):
    epsilon, length, M2, M10, M11, Mvirt2, weight_at_threshold, cubic, kappa = map(
        s.sympify,
        (epsilon, length, M2, M10, M11, Mvirt2, weight_at_threshold, cubic, kappa),
    )
    return real_threshold.remainder_majorant(
        epsilon, length, M2, M10, M11
    ) + epsilon * cubic**2 * s.Abs(weight_at_threshold) * Mvirt2 / (16 * s.pi * kappa)


@cache
def data():
    E, a, b = s.symbols("energy transverse_p longitudinal_p", positive=True)
    eta = s.diag(1, -1, -1, -1)
    p = s.Matrix([E, a, 0, b])
    r = s.Matrix([E, -a, 0, -b])
    P = p + r
    q = s.Matrix([1, 0, 0, 1])
    dot = lambda v, w: (v.T * eta * w)[0]
    J = p * p.T / dot(p, q) + r * r.T / dot(r, q) - P * P.T / dot(P, q)
    checks = {}
    checks["whole_three_leg_current_Ward"] = (J * eta * q).applyfunc(s.factor)
    sew = (
        sum(eta[i, i] * eta[j, j] * J[i, j] ** 2 for i in range(4) for j in range(4))
        - s.trace(eta * J) ** 2 / 2
    )
    checks["whole_conserved_current_matches_complete_TT"] = s.factor(
        sew - 2 * E * E * a**4 / (E * E - b * b) ** 2
    )
    B, x, nn = s.symbols("beta_squared angle heavy_mass_squared", positive=True)
    mm = nn * (1 - B) / 4
    Q = nn - 4 * mm
    V = nn * nn - 4 * mm * nn + 2 * mm * mm
    diaglight = (
        mm**2
        / 2
        * (
            1 / (s.sqrt(nn) / 2 - s.sqrt(nn * B) * x / 2) ** 2
            + 1 / (s.sqrt(nn) / 2 + s.sqrt(nn * B) * x / 2) ** 2
        )
    )
    diagheavy = nn / 2
    crosslight = 2 * V / (nn * (1 - B * x * x))
    crossheavy = -2 * (nn - 2 * mm) / (1 - B * x * x)
    kernel = (1 - x * x) ** 2 / (1 - B * x * x) ** 2
    checks["whole_covariant_radiation_sew_pointwise"] = s.factor(
        diaglight + diagheavy + crosslight + crossheavy - Q * Q * kernel / (2 * nn)
    )
    primitive = mm**2 * x / (nn / 4 * (1 - B * x * x))
    checks["whole_two_light_self_angle_primitive"] = s.factor(
        s.diff(primitive, x) - diaglight
    )
    checks["whole_two_light_self_angle_integral"] = s.factor(
        primitive.subs(x, 1) - primitive.subs(x, 0) - mm
    )
    I0 = s.atanh(s.sqrt(B)) / s.sqrt(B)
    F0 = (3 - B + (B - 1) * (B + 3) * I0) / (2 * B * B)
    R3 = mm + nn / 2 - 4 * mm * (1 - mm / nn) * I0
    checks["whole_invariant_soft_coefficient_vs_forward_kernel"] = s.factor(
        R3 - Q * Q * F0 / (2 * nn)
    )
    checks["whole_self_and_pair_coefficient_sum"] = s.factor(
        mm + nn / 2 + 2 * V * I0 / nn - 2 * (nn - 2 * mm) * I0 - R3
    )
    v, t, u = s.symbols("feynman_share cayley_share soft_radius", positive=True)
    vt = (1 - t) / (1 + t)
    Ah = nn * v + mm * (1 - v) ** 2
    checks["whole_heavy_light_cusp_cayley_denominator"] = s.factor(
        Ah.subs(v, vt) - nn * (1 - B * t * t) / (1 + t) ** 2
    )
    checks["whole_heavy_light_cusp_jacobian"] = s.factor(
        -s.diff(vt, t) - 2 / (1 + t) ** 2
    )
    checks["whole_heavy_light_cusp_transformed_integrand"] = s.factor(
        -s.diff(vt, t) / Ah.subs(v, vt) - 2 / (nn * (1 - B * t * t))
    )
    checks["whole_heavy_light_cusp_primitive"] = s.factor(
        s.diff(2 * s.atanh(s.sqrt(B) * t) / (nn * s.sqrt(B)), t)
        - 2 / (nn * (1 - B * t * t))
    )
    Al = mm - nn * v * (1 - v)
    checks["whole_light_light_cusp_to_frozen_massive_master"] = s.factor(
        Al.subs(v, (1 + t) / 2) - nn * (t * t - B) / 4
    )
    beta = s.Symbol("beta", positive=True)
    primitive = s.log((t - beta) / (t + beta)) / (2 * nn * beta)
    checks["whole_light_light_cusp_PV_primitive"] = s.factor(
        s.diff(primitive, t) - 1 / (nn * (t * t - beta * beta))
    )
    checks["whole_light_light_cusp_delta_jacobian"] = (
        s.diff(nn * (t * t - beta * beta), t).subs(t, beta) - 2 * nn * beta
    )
    e = s.Symbol("epsilon", positive=True)
    checks["whole_soft_radial_triangle_primitive"] = s.simplify(
        s.diff(u ** (2 * e) / (2 * e), u) - u ** (-1 + 2 * e)
    )
    mi, mj, z = s.symbols("first_mass_squared second_mass_squared dot", real=True)
    A = mi * v * v + mj * (1 - v) ** 2 - 2 * z * v * (1 - v)
    xx, yy = u * v, u * (1 - v)
    checks["whole_pair_triangle_radial_mass_homogeneity"] = s.expand(
        mi * xx * xx + mj * yy * yy - 2 * z * xx * yy - u * u * A
    )
    checks["whole_pair_triangle_simplex_jacobian"] = (
        s.det(
            s.Matrix([[s.diff(xx, u), s.diff(xx, v)], [s.diff(yy, u), s.diff(yy, v)]])
        )
        + u
    )
    checks["whole_literal_loop_phase"] = (s.I * (-s.I) ** 2 * s.I**3) * s.I / s.I + 1
    # Strip the common1/(16pi^2*kappa*epsilon).
    proper_ll = -2 * V * I0 / nn
    proper_hl = (nn - 2 * mm) * I0
    legs = -mm - nn / 2
    checks["whole_three_pair_plus_three_LSZ_real_IR"] = s.factor(
        proper_ll + 2 * proper_hl + legs + R3
    )
    gg, kk = s.symbols("cubic kappa", positive=True)
    K0 = gg * gg * Q * Q * F0 / (8 * s.pi * kk * nn)
    negative_delta = s.pi * gg * gg * (-2 * R3 / (16 * s.pi * s.pi * kk * e))
    checks["whole_formal_heavy_delta_and_real_cut_IR_pair"] = s.factor(
        K0 / (2 * e) + negative_delta
    )
    checks["whole_positive_cut_vs_invariant_radiation"] = s.factor(
        K0 / (s.pi * gg * gg) - R3 / (4 * s.pi * s.pi * kk)
    )
    rr, bb, ll, nu2, nnorm = s.symbols(
        "positive_R real_B1 angular_log_ratio scale_squared normalizing_mass_squared",
        positive=True,
    )
    kk0 = s.Symbol("positive_K0", positive=True)
    cut = kk0 * (s.EulerGamma + 1 + s.log(nnorm / (16 * s.pi * nu2)) + ll) / 2
    virtual = kk0 * (bb / rr - s.EulerGamma + s.log(4 * s.pi * nu2)) / 2
    checks["whole_known_finite_pair_scale_and_Gamma_cancellation"] = s.simplify(
        s.expand_log(
            cut + virtual - kk0 * (1 + s.log(nnorm / 4) + ll + bb / rr) / 2, force=True
        )
    )
    checks["whole_known_finite_pair_scale_derivative_zero"] = s.simplify(
        s.diff(cut + virtual, nu2)
    )
    rescale = s.Symbol("positive_mass_squared_rescaling", positive=True)
    checks["whole_known_finite_pair_mass_units_cancel"] = s.simplify(
        s.expand_log(
            s.log(rescale * nnorm / 4)
            + bb / rr
            - s.log(rescale)
            - (s.log(nnorm / 4) + bb / rr),
            force=True,
        )
    )
    A0, A1, A2 = s.symbols(
        "analytic_real_C0 analytic_real_C1 analytic_real_C2", real=True
    )
    checks["whole_virtual_Taylor_remainder_coefficient"] = s.expand(
        (A0 + e * A1 + e * e * A2 / 2) / e - A0 / e - A1 - e * A2 / 2
    )
    M2, M10, M11, Mvirt2, wn = s.symbols(
        "real_second real_y real_mixed virtual_second weight_at_threshold",
        positive=True,
    )
    return {
        "whole_conserved_three_leg_radiation_current": J,
        "whole_invariant_radiation_coefficient": proper.radiation_coefficient(),
        "whole_real_and_imaginary_transition_poles": proper.B0()
        / (16 * s.pi**2 * source.K * source.EP),
        "whole_known_finite_local_pair_coefficient": finite_local_coefficient(),
        "whole_known_paired_finite_functional": paired_finite_functional(),
        "whole_explicit_conditional_paired_remainder": paired_remainder(
            source.EP, L, M2, M10, M11, Mvirt2, wn
        ),
        "whole_formal_paired_definition": "D_W(e)=integral_n^{n(1+L)}rho_Hh,e(s)W(s)ds+pi*g^2*2Re[deltaT_minimal,proper+LSZ(e)/g]*W(n). This is a named formal distributional combination, not the complete amplitude or a chosen physical detector observable.",
        "whole_pole_pairing_derivation": "The complete real cut has K0=g^2 R3/(4pi*kappa). The three proper pair triangles and all three external residues give Re(deltaT/g)=-R3/(16pi^2*kappa*e). Hence the real K0/(2e) pole cancels the virtual heavy delta coefficient exactly. The computed pure-GR heavy mass shift is zero in D, so that known sector has no delta-prime. Additional mass/local matching remains separate. Imaginary Coulomb terms multiply principal values and other cuts; they have not been canceled.",
        "whole_finite_and_unit_boundary": "The full first e derivatives give the displayed finite functional. EulerGamma and nu cancel between its two pieces. Define Bhat_e=B_e/n^e; its real first derivative divided by R3 is Re(B1)/R3+ln n, so the local bracket is1-ln4+L_B/F0+Re(Bhat1)/R3 and is invariant under a common mass-squared rescaling. No finite distribution constant is guessed from a D0 residue.",
        "whole_conditional_remainder_proof": "Use the complete S291 real bound. Write C_e=Gamma(1-e)(4pi nu^2)^(-e)B_e, analytic at0, and Mvirt2=sup_0<=e<=1/8|Re d_e^2 C_e|. Taylor's theorem for C_e/e adds e*g^2*abs(W(n))*Mvirt2/(16pi*kappa). The subtracted physical-root moments and positive massive gaps prove finiteness of the virtual supremum. W is fixed C1 and0<L<=1/2. These are symbolic conditional bounds, not evaluated physical detector errors.",
        "whole_not_physical_IR_or_P8_closure": "H remains perturbatively unstable. H-metric and local/higher-EFT finite matching, other cuts and the imaginary Coulomb principal value, physical near-resonance/width, all-angle/dressed IR, fixed-transfer Regge and original V/G/B/P8 remain open. A finite named subset is not a positive physical spectral measure or full-source b20 verdict.",
        "checks": checks,
        "gates": {
            "whole_three_leg_current_and_all_self_cross_terms": True,
            "entire_real_and_virtual_soft_delta_poles_cancel": True,
            "full_evanescent_finite_pair_not_D0_guess": True,
            "loop_scale_Gamma_and_mass_unit_consistency": True,
            "explicit_compact_real_plus_virtual_remainder": True,
            "Coulomb_metric_matching_width_and_P8_boundaries_retained": True,
        },
    }
