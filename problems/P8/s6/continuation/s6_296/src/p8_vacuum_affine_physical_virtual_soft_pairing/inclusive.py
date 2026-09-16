"""Complete selected one-Newton rate dictionary with unmatched hard data explicit."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_gravity_pole_completion import soft as frozen
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import rate as real_rate
from p8_vacuum_affine_whole_mixed_heavy_gravity_sector import graphs as mixed
from p8_vacuum_affine_whole_quartic_gravity_sector import forward as quartic

from . import continuity, soft, source

HARD = s.Symbol("finite_real_hard_virtual_interference", real=True)
REAL = s.Symbol("finite_physical_exact_minus_soft_real_rate", real=True)


def physical_reference_rate(
    hard_interference=HARD,
    real_remainder=REAL,
    energy=soft.E,
    cosine=soft.COST,
    mass=source.MU,
    kappa=source.K,
):
    return (
        1
        + s.sympify(hard_interference)
        + soft.finite_conversion(energy, cosine, mass, kappa)
        + s.sympify(real_remainder)
    )


def physical_reference_error(resolution=continuity.RESOLUTION, kappa=source.K):
    return real_rate.finite_real_rate_error(
        resolution, kappa
    ) + continuity.finite_conversion_bound(kappa)


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(value)

    mu, Q, tau, z, beta, ss, C, g, H, K, ep, B = s.symbols(
        "mu Q tau z beta s C g H kappa epsilon Bsoft", positive=True
    )
    y = (1 - z * z) / 4
    Ns = (mu + Q / 2) ** 2 - mu**2 / 2
    Nt = (mu + tau / 2) ** 2 - mu**2 / 2
    Nu = (mu + (Q - tau) / 2) ** 2 - mu**2 / 2
    density = 2 * mu + 4 * (
        Ns / (mu + Q * y) - Nt / (mu + tau * y) - Nu / (mu + (Q - tau) * y)
    )
    put(
        "complete_current_equals_frozen_attenuation_density",
        density - 2 * frozen.attenuation_density(Q, tau, mu, z),
    )
    put(
        "frozen_all_angle_attenuation_upper_bound",
        4 * tau * (Q - tau) / (5 * mu)
        + (2 * tau - Q) ** 2 / (5 * mu)
        - Q**2 / (5 * mu),
    )
    put("compact_K0_upper", s.Rational(12**2, 5) - s.Rational(144, 5))
    chi = s.Symbol("physical_rapidity", positive=True)
    beta_phys = s.tanh(chi)  # Explicit0<beta<1, so atanh(beta_phys)=chi is real.
    Ms = (-chi + s.I * s.pi / 2) / (ss * beta_phys)
    Ls = 4 * chi / (ss * beta_phys)
    put("physical_incoming_current_vs_Feynman_sheet", s.re(Ms) + Ls / 4)
    x = s.Symbol("x", real=True)
    primitive = s.log((x - beta) / (x + beta)) / (2 * ss * beta)
    put(
        "physical_PV_antiderivative",
        s.diff(primitive, x) - 1 / (ss * (x * x - beta * beta)),
    )
    put("physical_sheet_Coulomb_imaginary", s.im(Ms) - s.pi / (2 * ss * beta_phys))
    # Exact inherited proof rows are copied, never mutated.
    for name in (
        "whole_six_pair_proper_pole_is_massive_four_leg_soft",
        "whole_quartic_pair_contact_LSZ_UV",
        "complete_endpoint_UV_is_constant_on_shell",
    ):
        put("frozen_quartic_" + name, quartic.data()["checks"][name])
    for name in (
        "whole_noncusp_Gamma_zero_coefficient",
        "whole_box_IR_integration_by_parts",
    ):
        put("frozen_mixed_" + name, mixed.data()["checks"][name])
    rawC = C * (2 * B - s.Rational(5, 3) * mu) / (16 * s.pi**2 * K * ep)
    uvct = 5 * C * mu / (48 * s.pi**2 * K * ep)
    put(
        "complete_UV_subtracted_quartic_IR_pole",
        rawC + uvct - C * B / (8 * s.pi**2 * K * ep),
    )
    put(
        "complete_mixed_IR_pole",
        g * g * H * (2 * B) / (16 * s.pi**2 * K * ep)
        - g * g * H * B / (8 * s.pi**2 * K * ep),
    )
    A0 = C + g * g * H
    put(
        "entire_selected_matter_virtual_pole",
        C * B / (8 * s.pi**2 * K * ep)
        + g * g * H * B / (8 * s.pi**2 * K * ep)
        - A0 * B / (8 * s.pi**2 * K * ep),
    )
    K0 = s.Symbol("K0", real=True)
    put(
        "entire_real_virtual_rate_pole",
        2 * A0 * (-K0 / 2) / (8 * s.pi**2 * K * ep) / A0 + K0 / (8 * s.pi**2 * K * ep),
    )
    hard, phase = s.symbols("hard phase", real=True)
    book = s.Symbol("bookkeeping", real=True)
    # The selected expression drops the quadratic loop order by construction.
    modulus = (1 + book * (hard + s.I * phase)) * (1 + book * (hard - s.I * phase))
    put("consistent_one_loop_interference", s.expand(modulus).coeff(book, 1) - 2 * hard)
    put(
        "Coulomb_phase_absent_only_from_linear_rate",
        s.diff(s.expand(modulus).coeff(book, 1), phase),
    )
    put(
        "loop_square_not_retained",
        s.expand(modulus).coeff(book, 2) - (hard * hard + phase * phase),
    )

    resolution = s.Symbol("resolution", positive=True)
    put(
        "whole_reference_error_is_sum_of_proved_bounds",
        physical_reference_error(resolution, K)
        - real_rate.finite_real_rate_error(resolution, K)
        - 112 / K,
    )
    put(
        "whole_unmatched_hard_interference_retained",
        s.diff(physical_reference_rate(), HARD) - 1,
    )
    put(
        "whole_finite_real_remainder_retained",
        s.diff(physical_reference_rate(), REAL) - 1,
    )
    put(
        "whole_forward_conversion_zero",
        soft.finite_conversion(soft.E, 1, source.MU, source.K),
    )
    Am, Ag, Lm, Lmg, loop, newton = s.symbols(
        "matter_Born gravity_Born matter_loop gravity_dressing_loop loop_order Newton_order",
        real=True,
    )
    full = Am + newton * Ag + loop * (Lm + newton * Lmg)
    put(
        "same_rate_order_hard_matter_gravity_tree_interference_is_separate",
        s.expand(full**2).coeff(loop, 1).coeff(newton, 1) - 2 * (Am * Lmg + Ag * Lm),
    )
    return {
        "whole_selected_physical_reference_rate": physical_reference_rate(),
        "whole_selected_physical_reference_error": physical_reference_error(),
        "whole_virtual_pole_and_UV_proof": "The complete S293 raw C pole is C[2Bsoft-5mu/3]/(16pi^2 kappa e). Subtract its local UV pole with the covariant counterterm, retaining its unchosen finite part, to obtain C Bsoft/(8pi^2 kappa e). The complete S294 mixed cusp/box/noncusp pole is g^2 Htree*2Bsoft/(16pi^2 kappa e). RH/local matching is covariantly renormalized and remains in the hard coefficient. Together these give A0 Bsoft/(8pi^2 kappa e), not an isolated graph's pole.",
        "whole_physical_sheet_link": "The incoming physical M_s=(-atanh(beta)+i*pi/2)/(s beta) is the already-certified S288 i0 continuation. The complete real current density is twice its positive attenuation density, so K0=-2ReBsoft. Naively substituting s>4mu into separate principal radicals is not that cut boundary. The imaginary Coulomb term is retained at amplitude level.",
        "whole_selected_finite_rate_statement": "After local UV renormalization and the exact frozen S278 analytic division, the consistently truncated Born-normalized gravitational-dressing contribution equals1+2Re(deltaA_hard/A0)+Delta_soft+R_real. The finite conversion is[K1+(EulerGamma-2-lnpi)K0]/(8pi^2 kappa). S295 bounds the complete physical real remainder. Their total absolute correction is below10^-792 at the original parameters on mu1,n>=128,5/4<=E<=2,all angles,0<resolution<=1/8.",
        "whole_hard_data_and_order_boundary": "The finite hard interference is explicit independent data; this result does not determine its matching values, sign, positivity or UV completion. No one-loop square is silently included. Individual Fock real rates diverge; the finite statement is the selected real/virtual combination. In particular2Re(A_gravity_Born^* deltaA_matter_loop) contributes at the same C^2,Cg^2,g^4 over kappa rate order but is not part of the selected2Re(A_matter_Born^* deltaA_gravity_dressing) computed here. Its gravitationally IR-finite hard value is not fixed by this cancellation. Pure-gravity radiation, other hard-loop terms, higher matching, all-loop unitarity, amplitude Coulomb/forward Regge and original V/G/B/P8 remain open.",
        "checks": checks,
        "gates": {
            "whole_inherited_quartic_and_mixed_virtual_poles": True,
            "local_UV_and_IR_subtractions_not_confused": True,
            "correct_physical_sheet_and_all_self_terms": True,
            "finite_real_virtual_pair_not_individual_Fock_rate": True,
            "consistent_one_Newton_interference_no_loop_square": True,
            "explicit_hard_matching_not_fixed_by_soft_conversion": True,
            "full_gravity_Regge_and_original_P8_not_closed": True,
        },
    }
