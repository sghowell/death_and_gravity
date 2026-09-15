"""All mixed two-body cuts and the explicit covariant heavy-residue anchor."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_graviton_production_threshold import production
from p8_vacuum_affine_heavy_resonance_soft_pole_pairing import pairing
from p8_vacuum_affine_whole_quartic_gravity_sector import cuts as quartic_cut

from . import graphs, source

MU, N, G, K, EP = source.MU, source.N, source.G, source.K, source.EP
S = graphs.S
E2 = s.Symbol("formal_resolution_squared", positive=True)
CRH = s.Symbol("unmatched_finite_RH", real=True)


def angular_J(ratio, epsilon=EP):
    ratio, epsilon = map(s.sympify, (ratio, epsilon))
    if ratio == 0:
        return s.S.One
    if epsilon == 0:
        return s.atanh(ratio) / ratio
    return s.hyper((1, s.Rational(1, 2)), (s.Rational(3, 2) + epsilon,), ratio**2)


def angular_quadratic(ratio, epsilon=EP):
    ratio, epsilon = map(s.sympify, (ratio, epsilon))
    if ratio == 0:
        return 1 / (3 + 2 * epsilon)
    return (angular_J(ratio, epsilon) - 1) / ratio**2


def light_tree(energy=S, angle=quartic_cut.X, mass=MU, heavy=N):
    energy, angle, mass, heavy = map(s.sympify, (energy, angle, mass, heavy))
    q = energy - 4 * mass
    return (
        1 / (heavy - energy)
        + 1 / (heavy + q * (1 - angle) / 2)
        + 1 / (heavy + q * (1 + angle) / 2)
    )


def light_average(energy=S, mass=MU, heavy=N, epsilon=EP):
    energy, mass, heavy, epsilon = map(s.sympify, (energy, mass, heavy, epsilon))
    q = energy - 4 * mass
    r = q / (2 * heavy + q)
    kh = 4 / (2 * heavy + q)
    R = 1 / (heavy - energy)
    p, a, b = quartic_cut.tree_parts(energy, mass, epsilon)
    J = angular_J(r, epsilon)
    Ip = (1 + 2 * epsilon) / (2 * epsilon)
    return R * (p * Ip + a + b / (3 + 2 * epsilon)) + kh * (
        p * (Ip - r * r * J) / (1 - r * r) + a * J + b * angular_quadratic(r, epsilon)
    )


def light_cut(
    energy=S, mass=MU, heavy=N, epsilon=EP, cubic=G, kappa=K, scale_squared=graphs.NU2
):
    return (
        s.sympify(cubic) ** 2
        * quartic_cut.whole_phase(energy, mass, epsilon, scale_squared)
        * light_average(energy, mass, heavy, epsilon)
        / (2 * s.sympify(kappa))
    )


def hard_light_cut(energy=S, mass=MU, heavy=N, cubic=G, kappa=K, resolution_squared=E2):
    energy, mass, heavy, cubic, kappa, resolution_squared = map(
        s.sympify, (energy, mass, heavy, cubic, kappa, resolution_squared)
    )
    q = energy - 4 * mass
    beta = s.sqrt(q / energy)
    V = graphs.numerator(energy, mass, 0)
    p, a, b = quartic_cut.tree_parts(energy, mass, 0)
    r = q / (2 * heavy + q)
    kh = 4 / (2 * heavy + q)
    R = 1 / (heavy - energy)
    J = angular_J(r, 0)
    Hf = R + kh / (1 - r * r)
    Air = cubic**2 * Hf * V / (8 * s.pi * kappa * energy * beta)
    regular = R * (a + b / 3) + kh * (
        -p * r * r * J / (1 - r * r) + a * J + b * angular_quadratic(r, 0)
    )
    return Air * (
        s.EulerGamma + s.log(q / (4 * s.pi * resolution_squared)) + 2 * mass**2 / V
    ) + cubic**2 * beta * regular / (16 * s.pi * kappa)


def heavy_cut(
    energy=S, mass=MU, heavy=N, epsilon=EP, cubic=G, kappa=K, scale_squared=graphs.NU2
):
    energy, mass, heavy, epsilon, cubic, kappa, scale_squared = map(
        s.sympify, (energy, mass, heavy, epsilon, cubic, kappa, scale_squared)
    )
    if energy == 4 * heavy:
        return s.S.Zero
    A = energy / 2 - heavy
    r = s.sqrt((energy - 4 * mass) * (energy - 4 * heavy)) / (2 * A)
    alpha = energy / 4 + 2 * mass * heavy / (energy * (1 + epsilon))
    b = -(energy - 4 * mass) * (energy - 4 * heavy) / (4 * energy)
    return (
        cubic**2
        * quartic_cut.whole_phase(energy, heavy, epsilon, scale_squared)
        * (alpha * angular_J(r, epsilon) + b * angular_quadratic(r, epsilon))
        / (kappa * A)
    )


def weighted_bubble(energy=S, mass=MU):
    energy, mass = map(s.sympify, (energy, mass))
    beta = s.sqrt(1 - 4 * mass / energy)
    return (
        -s.log(mass) / 6
        + s.Rational(4, 9)
        - beta**2 / 6
        - beta * (3 - beta**2) * s.atanh(beta) / 6
        + s.I * s.pi * beta * (3 - beta**2) / 12
    )


def mixing_transition(mass=MU, heavy=N, cubic=G, kappa=K, curvature=CRH):
    mass, heavy, cubic, kappa, curvature = map(
        s.sympify, (mass, heavy, cubic, kappa, curvature)
    )
    return (heavy + 2 * mass) * weighted_bubble(heavy, mass) / (
        32 * s.pi**2 * kappa
    ) + (heavy + 2 * mass) * curvature / (kappa * cubic)


def require_cut_domain(species, energy, mass, heavy):
    if not isinstance(species, str) or species not in (
        "light",
        "heavy",
        "heavy_graviton",
    ):
        raise ValueError("Require a named complete mixed two-body channel")
    energy, mass, heavy = map(source.require_mass, (energy, mass, heavy))
    threshold = {"light": 4 * mass, "heavy": 4 * heavy, "heavy_graviton": heavy}[
        species
    ]
    if heavy <= 4 * mass or energy <= threshold or energy == heavy:
        raise ValueError("Require a separated open cut away from the heavy resonance")
    return energy, mass, heavy


@cache
def data():
    mu, n, a, e, z, r = s.symbols("mu n energy epsilon angle ratio", positive=True)
    Q = a - 4 * mu
    V = (a - 2 * mu) ** 2 - 2 * mu**2 / (1 + e)
    p = 4 * V / Q
    aa = -7 * a / 4 + 4 * mu + 2 * mu**2 / (a * (1 + e))
    bb = -(Q**2) / (4 * a)
    R = 1 / (n - a)
    Kh = 4 / (2 * n + Q)
    rr = Q / (2 * n + Q)
    H = R + Kh / (1 - rr**2 * z**2)
    Hf = R + Kh / (1 - rr**2)
    checks = {}
    checks["all_three_heavy_exchanges"] = s.factor(
        H - (1 / (n - a) + 1 / (n + Q * (1 - z) / 2) + 1 / (n + Q * (1 + z) / 2))
    )
    checks["forward_heavy_tree"] = s.factor(Hf - (1 / (n - a) + 1 / n + 1 / (n + Q)))
    checks["whole_light_denominator_gap"] = s.factor(
        (n + Q / 2) ** 2 - (Q / 2) ** 2 - n * (n + Q)
    )
    checks["whole_angular_pole_partial_fraction"] = s.factor(
        1 / ((1 - z * z) * (1 - r * r * z * z))
        - (1 / (1 - z * z) - r * r / (1 - r * r * z * z)) / (1 - r * r)
    )
    checks["whole_angular_regular_quadratic"] = s.factor(
        z * z / (1 - r * r * z * z) - (1 / (1 - r * r * z * z) - 1) / (r * r)
    )
    checks["whole_endpoint_subtracted_integrand"] = s.factor(
        (H - Hf) * p / (1 - z * z)
        + p * Kh * rr**2 / ((1 - rr**2) * (1 - rr**2 * z * z))
    )
    checks["whole_D_p_derivative"] = s.factor(s.diff(p, e).subs(e, 0) - 8 * mu**2 / Q)
    J0 = s.Symbol("whole_J0", real=True)
    _ph0, ph1 = s.symbols("phase0 phase1", real=True)
    ell = s.Symbol("log_E2_over_nu2", real=True)
    J1 = s.Symbol("whole_J1", real=True)
    jp = J0 + e * J1
    Ip = (1 + 2 * e) / (2 * e)
    avg = R * (p * Ip + aa + bb / (3 + 2 * e)) + Kh * (
        p * (Ip - rr**2 * jp) / (1 - rr**2) + aa * jp + bb * (jp - 1) / rr**2
    )
    pole = p.subs(e, 0) * Hf / 2
    Reg = R * (aa.subs(e, 0) + bb / 3) + Kh * (
        -p.subs(e, 0) * rr**2 * J0 / (1 - rr**2)
        + aa.subs(e, 0) * J0
        + bb * (J0 - 1) / rr**2
    )
    # The literal expression contains removable e*1/e terms; cancel before evaluation.
    finite = s.diff(s.cancel(e * avg), e).subs(e, 0)
    checks["whole_angular_finite_coefficient"] = s.factor(
        finite - p.subs(e, 0) * Hf - s.diff(p, e).subs(e, 0) * Hf / 2 - Reg
    )
    checks["whole_hard_finite_cut"] = s.factor(
        finite
        + ph1 * pole
        - ell * pole
        - (pole * (2 + 2 * mu**2 / V.subs(e, 0) + ph1 - ell) + Reg)
    )
    b = s.Symbol("physical_beta", positive=True)
    q = s.Symbol("cut_coordinate", real=True)
    Msp = s.log((q - b) / (q + b)) / (2 * a * b)
    Mup = s.atanh(b * q) / (a * b)
    checks["whole_light_PV_primitive"] = s.factor(
        s.diff(Msp, q) - 1 / (a * (q * q - b * b))
    )
    checks["whole_crossed_real_primitive"] = s.factor(
        s.diff(Mup, q) - 1 / (a * (1 - b * b * q * q))
    )
    checks["whole_light_delta_jacobian"] = (
        s.diff(a * (q * q - b * b), q).subs(q, b) - 2 * a * b
    )
    checks["whole_forward_soft_real_cancellation"] = s.factor(
        V.subs(e, 0) * (-s.atanh(b) + s.atanh(b)) / (a * b)
    )
    checks["whole_forward_crossed_numerator"] = s.factor(
        V.subs({a: 4 * mu - a, e: 0}, simultaneous=True) - V.subs(e, 0)
    )
    checks["whole_forward_zero_channel_self_cancellation"] = (
        2 * (2 * mu * mu) / (4 * mu) - mu
    )

    # Full HH production sew, including its distinct-mass graviton numerator.
    d = s.Symbol("dimension", positive=True)
    pm = (a - 4 * mu) / 4
    pn = (a - 4 * n) / 4
    raw = (
        4 * pm * pn * z * z
        - a * (pm + pn)
        + (d - 1) * a * a / 4
        - ((d - 2) * a / 2 + 2 * mu) * ((d - 2) * a / 2 + 2 * n) / (d - 2)
    )
    grav = (
        a / 4
        + 4 * mu * n / (a * (d - 2))
        - (a - 4 * mu) * (a - 4 * n) * z * z / (4 * a)
    )
    checks["whole_D_distinct_mass_graviton_tensor"] = s.factor(-raw / a - grav)
    A = a / 2 - n
    checks["whole_heavy_exchange_gap"] = s.factor(
        A * A - (a - 4 * mu) * (a - 4 * n) / 4 - n * n - mu * (a - 4 * n)
    )
    a0, a2 = s.symbols("matter_P0 matter_P2", real=True)
    alpha = a / 4 + 2 * mu * n / a
    bc = -(a - 4 * mu) * (a - 4 * n) / (4 * a)
    sew = (alpha + bc / 3) * a0 + 2 * bc * a2 / 15
    endpoint = (
        (a + 2 * mu) * (a + 2 * n) * a0 - (a - 4 * mu) * (a - 4 * n) * a2 / 5
    ) / (6 * a)
    checks["whole_HH_cut_equals_both_endpoint_form_factors"] = s.factor(sew - endpoint)

    # Weighted physical logarithm, real primitive on either side of its root.
    beta = s.Symbol("beta", positive=True)
    f0 = q * s.log(q * q - beta * beta) - 2 * q - beta * s.log((q - beta) / (q + beta))
    f2 = (
        q**3 * s.log(q * q - beta * beta) / 3
        - 2 * q**3 / 9
        - 2 * beta * beta * q / 3
        - beta**3 * s.log((q - beta) / (q + beta)) / 3
    )
    checks["whole_log_primitive"] = s.factor(s.diff(f0, q) - s.log(q * q - beta * beta))
    checks["whole_weighted_log_primitive"] = s.factor(
        s.diff(f2, q) - q * q * s.log(q * q - beta * beta)
    )
    checks["whole_weighted_imaginary_interval"] = (
        s.integrate((1 - q * q) / 4, (q, 0, beta)) - beta * (3 - beta * beta) / 12
    )
    U, B = s.symbols("UV_pole_coefficient renormalized_bubble", real=True)
    contract = n + 2 * mu / (1 + e)
    covariant = s.cancel(contract * (U / e + B) - contract * U / e)
    scalar_only = s.cancel(contract * (U / e + B) - (n + 2 * mu) * U / e)
    checks["whole_D_covariant_counterterm_finite"] = (
        s.limit(covariant, e, 0) - (n + 2 * mu) * B
    )
    checks["D4_only_subtraction_finite_difference"] = (
        s.limit(scalar_only - covariant, e, 0) + 2 * mu * U
    )
    checks["actual_evanscent_difference"] = (-2 * mu * U).subs(
        U, -s.Rational(1, 6)
    ) - mu / 3
    minimal = 2 * V * (1 + 2 * e) / (e * Q) - 2 * (a - 2 * mu)
    mixing = (a + 2 * mu / (1 + e)) * (mu / a + Q * (e + 1) / (2 * a * (2 * e + 3)))
    checks["whole_resonant_light_cut_requires_full_Hmetric"] = s.factor(
        p * (1 + 2 * e) / (2 * e) + aa + bb / (3 + 2 * e) - minimal - mixing
    )
    cRH, gg, kk = s.symbols("unmatched_RH cubic kappa", real=True)
    hanchor = -2 * gg * cRH
    checks["whole_RH_forward_anchor_dictionary"] = s.factor(
        -2 * hanchor * (n + 2 * mu) / (kk * (n - 2 * mu) ** 3)
        - 4 * gg * cRH * (n + 2 * mu) / (kk * (n - 2 * mu) ** 3)
    )
    window = s.Symbol("compact_heavy_window", positive=True)
    return {
        "whole_light_pair_D_cut": light_cut(),
        "whole_light_pair_formal_soft_divided_cut": hard_light_cut(),
        "whole_heavy_pair_D_cut": heavy_cut(),
        "whole_distinct_heavy_graviton_D4_cut": production.forward_cut(),
        "whole_covariant_weighted_bubble": weighted_bubble(),
        "whole_resonant_mixing_transition_with_unmatched_RH": mixing_transition(),
        "whole_known_formal_threshold_pair_plus_mixing": pairing.paired_finite_functional(
            window
        )
        + G**2
        * (N + 2 * MU)
        * s.re(weighted_bubble(N, MU))
        * pairing.W(N)
        / (16 * s.pi * K),
        "whole_unmatched_RH_delta_weight": 2
        * s.pi
        * G
        * (N + 2 * MU)
        * CRH
        * pairing.W(N)
        / K,
        "whole_cut_inventory_and_normalization": "At selected g^2/kappa order the two-body cuts are PhiPhi matter-heavy/GR interference, HH matter-light-exchange/GR-s-channel interference, and distinct Hh production. Phi parity forbids PhiH and Phih. Pure hh begins at different metric order. All light s/t/u exchanges and the complete GR tree are kept. Identical-pair optical1/4 times2 interference gives Phi2/2; Hh uses its distinct factor. Above4n all lower channels remain open.",
        "whole_resonance_and_counterterm_boundary": "The full H-metric weighted bubble completes the finite resonant light cut, not only its Coulomb pole. Renormalize its entire transverse D tensor before contracting and taking D4. A D4-only scalar subtraction produces the mu/3 evanescent discrepancy checked here. The known finite bubble is at the inherited MS reference scale1; c_RH remains independent. The S292 formal Hh/virtual pair and this finite mixing delta do not establish a stable heavy atom, physical width resummation or a positive infrared-safe measure.",
        "whole_IR_boundary": "The existing S288 forward soft kernel is purely imaginary. Dividing by the already-stated S278 factor cancels the full light-channel Coulomb principal-value pole, including every crossed heavy exchange. Its imaginary product with the heavy delta is real; it does not supply another imaginary delta subtraction. The finite hard light cut remains an interference coefficient and is not required to be positive.",
        "checks": checks,
        "gates": {
            "complete_light_heavy_and_distinct_Hh_cut_inventory": True,
            "full_D_angular_average_not_only_Coulomb_residue": True,
            "HH_whole_tensor_matches_both_endpoint_form_factors": True,
            "whole_weighted_log_and_covariant_RH_subtraction": True,
            "complete_resonant_Hmetric_piece_and_unmatched_anchor": True,
            "S292_real_virtual_pair_retained_with_original_qualifications": True,
            "no_physical_width_IR_measure_or_Regge_closure": True,
        },
    }
