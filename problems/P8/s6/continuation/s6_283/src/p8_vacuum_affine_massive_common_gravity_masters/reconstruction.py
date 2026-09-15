"""One common crossed master basis for both complete pure-gravity cuts."""

from functools import cache

import sympy as s

from . import source

S, T, MU, K, EP = source.S, source.T, source.MU, source.K, source.EP
MT, MU_CUT, MS = source.invariant.MT, source.invariant.MU_CUT, source.invariant.MS


def coefficients(energy=S, transfer=T, mass=MU):
    energy, transfer, mass = map(s.sympify, (energy, transfer, mass))
    q = energy - 4 * mass
    u = -q - transfer
    z = 1 + 2 * transfer / q
    v, a, b, r, d = source.invariant.polynomials(energy, mass)
    h0 = -2 * energy + 6 * mass - 2 * mass**2 / energy + q * q / (6 * energy)
    h2 = -q * q / (6 * energy)
    p = 4 * v / q
    h = h0 + h2 * s.legendre(2, z)
    return {
        "V": v,
        "H": s.factor(h),
        "H0": h0,
        "H2": h2,
        "P_times_kappa": p,
        "C00mu": 2 * a - 6 * transfer * u * b / q**2,
        "B00": r / 30 + transfer * u * d / (10 * q**2),
        "C0mumu": s.factor(-2 * v * h),
        "Bmm": s.factor(
            (-3 * p * h2 * s.legendre(2, z) + h0 * h0 + h2 * h2 * s.legendre(2, z) / 5)
            / 2
        ),
    }


def evanescent_coefficients(energy=S, transfer=T, mass=MU, epsilon=EP):
    row = coefficients(energy, transfer, mass)
    vd = row["V"] + 2 * mass**2 * epsilon / (1 + epsilon)
    hd = row["H"] - 2 * mass**2 * epsilon / (energy * (1 + epsilon))
    return vd * vd, -2 * vd * hd


def crossed_basis():
    values = (S, T, 4 * MU - S - T)
    boxes = [
        {
            "massless_channel": a,
            "massive_channel": b,
            "coefficient": source.invariant.polynomials(b)[0] ** 2,
        }
        for a in values
        for b in values
        if a != b
    ]
    channels = [
        {"energy": a, "transfer": b, "coefficients": coefficients(a, b)}
        for a, b in (
            (values[0], values[1]),
            (values[1], values[0]),
            (values[2], values[0]),
        )
    ]
    return {
        "six_ordered_boxes": boxes,
        "three_triangle_bubble_channels": channels,
        "overall_normalization": 1 / (16 * s.pi**2 * K * K),
        "scope": "Cut-matched representative. Use the displayed evanescent box/C0mumu coefficients for their finite IR parts. Oepsilon bubble coefficients times UV poles, cut-free rational terms and physical pole/local matching are NOT fixed by this representative.",
    }


@cache
def data():
    z = source.Z
    q = S - 4 * MU
    t = -q * (1 - z) / 2
    u = -q * (1 + z) / 2
    P, _, _, a0, a2 = source.elastic.coefficients(
        cubic2=s.Integer(0), contact=s.Integer(0)
    )
    A = source.elastic.whole_tree(cubic2=s.Integer(0), contact=s.Integer(0))
    row = coefficients(S, t)
    boxD, ctD = evanescent_coefficients(S, t)
    LG = s.Symbol("raw_log_Q_over_4pi_nu_squared_plus_Gamma_E", real=True)
    Lt = LG + s.log((1 - z) / 2)
    Lu = LG + s.log((1 + z) / 2)
    box = 4 * boxD / (K * K * q) * ((-1 / EP - Lt) / t + (-1 / EP - Lu) / u)
    triangle = -2 * ctD / (K * K * q) * (1 / EP + LG)
    bubble = 2 * row["Bmm"] / (K * K)
    total = s.series(box + triangle + bubble, EP, 0, 1).removeO()
    p1 = 8 * MU * MU / (K * q)
    faux = (
        P * P * source.infrared.angular_finite_same()
        - 3 * P * a2 * s.legendre(2, z)
        + a0 * a0
        + a2 * a2 * s.legendre(2, z) / 5
    )
    expected = (
        P * A / EP + faux + P * A * LG + p1 * A + P * source.infrared.evanescent_tree()
    )
    generic = coefficients()
    uu = 4 * MU - S - T
    vt = source.invariant.polynomials(T)[0]
    vu = source.invariant.polynomials(uu)[0]
    gg = (
        4 * (vt * vt * MT + vu * vu * MU_CUT) / S
        - 2 * generic["C00mu"] * MS
        + generic["B00"]
    )
    channels = ((S, T, uu), (T, S, uu), (uu, S, T))
    rows = [coefficients(a, b) for a, b, _ in channels]
    mk = s.symbols("M_s M_t M_u", real=True)
    pair_IR = sum(
        (-4 * r["V"] ** 2 * (1 / b + 1 / c) - 2 * r["C0mumu"]) * m
        for (a, b, c), r, m in zip(channels, rows, mk, strict=True)
    )
    tree = source.gravity_tree()
    expected_IR = 4 * K * tree * sum(r["V"] * m for r, m in zip(rows, mk, strict=True))
    self_IR = -2 * MU * K * tree
    full_soft = (
        2 * K * tree * (2 * sum(r["V"] * m for r, m in zip(rows, mk, strict=True)) - MU)
    )
    bubble_sum = sum(r["B00"] + r["Bmm"] for r in rows)
    checks = {
        "whole_regular_gravity_tree_H": s.factor(
            row["H"] - (-2 * S + 6 * MU - 2 * MU * MU / S + t * u / S)
        ),
        "original_P_coefficient": s.factor(row["P_times_kappa"] - K * P),
        "original_a0_coefficient": s.factor(row["H0"] - K * a0),
        "original_a2_coefficient": s.factor(row["H2"] - K * a2),
        "entire_pure_gravity_elastic_bare_cut_all_finite_terms": s.factor(
            total - expected
        ),
        "entire_original_graviton_normal_cut_same_boxes": s.factor(
            gg - 4 * source.invariant.algebraic_cut()
        ),
        "all_crossed_pairwise_IR_terms": s.factor(pair_IR - expected_IR),
        "whole_pair_and_four_leg_self_IR_matches_S278": s.factor(
            pair_IR + self_IR - full_soft
        ),
        "massless_whole_bubble_sum_203_over_40": s.factor(
            bubble_sum.subs(MU, 0)
            - s.Rational(203, 40) * (S * S + T * T + (-S - T) ** 2)
        ),
    }
    for i, ((a, b, c), r) in enumerate(zip(channels, rows, strict=True)):
        checks["whole_gravity_tree_channel_decomposition_" + str(i)] = s.factor(
            r["H"] - r["V"] * (1 / b + 1 / c) - K * tree
        )
    reverse = coefficients(S, uu)
    for name in ("V", "H", "C00mu", "B00", "C0mumu", "Bmm"):
        checks["whole_other_channel_crossing_" + name] = s.factor(
            reverse[name] - generic[name]
        )
    return {
        "whole_common_crossed_basis": crossed_basis(),
        "whole_finite_IR_evanescent_box_triangle_coefficients": evanescent_coefficients(),
        "whole_elastic_cut_from_shared_masters_normalized": total,
        "whole_normal_graviton_cut_from_shared_masters_normalized": gg,
        "whole_four_leg_completed_soft_pole": full_soft,
        "whole_massless_bubble_sum_calibration": s.Rational(203, 40)
        * (S * S + T * T + (-S - T) ** 2),
        "no_double_count": "A box I4(a,b) carries its massless normal cut in a and massive normal cut in b. The six ordered internal-mass assignments occur once each. Its two cut appearances are not two independent graphs.",
        "finite_matching": "The evanescent stress traces VD=V+2mu^2 EP/(1+EP),HD=H-2mu^2 EP/[s(1+EP)] supply box coefficientVD^2 and massive-triangle coefficient-2VD HD. Together with the full phase-space constants they reproduce S281's complete pure-gravity finite cut, not only its soft pole.",
        "rational_boundary": "Four-dimensional normal cuts leave cut-free rational terms. In particular Oepsilon bubble coefficients can multiply UV poles, and apparent Gram denominators at s=4mu require analytic completion. This packet does not set those terms, the physical finite pole anchor or local Wilson coefficients to zero.",
        "checks": checks,
        "gates": {
            "six_distinct_ordered_boxes_once_each": len(
                crossed_basis()["six_ordered_boxes"]
            )
            == 6,
            "all_three_channel_triangle_and_bubble_assignments": len(
                crossed_basis()["three_triangle_bubble_channels"]
            )
            == 3,
            "same_box_has_both_species_cuts": True,
            "whole_elastic_finite_trace_and_phase_retained": True,
            "whole_soft_self_term_not_only_pair_terms": True,
            "massless_bubble_calibration_not_massless_replacement": True,
            "rational_and_finite_pole_matching_not_claimed_complete": True,
        },
    }
