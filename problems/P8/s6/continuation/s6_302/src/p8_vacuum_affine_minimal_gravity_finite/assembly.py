"""Full finite minimal-gravity representative with three symbolic anchors."""

from functools import cache

import sympy as s

from . import masters, source

dimensional, poles, threshold = source.dimensional, source.poles, source.threshold
A, B, MU, ELL, EPS = s.symbols("energy transfer positive_mass raw_log epsilon")
ALPHA, BETA, DK = s.symbols("finite_alpha finite_beta finite_delta_kappa")


def channels(energy, transfer, mass):
    a, b, mu = map(s.sympify, (energy, transfer, mass))
    c = 4 * mu - a - b
    return ((a, b, c), (b, a, c), (c, a, b))


def eikonal(channel, mass):
    return poles.eikonal(channel, mass, 4)


def soft_kernel(energy, transfer, mass, rows):
    return (
        sum(
            eikonal(a, mass) * rows[i]["J"]
            for i, (a, _, _) in enumerate(channels(energy, transfer, mass))
        )
        / 2
        - mass
    )


def tree_jets(energy, transfer, mass):
    td = poles.newton_shape(energy, transfer, mass, 4 + 2 * EPS)
    return td.subs(EPS, 0), s.diff(td, EPS).subs(EPS, 0)


def completion_jets(energy, transfer, mass):
    rd = poles.completion(energy, transfer, mass, 4 + 2 * EPS)
    return rd.subs(EPS, 0), s.diff(rd, EPS).subs(EPS, 0)


def uv_polynomial(energy, transfer, mass):
    vals = [row[0] for row in channels(energy, transfer, mass)]
    return (
        s.Rational(203, 40) * sum(a * a for a in vals)
        - s.Rational(169, 3) * mass * mass
    )


def known_components(energy, transfer, mass, ell, rows):
    mu = s.sympify(mass)
    cs = channels(energy, transfer, mu)
    vals = [row[0] for row in cs]
    vs = [eikonal(a, mu) for a in vals]
    result = {
        "ordered_boxes": s.Integer(0),
        "massless_triangles": s.Integer(0),
        "massive_triangles": s.Integer(0),
        "bubbles": s.Integer(0),
        "old_Gram_finite_subtraction": s.Integer(0),
    }
    for i, a in enumerate(vals):
        for j, _ in enumerate(vals):
            if i != j:
                result["ordered_boxes"] += (
                    vs[j] ** 2 * rows[j]["J"] * rows[i]["la"]
                    - 4 * mu * mu * vs[j] * rows[j]["J"]
                ) / a
    for i, (a, b, _) in enumerate(cs):
        row = dimensional.master_coefficients(a, b, mu, EPS)
        finite = masters.finite_dictionary(rows[i], ell)
        result["massless_triangles"] += row["C00mu"].subs(EPS, 0) * finite["C00mu"]
        result["massive_triangles"] += (
            row["C0mumu"].subs(EPS, 0) * finite["C0mumu"]
            - s.diff(row["C0mumu"], EPS).subs(EPS, 0) * rows[i]["J"] / 2
        )
        result["bubbles"] += sum(
            row[name].subs(EPS, 0) * finite[name] for name in ("B00", "Bmm")
        )
        result["old_Gram_finite_subtraction"] -= threshold.old_gram_choice(
            a, b, mu, bubble=ell + 2
        )
    t0, t1 = tree_jets(energy, transfer, mu)
    r0, r1 = completion_jets(energy, transfer, mu)
    result["compact_evanescent_bubbles_and_Gram"] = dimensional.compact_crossed_evans(
        energy, transfer, mu
    )
    result["four_external_LSZ_finite"] = 6 * mu * ell * t0 + 14 * mu * t0 - 6 * mu * t1
    result["physical_pole_completion_finite"] = -ell * r0 + r1
    result["analytic_soft_division_tree_derivative"] = (
        2 * t1 * soft_kernel(energy, transfer, mu, rows)
    )
    return result


def known_finite(energy, transfer, mass, ell, rows):
    return sum(known_components(energy, transfer, mass, ell, rows).values())


def complete_finite(energy, transfer, mass, ell, rows, kappa, alpha, beta, delta_kappa):
    """Matching coordinates are required arguments; no implicit zero assignment."""
    t0, _ = tree_jets(energy, transfer, mass)
    local = (
        alpha * sum(a * a for a, _, _ in channels(energy, transfer, mass))
        + beta * mass * mass
    )
    return (
        known_finite(energy, transfer, mass, ell, rows)
        + local
        + 16 * s.pi**2 * delta_kappa * t0
    ) / (16 * s.pi**2 * kappa * kappa)


def raw_laurent_coefficients(energy, transfer, mass, ell, rows):
    """Independent full-D coefficient times master pole/finite decomposition."""
    mu = s.sympify(mass)
    cs = channels(energy, transfer, mu)
    vals = [row[0] for row in cs]
    vd = [poles.eikonal(a, mu, 4 + 2 * EPS) for a in vals]
    principal, finite = s.Integer(0), s.Integer(0)
    for i, a in enumerate(vals):
        for j, _ in enumerate(vals):
            if i != j:
                principal -= vd[j] ** 2 * rows[j]["J"] / a
                finite += vd[j] ** 2 * rows[j]["J"] * rows[i]["la"] / a
    for i, (a, b, _) in enumerate(cs):
        coeff = dimensional.master_coefficients(a, b, mu, EPS)
        ms = masters.finite_dictionary(rows[i], ell)
        principal -= coeff["C0mumu"] * rows[i]["J"] / 2 + coeff["B00"] + coeff["Bmm"]
        finite += sum(coeff[name] * ms[name] for name in ms)
        principal -= threshold.full_gram_choice(a, b, mu, bubble=-1)
        finite -= threshold.full_gram_choice(a, b, mu, bubble=ell + 2)
        # full_gram_choice contains a fixed rational part; that belongs only to finite.
        principal += threshold.full_gram_choice(a, b, mu, bubble=0)
    td = poles.newton_shape(energy, transfer, mu, 4 + 2 * EPS)
    rd = poles.completion(energy, transfer, mu, 4 + 2 * EPS)
    residue = 2 * mu * poles.residue_ratio(4 + 2 * EPS) * td - rd
    principal -= residue
    finite += ell * residue
    return principal, finite


@cache
def data():
    a, b, mu, ell = A, B, MU, ELL
    cs = channels(a, b, mu)
    rows = [masters.symbols("channel" + str(i)) for i in range(3)]
    vs = [eikonal(x, mu) for x, _, _ in cs]
    coeffs = [dimensional.master_coefficients(x, y, mu, EPS) for x, y, _ in cs]
    t0, t1 = tree_jets(a, b, mu)
    r0, r1 = completion_jets(a, b, mu)
    bsoft = soft_kernel(a, b, mu, rows)
    uv = uv_polynomial(a, b, mu)
    infrared = -sum(
        vs[j] ** 2 * rows[j]["J"] / cs[i][0]
        for i in range(3)
        for j in range(3)
        if i != j
    )
    infrared -= sum(
        coeffs[i]["C0mumu"].subs(EPS, 0) * rows[i]["J"] / 2 for i in range(3)
    )
    infrared += 2 * mu * t0
    principal, finite = raw_laurent_coefficients(a, b, mu, ell, rows)
    ultraviolet = s.factor(principal.subs(EPS, 0) - infrared)
    known = known_finite(a, b, mu, ell, rows)
    all_finite = (
        s.diff(principal, EPS).subs(EPS, 0) + finite.subs(EPS, 0) + 2 * t1 * bsoft
    )
    c = cs[2][0]
    checks = {
        "whole_tagged_IR_residue_including_four_legs": s.factor(
            infrared + 2 * t0 * bsoft
        ),
        "whole_tagged_UV_residue_is_local": s.factor(ultraviolet + uv),
        "whole_UV_then_IR_pole_subtraction": s.factor(
            principal.subs(EPS, 0) + uv + 2 * t0 * bsoft
        ),
        "whole_finite_all_D_Laurent_assembly": s.factor(all_finite - known),
        "whole_tree_epsilon_coefficient": s.factor(
            t1 - 2 * mu * mu * sum(1 / x for x, _, _ in cs)
        ),
        "whole_tree_D4_regular_constant": s.factor(
            t0
            - (
                2 * mu * mu * sum(1 / x for x, _, _ in cs)
                - sum(y * z / x for x, y, z in cs)
                - 6 * mu
            )
        ),
        "LSZ_residue_ratio_zero": poles.residue_ratio(4) - 3,
        "LSZ_residue_ratio_first": s.diff(poles.residue_ratio(4 + 2 * EPS), EPS).subs(
            EPS, 0
        )
        + 7,
        "raw_D_tree_soft_finite_subtraction": s.diff(
            2 * poles.newton_shape(a, b, mu, 4 + 2 * EPS) * bsoft, EPS
        ).subs(EPS, 0)
        - 2 * t1 * bsoft,
        "resolution_change_at_fixed_reference": s.diff(
            2
            * poles.newton_shape(a, b, mu, 4 + 2 * EPS)
            * bsoft
            * s.exp(2 * EPS * s.Symbol("log_resolution")),
            EPS,
        ).subs(EPS, 0)
        - 2 * t1 * bsoft
        - 4 * t0 * bsoft * s.Symbol("log_resolution"),
        "old_and_evanescent_Gram_both_retained": s.factor(
            sum(
                -s.diff(coeffs[i]["B00"] + coeffs[i]["Bmm"], EPS).subs(EPS, 0)
                - threshold.evanescent_gram_choice(x, y, mu)
                for i, (x, y, _) in enumerate(cs)
            )
            - dimensional.compact_crossed_evans(a, b, mu)
        ),
        "full_finite_Newton_shift_sign": s.factor(
            s.diff(-t0 / (s.Symbol("kappa") + DK), DK).subs(DK, 0)
            - t0 / s.Symbol("kappa") ** 2
        ),
    }
    for i, (x, y, z) in enumerate(cs):
        vd = poles.eikonal(x, mu, 4 + 2 * EPS)
        h = -2 * x + 6 * mu - 2 * mu * mu / x + y * z / x
        cm1 = 4 * mu * mu * (eikonal(x, mu) / x - h)
        checks["box_full_D_derivative_" + str(i)] = s.factor(
            s.diff(vd**2, EPS).subs(EPS, 0) - 4 * mu * mu * eikonal(x, mu)
        )
        checks["massive_triangle_compact_D4_" + str(i)] = s.factor(
            coeffs[i]["C0mumu"].subs(EPS, 0) + 2 * eikonal(x, mu) * h
        )
        checks["massive_triangle_first_derivative_" + str(i)] = s.factor(
            s.diff(coeffs[i]["C0mumu"], EPS).subs(EPS, 0) - cm1
        )
    coefficients0 = tuple(v.subs(poles.D, 4) for v in poles.completion_coefficients())
    coefficients1 = tuple(
        2 * s.diff(v, poles.D).subs(poles.D, 4) for v in poles.completion_coefficients()
    )
    expected0 = (
        s.Rational(164, 15),
        -s.Rational(2, 15),
        -s.Rational(73, 15),
        s.Rational(1, 15),
    )
    expected1 = (
        -s.Rational(4879, 225),
        -s.Rational(28, 225),
        s.Rational(3128, 225),
        -s.Rational(16, 225),
    )
    for i in range(4):
        checks["physical_completion_constant_" + str(i)] = (
            coefficients0[i] - expected0[i]
        )
        checks["physical_completion_first_" + str(i)] = coefficients1[i] - expected1[i]
    for label, order in (
        ("swap_st", (1, 0, 2)),
        ("swap_su", (2, 1, 0)),
        ("swap_tu", (0, 2, 1)),
    ):
        vals = [cs[i][0] for i in order]
        permuted = [rows[i] for i in order]
        checks["whole_finite_crossing_" + label] = s.factor(
            known_finite(vals[0], vals[1], mu, ell, permuted) - known
        )
    kapp = s.Symbol("positive_kappa", positive=True)
    matched = complete_finite(a, b, mu, ell, rows, kapp, ALPHA, BETA, DK)
    for key, variable, expected in (
        ("alpha", ALPHA, (a * a + b * b + c * c) / (16 * s.pi**2 * kapp * kapp)),
        ("beta", BETA, mu * mu / (16 * s.pi**2 * kapp * kapp)),
        ("Newton", DK, t0 / (kapp * kapp)),
    ):
        checks["independent_finite_anchor_" + key] = s.factor(
            s.diff(matched, variable) - expected
        )
    lambdas = s.symbols("lambda0:3")
    scale_known = known.subs({rows[i]["la"]: ell - lambdas[i] for i in range(3)})
    checks["whole_common_log_derivative"] = s.factor(
        s.diff(scale_known, ell) - 2 * t0 * bsoft - uv
    )
    return {
        "whole_ordered_master_symbols": rows,
        "whole_D_master_coefficients": coeffs,
        "whole_tagged_UV_residue": ultraviolet,
        "whole_tagged_IR_residue": infrared,
        "whole_local_UV_polynomial": uv,
        "whole_tree_and_physical_completion_jets": (t0, t1, r0, r1),
        "whole_known_finite_components": known_components(a, b, mu, ell, rows),
        "whole_matched_finite_representative": matched,
        "raw_convention_and_reference": "D=4+2epsilon with unchanged S283 raw factors. "
        "Subtract the local UV pole before dividing the analytic soft factor "
        "at resolution/nu=1. The D-tree derivative leaves +2*T1*Bsoft. "
        "This fixes coordinates, not numerical matching values.",
        "finite_matching": "alpha*(s^2+t^2+u^2)+beta*mu^2+16*pi^2*delta_kappa*T0 "
        "inside the overall1/(16*pi^2*kappa^2); all three anchors unassigned.",
        "checks": checks,
        "gates": {
            "all_six_ordered_boxes_and_three_crossed_channels": True,
            "both_triangle_species_and_both_bubbles": True,
            "all_D_coefficients_and_finite_evanescent_terms_retained": True,
            "old_Gram_finite_not_confused_with_evanescent_Gram": True,
            "four_LSZ_and_physical_completion_finite_terms_retained": True,
            "UV_tag_removed_before_IR_division": True,
            "D_tree_soft_derivative_not_dropped": True,
            "three_physical_matching_coordinates_required_not_chosen": True,
        },
    }
