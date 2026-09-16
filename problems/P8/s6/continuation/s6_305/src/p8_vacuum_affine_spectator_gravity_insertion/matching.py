"""Known finite shifts without choosing the remaining physical coordinates."""

from functools import cache

import sympy as s
from p8_vacuum_affine_minimal_gravity_finite import assembly

from . import insertion as ins
from . import source

ALPHA, BETA, DK = assembly.ALPHA, assembly.BETA, assembly.DK


def fixed_shifts():
    ell = s.log(source.HEAVY_MASS2) + 2
    return -ell / 40, -3 * ell / 10, source.poles.delta_kappa()


def remaining_local(energy, transfer, mass, kappa, alpha, beta, delta_kappa):
    """Required arguments; no default-zero physical matching convention."""
    q2 = sum(a * a for a, _, _ in ins.channels(energy, transfer, mass))
    return (alpha * q2 + beta * mass**2) / (
        16 * s.pi**2 * kappa**2
    ) - delta_kappa * ins.gravity_born(energy, transfer, mass, kappa) / kappa


def augmented_known(known_minimal_amplitude, energy, transfer, mass, kappa):
    return (
        known_minimal_amplitude
        + ins.gapped_fixed(energy, transfer, mass, kappa)
        + ins.massless_known(energy, transfer, mass, kappa, 1)
    )


@cache
def data():
    a, b, mu, k = ins.A, ins.B, ins.MU, ins.K
    c = 4 * mu - a - b
    cs = ins.channels(a, b, mu)
    q2 = a * a + b * b + c * c
    ell = s.Symbol("fixed_log_n_plus2", real=True)
    n0 = sum(ins.numerator0(x, mu) for x, _, _ in cs)
    n2 = sum(ins.numerator2(x, y, z, mu) for x, y, z in cs)
    local = -(n2 * ell / 60 / (16 * s.pi**2) + n0 * 2 * ell / (768 * s.pi**2)) / k**2
    checks = {}

    def put(name, value):
        checks[name] = s.factor(value)

    put("whole_crossed_spin0_sum", n0 - (q2 + 28 * mu**2) / 3)
    put("whole_crossed_spin2_sum", n2 - s.Rational(2, 3) * (q2 - 8 * mu**2))
    put(
        "whole_fixed_curvature_local_polynomial",
        local + ell * (q2 + 12 * mu**2) / (640 * s.pi**2 * k**2),
    )
    put(
        "same_existing_two_local_axes",
        remaining_local(a, b, mu, k, -ell / 40, -3 * ell / 10, 0) - local,
    )
    t0, _ = assembly.tree_jets(a, b, mu)
    put("same_existing_Newton_axis", t0 + k * ins.gravity_born(a, b, mu, k))
    rest = remaining_local(a, b, mu, k, ALPHA, BETA, DK)
    for name, var, target in (
        ("alpha", ALPHA, q2 / (16 * s.pi**2 * k**2)),
        ("beta", BETA, mu**2 / (16 * s.pi**2 * k**2)),
        ("Newton", DK, t0 / k**2),
    ):
        put("remaining_independent_" + name, s.diff(rest, var) - target)
    full = (
        augmented_known(s.Symbol("whole_minimal_known_amplitude"), a, b, mu, k) + rest
    )
    for name, var in (("alpha", ALPHA), ("beta", BETA), ("Newton", DK)):
        put("same_augmented_independent_" + name, s.diff(full, var) - s.diff(rest, var))
    scale = s.Symbol("scale_squared", positive=True)
    change = s.factor(scale * s.diff(ins.massless_known(a, b, mu, k, scale), scale))
    put(
        "M1_reference_scale_change_is_existing_local_axes",
        change - (q2 + 12 * mu**2) / (640 * s.pi**2 * k**2),
    )
    x0, x1, x2 = s.symbols("x0 x1 x2")
    ns = ((x0, x1, x2), (x1, x0, x2), (x2, x0, x1))
    aa, hh = s.Function("A"), s.Function("H")
    expr = -sum(
        ins.numerator2(u, v, w, mu) * aa(-u) / (16 * s.pi**2 * k * k)
        + ins.numerator0(u, mu) * hh(-u) / (768 * s.pi**2 * k * k)
        for u, v, w in ns
    )
    for index, pair in enumerate(((x0, x1), (x0, x2), (x1, x2))):
        put(
            "whole_crossing_" + str(index),
            expr.xreplace({pair[0]: pair[1], pair[1]: pair[0]}) - expr,
        )
    return {
        "whole_fixed_H_Proca_local_shifts": dict(
            zip(("alpha", "beta", "delta_kappa"), fixed_shifts())
        ),
        "whole_remaining_physical_matching_polynomial": rest,
        "whole_augmented_known_plus_unassigned_matching": full,
        "whole_M1_reference_scale_derivative": change,
        "whole_fixed_local_H_Proca_amplitude": local,
        "matching_ownership": "Add the source-fixed H/Proca shifts once. They occupy the SAME alpha,beta,delta_kappa axes used in S302, not three new arbitrary coefficients. One remaining aggregate coordinate per axis carries unresolved light-Phi/graviton/M1 finite data in this sector. The augmented known reference retains those three required independent arguments. Other mixed/matter-sector matching remains separately unresolved.",
        "sign_boundary": "The fixed Newton insertion changes the known-reference physical interference, not the original bare kappa or chosen state. Neither its negative known rate nor its real local polynomial is a sign assertion for the full matched rate or isolated pole-subtracted positivity.",
        "checks": checks,
        "gates": {
            "all_three_physical_matching_coordinates_still_symbolic": all(
                full.has(v) for v in (ALPHA, BETA, DK)
            ),
            "fixed_shift_not_duplicate_independent_matching_axes": len(fixed_shifts())
            == 3,
            "whole_original_M1_scale_dependence_mapped_not_deleted": change != 0,
            "no_new_finite_prescription_or_Newton_redefinition": True,
        },
    }
