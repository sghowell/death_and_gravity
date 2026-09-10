"""Canonical field factors, fixed light LSZ and the even-vacuum connected kernel."""

from functools import cache

import sympy as sp
from p8_vacuum_two_loop_light_pole import calibration as pole
from p8_vacuum_two_loop_light_pole import subtraction as pole_subtraction


def truncate(expr, parameter):
    return sp.series(expr, parameter, 0, 3).removeO().expand()


@cache
def data():
    hbar = sp.Symbol("formal_loop_order")
    z1, z2 = sp.symbols("bare_field_factor_first bare_field_factor_second")
    Z = 1 + hbar * z1 + hbar**2 * z2
    A0, A1, A2 = sp.symbols(
        "canonical_tree_amplitude canonical_one_loop_amplitude canonical_two_loop_amplitude"
    )
    ren = A0 + hbar * A1 + hbar**2 * A2
    bare_amputated = truncate(ren / Z**2, hbar)
    bare_lsz = truncate(Z**2, hbar)
    assembled = truncate(bare_lsz * bare_amputated, hbar)
    L, G, M = sp.symbols("canonical_L canonical_G canonical_M")
    L1, L2, G1, G2 = sp.symbols(
        "total_canonical_L1 total_canonical_L2 total_canonical_G1 total_canonical_G2"
    )
    Ltotal = L + hbar * L1 + hbar**2 * L2
    Gtotal = G + hbar * G1 + hbar**2 * G2
    Lbare = truncate(Ltotal / Z**2, hbar)
    Gbare = truncate(Gtotal / Z, hbar)
    gbare = truncate(Gbare**2, hbar)
    gtotal = truncate(Gtotal**2, hbar)
    phi, H = sp.symbols("light_field heavy_field")
    V = phi**2 / 2 + M * H**2 / 2 + G * H * phi**2 / 2 + L * phi**4 / 24
    ct_mass, ct_kinetic, ct_heavy_source = sp.symbols(
        "local_light_mass_CT local_light_kinetic_CT heavy_one_point_CT"
    )
    symbolic_counterterms = (
        ct_mass * phi**2 / 2
        + ct_kinetic * phi**2
        + ct_heavy_source * H
        + G1 * H * phi**2 / 2
        + L1 * phi**4 / 24
    )
    J, K, c3, c4 = sp.symbols(
        "diagnostic_source positive_kernel diagnostic_cubic diagnostic_quartic"
    )
    inverse_source = (
        J / K - c3 * J**2 / (2 * K**3) + (c3**2 / (2 * K**5) - c4 / (6 * K**4)) * J**3
    )
    source_equation = (
        K * inverse_source + c3 * inverse_source**2 / 2 + c4 * inverse_source**3 / 6
    )
    connected_fourth = sp.diff(inverse_source, J, 3).subs(J, 0)
    inherited_inverse = pole_subtraction.data()["two_loop_inverse_propagator"]
    inherited_v = next(
        s for s in inherited_inverse.free_symbols if s.name == "light_shell_distance"
    )
    checks = {
        "bare_amputation_first_order": sp.expand(
            bare_amputated.coeff(hbar, 1) - (A1 - 2 * z1 * A0)
        ),
        "bare_amputation_second_order": sp.expand(
            bare_amputated.coeff(hbar, 2)
            - (A2 - 2 * z1 * A1 + (3 * z1**2 - 2 * z2) * A0)
        ),
        "bare_LSZ_first_order": sp.expand(bare_lsz.coeff(hbar, 1) - 2 * z1),
        "bare_LSZ_second_order": sp.expand(bare_lsz.coeff(hbar, 2) - z1**2 - 2 * z2),
        "complete_bare_to_canonical_amplitude": sp.expand(assembled - ren),
        "literal_quartic_canonical_field_relation": truncate(
            Lbare * Z**2 - Ltotal, hbar
        ),
        "literal_cubic_canonical_field_relation": truncate(Gbare * Z - Gtotal, hbar),
        "derived_squared_coupling_field_relation": truncate(
            gbare * Z**2 - gtotal, hbar
        ),
        "derived_squared_coupling_second_order_keeps_product": sp.expand(
            gtotal.coeff(hbar, 2) - 2 * G * G2 - G1**2
        ),
        "literal_even_polynomial_vacuum_action": sp.expand(V.subs(phi, -phi) - V),
        "all_selected_local_counterterms_even_in_light_field": sp.expand(
            symbolic_counterterms.subs(phi, -phi) - symbolic_counterterms
        ),
        "diagnostic_inverse_Legendre_source_equation": sp.series(
            source_equation - J, J, 0, 4
        )
        .removeO()
        .expand(),
        "diagnostic_amputated_connected_four_point_with_cubic_exchange": sp.factor(
            K**4 * connected_fourth - (3 * c3**2 / K - c4)
        ),
        "diagnostic_even_vacuum_has_no_three_point_exchange": sp.factor(
            (K**4 * connected_fourth).subs(c3, 0) + c4
        ),
        "literal_inherited_light_pole_residue": sp.limit(
            inherited_v / inherited_inverse, inherited_v, 0
        )
        - 1,
        "actual_inherited_unit_disc_residue": pole.point(1)[
            "light_pole_residue_through_two_loops"
        ]
        - 1,
    }
    return {
        "formal_loop_parameter": hbar,
        "bare_to_canonical_light_field_factor": Z,
        "canonical_amputated_amplitude_through_two_loops": ren,
        "bare_amputated_amplitude_through_two_loops": bare_amputated,
        "bare_four_external_leg_LSZ_factor": bare_lsz,
        "assembled_bare_amplitude_through_two_loops": assembled,
        "bare_quartic_parameter": Lbare,
        "bare_cubic_parameter": Gbare,
        "derived_bare_squared_cubic_parameter": gbare,
        "total_canonical_squared_cubic_parameter": gtotal,
        "diagnostic_nonzero_extra_LSZ_error": truncate(bare_lsz * ren - ren, hbar),
        "literal_even_vacuum_potential": V,
        "diagnostic_connected_four_point_Legendre_identity": connected_fourth,
        "actual_inherited_light_pole": pole.point(1),
        "scope": "The symmetric perturbative Phi vacuum has no odd light vertices. After exact Gaussian H elimination, connected four-point light-reducible contributions are external two-point chains or vanish through the odd three-point kernel. The fixed OS canonical light residue is one through two loops. Bare field Z is not that canonical pole residue; its bare LSZ and amputation factors cancel. The heavy field retains its unchanged canonical reference, with no new heavy on-shell or kinetic counterterm condition. The Legendre polynomial check is a diagnostic of the functional derivative coefficients, not a substitute for the labelled momentum-space graph inventory. No two-loop derivative-coordinate/source map is claimed.",
        "checks": checks,
    }
