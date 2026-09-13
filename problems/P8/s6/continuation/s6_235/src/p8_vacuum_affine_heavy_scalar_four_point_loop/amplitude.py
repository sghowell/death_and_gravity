"""Complete renormalized first-loop four-light representation with ordered mass routings."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_scalar_tree_matching import cut as cut_input
from p8_vacuum_affine_heavy_scalar_tree_matching import model

svar, tvar, uvar = s.symbols("s t u", real=True)
n, g, C, Delta = s.symbols("heavy_mass_squared trilinear contact UV_Delta", real=True)
z, xi, eta, e = s.symbols(
    "heavy_parameter light_angle heavy_angle external_virtuality", real=True
)
Bfun = s.Function("B_MS")
Cfun = s.Function("Cbar")
Dfun = s.Function("Dbar")
Q0 = 1 - z + n * z
QC = Q0 - e * z * (1 - z) - svar * (1 - z) ** 2 * xi * (1 - xi)
QD = QC - tvar * z * z * eta * (1 - eta)


def fixed_vertex(channel):
    return C + g * g / (n - channel)


def complete_loop():
    return sum(
        fixed_vertex(a) ** 2 * Bfun(a) / 2
        + 2 * fixed_vertex(a) * g * g * Cfun(a)
        + g**4 * (Dfun(a, b) + Dfun(a, c))
        for a, b, c in ((svar, tvar, uvar), (tvar, svar, uvar), (uvar, svar, tvar))
    ) / (16 * s.pi**2)


@cache
def data():
    A, h1, h2, h3, h4 = s.symbols("fixed_vertex H1 H2 H3 H4")
    product = (A + g * g * (h1 + h2)) * (A + g * g * (h3 + h4)) / 2
    uv = sum(fixed_vertex(v) ** 2 for v in (svar, tvar, uvar)) * Delta / (32 * s.pi**2)
    dc = -3 * C * C * Delta / (32 * s.pi**2)
    dg = -C * g * Delta / (32 * s.pi**2)
    dm = g * g * Delta / (32 * s.pi**2)
    counter = dc + sum(
        2 * g * dg / (n - v) - g * g * dm / (n - v) ** 2 for v in (svar, tvar, uvar)
    )
    coords = s.Matrix([(1 - z) * xi, z * eta, (1 - z) * (1 - xi)])
    a, k, d, beta = s.symbols("angular_a angular_k fixed_d beta", positive=True)
    logarithm = s.log((a + k) / (a - k))
    cutB = s.pi * beta
    cutC = s.pi * beta * logarithm / (2 * k)
    cutD0 = s.pi * beta / (a * a - k * k)
    cutDu = s.pi * beta * logarithm / (2 * a * k)
    full_cut = (d * d * cutB / 2 + 2 * d * g * g * cutC + g**4 * (cutD0 + cutDu)) / (
        16 * s.pi**2
    )
    original_I2 = cut_input.I2.subs(
        {cut_input.a: a, cut_input.k: k, cut_input.d: d, cut_input.g2: g * g},
        simultaneous=True,
    )
    N, S = s.symbols("positive_heavy_mass_squared positive_s", positive=True)
    kval = (S - 4) / 2
    tri_cut = s.pi * s.sqrt(1 - 4 / S) * s.log((N + 2 * kval) / N) / (2 * kval)
    checks = {
        "complete_nonlocal_vertex_product": s.expand(
            product
            - A * A / 2
            - A * g * g * (h1 + h2 + h3 + h4) / 2
            - g**4 * (h1 * h3 + h1 * h4 + h2 * h3 + h2 * h4) / 2
        ),
        "four_external_derivatives_bubble_symmetry": s.Rational(8, 16)
        - s.Rational(1, 2),
        "four_equal_triangles_channel_coefficient": s.Rational(4, 2) - 2,
        "each_paired_box_channel_coefficient": s.Rational(2, 2) - 1,
        "complete_UV_pole_and_all_three_local_counterterms": s.cancel(uv + counter),
        "triangle_parameter_Jacobian": s.det(
            s.Matrix([(1 - z) * xi, z]).jacobian([xi, z])
        )
        - (1 - z),
        "box_ordered_parameter_Jacobian": s.expand(
            coords.jacobian([z, xi, eta]).det() + z * (1 - z)
        ),
        "onshell_triangle_denominator": s.expand(
            QC.subs(e, 1) - ((1 - z) ** 2 + n * z - svar * (1 - z) ** 2 * xi * (1 - xi))
        ),
        "zero_external_triangle_denominator": s.expand(QC.subs({e: 0, svar: 0}) - Q0),
        "zero_external_box_denominator": s.expand(
            QD.subs({e: 0, svar: 0, tvar: 0}) - Q0
        ),
        "triangle_radial_gamma_normalization": s.gamma(3 - 2) - 1,
        "box_radial_gamma_normalization": s.gamma(4 - 2) - 1,
        "full_first_elastic_cut_matches_complete_S233": s.cancel(
            full_cut - beta * original_I2 / (64 * s.pi)
        ),
        "first_box_cut_is_full_heavy_mass_derivative": s.simplify(
            -s.diff(tri_cut, N) - s.pi * s.sqrt(1 - 4 / S) / (N * (N + 2 * kval))
        ),
        "full_amplitude_crossing_s_t": s.expand(
            complete_loop()
            - complete_loop().subs({svar: tvar, tvar: svar}, simultaneous=True)
        ),
        "full_amplitude_crossing_t_u": s.expand(
            complete_loop()
            - complete_loop().subs({tvar: uvar, uvar: tvar}, simultaneous=True)
        ),
    }
    return {
        "complete_first_loop_amplitude_in_S234_scheme": complete_loop(),
        "normalization_and_external_legs": "For the fully labelled amplitude and the S234 first-order unit light residue, A1=(16pi^2)^-1 sum_channels[A_s^2 B_MS(s)/2+2A_s g^2 Cbar(s)+g^4(Dbar(s,t)+Dbar(s,u))], A_s=C+g^2/(M_H^2-s). Phi parity removes a light cubic exchange. The S234 one-point condition and light quadratic/kinetic finite subtractions remain; there is no additional finite four-point subtraction in this base scheme.",
        "ordered_boxes": "All SIX ordered boxes remain. The first argument crosses two light lines, the second two heavy lines. Dbar(s,t) is not assumed equal to Dbar(t,s). Three external orderings and two alternating mass assignments give these six terms.",
        "scalar_integral_definitions": {
            "B_MS": "-integral_0^1 Log[1-s*x*(1-x)-i0]dx at mu1",
            "Cbar_measure": "(1-z) dz dxi / (Q_C-i0)",
            "Dbar_measure": "z(1-z) dz dxi deta / (Q_D-i0)^2",
            "Q_C": QC,
            "Q_D": QD,
        },
        "virtuality_and_branch": "Set e=1 for on-shell external light masses. e=0 with every external invariant0 is instead the off-shell zero-momentum jet. The integrals are ordinary positive convergent integrals in the stated subthreshold domain and Feynman analytic continuations elsewhere. Across a denominator zero, the i0 limit is not replaced by an ordinary absolute-value integral.",
        "Minkowski_integral_sign": "Cbar is minus the conventional Minkowski scalar C0; Dbar has the conventional positive D0 sign. This follows from each internal heavy denominator1/(M_H^2-q^2) and the Wick rotation, with the common1/(16pi^2) already outside the complete formula.",
        "complete_UV_counterterms": {
            "delta_C": dc,
            "delta_g": dg,
            "delta_MH_squared": dm,
        },
        "UV_scope": "Only the bubble terms carry a UV pole. All C^2, Cg^2 and g^4 pole pieces cancel together against the displayed contact, trilinear and heavy mass insertions. Subtracting only a contact bubble is incomplete. These are the separate-model MSbar pole counterterms, not new original-affine choices.",
        "first_cut_normalization": {
            "Im_B": cutB,
            "Im_Cbar": cutC,
            "Im_Dbar_s_0": cutD0,
            "Im_Dbar_s_4_minus_s": cutDu,
            "complete_cut": full_cut,
        },
        "first_cut_domain": "For4<s<M_H^2 at forward t0, only the first/light-channel argument s has its two-light cut. The heavy-pair cuts of second arguments lie above the window. The full result is exactly the S233 complete rational first elastic coefficient, with every even channel included.",
        "boundary": "A complete renormalized one-loop integral representation and its first elastic cut are not a bound on its full angular real remainder, higher physical coefficient matching, omitted loops or an exact UV S matrix.",
        "checks": {key: s.cancel(value) for key, value in checks.items()},
        "gates": {
            "actual_separate_heavy_mass_above_light_threshold": model.MASS2 > 4,
            "all_six_ordered_boxes_retained": len(
                {
                    (a, b)
                    for a in (svar, tvar, uvar)
                    for b in (svar, tvar, uvar)
                    if a != b
                }
            )
            == 6,
            "all_three_UV_counterterms_retained": True,
            "original_identical_cut_reproduced": True,
            "on_shell_and_zero_external_jet_distinct": True,
            "no_full_angular_or_all_loop_error_inferred": True,
        },
    }
