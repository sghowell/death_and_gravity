"""Independent bare-coupling and local-counterterm scheme conversion."""

from functools import cache

import sympy as sp
from p8_vacuum_gauge_yukawa_screen import flow


@cache
def data():
    L, G, M, Y = sp.symbols(
        "MS_quartic MS_cubic MS_heavy_mass_squared MS_Yukawa_squared", positive=True
    )
    h, Ibar, ell, sigma, r, fp = sp.symbols(
        "loop_order MS_dimensional_reference scale_logarithm fixed_scalar_finite_contact scalar_residue_slope fermion_residue_slope",
        real=True,
    )
    N, Q = sp.Integer(6), 16 * sp.pi**2
    g = G * G
    zF = -2 * N * Y * Ibar / Q
    Zms = 1 + h * zF
    Zstar = 1 + h * (zF - r + fp)
    dG = -L * G * ell / (2 * Q) + (fp - r) * G
    dM = -g * ell / (2 * Q)
    dL = -3 * L * L * ell / (2 * Q) - sigma + 2 * (fp - r) * L
    ms_ctG = L * G * Ibar / (2 * Q)
    ms_ctM = g * Ibar / (2 * Q)
    fermion_ctL = -24 * N * Y * Y * Ibar / Q
    ms_ctL = 3 * L * L * Ibar / (2 * Q) + fermion_ctL
    star_ctG = L * G * (Ibar + ell) / (2 * Q)
    star_ctM = g * (Ibar + ell) / (2 * Q)
    star_ctL = 3 * L * L * (Ibar + ell) / (2 * Q) + sigma + fermion_ctL
    bareGms = (G + h * ms_ctG) / Zms
    bareGstar = (G + h * dG + h * star_ctG) / Zstar
    bareLms = (L + h * ms_ctL) / Zms**2
    bareLstar = (L + h * dL + h * star_ctL) / Zstar**2
    bareMms = M + h * ms_ctM
    bareMstar = M + h * dM + h * star_ctM
    hs = sp.symbols("heavy_exchange_h_s heavy_exchange_h_t heavy_exchange_h_u")
    A0 = -L + g * sum(hs)
    square = sum((-L + g * v) ** 2 for v in hs)
    tree_shift = -dL + 2 * G * dG * sum(hs) - g * dM * sum(v * v for v in hs)
    expected = ell * square / (2 * Q) + sigma + 2 * (fp - r) * A0
    ms_amp_ct = -ms_ctL + 2 * G * ms_ctG * sum(hs) - g * ms_ctM * sum(v * v for v in hs)
    star_amp_ct = (
        -star_ctL + 2 * G * star_ctG * sum(hs) - g * star_ctM * sum(v * v for v in hs)
    )
    bareLpole = sp.factor(ms_ctL - 2 * zF * L)
    yflow, Lflow = sp.symbols("positive_Yukawa quartic_L", positive=True)
    inherited_beta = flow.data()["beta_quartic_times_loop_denominator"].subs(
        {yflow: sp.sqrt(Y), Lflow: L}
    )
    first = lambda value: sp.factor(sp.diff(value, h).subs(h, 0))
    checks = {
        "bare_cubic_equal_through_one_loop": first(bareGstar - bareGms),
        "bare_quartic_equal_through_one_loop": first(bareLstar - bareLms),
        "bare_heavy_mass_equal_through_one_loop": first(bareMstar - bareMms),
        "MS_bare_quartic_pole_and_wavefunction_reproduce_inherited_flow": sp.factor(
            2 * Q * bareLpole / Ibar - inherited_beta
        ),
        "MS_bare_cubic_pole_includes_wavefunction": sp.factor(
            ms_ctG - zF * G - (L * G / 2 + 2 * N * Y * G) * Ibar / Q
        ),
        "canonical_field_ratio_first_order": first(Zms / (1 + h * (r - fp)) - Zstar),
        "tree_variation_matches_complete_reference_conversion": sp.factor(
            tree_shift - expected
        ),
        "independent_total_counterterm_difference": sp.factor(
            ms_amp_ct - star_amp_ct - ell * square / (2 * Q) - sigma
        ),
        "counterterm_and_field_conversion_agree": sp.factor(
            ms_amp_ct - star_amp_ct + 2 * (fp - r) * A0 - tree_shift
        ),
        "cubic_square_conversion_uses_two_G": first((G + h * dG) ** 2 - g) - 2 * G * dG,
        "whole_dimensional_reference_cancels_in_scheme_change": sp.diff(expected, Ibar),
        "finite_contact_not_a_second_derivative_fit": sp.diff(sigma, hs[0]),
    }
    return {
        "parameters": {"L": L, "G": G, "M": M, "Y": Y},
        "reference_symbols": {
            "h": h,
            "Ibar": Ibar,
            "ell": ell,
            "sigma": sigma,
            "r": r,
            "fp": fp,
        },
        "minimal_scalar_field_UV_increment": zF,
        "MS_field_factor": Zms,
        "canonical_star_field_factor": Zstar,
        "MS_total_vertex_counterterms": {"L": ms_ctL, "G": ms_ctG, "M": ms_ctM},
        "canonical_star_total_vertex_counterterms": {
            "L": star_ctL,
            "G": star_ctG,
            "M": star_ctM,
        },
        "canonical_star_parameter_increments": {"L": dL, "G": dG, "M": dM},
        "MS_bare_quartic_pole": bareLpole,
        "bare_coupling_identities": {
            "MS_G": bareGms,
            "star_G": bareGstar,
            "MS_L": bareLms,
            "star_L": bareLstar,
            "MS_M": bareMms,
            "star_M": bareMstar,
        },
        "exact_first_order_tree_conversion": expected,
        "heavy_exchange_variables": hs,
        "scope": "Interaction parameters are MS at scale mF; scalar pole mass, H one-point function and vacuum energy have explicit physical local references. The star parameters are those used in the canonical old-scalar-subtraction plus new-fermion-box representation. Its unit-residue amplitude gets no second LSZ factor. The entire finite logarithm and the old fixed contact are retained, not fitted.",
        "checks": checks,
    }
