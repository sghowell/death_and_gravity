"""Whole unchanged vacuum jets and all leading two-body species below H threshold."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_graviton_cut_isolated_replay import source as previous
from p8_vacuum_affine_proca_gaussian import bridge as proca_parent

current = previous.current
S, MU, K = previous.S, previous.MU, previous.K
Z = s.Symbol("physical_scattering_cosine", real=True)
VECTOR_MASS = proca_parent.MASS
VECTOR_MASS2 = VECTOR_MASS**2
HEAVY_MASS2 = current.heavy.MASS2
CUBIC = current.heavy.G
CONTACT = current.heavy.CONTACT


@cache
def data():
    inherited = previous.data()
    full = current.fixed_functions()
    u, X = current.u, current.X
    zero = {u: 0, X: 0}
    F, R = full["F_full"], full["R_full"]
    j = current.heavy.coefficients()["normalized_heavy_source"]
    checks = {key: value for key, value in inherited["checks"].items()}
    checks.update(
        {
            "entire_F_XX_no_vacuum_four_derivative_contact": s.cancel(
                s.diff(F, X, 2).subs(zero)
            ),
            "entire_R_XX_no_vacuum_Ia_quartic": s.cancel(s.diff(R, X, 2).subs(zero)),
            "entire_F_uuX_no_Phi_squared_Y": s.cancel(s.diff(F, u, 2, X).subs(zero)),
            "entire_F_uuu_no_light_cubic": s.cancel(s.diff(F, u, 3).subs(zero)),
            "entire_F_uX_no_derivative_cubic": s.cancel(s.diff(F, u, X).subs(zero)),
            "entire_F_uuuu_retained_contact": s.cancel(
                s.diff(F, u, 4).subs(zero) / previous.KAPPA - CONTACT
            ),
            "entire_J_uu_retained_heavy_cubic": s.cancel(
                s.diff(j, u, 2).subs(zero) / s.sqrt(previous.KAPPA) - CUBIC
            ),
            "literal_retained_Proca_mass": VECTOR_MASS - 1000,
            "literal_retained_Proca_zeta": 1 / proca_parent.ZETA - VECTOR_MASS2,
            "actual_heavy_not_rounded_mass": HEAVY_MASS2
            - (s.Integer(10) ** 200 / 512 + 2),
        }
    )
    a, au, ax, H, box, zz = s.symbols("a a_u a_X Href Box_u Z_u", real=True)
    rg = 1 + X**2 * a
    rx = 2 * X * a + X**2 * ax
    ru = X**2 * au
    literal = (rg - 1) * (
        -3 * H + 3 * ru / (4 * rg) + box / X + (-1 / X**2 + 3 * rx / (2 * rg * X)) * zz
    )
    regular = (
        -3 * H * X**2 * a
        + 3 * X**4 * a * au / (4 * rg)
        + X * a * box
        + (-a + 3 * X * a * rx / (2 * rg)) * zz
    )
    checks["entire_regular_retained_Proca_source"] = s.cancel(literal - regular)
    f = s.Symbol("light_field_degree")
    marked = f * regular.subs(
        {X: f * f * X, box: f * box, zz: f**3 * zz}, simultaneous=True
    )
    numerator, denominator = s.fraction(s.cancel(marked))
    for degree in range(4):
        checks["Proca_source_has_no_degree_" + str(degree)] = s.expand(numerator).coeff(
            f, degree
        )
    return {
        "entire_original_R_F": (R, F),
        "entire_retained_heavy_source": j,
        "entire_formal_loop_marker": inherited[
            "entire_formal_loop_marked_scalar_coefficient"
        ],
        "all_three_vacuum_constants": inherited["all_three_retained_vacuum_constants"],
        "parameters": {
            "kappa": previous.KAPPA,
            "external_mass_squared": previous.MASS2,
            "vector_mass_squared": VECTOR_MASS2,
            "heavy_mass_squared": HEAVY_MASS2,
            "cubic": CUBIC,
            "contact": CONTACT,
        },
        "entire_regular_Proca_source_one_form_scalar": regular,
        "leading_four_point_action": "Canonical massive Phi, massless original M1, three-polarization mass1000 Proca, TT Einstein gravity, and heavy H with H Phi^2 cubic and unchanged Phi^4 contact. The full four-Phi higher derivative jets cancel; no deleted bare operator is reinstated.",
        "leg_count": "The retained Proca source one-form begins at least four Phi legs; its square begins eight. There is no Phi Phi A, Phi Phi A A contact beyond minimal gravity, H A A, H h h, or light cubic. R-1 begins at least four fields; the regular Ia operators cannot create a lower-arity vertex. Minimal physical-metric vector mass and kinetic terms are independent of Phi.",
        "complete_first_loop_species_below_heavy_threshold": "For 4m^2<s<M_H^2, the physical two-body cuts are Phi Phi, original M1 M1, TT graviton graviton, and Proca Proca when s>=4M_A^2. H plus a massless state first opens at M_H^2. Mixed light-species production is absent at this tree order by the listed vertices. Three-body cuts begin at a higher formal loop order.",
        "dimensional_continuation": "Only the actual leading vacuum vertices are continued minimally to D=4+2epsilon for the elastic soft calculation. This is an explicitly specified regulator convention, not a unique full off-background DHOST continuation or new physical counterterm.",
        "checks": checks,
        "gates": {
            "entire_parent_fixed_functions_retained": R
            == current.heavy.coefficients()["R"],
            "no_vacuum_retuning_or_state_repreparation": True,
            "same_original_positive_vector_mass_and_zeta": VECTOR_MASS == 1000
            and proca_parent.ZETA == s.Rational(1, 10**6),
            "source_degree_denominator_regular": denominator.subs(f, 0) != 0,
            "formal_named_window_below_actual_heavy_mass": bool(
                4 * VECTOR_MASS2 < s.Integer(10) ** 196 < HEAVY_MASS2
            ),
            "formal_loop_scope_not_exact_quantum_vacuum": True,
        },
    }
