"""Full scalar, Yukawa and counterterm substitution with a fixed map."""

from functools import cache

import sympy as s
from p8_vacuum_quantum_map import interactions as first


@cache
def data():
    x, y, R, S, t, h = s.symbols("site_x site_y cubic_R_x cubic_R_y field_degree h")
    old = first.data()
    Q = old["two_site_free_action"]
    V = old["two_site_full_nonlocal_quartic_action"]
    F = {x: t * x + t**3 * R, y: t * y + t**3 * S}
    # Match the frozen parent's symbol assumptions and names exactly.
    oldvars = {str(z): z for z in Q.free_symbols | V.free_symbols}
    xx, yy = oldvars["site_x"], oldvars["site_y"]
    F = {xx: t * x + t**3 * R, yy: t * y + t**3 * S}
    scalar = s.expand((Q + V).subs(F, simultaneous=True))
    a1, a2, b1, b2, c1, c2, Y, m = s.symbols("a1 a2 b1 b2 c1 c2 Yukawa mass")
    eta, zeta = s.symbols("fermion_bilinear_x fermion_bilinear_y")
    fermion = m * (eta + zeta) + Y * (xx * eta + yy * zeta)
    fermion_mapped = s.expand(fermion.subs(F, simultaneous=True))
    ct = (
        (h * a1 + h * h * a2) * (xx**2 + yy**2) / 2
        + (h * b1 + h * h * b2) * (xx**4 + yy**4) / 24
        + (h * c1 + h * h * c2) * (xx * eta + yy * zeta)
    )
    ct_mapped = s.expand(ct.subs(F, simultaneous=True))
    G, M, H, J = s.symbols("G M H J")
    heavy = M * H * H / 2 + H * (J + G * (x + R) ** 2 / 2)
    eliminated = s.expand(heavy.subs(H, -(J + G * (x + R) ** 2 / 2) / M))
    parent_eliminated = -((J + G * x * x / 2) ** 2) / (2 * M)
    return {
        "all_scalar_generated_degrees": {
            n: scalar.coeff(t, n) for n in (2, 4, 6, 8, 10, 12)
        },
        "full_generated_Yukawa_action": fermion_mapped,
        "full_transformed_parent_counterterms": ct_mapped,
        "fixed_map_prescription": "The literal S6.111 R coefficients are fixed numerical values at the named canonical interaction boundary. F has no extra h dependence. Every parent coefficient/counterterm is transformed with that same F; a different coupling-dependent F prescription would require its own coordinate variations.",
        "checks": {
            "all_scalar_degrees_retained": s.expand(
                scalar - sum(t**n * scalar.coeff(t, n) for n in (2, 4, 6, 8, 10, 12))
            ),
            "ordinary_Yukawa_vertex": fermion_mapped.coeff(t, 1)
            - Y * (x * eta + y * zeta),
            "generated_three_Phi_Yukawa_vertex": fermion_mapped.coeff(t, 3)
            - Y * (R * eta + S * zeta),
            "fermion_mass_unchanged": fermion_mapped.coeff(t, 0) - m * (eta + zeta),
            "first_counterterm_sextic": ct_mapped.coeff(h, 1).coeff(t, 6)
            - a1 * (R * R + S * S) / 2
            - b1 * (x**3 * R + y**3 * S) / 6,
            "second_quadratic_counterterm_quartic": ct_mapped.coeff(h, 2).coeff(t, 4)
            - a2 * (x * R + y * S)
            - b2 * (x**4 + y**4) / 24,
            "first_generated_Yukawa_counterterm": ct_mapped.coeff(h, 1).coeff(t, 3)
            - c1 * (R * eta + S * zeta),
            "heavy_elimination_commutes_with_fixed_F": s.expand(
                eliminated - parent_eliminated.subs(x, x + R)
            ),
            "source_H_and_vacuum_terms_kept_under_F": s.diff(eliminated, J, 2) + 1 / M,
        },
        "scope": "Bilinear symbols are even fermion spurions for an exact substitution identity, not commuting replacements inside fermion loops. Heavy inverses are full kernels. No generated higher derivative is resummed into a new free propagator.",
    }
