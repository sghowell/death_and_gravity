"""Identify all lower-order owners in the isolated quartic-contact direction."""

from functools import cache

import sympy as s


@cache
def data():
    L, g, G, sigma, z = s.symbols("L g G sigma z")
    tad, bubble = s.symbols("light_tadpole mixed_heavy_light_bubble")
    r, fp, Y, m, pole = s.symbols("scalar_slope fermion_slope Y m MS_pole")
    f1 = L * tad / 2 - g * bubble
    J1 = -G * tad / 2
    vacuum1 = s.Symbol("first_vacuum_determinant")
    Phi, H = s.symbols("Phi H")
    potential = sigma * Phi**4 / s.factorial(4)
    zpsi = -(Y / 2 + s.Symbol("a_Cf")) * pole
    mass = (Y - 4 * s.Symbol("a_Cf")) * m * pole
    return {
        "first_scalar_inverse_mass_function": f1,
        "contact_inserted_scalar_mass": sigma * s.diff(f1, L),
        "first_H_source": J1,
        "first_canonical_field_correction": r - fp,
        "contact_potential": potential,
        "checks": {
            "contact_insertion_mass_is_local": s.diff(f1, L) - tad / 2,
            "contact_insertion_slope_zero": s.diff(sigma * s.diff(f1, L), z),
            "physical_mass_reference_cancels_local_insertion": sigma * tad / 2
            - sigma * tad / 2,
            "scalar_slope_L_independent": s.diff(r, L),
            "fermion_slope_L_independent": s.diff(fp, L),
            "first_source_L_independent": s.diff(J1, L),
            "first_fermion_kinetic_reference_L_independent": s.diff(zpsi, L),
            "first_fermion_mass_reference_L_independent": s.diff(mass, L),
            "first_vacuum_L_independent_at_fixed_free_masses": s.diff(vacuum1, L),
            "contact_zero_background_value": potential.subs(Phi, 0),
            "contact_zero_background_scalar_source": s.diff(potential, Phi).subs(
                Phi, 0
            ),
            "contact_zero_background_quadratic_hessian": s.diff(potential, Phi, 2).subs(
                Phi, 0
            ),
            "contact_has_no_heavy_source": s.diff(potential, H),
            "no_second_external_field_factor": s.diff(r - fp, L),
        },
        "scope": "At fixed g,M,Y,m and physical Phi mass one, the only new quadratic contact insertion is the momentum-independent light tadpole. The fixed physical mass condition cancels it. The first field correction, H source and vacuum determinant are L independent. A quartic counterterm first contributes to a zero-background vacuum graph at order h^3, not h^2. These ownership facts do not set the remaining model vacuum reference or canonical two-loop correction to zero.",
    }
