"""Both active flavor masses and all fourteen vacuum-energy constants."""

from functools import cache

import sympy as sp


@cache
def data():
    phi = sp.symbols("reference_scalar_Phi", real=True)
    m, y, nu = sp.symbols(
        "positive_fermion_mass positive_Yukawa positive_reference_scale", positive=True
    )
    r = sp.symbols("dimensionless_field_ratio", real=True)
    Nc, N, Ntotal = sp.Integer(3), sp.Integer(6), sp.Integer(42)
    q = 16 * sp.pi**2
    f = lambda M: M**4 * (sp.log(M * M / (nu * nu)) - sp.Rational(3, 2))
    active = -Nc * (f(m + y * phi) + f(m - y * phi)) / q
    inactive = -12 * Nc * f(m) / q
    V = active + inactive
    Ibar = sp.symbols("entire_dimensional_reference_Ibar")
    UV = (Nc * ((m + y * phi) ** 4 + (m - y * phi) ** 4) + 12 * Nc * m**4) * Ibar / q
    f0 = 4 * N * y * y * m * m / q
    v4 = -64 * N * y**4 / q
    vacuum = sp.Rational(3, 2) * Ntotal * m**4 / q
    k = sp.symbols("even_field_degree", integer=True, positive=True)
    high_coefficient = 96 * Nc / (k * (k - 1) * (k - 2) * (k - 3) * (k - 4))
    normalized = -Nc * (
        (1 + r) ** 4 * (2 * sp.log(1 + r) - sp.Rational(3, 2))
        + (1 - r) ** 4 * (2 * sp.log(1 - r) - sp.Rational(3, 2))
    )
    checks = {
        "all_fourteen_color_flavor_count": Ntotal - 14 * Nc,
        "active_color_flavor_count": N - 2 * Nc,
        "inert_color_flavor_count": Ntotal - N - 12 * Nc,
        "full_vacuum_energy_constant": sp.simplify(V.subs({phi: 0, nu: m}) - vacuum),
        "all_flavor_UV_vacuum_reference": UV.subs(phi, 0) - Ntotal * m**4 * Ibar / q,
        "UV_mass_reference_matches_two_point": sp.diff(UV, phi, 2).subs(phi, 0)
        - 12 * N * y * y * m * m * Ibar / q,
        "UV_quartic_reference": sp.diff(UV, phi, 4).subs(phi, 0)
        - 24 * N * y**4 * Ibar / q,
        "mass_threshold_from_potential": sp.simplify(
            sp.diff(V, phi, 2).subs({phi: 0, nu: m}) - f0
        ),
        "quartic_threshold_from_potential": sp.simplify(
            sp.diff(V, phi, 4).subs({phi: 0, nu: m}) - v4
        ),
        "odd_linear_threshold_zero": sp.simplify(sp.diff(V, phi).subs(phi, 0)),
        "odd_cubic_threshold_zero": sp.simplify(sp.diff(V, phi, 3).subs(phi, 0)),
        "fifth_derivative_seed": sp.simplify(sp.diff(f(m), m, 5) - 48 / m),
        "positive_sixth_remainder_coefficient": high_coefficient.subs(k, 6)
        - sp.Rational(2, 5),
    }
    rows = []
    for degree in range(0, 15, 2):
        coefficient = sp.simplify(
            sp.diff(normalized, r, degree).subs(r, 0) / sp.factorial(degree)
        )
        expected = {0: 9, 2: 12, 4: -16}.get(
            degree, high_coefficient.subs(k, degree) if degree >= 6 else 0
        )
        checks[f"independent_even_field_coefficient_{degree}"] = coefficient - expected
        rows.append({"degree": degree, "coefficient": coefficient})
    return {
        "reference_field": phi,
        "fermion_mass": m,
        "Yukawa": y,
        "active_one_loop_MSbar_potential": active,
        "inert_one_loop_MSbar_potential": inactive,
        "all_flavor_vacuum_energy_at_scale_m": vacuum,
        "entire_UV_potential_reference": sp.expand(UV),
        "zero_momentum_mass_threshold_at_scale_m": f0,
        "quartic_threshold_at_scale_m": v4,
        "positive_higher_even_field_coefficient": high_coefficient,
        "independent_field_coefficients": rows,
        "field_domain": "|y Phi_reference/mF|<1, with a uniform remainder bound on the closed half-domain",
        "higher_field_remainder": "All even degrees k>=6 have coefficient 96 Nc/[k(k-1)(k-2)(k-3)(k-4)]>0. On |r|<=1/2 their sum is at most (2/5)r^6/(1-r^2)<=8r^6/15.",
        "scope": "Complete one-loop fermion determinant only, at reference scale nu=mF. The inert flavors contribute to the fixed vacuum energy even though their Phi derivatives vanish. Old scalar loops and later gauge/fermion loops are not included in this functional.",
        "checks": checks,
    }
