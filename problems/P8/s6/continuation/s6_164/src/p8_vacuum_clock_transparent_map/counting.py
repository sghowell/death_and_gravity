"""Vacuum field grading and clock-jet transparency."""

from functools import cache

import sympy as s

from . import gate


@cache
def data():
    eps, a, b, R0 = s.symbols("variation a b R0", real=True)
    t, r, kappa, Y = s.symbols("t R kappa Y", real=True)
    q = gate.data()["gate"]
    # The exact quotient proves all lower variations vanish without
    # expanding the full degree-thirty variation polynomial.
    quotient = gate.data()["clock_eighth_factor_quotient"]
    factorized = (
        (eps * a + eps**2 * b) ** 8
        * quotient.subs(gate.X, 1 + eps * a + eps**2 * b)
        * R0
    )
    map_difference = s.expand((q.subs(gate.X, t * t * Y / kappa) - 1) * t**3 * r)
    degrees = sorted(s.Poly(map_difference, t).monoms())
    field_degrees = [d[0] for d in degrees]
    checks = {
        "first_new_vacuum_map_field_degree": min(field_degrees) - 19,
        "last_new_vacuum_map_field_degree": max(field_degrees) - 33,
        "new_map_scalar_degrees_are_nineteen_through_thirty_three": sum(
            v != 19 + 2 * j for j, v in enumerate(field_degrees)
        ),
        "vacuum_map_first_changed_coefficient": s.expand(map_difference).coeff(t, 19)
        + 6435 * Y**8 * r / kappa**8,
        "clock_map_eighth_variation_coefficient": s.limit(factorized / eps**8, eps, 0)
        - 6435 * a**8 * R0,
        "clock_map_has_exact_eighth_variation_factor": s.expand(
            q - (gate.X - 1) ** 8 * quotient
        ),
        "gate_unchanged_when_clock_X_is_identically_one": q.subs(gate.X, 1),
        "new_field_support_radius": max(field_degrees) - 33,
        "new_quadratic_source_support_radius": 2 * max(field_degrees) - 66,
        "full_source_spectral_square_radius": 66**2 - 4356,
    }
    return {
        "vacuum_map_difference_scalar_field_degrees": field_degrees,
        "clock_variations_identical_through_order": 7,
        "first_possible_nonidentity_clock_variation_order": 8,
        "new_mapped_field_Fourier_radius": 33,
        "new_mapped_source_Fourier_radius": 66,
        "clock_identity": "On every canonical unit clock Psi=sqrt(kappa)t with physical X=1, q8(X) is identically zero and Fhat=Psi, including its derivatives along the trajectory. This holds independently of the flat/curved coordinate expression of the old cubic R, provided those jets are defined. The gate's first seven field/metric variations vanish there.",
        "vacuum_source_identity": "Fhat-Fold starts at scalar field degree nineteen. Hence all action/source terms that can enter the already named two-loop vacuum observables coincide. Use the complete transformed action, counterterms and source J Fhat; do not discard the generated higher terms from its definition.",
        "regulator_prescription": "In the same fixed light-mass units use X_D=mu^(2epsilon)(partial Psi_D)^2/(kappa m_Phi^4), with m_Phi=1 in the numerical formulas and dimensionless kappa. The q8 monomial X_D^j supplies mu^(2j epsilon), together with the same regulated cubic R_D and restored light-mass factors of S6.160. Apply this common lift to action, counterterms, Jacobian and physical source before taking finite parts.",
        "not_a_background_solution": "Clock transparency is a map identity. It does not make the polynomial SAT8 parent support the target bounce or prove a global inverse, common-parent scalar matching or a rolling quantum remainder.",
        "checks": checks,
    }
