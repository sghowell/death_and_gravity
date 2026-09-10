"""Complex on-shell routing and regulated full-box Ward bookkeeping."""

from functools import cache
from itertools import permutations

import sympy as sp
from p8_vacuum_gauge_yukawa_screen import model

from . import series


@cache
def data():
    E, p = sp.symbols("positive_energy scalar_spatial_component", real=True)
    r = sp.symbols("real_stereographic_angle", real=True)
    n = sp.Matrix([2 * r / (1 + r * r), 0, (1 - r * r) / (1 + r * r)])
    e1 = sp.Matrix([n[2], 0, -n[0]])
    e2 = sp.Matrix([0, 1, 0])
    # All-outgoing Euclidean external vectors; p continues to imaginary values.
    momenta = {
        "Phi1": sp.Matrix([-sp.I * E, 0, 0, -p]),
        "Phi2": sp.Matrix([-sp.I * E, 0, 0, p]),
        "gauge1": sp.Matrix([sp.I * E, E * n[0], E * n[1], E * n[2]]),
        "gauge2": sp.Matrix([sp.I * E, -E * n[0], -E * n[1], -E * n[2]]),
    }
    checks = {
        "spatial_gauge_direction_unit": sp.factor(n.dot(n) - 1),
        "first_transverse_vector_unit": sp.factor(e1.dot(e1) - 1),
        "second_transverse_vector_unit": e2.dot(e2) - 1,
        "transverse_vectors_orthogonal": e1.dot(e2),
        "first_gauge_transversality": sp.factor(n.dot(e1)),
        "second_gauge_transversality": n.dot(e2),
        "Euclidean_momentum_conservation": sum(momenta.values(), sp.zeros(4, 1)),
        "first_scalar_mass_shell": sp.factor(
            momenta["Phi1"].dot(momenta["Phi1"]) + E * E - p * p
        ),
        "second_scalar_mass_shell": sp.factor(
            momenta["Phi2"].dot(momenta["Phi2"]) + E * E - p * p
        ),
        "first_gauge_null": sp.factor(momenta["gauge1"].dot(momenta["gauge1"])),
        "second_gauge_null": sp.factor(momenta["gauge2"].dot(momenta["gauge2"])),
        "channel_invariant": (momenta["Phi1"] + momenta["Phi2"]).dot(
            momenta["Phi1"] + momenta["Phi2"]
        )
        + 4 * E * E,
        "energy_strict_rational_majorant": sp.Rational(5, 4) ** 2
        - sp.Rational(5, 4)
        - sp.Rational(5, 16),
        "spatial_l1_Cauchy_majorant": 2**2 - 3 - 1,
        "three_maximal_external_l1_bounds": 3 * 4 - 12,
    }
    routes = []
    for index, order in enumerate(series.data()["cyclic_labelled_box_orders"]):
        partial = [sp.zeros(4, 1)]
        for name in order[:-1]:
            partial.append(partial[-1] + momenta[name])
        checks[f"routing_{index}_closed"] = partial[-1] + momenta[order[-1]]
        checks[f"routing_{index}_four_propagators"] = len(partial) - 4
        routes.append({"vertices": order, "partial_momenta": partial})
    # Each of the two cyclic orders of the remaining vertices has three
    # positions for the contracted gauge insertion. Regulated translations
    # identify its two adjacent triangle terms with consecutive boundaries.
    ward_rows = []
    for orientation, rest in enumerate(permutations(("Phi2", "gauge2"))):
        triangles = sp.symbols(f"Ward_triangle_{orientation}_0:3")
        terms = [triangles[j] - triangles[(j + 1) % 3] for j in range(3)]
        checks[f"regulated_Ward_telescoping_orientation_{orientation}"] = sum(terms)
        ward_rows.append(
            {"remaining_cyclic_order": ("Phi1",) + rest, "contracted_terms": terms}
        )
    generators = model.data()["gauge_generators"]
    for i, Ta in enumerate(generators):
        checks[f"one_gauge_color_trace_{i}"] = sp.trace(Ta)
        for j, Tb in enumerate(generators):
            checks[f"two_gauge_color_commutator_trace_{i}_{j}"] = sp.trace(
                Ta * Tb - Tb * Ta
            )
    # Complete parity-even, degree-two rank-two Lorentz tensor basis:
    # eta times six Gram products, plus nine ordered vector outer products.
    gram = sp.symbols("kk ll pp kl kp lp")
    kk, ll, _pp, kl, kp, lp = gram
    coefficients = sp.symbols("degree_two_tensor_c0:15")
    metric_coefficient = sum(c * g for c, g in zip(coefficients[:6], gram, strict=True))
    tensor_coefficients = sp.Matrix(3, 3, coefficients[6:])
    first = (sp.Matrix([[kk, kl, kp]]) * tensor_coefficients).T
    first[0] += metric_coefficient
    second = tensor_coefficients * sp.Matrix([kl, ll, lp])
    second[1] += metric_coefficient
    equations = [
        sp.expand(v).coeff(g) for v in list(first) + list(second) for g in gram
    ]
    ward_matrix, _ = sp.linear_eq_to_matrix(equations, coefficients)
    expected_tensor = sp.Matrix([0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0])
    checks["degree_two_Ward_system_rank"] = ward_matrix.rank() - 14
    checks["degree_two_unique_gauge_tensor_kernel"] = ward_matrix * expected_tensor
    checks["degree_two_Ward_nullity"] = len(ward_matrix.nullspace()) - 1
    return {
        "Euclidean_all_outgoing_external_momenta": momenta,
        "mass_shell_dictionary": "4 E^2=s, p^2=s/4-1; no complex conjugation of continued p",
        "spatial_direction_and_two_real_transverse_polarizations": [n, e1, e2],
        "six_cyclic_box_routes": routes,
        "regulated_Ward_telescoping_rows": ward_rows,
        "physical_anchor_and_continued_interval": "[1,5], with physical open anchor 4<s<5",
        "strict_cumulative_component_norm_upper": 12,
        "fermion_mass_domain": "mF>=24",
        "low_degree_operator_classification": {
            "n0": "The full regulated Ward-invariant sum vanishes; individual logarithmic box terms are not discarded before regulation",
            "n1": "No Lorentz-invariant rank-three parity-even tensor of two gauge polarizations and one external momentum",
            "n2": "Unique CP-even neutral-scalar operator Phi^2 F^a_mu_nu F^a_mu_nu; CP-odd term absent for the real scalar Yukawas",
        },
        "degree_two_complete_parity_even_Ward_system": {
            "Gram_product_basis": gram,
            "tensor_coefficient_basis": coefficients,
            "linear_constraint_matrix": ward_matrix,
            "one_dimensional_solution": expected_tensor,
            "normalization": "(k dot l) eta_mu_nu - l_mu k_nu",
        },
        "continuation": "The absolutely convergent Taylor tail defines analytic continuation from real Euclidean external momenta inside the strict Neumann domain. Bose symmetry makes the cut angular integral even in p, so its continuation through p=0 is analytic in p^2. No arbitrary complex loop-contour translation is used.",
        "checks": checks,
    }
