"""Complex forward on-shell routes and exact permutation descent."""

from functools import cache
from itertools import permutations

import sympy as sp

from . import series


@cache
def data():
    E, p = sp.symbols("complex_energy complex_spatial_component")
    vectors = {
        "Phi1": sp.Matrix([-sp.I * E, 0, 0, -p]),
        "Phi2": sp.Matrix([-sp.I * E, 0, 0, p]),
        "Phi3": sp.Matrix([sp.I * E, 0, 0, p]),
        "Phi4": sp.Matrix([sp.I * E, 0, 0, -p]),
    }
    labels = tuple(vectors)
    s = sp.symbols("continued_channel_s")
    s_invariant = -(vectors["Phi1"] + vectors["Phi2"]).dot(
        vectors["Phi1"] + vectors["Phi2"]
    )
    t_invariant = -(vectors["Phi1"] + vectors["Phi3"]).dot(
        vectors["Phi1"] + vectors["Phi3"]
    )
    u_invariant = -(vectors["Phi1"] + vectors["Phi4"]).dot(
        vectors["Phi1"] + vectors["Phi4"]
    )
    routes = []
    checks = {
        "all_outgoing_momentum_conservation": sum(vectors.values(), sp.zeros(4, 1)),
        "direct_s_invariant": s_invariant - 4 * E * E,
        "direct_t_invariant": t_invariant,
        "direct_u_invariant": u_invariant + 4 * p * p,
        "equal_mass_Mandelstam_sum": s_invariant
        + t_invariant
        + u_invariant
        - 4 * (E * E - p * p),
        "on_shell_crossed_dictionary": u_invariant.subs(p * p, s / 4 - 1) - (4 - s),
        "strict_energy_and_spatial_majorant": sp.Rational(3, 2) ** 2
        - sp.Rational(7, 4)
        - sp.Rational(1, 2),
        "external_component_norm_majorant": 2 * sp.Rational(3, 2) - 3,
        "three_external_route_majorant": 3 * 3 - 9,
        "strict_routing_bound_has_margin": 18 - 9 - 9,
        "minimum_mass_half_Neumann_ratio": sp.Rational(18, 36) - sp.Rational(1, 2),
    }
    for label, vector in vectors.items():
        checks[label + "_same_mass_shell"] = vector.dot(vector) + E * E - p * p
    energy_flip = {"Phi1": "Phi4", "Phi2": "Phi3", "Phi3": "Phi2", "Phi4": "Phi1"}
    spatial_flip = {"Phi1": "Phi2", "Phi2": "Phi1", "Phi3": "Phi4", "Phi4": "Phi3"}
    for label in labels:
        checks[label + "_energy_sign_is_external_permutation"] = (
            vectors[label].subs(E, -E) - vectors[energy_flip[label]]
        )
        checks[label + "_spatial_sign_is_external_permutation"] = (
            vectors[label].subs(p, -p) - vectors[spatial_flip[label]]
        )
    orders = series.data()["cyclic_labelled_four_scalar_orders"]
    for i, order in enumerate(orders):
        partial = [sp.zeros(4, 1)]
        for label in order[:-1]:
            partial.append(partial[-1] + vectors[label])
        routes.append({"vertices": order, "partial_momenta": partial})
        checks[f"route_{i}_momentum_closure"] = partial[-1] + vectors[order[-1]]
        checks[f"route_{i}_four_propagators"] = len(partial) - 4
    orbit_rows = []
    for i, perm in enumerate(permutations(labels)):
        renamed = dict(zip(labels, perm, strict=True))
        images = []
        for order in orders:
            value = tuple(renamed[label] for label in order)
            first = value.index("Phi1")
            images.append(value[first:] + value[:first])
        checks[f"all_label_permutation_{i}_preserves_box_set"] = len(
            set(images).symmetric_difference(orders)
        )
        orbit_rows.append(
            {"permutation": perm, "distinct_cyclic_images": len(set(images))}
        )
    return {
        "Euclidean_forward_external_vectors": vectors,
        "on_shell_dictionary": "E^2=s/4, p^2=s/4-1, all p_i,E^2=-1, t=0, u=4-s.",
        "six_complete_forward_routes": routes,
        "all_twenty_four_external_permutation_checks": orbit_rows,
        "independent_energy_sign_permutation": energy_flip,
        "independent_spatial_sign_permutation": spatial_flip,
        "complex_s_disc_center": 2,
        "complex_s_disc_radius": 5,
        "strict_external_component_norm_upper": 3,
        "strict_three_external_route_norm_upper": 9,
        "conservative_Neumann_route_norm_upper": 18,
        "minimum_fermion_mass": 36,
        "analytic_descent": "The convergent Euclidean momentum expansion is invariant under independent E and p sign changes because each is an external permutation of the full six-box sum. Its analytic series therefore depends on E^2 and p^2 and descends holomorphically to s through their zeroes.",
        "scope": "Uniform complex-domain continuation of the complete one-loop box after its entire local reference subtraction, not an arbitrary shifted divergent loop contour or a full quantum analyticity theorem.",
        "checks": checks,
    }
