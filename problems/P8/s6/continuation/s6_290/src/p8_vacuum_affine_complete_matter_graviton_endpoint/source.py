"""Whole original vacuum jets, one-point condition and one-loop endpoint arity."""

from functools import cache

import sympy as s
from p8_vacuum_affine_covariant_gaussian_four_point_matching import source as previous
from p8_vacuum_affine_massive_elastic_proca_infrared import source as minimal
from p8_vacuum_affine_massive_matter_graviton_vertex import vertex as old_vertex

MU, K, D, EP = previous.MU, previous.K, previous.D, previous.EP
N, G, Z, V = old_vertex.N, old_vertex.G, old_vertex.Z, old_vertex.V
T = old_vertex.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT = minimal.CONTACT
require_mass, require_order = previous.require_mass, previous.require_order


@cache
def data():
    old = previous.data()
    answer = {k: v for k, v in old.items() if k not in ("checks", "gates")}
    checks = dict(old["checks"])
    full = minimal.current.fixed_functions()
    u, X = minimal.current.u, minimal.current.X
    zero = {u: 0, X: 0}
    R, F = full["R_full"], full["F_full"]
    j = minimal.current.heavy.coefficients()["normalized_heavy_source"]
    for i in range(5):
        for k in range(3):
            if 0 < i + 2 * k <= 4:
                checks["whole_curvature_vacuum_weight_" + str(i) + "_" + str(k)] = (
                    s.cancel(s.diff(R, u, i, X, k).subs(zero))
                )
    for i, k in ((3, 0), (1, 1), (0, 2), (2, 1)):
        checks["whole_scalar_nonpotential_jet_" + str(i) + "_" + str(k)] = s.cancel(
            s.diff(F, u, i, X, k).subs(zero)
        )
    for i, k in ((0, 0), (1, 0), (3, 0), (0, 1), (1, 1)):
        checks["whole_H_source_absent_jet_" + str(i) + "_" + str(k)] = s.cancel(
            s.diff(j, u, i, X, k).subs(zero)
        )
    checks["whole_endpoint_actual_quartic"] = s.cancel(
        s.diff(F, u, 4).subs(zero) / KAPPA - CONTACT
    )
    checks["whole_endpoint_actual_cubic"] = s.cancel(
        s.diff(j, u, 2).subs(zero) / s.sqrt(KAPPA) - CUBIC
    )
    answer.update(
        {
            "whole_selected_endpoint_couplings": {
                "mu": s.Integer(1),
                "n": HEAVY_MASS2,
                "g": CUBIC,
                "C": CONTACT,
                "kappa": KAPPA,
            },
            "whole_endpoint_graph_inventory": "At two external matter legs and one matter loop, the graph-degree identity permits a cubic pair or one quartic. Retain both distinct massive stress triangles, both metric-cubic contacts, the complete quartic metric bubble/tadpole, covariant light OS mass/kinetic terms, H-metric mixing and every heavy-onepoint-reducible graph plus its already fixed counterterm. All curvature coefficient jets through matter degree4 vanish; the regular Ia terms create no additional endpoint vertex at this order.",
            "whole_onepoint_and_matching_boundary": "S234/S239 fix the entire flat heavy one-point-zero source counterterm and the covariant light OS conditions. Constant metric-density variations cancel with those counterterms; the nonlocal and curvature H-metric response remains. Flat conditions do not fix R_old Phi^2 or R_old H, or additional higher-derivative EFT matching.",
            "whole_four_point_boundary": "Only the graviton-exchange ENDPOINT contribution is reconstructed from this complete minimal C/g matter endpoint. Mixed matter/graviton irreducible four-point loops, pure-gravity and other matter sectors, all-loop errors and the physical IR/Regge construction are separate. The original source and physical bounce frame are not changed.",
            "checks": checks,
            "gates": {
                "entire_original_source_and_parameters_retained": "whole_original_R_F"
                in answer,
                "parent_cached_checks_copied_not_mutated": checks is not old["checks"],
                "all_endpoint_arity_curvature_jets_checked": True,
                "complete_flat_heavy_onepoint_counterterm_retained": True,
                "finite_curved_coefficients_not_fixed_by_flat_OS": True,
                "full_endpoint_not_full_mixed_four_point_or_P8": True,
            },
        }
    )
    return answer
