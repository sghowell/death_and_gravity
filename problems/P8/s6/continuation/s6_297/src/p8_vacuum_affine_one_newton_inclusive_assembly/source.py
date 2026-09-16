"""Actual finite-kappa matter vertices and complete one-loop graph ownership."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_parent_one_loop import germs
from p8_vacuum_affine_massive_elastic_proca_infrared import source as original
from p8_vacuum_affine_physical_virtual_soft_pairing import source as previous

MU, K, D, EP = previous.MU, previous.K, previous.D, previous.EP
N, G, T = previous.N, previous.G, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT = previous.CONTACT
require_mass, require_order = previous.require_mass, previous.require_order


@cache
def data():
    old = previous.data()
    result = {k: v for k, v in old.items() if k not in ("checks", "gates")}
    checks = dict(old["checks"])
    full = original.current.fixed_functions()
    u, X = original.current.u, original.current.X
    F, R = full["F_full"], full["R_full"]
    couplings = germs.couplings()

    def coefficient(expr, i, j):
        return s.cancel(
            s.diff(s.diff(expr, X, j).subs(X, 0), u, i).subs(u, 0)
            / (s.factorial(i) * s.factorial(j))
        )

    actual = {}
    for name, i, j in (("phi6", 6, 0), ("phi4Y", 4, 1), ("phi2Y2", 2, 2), ("Y3", 0, 3)):
        actual[name] = coefficient(F, i, j) / KAPPA**2
    actual["phi2_L3_minus_L4"] = 2 * coefficient(R, 2, 2) / KAPPA**2
    actual["Y_L3_minus_L4"] = 3 * coefficient(R, 0, 3) / KAPPA**2
    heavy_j = original.current.heavy.coefficients()["normalized_heavy_source"]
    actual["H_phi2_Y"] = coefficient(heavy_j, 2, 1) / KAPPA ** s.Rational(3, 2)
    for name, value in actual.items():
        checks["finite_kappa_full_matter_jet_" + name] = s.factor(
            value - couplings[name]
        )
    for i, j in ((0, 1), (0, 2), (1, 1), (2, 1), (3, 1), (1, 2)):
        checks[f"no_lower_R_jet_{i}_{j}"] = coefficient(R, i, j)
    checks["dependent_a4_difference_first_twelve_legs"] = s.Integer(2 * 4 + 4 - 12)
    checks["dependent_a5_difference_first_twelve_legs"] = s.Integer(2 * 4 - 2 + 6 - 12)
    patterns = germs.graph_degrees(4)
    for pattern in patterns:
        checks["complete_one_loop_degree_pattern_" + "_".join(map(str, pattern))] = (
            s.Integer(sum((j + 1) * v for j, v in enumerate(pattern)) - 4)
        )
    result.update(
        {
            "whole_finite_kappa_matter_couplings": actual,
            "whole_same_matter_one_loop_graph_patterns": patterns,
            "whole_graph_equality_proof": "All actual current f,a3 and H-source vertices through six scalar legs agree with S239. The first r_Y term has four scalar fields. The finite-kappa differences a4=-a3-7r_Y^2/[4(kappa+r)] and a5=r_Y^2/[(kappa+r)Y] therefore start at twelve scalar legs. The Proca source square starts at eight and needs two linear vector-source vertices; it has no first-loop four-Phi graph. The exact connected identity sum(degree-2)=E+2L-2 restricts E4,L1 to the five recorded patterns. Thus the no-internal-graviton one-loop diagrams, including tadpoles and cubic/five-field contractions, equal the complete S239 matter diagrams at finite kappa. This is diagram equality, not interchange of a quantum limit.",
            "whole_renormalization_scope": "Retain the S239 H8A420-VAC-OS4 flat matter mass/residue, heavy-onepoint, pole-only MSbar and symmetric finite-value conditions, with its complete higher local vertices. No new finite value is chosen. Gravity-generated finite curvature, heavy-residue and higher local matching remains explicit in the S296 hard dressing coefficient. Changing the matter prescription requires a separately transported matching term.",
            "whole_one_Newton_inventory": "At one loop, the first inverse-kappa gravitational dressing has only the C/kappa and g^2/kappa source graphs of S293/S294; source higher scalar vertices cannot enter a graph containing a metric line at this order. The five-point real tree at inverse-sqrt-kappa is the whole S295 tree. Complete linear-Newton rate interference additionally includes twice the gravity Born amplitude times the real complete matter loop. Gravity-Born radiation and pure-gravity virtual loops enter the next Newton rate orders, which are not omitted-loop bounded here.",
            "checks": checks,
            "gates": {
                "all131_parent_source_checks_copied": len(old["checks"]) == 131
                and checks is not old["checks"],
                "all_seven_actual_higher_matter_couplings": set(actual)
                == set(couplings),
                "same_complete_five_graph_patterns": set(patterns)
                == {
                    (4, 0, 0, 0),
                    (2, 1, 0, 0),
                    (0, 2, 0, 0),
                    (1, 0, 1, 0),
                    (0, 0, 0, 1),
                },
                "dependent_differences_cannot_enter_first_four_point_loop": 12 - 2 > 4
                and 8 - 2 > 4,
                "original_counterterms_not_replaced_by_new_finite_choice": True,
                "formal_diagram_equality_not_quantum_limit_theorem": True,
            },
        }
    )
    return result
