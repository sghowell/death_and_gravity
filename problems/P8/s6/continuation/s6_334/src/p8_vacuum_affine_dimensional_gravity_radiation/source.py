"""Unchanged source and credit for pre-existing general-D matter results."""

from functools import cache

import sympy as s
from p8_vacuum_affine_joint_soft_regulator import source as previous
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import source as matter

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_energies, require_multiplicity = (
    previous.require_energies,
    previous.require_multiplicity,
)
original_parameters = previous.original_parameters


@cache
def data():
    checks = dict(previous.data()["checks"])
    for name in ("KAPPA", "HEAVY_MASS2", "CUBIC", "CONTACT"):
        checks["same_original_matter_" + name] = s.simplify(
            globals()[name] - getattr(matter, name)
        )
    checks["selected_all47_graph_count"] = s.Integer(26 + 3 * (4 + 2 + 1) - 47)
    e = s.Symbol("epsilon", nonnegative=True)
    checks["original_positive_dimension_ratio"] = s.cancel(
        ((s.Symbol("D") - 4) / (s.Symbol("D") - 2)).subs(s.Symbol("D"), 4 + 2 * e)
        - e / (1 + e)
    )
    return {
        "checks": checks,
        "gates": {
            "original_action_mass_couplings_and_matching_unchanged": True,
            "S295_already_general_D_matter_Ward_and_complete26": True,
            "S296_already_rank_two_sew_and_fixed_disk_matter_measure": True,
            "S304_21_gravity_graphs_extended_in_both_dimension_occurrences": True,
            "auxiliary_massless_scalar_is_identity_not_new_original_particle": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_source_boundary": "Same original canonical minimal action, physical mu=1, fixed couplings and S296 positive D=4+2epsilon convention. New work supplies the gravity21 graph continuation and exact all47 tree sew; earlier S295 matter26 and S296 rank-two measure results are credited. It does not supply finite radiative hard loops, integrated gravity real-minus-soft pairing or original P8 closure.",
    }
