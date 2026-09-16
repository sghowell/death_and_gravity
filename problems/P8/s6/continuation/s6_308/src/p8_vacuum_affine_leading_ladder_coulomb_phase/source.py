"""Unchanged source and explicit leading-ladder ownership."""

from functools import cache

from p8_vacuum_affine_finite_forward_phase import phase as old_phase
from p8_vacuum_affine_uniform_radiation_soft_limit import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
born = previous.born
__all__ = [
    "CONTACT",
    "CUBIC",
    "EP",
    "HEAVY_MASS2",
    "KAPPA",
    "MU",
    "G",
    "K",
    "N",
    "T",
    "born",
    "data",
    "old_phase",
    "require_mass",
    "require_order",
]


@cache
def data():
    old = previous.data()
    result = {
        key: value for key, value in old.items() if key not in ("checks", "gates")
    }
    checks = dict(old["checks"])
    checks.update(
        {
            "same_original_Coulomb_Newton_parameter": KAPPA - old_phase.source.KAPPA,
            "same_original_full_Born_Newton_parameter": KAPPA - born.source.KAPPA,
            "same_original_physical_light_mass": MU - old_phase.source.MU,
        }
    )
    result.update(
        {
            "whole_selected_graph_ownership": "The leading straight-worldline/free-graviton ladder and crossed-ladder class is constructed from the same original full-D pure-Einstein massive-scalar Born pole. It is a named leading-forward graph class, not a replacement for the original nonlinear action or a complete interacting quantum amplitude. S303 owns the independent entire one-loop principal Coulomb phase used for calibration.",
            "whole_other_sectors_not_deleted": "Finite Newton matching and spectator-induced Newton shifts, recoil, nonlinear graviton vertices, hard loops, real radiation and short-distance operators remain outside this selected leading class. The S307 single-real uniform bound and S299 all-N leading-soft rate retain their distinct ownership and do not close these omissions.",
            "whole_reference_convention": "Original mu=nu=1, Ddim=4+2e, kappa=1e800 and25/4<=s<=16. The same analytic IR phase is divided before taking e->0 at each fixed perturbative order. The resulting leading-class coefficients are summed for|eta|<1. No massless-scalar or trans-Planckian classical hierarchy is imposed on the original parameters.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "new_source_checks_not_parent_alias": checks is not old["checks"],
                "named_leading_ladder_class_not_full_QFT": True,
                "original_nonlinear_action_and_matching_not_replaced": True,
                "full_one_loop_S303_phase_is_independent_calibration": True,
            },
        }
    )
    return result
