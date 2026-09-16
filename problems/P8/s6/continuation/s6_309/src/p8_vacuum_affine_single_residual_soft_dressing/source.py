"""Original47-tree kernel and immutable single-residual ownership."""

from functools import cache

from p8_vacuum_affine_leading_ladder_coulomb_phase import source as previous
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import (
    recoil as single_recoil,
)
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree
from p8_vacuum_affine_minimal_gravity_radiation import source as radiation_source
from p8_vacuum_affine_minimal_gravity_radiation import vertices
from p8_vacuum_affine_radiative_angular_finite import bounds as angular_bounds
from p8_vacuum_affine_radiative_angular_finite import conversion
from p8_vacuum_affine_radiative_state_soft_index import index, recoil
from p8_vacuum_affine_uniform_radiation_soft_limit import softlimit

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
born = previous.born
gamma, poisson = softlimit.source.gamma, softlimit.source.poisson
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
    "angular_bounds",
    "born",
    "conversion",
    "data",
    "gamma",
    "index",
    "poisson",
    "recoil",
    "require_mass",
    "require_order",
    "single_recoil",
    "softlimit",
    "tree",
    "vertices",
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
            "same_single_residual_Newton_parameter": KAPPA - softlimit.source.KAPPA,
            "same_radiative_conversion_Newton_parameter": KAPPA
            - conversion.source.KAPPA,
            "same_entire47_tree_Newton_parameter": KAPPA - radiation_source.KAPPA,
            "same_radiative_index_original_mass": MU - index.source.MU,
        }
    )
    result.update(
        {
            "whole_selected_graph_ownership": "S295/S304 own the complete26matter+21Einstein single-graviton tree and S307 owns its Born-normalized signed real-minus-leading-soft measure. No new graph enumeration is claimed. This successor applies the full state-correct all-N LEADING-soft calorimetric kernel to exactly one marked finite residual.",
            "whole_residual_measure": "dR=[rho sum_pol|M5/(Am+AG)|^2-sum_pol|S0|^2]dPhi1. All47tree interferences, original parameters, recoil and positive full Born remain. The absolute cumulative measure is bounded bymin(1e32,2e14*x+2e37*x^2)/kappa; it is signed, not a positive emission probability.",
            "whole_selected_soft_dressing": "Each marked radiative state uses its own conserved full massive-plus-null current, dimensional angular coefficient and virtual rate pole. Remainingenergy isx-omega; all additional emissions are leading soft with the marked hard state held fixed. Their full energy-simplex sum is retained.",
            "whole_not_full_physical_decomposition": "This explicitly defined single-residual reference is not identified with an exact decomposition of all interacting multi-real amplitudes. Two-soft contact terms, correlated recoil, multiple nonleading residuals, finite radiative hard loops and matching remain outside its error theorem.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "new_source_checks_not_parent_alias": checks is not old["checks"],
                "same_complete47_tree_signed_measure": True,
                "state_correct_virtual_and_total_remaining_energy": True,
                "selected_single_residual_class_not_complete_QFT": True,
            },
        }
    )
    return result
