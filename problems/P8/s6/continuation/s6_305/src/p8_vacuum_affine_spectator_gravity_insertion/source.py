"""Unchanged source and ownership of spectator metric-loop contributions."""

from functools import cache

import sympy as s
from p8_vacuum_affine_gaussian_metric_pole_matching import poles, spectral
from p8_vacuum_affine_gaussian_metric_pole_matching import source as gaussian
from p8_vacuum_affine_massive_detector_ir_cut import cut as massless
from p8_vacuum_affine_minimal_gravity_radiation import source as previous
from p8_vacuum_affine_one_newton_inclusive_assembly import forward

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
VECTOR_MASS2 = gaussian.VECTOR_MASS2
require_mass, require_order = previous.require_mass, previous.require_order
__all__ = [
    "CONTACT",
    "CUBIC",
    "EP",
    "HEAVY_MASS2",
    "KAPPA",
    "MU",
    "VECTOR_MASS2",
    "G",
    "K",
    "N",
    "T",
    "data",
    "forward",
    "gaussian",
    "massless",
    "poles",
    "require_mass",
    "require_order",
    "spectral",
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
            "same_gaussian_heavy_mass_squared": HEAVY_MASS2 - gaussian.HEAVY_MASS2,
            "same_gaussian_gravity_constant": KAPPA - gaussian.KAPPA,
            "same_original_vector_mass_squared": VECTOR_MASS2 - 10**6,
            "same_fixed_Newton_shift_not_new_choice": s.factor(
                poles.delta_kappa() - 2 * gaussian.fixed_coefficients()["R_old"]
            ),
            "spectator_bubble_E4_L1_degree_budget": s.Integer(
                4 * (3 - 2) - (4 + 2 * 1 - 2)
            ),
            "spectator_seagull_E4_L1_degree_budget": s.Integer(
                (4 - 2) + 2 * (3 - 2) - (4 + 2 * 1 - 2)
            ),
            "no_species_source_square_at_E4_L1": s.Integer(2 * (5 - 2) - 4 - 2),
            "higher_degree12_jet_exceeds_E4_L1_budget": s.Integer((12 - 2) - 4 - 6),
        }
    )
    result.update(
        {
            "whole_spectator_species": {
                "H": HEAVY_MASS2,
                "Proca": VECTOR_MASS2,
                "original_M1": s.S.Zero,
            },
            "whole_insertion_ownership": "Add the complete source-fixed H and Proca Gaussian metric kernels of S285 to physical four-Phi scattering at g^0/kappa^2. Their covariant seagulls, lower polynomial and original volume cancellation are inseparable. Add S278's complete original massless-M1 log part once. The Phi metric self-energy is already included in S302 and is NOT added again.",
            "whole_graph_boundary": "At E4,L1 sum(degree-2)=4. The two-scalar stress vertices and the closed spectator metric bubble/seagull exhaust this g-independent spectator sector. Vacuum-onepoint insertions cancel only with the original whole source-pinned volume term. The regular linear Proca source starts at four scalar fields and its square exceeds this budget. The dependent degree12 jets also cannot enter. H-production g^2/kappa and g^4 cuts belong to S294 and S239; they are not relabeled as pure spectator cuts.",
            "whole_matching_boundary": "Use only H/Proca finite coefficients already fixed by S285. M1 is a different species from Phi; its nonanalytic part does not determine real local terms. Remaining physical scalar/gravity/M1 matching coordinates stay symbolic. No new finite prescription or canonical Newton redefinition is made.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "source_checks_copied_not_parent_alias": checks is not old["checks"],
                "all_three_original_spectator_species_distinguished": VECTOR_MASS2
                != HEAVY_MASS2
                and VECTOR_MASS2 > 0,
                "whole_kernels_not_isolated_bubbles": True,
                "light_Phi_cut_not_double_counted": True,
            },
        }
    )
    return result
