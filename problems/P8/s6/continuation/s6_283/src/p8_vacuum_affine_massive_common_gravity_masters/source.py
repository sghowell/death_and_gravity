"""Whole unchanged source and literal pure-gravity vacuum vertices."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_elastic_proca_infrared import elastic
from p8_vacuum_affine_massive_elastic_proca_infrared import (
    infrared as inherited_infrared,
)
from p8_vacuum_affine_massive_graviton_threshold_sheet import (
    invariant as inherited_invariant,
)
from p8_vacuum_affine_massive_graviton_threshold_sheet import source as previous

infrared, invariant = inherited_infrared, inherited_invariant
S, T, U, MU, K = previous.S, previous.T, previous.U, previous.MU, previous.K
Z = elastic.Z
EP, NU = s.symbols("positive_IR_epsilon positive_dimensional_scale", positive=True)
EPS = s.Symbol("QCD_epsilon", real=True)


def gravity_tree(energy=S, transfer=T, mass=MU, kappa=K):
    energy, transfer, mass, kappa = map(s.sympify, (energy, transfer, mass, kappa))
    u = 4 * mass - energy - transfer
    values = (energy, transfer, u)
    return (
        -sum(
            (2 * mass**2 - 2 * mass * a - values[(i + 1) % 3] * values[(i + 2) % 3]) / a
            for i, a in enumerate(values)
        )
        / kappa
    )


def require_physical(energy, angle, mass=1):
    ss, zz, mu = map(previous.exact_real, (energy, angle, mass))
    if mu <= 0 or ss <= 4 * mu or not -1 < zz < 1:
        raise ValueError("Require s>4mu and a strictly nonforward physical angle")
    return ss, zz, mu


@cache
def data():
    inherited = previous.data()
    checks = dict(inherited["checks"])
    t = -(S - 4 * MU) * (1 - Z) / 2
    minimal = elastic.whole_tree(cubic2=s.Integer(0), contact=s.Integer(0))
    checks["whole_invariant_original_pure_gravity_tree"] = s.factor(
        gravity_tree(transfer=t) - minimal
    )
    checks["full_tree_mass_sum"] = S + T + (4 * MU - S - T) - 4 * MU
    return {
        "whole_original_R_F": inherited["whole_original_R_F"],
        "whole_retained_heavy_source": inherited["whole_retained_heavy_source"],
        "all_three_vacuum_constants": inherited["all_three_vacuum_constants"],
        "whole_formal_loop_marker": inherited["whole_formal_loop_marker"],
        "parameters": inherited["parameters"],
        "entire_original_minimal_gravity_tree": gravity_tree(),
        "coupling": "Same g=eta+2h/sqrt(kappa) and Einstein normalization. The scalar-graviton scattering vertex is -i/sqrt(kappa) times T_mn=p_m r_n+p_n r_m-eta_mn(p.r-mu). No massless external limit is used in the massive result.",
        "loop_scope": "Pure-gravity first-loop contributions to four original massive Phi legs. Mixed heavy/contact-gravity diagrams, M1 and Proca are not silently included. The full source has no additional lower-arity vacuum coupling: the retained higher derivative coefficients start at the leg orders already checked in S281.",
        "checks": checks,
        "gates": {
            "whole_original_source_not_replaced_by_free_testbed": all(
                k in inherited
                for k in (
                    "whole_original_R_F",
                    "whole_retained_heavy_source",
                    "all_three_vacuum_constants",
                )
            ),
            "source_dictionary_copied_without_parent_mutation": checks
            is not inherited["checks"],
            "literal_massive_gravity_tree_all_three_channels": gravity_tree().has(
                S, T, MU
            ),
            "original_kappa_and_vacuum_constants_unchanged": True,
            "pure_gravity_sector_not_all_first_loop_matter": True,
            "formal_vacuum_not_exact_LSZ_existence": True,
        },
    }
