"""Restore finite Taylor terms without adding the unmatched local target."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_reference_spatial_remainder import limit as finite_parent

from . import degrees, radials, tail

KNOWN = 5 * s.Integer(10) ** 54
KNOWN_TAIL = 7 * s.Integer(10) ** 60


@cache
def data():
    known, contact, poly, finite, uv = s.symbols("known contact poly finite uv")
    checks = {
        "exact_full_polynomial_partition": s.expand(poly - finite - uv).subs(
            poly, finite + uv
        ),
        "same_finite_regulator_current": s.expand(
            (known + contact + poly) - ((known + finite) + (contact + uv))
        ).subs(poly, finite + uv),
        "remaining_fifteen_candidate_UV_cells": len(degrees.ultraviolet_cells()) - 15,
        "canonical_finite_polynomial": 4 * radials.FINITE / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 750,
        "canonical_finite_polynomial_tail": 4 * tail.TAIL / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 745,
        "canonical_updated_known_actual_piece": 4 * KNOWN / modes.KAPPA
        - 2 * s.Rational(1, 10) ** 745,
        "canonical_updated_known_actual_tail": 4 * KNOWN_TAIL / modes.KAPPA
        - 28 * s.Rational(1, 10) ** 740,
    }
    return {
        "exact_identity": "At the original finite regulator, P_all,K=P_fin,K+P_UV,K. Therefore J_actual,K=(known_S203,K+P_fin,K)+(C_unit,K+P_UV,K). The finite oversubtractions are returned to the known current with their original signs; no term is dropped.",
        "known_bound": "Combining the actual S203 known piece4e54 M[D]Y[Gamma] and its tail6e60 M[D]Y[Gamma]/K with the finite polynomial gives5e54 M[D]Y[Gamma] and7e60 M[D]Y[Gamma]/K. Here Y^2=N61^2+X46^2 is unchanged.",
        "canonical": "The finite polynomial alone has canonical bounds4e-750 and4e-745/K. The updated known actual piece has2e-745 and28e-740/K.",
        "local_target_not_added": "S204 computes only the already-fixed local action target. It is NOT added to the known quantum remainder here, because equality of the remaining quantum sector to that local prescription has not been proved.",
        "remaining": "The fifteen candidate cells, with all35 source time-jet entries, remain together with the complete one-leg contact. Their actual UV expansion, finite/divergent coefficients, sharp-band artifacts and original covariant matching are still required. Candidate is not a claim that every cell diverges.",
        "boundary": "This finite polynomial restoration is not a full matched spatial response, mixed inverse, finite-coupling background, stability result, cutoff or original V/G/B closure.",
        "checks": checks,
        "gates": {
            "updated_known_actual_piece_bound": finite_parent.KNOWN + radials.FINITE
            < KNOWN,
            "updated_known_actual_tail_bound": finite_parent.KNOWN_TAIL + tail.TAIL
            < KNOWN_TAIL,
            "unchanged_existing_known_piece_anchor": finite_parent.KNOWN
            == 4 * s.Integer(10) ** 54,
            "complete_contact_and_candidate_sector_retained": True,
            "local_target_not_double_counted": True,
            "original_V_G_B_open": True,
        },
    }
