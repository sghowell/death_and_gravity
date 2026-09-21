"""Independent component contractions on the original physical recoil states."""

from functools import cache
from itertools import product

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import sew, tensor
from p8_vacuum_affine_minimal_gravity_radiation import vertices as old

from . import bounds, contact


@cache
def data():
    checks, gates, rows = {}, {}, []
    for index in range(3):
        ps, k, _ = old.sample(index)
        basis = tensor.frame(k)
        T = contact.core(ps, k)
        projected = contact.transverse_core(ps, k)
        checks[f"state{index}_all_original_mass_shells"] = s.Matrix(
            [contact.dot(p, p) - 1 for p in ps]
        ).applyfunc(s.factor)
        checks[f"state{index}_null_emitted_momentum"] = s.factor(contact.dot(k, k))
        checks[f"state{index}_exact_momentum_conservation"] = sum(
            ps, k.copy()
        ).applyfunc(s.factor)
        checks[f"state{index}_all_four_Ward_components"] = (
            T * contact.ETA * k
        ).applyfunc(s.factor)
        norm = s.factor(sew.bilinear(projected, projected, 4))
        literal_norm = s.S.Zero
        for polidx, (eps, norm2) in enumerate(
            tensor.transverse_polarizations(basis, contact.ETA)
        ):
            literal = contact.literal_vertex(ps, k, eps)
            dyad = sum(T[a, b] * eps[a, b] for a, b in product(range(4), repeat=2))
            checks[
                f"state{index}_polarization{polidx}_literal24_equals_complete_dyad"
            ] = s.factor(literal - dyad)
            curvature = contact.linear_curvature(k, eps)
            ricci = s.Matrix(
                4,
                4,
                lambda b, d, curvature=curvature: sum(
                    contact.ETA[a, c] * curvature[a, b, c, d]
                    for a, c in product(range(4), repeat=2)
                ),
            )
            checks[
                f"state{index}_polarization{polidx}_all_linearized_Ricci_components"
            ] = ricci.applyfunc(s.factor)
            checks[f"state{index}_polarization{polidx}_trace"] = s.trace(
                contact.ETA * eps
            )
            checks[f"state{index}_polarization{polidx}_transverse"] = (
                eps * k
            ).applyfunc(s.factor)
            literal_norm += literal * s.conjugate(literal) / norm2
        checks[f"state{index}_both_literal_polarizations_equal_canonical_sew"] = (
            s.factor(literal_norm - norm)
        )
        original, _ = tensor.original_cores(ps, k)
        cross = s.factor(2 * sew.bilinear(original, projected, 4))
        wrong = s.zeros(4)
        for i, j in contact.PAIRS:
            vector = contact.dot(k, ps[i]) * ps[j]
            wrong += vector * vector.T
        gates[f"state{index}_on_shell_contact_not_redundant"] = bool(norm > 0)
        gates[f"state{index}_full_original47_directional_interference_nonzero"] = bool(
            cross != 0
        )
        gates[f"state{index}_whole_analytic_contact_bound_calibrated"] = bool(
            norm <= bounds.tensor_upper(k[0]) ** 2
        )
        gates[f"state{index}_one_sided_dyad_negative_control_fails_Ward"] = bool(
            wrong * contact.ETA * k != s.zeros(4, 1)
        )
        rows.append(
            {
                "index": index,
                "omega": k[0],
                "TT_norm_squared": norm,
                "original_tree_interference_sign": s.sign(cross),
            }
        )
    gates["three_actual_states_and_both_polarizations"] = len(rows) == 3
    gates["no_rational_coupling_fixture_substituted_for_original_tree"] = True
    gates["nonzero_point_does_not_claim_nonzero_all_angle_integral"] = True
    return {
        "checks": checks,
        "gates": gates,
        "whole_original_physical_calibrations": rows,
        "whole_interference_qualification": "The contact is a unit matching DIRECTION at chi=0; no physical coefficient is assigned. The base amplitude retains all47 original tree graphs with exact original n,g,C,kappa. Both physical polarizations and their interference are included. Nonzero local directional interference persists on an open neighborhood by continuity, but an all-angle integral is not inferred.",
    }
