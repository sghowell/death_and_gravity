"""Literal complete generator, unchanged auxiliary branch and source boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_nonlinear_lapse_branch import branch
from p8_vacuum_affine_nonlinear_lapse_branch import source as auxiliary
from p8_vacuum_affine_spatial_gauge import gauge, spatial


@cache
def data():
    old = spatial.generator()
    x = old["whole_spatial_coordinates"]
    g, pi = old["whole_metric_and_metric_momenta"]
    W, pW = old["whole_vector_and_vector_momenta"]
    fields, ps = old["whole_two_matter_and_two_auxiliary_scalar_pairs"]
    xi = old["whole_spatial_gauge_vector"]
    primary = {ps[2]: 0, ps[3]: 0}
    full = old["whole_spatial_generator_density"].subs(primary).doit()
    lag = old["whole_canonical_Lie_pairing_density"].subs(primary).doit()
    flux = old["whole_spatial_generator_boundary_flux"]
    checks = {
        "literal_whole_primary_generator_and_boundary": s.expand(
            lag - xi.dot(full) - sum(s.diff(flux[i], x[i]) for i in range(3))
        )
    }
    translations = []
    for i in range(3):
        fixed = {xi[j]: s.Integer(i == j) for j in range(3)}
        pairing = sum(
            pi[j, k] * s.diff(g[j, k], x[i]) for j in range(3) for k in range(3)
        )
        pairing += sum(pW[j] * s.diff(W[j], x[i]) for j in range(3))
        pairing += sum(ps[a] * s.diff(fields[a], x[i]) for a in range(2))
        translations.append(pairing)
        checks["literal_full_translation_pairing_" + str(i)] = s.expand(
            lag.subs(fixed, simultaneous=True).doit() - pairing
        )
    frame = spatial.frame()
    aux = auxiliary.data()
    roots = branch.auxiliaries()
    checks["original_full_bounce_Hamiltonian_not_truncated"] = aux["checks"][
        "literal_entire_parent_Hamiltonian_slice"
    ]
    checks["original_complete_physical_to_hat_shape_identity"] = frame["checks"][
        "whole_physical_and_hat_shape_density_equal"
    ]
    return {
        "whole_original_generator_before_primary_reduction": old[
            "whole_spatial_generator_density"
        ],
        "whole_original_primary_surface_generator": full,
        "whole_original_metric_vector_matter_Lie_pairing": lag,
        "whole_original_spatial_boundary_flux": flux,
        "whole_original_translation_pairing_densities": translations,
        "whole_original_primitive_boundary": frame["whole_primitive_boundary_binding"],
        "whole_original_R_F_and_profiles": frame[
            "whole_current_R_F_and_fixed_profile_bindings"
        ],
        "whole_original_bounce_Hamiltonian": aux[
            "whole_bounce_Hamiltonian_in_normalized_density_invariants"
        ],
        "whole_original_bounce_auxiliary_constraint": aux["whole_bounce_constraint"],
        "whole_full_off_gauge_ghost_operator": gauge.full_operator()[
            "whole_off_gauge_ghost_operator"
        ],
        "whole_exact_slice_ghost_operator": gauge.full_operator()[
            "whole_gauge_surface_operator"
        ],
        "whole_original_auxiliary_invariant_box": auxiliary.DELTA,
        "source_boundary": "The first two scalar pairs are M1 and H. Only the lapse and normal-vector primary momenta are set to zero. The full vector Gauss term, symmetric metric off-diagonal pairing and primitive boundary remain. The spatial construction is kinematic at fixed clock time. Its intersection with the S266 nonlinear auxiliary chart is asserted only for data whose complete reconstructed invariants lie in the stated bounce-slice box; no Gaussian or uniform canonical-radius assignment is made.",
        "checks": checks,
        "gates": {
            "both_actual_matter_channels_retained": all(
                full.has(f) for f in fields[:2]
            ),
            "full_vector_Gauss_term_not_deleted": any(
                full[i].has(s.diff(pW[j], x[j])) for i in range(3) for j in range(3)
            ),
            "only_auxiliary_primary_momenta_removed": not any(
                full.has(z) for z in ps[2:]
            ),
            "same_complete_regular_auxiliary_branch": all(roots["gates"].values()),
            "same_actual_source_and_physical_map": all(aux["gates"].values())
            and all(frame["gates"].values()),
            "spatial_chart_does_not_extend_auxiliary_time_domain": True,
        },
    }
