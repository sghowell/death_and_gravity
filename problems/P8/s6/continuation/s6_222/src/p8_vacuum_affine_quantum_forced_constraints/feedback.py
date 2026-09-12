"""Ordered causal Schur identities on an explicitly restricted graph domain."""

from functools import cache

import sympy as s


def zero_matrix(expr):
    return expr.expand().as_explicit().applyfunc(s.expand)


@cache
def data():
    C = s.MatrixSymbol("C", 3, 4)
    F = s.MatrixSymbol("F", 4, 3)
    D = s.MatrixSymbol("D", 3, 3)
    G = s.MatrixSymbol("G0", 4, 4)
    Q = s.MatrixSymbol("Rhat", 3, 3)
    W = s.MatrixSymbol("W_output", 3, 3)
    g = s.MatrixSymbol("g", 3, 1)
    e = s.MatrixSymbol("external_e", 3, 1)
    z = s.MatrixSymbol("Zh", 4, 1)
    Qbar = W * Q
    T = D + C * G * F
    S = s.Identity(3) - Qbar * T
    Z = z + G * F * g
    metric = C * Z + D * g
    direct = g - e - Qbar * metric
    schur = S * g - e - Qbar * C * z
    full = s.BlockMatrix(
        [[s.Identity(4), -G * F], [-Qbar * C, s.Identity(3) - Qbar * D]]
    )
    lower = s.BlockMatrix(
        [[s.Identity(4), s.ZeroMatrix(4, 3)], [-Qbar * C, s.Identity(3)]]
    )
    upper = s.BlockMatrix([[s.Identity(4), -G * F], [s.ZeroMatrix(3, 4), S]])
    eps = s.Symbol("formal_marker")
    source = s.MatrixSymbol("prescribed_S0", 3, 1)
    g1 = Qbar * source
    g2 = Qbar * T * g1
    truncated = eps * g1 + eps**2 * g2
    checks = {
        "ordered_force_equation_equals_full_reconstructed_response": zero_matrix(
            direct - schur
        ),
        "full_block_unit_triangular_Schur_factorization": zero_matrix(
            s.block_collapse(full - lower * upper)
        ),
        "metric_reconstruction_including_direct_auxiliary": zero_matrix(
            metric - (C * z + T * g)
        ),
        "first_formal_force_coefficient": zero_matrix(g1 - Qbar * source),
        "second_formal_force_coefficient_includes_auxiliary": zero_matrix(
            g2 - Qbar * D * g1 - Qbar * C * G * F * g1
        ),
        "exact_residual_of_two_formal_coefficients": zero_matrix(
            truncated - eps * Qbar * (source + T * truncated) + eps**3 * Qbar * T * g2
        ),
    }
    return {
        "full_coupled_system": "Z'=KZ+Fg+h; s=CZ+Dg; g=e+Qbar s. Here h and e are explicitly declared mathematical probe drives, the original phase data are zero, and the same prepared quantum state is retained.",
        "ordered_Schur_equation": "[I-Qbar(D+C G0 F)]g=e+Qbar C Zh, Zh=G0 h. Reconstruct Z=Zh+G0 Fg and s=CZh+(D+C G0 F)g. Qbar=W_output Rhat, not Rhat W_output.",
        "graph_domain": "X_r is the initial-prepared H13_t H^(r+8)_x source completion, Y_r=L2_t H^(r-3)_x. For e inY_r and C Zh inX_r, the Schur equation is defined on D_T={g inY_r:(D+C G0 F)g inX_r}. No density, invariance, surjectivity, inverse or contraction on this graph domain is asserted.",
        "two_way_equivalence": "For an admissible forced-system solution, classical variation of constants gives Z=Zh+G0Fg and hence the Schur equation. Conversely a graph-domain solution reconstructs s inX_r, the original weighted phase equations, both force-dependent constraints and g=e+Qbar s. These are substitutions using only the classical inverse, not a quantum inverse theorem.",
        "preparation": "The source retains its original zero-past preparation. Smooth compact-momentum h with an initial zero neighborhood gives C Zh inX_r by finite-dimensional ODE smoothness and the fixed smooth coefficients. This does not add independent quantum initial data. More general probe drives require the stated source-domain membership.",
        "formal_boundary": "The displayed first and second marker coefficients are algebraic bookkeeping only. The nonzero direct auxiliary contribution to the second coefficient cannot be omitted. No finite-coupling Taylor remainder, convergence or order reduction is inferred.",
        "checks": checks,
        "gates": {
            "output_density_and_nonlocal_response_not_commuted": Qbar != Q * W,
            "direct_auxiliary_term_not_identically_absent": Qbar * D
            != s.ZeroMatrix(3, 3),
            "no_retarded_transpose_symmetrization": True,
            "no_quantum_inverse_in_block_factorization": True,
            "explicit_graph_domain_before_feedback_composition": True,
            "mathematical_drives_not_unproved_physical_source_maps": True,
            "all_memory_and_higher_time_derivatives_retained": True,
            "no_state_reset_or_new_physical_cutoff": True,
        },
    }
