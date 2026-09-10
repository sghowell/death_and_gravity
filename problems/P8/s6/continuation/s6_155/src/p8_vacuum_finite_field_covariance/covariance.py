"""Bare vertex/line counting, complete renormalized amplitudes and pole products."""

from functools import cache

import sympy as s


@cache
def data():
    K, h, e, I = s.symbols("K h epsilon internal_Phi_lines", positive=True)
    k0, k1, t0, A0, A1, A2, pole, finite = s.symbols("k0 k1 t0 A0 A1 A2 pole finite")
    kD = k0 + e * k1
    raw = pole / e + finite
    ct = -pole / e
    graph_vertex = -(I + 2) * kD * raw
    graph_lines = I * kD * raw
    combined = graph_vertex + graph_lines
    amp = (1 - 2 * h * k0 + h**2 * (3 * k0**2 - 2 * t0)) * (A0 + h * A1 + h**2 * A2)
    first = s.expand(amp).coeff(h, 1)
    second = s.expand(amp).coeff(h, 2)
    checks = {
        "all_Phi_half_edges_accounted": 2 * I + 4 - 2 * (I + 2),
        "vertex_plus_internal_line_covariance": s.expand(combined + 2 * kD * raw),
        "finite_renormalized_product_before_epsilon_zero": s.expand(
            s.limit(-2 * kD * (raw + ct), e, 0) + 2 * k0 * finite
        ),
        "raw_finite_part_contains_epsilon_pole_product": s.expand(
            s.expand(-2 * kD * raw).coeff(e, 0) + 2 * k0 * finite + 2 * k1 * pole
        ),
        "counterterm_finite_part_cancels_pole_product": s.expand(
            s.expand(-2 * kD * ct).coeff(e, 0) - 2 * k1 * pole
        ),
        "canonical_amplitude_first_coefficient": s.expand(first - A1 + 2 * k0 * A0),
        "canonical_amplitude_second_coefficient": s.expand(
            second - A2 + 2 * k0 * A1 - (3 * k0**2 - 2 * t0) * A0
        ),
        "single_canonical_field_factor_second_coefficient": s.expand(
            second - (A2 - 2 * k0 * A1 + (3 * k0**2 - 2 * t0) * A0)
        ),
    }
    L, g, Y, F1, F2, F3, Ff = s.symbols("L g Y kernel1 kernel2 kernel3 fermion_kernel")
    # Full heavy propagators and external invariants are in the kernels:
    # field rescaling changes none of them.
    scalar = L**2 * F1 + L * g * F2 + g**2 * F3
    fermion = Y**2 * Ff
    Euler = lambda f: 2 * L * s.diff(f, L) + 2 * g * s.diff(f, g) + Y * s.diff(f, Y)
    checks["complete_scalar_one_loop_vertex_weight"] = s.expand(
        Euler(scalar) - 4 * scalar
    )
    checks["complete_fermion_box_vertex_weight"] = s.expand(
        Euler(fermion) - 2 * fermion
    )
    M = s.Symbol("M")
    margin = L - 3 * g / M
    checks["common_field_factor_preserves_completed_square_margin"] = s.factor(
        L / K**2 - 3 * (g / K**2) / M - margin / K**2
    )
    return {
        "one_loop_scalar_Phi_line_count": 2,
        "one_loop_fermion_box_Phi_line_count": 0,
        "general_regulated_vertex_variation": graph_vertex,
        "general_regulated_line_variation": graph_lines,
        "canonical_amplitude_first_coefficient": first,
        "canonical_amplitude_second_coefficient": second,
        "first_normalization_only_second_amplitude": 3 * k0**2 * A0 - 2 * k0 * A1,
        "new_second_normalization_amplitude_kept_separate": -2 * t0 * A0,
        "checks": checks,
        "scope": "At the common regulator a Phi^4 graph has sum valences=2 I_Phi+4. Vertices give K^-(I_Phi+2), internal propagators K^I_Phi, so the complete amputated graph gives K^-2. Assigned counterterms are transformed with the same regulator. After the full one-loop sum is finite, its physical amplitude uses k0 only. Do not add the equivalent finite coefficient-map commutator again.",
    }
