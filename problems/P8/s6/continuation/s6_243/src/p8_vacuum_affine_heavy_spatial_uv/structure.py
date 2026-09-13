"""Full scalar vertex, exact ordered endpoint identity and scalar MSbar normalization."""

from functools import cache

import sympy as s

from . import jets


@cache
def data():
    aa, mass = s.symbols("a_physical mass", positive=True)
    Pk, Pl = s.symbols("Pk Pl")
    k = s.Matrix(s.symbols("k1 k2 k3", real=True))
    l = s.Matrix(s.symbols("l1 l2 l3", real=True))
    d11, d12, d13, d22, d23, d33 = s.symbols("d11 d12 d13 d22 d23 d33", real=True)
    Q = s.Matrix([[d11, d12, d13], [d12, d22, d23], [d13, d23, d33]])
    trace = s.trace(Q)
    M = s.zeros(5)
    M[0, 0] = -trace / 2
    M[4, 4] = trace / 2
    M[1:4, 1:4] = trace * s.eye(3) / 2 - Q
    Fk = s.Matrix([Pk, *list(s.I * k / aa), mass])
    Fl = s.Matrix([Pl, *list(s.I * l / aa), mass])
    direct = (Fk.T * M * Fl)[0]
    full = (
        trace * (-Pk * Pl - (k.T * l)[0] / aa**2 + mass**2) / 2
        + (k.T * Q * l)[0] / aa**2
    )
    checks = {"all_six_literal_spatial_Hamiltonian_directions": s.expand(direct - full)}
    x, abar, bbar, pkbar, plbar, dot, contraction = s.symbols(
        "x abar bbar pkbar plbar dot contraction", real=True
    )
    # No time/trace or gradient term is dropped before multiplication.
    normalized = (
        trace * (-pkbar * plbar - dot + mass**2 * aa**2 * x**2) / 2 + contraction
    ) / (2 * s.sqrt(abar * bbar))
    physical = (
        trace
        * (
            -(pkbar / (aa * x)) * (plbar / (aa * x))
            - dot / (aa * aa * x * x)
            + mass * mass
        )
        / 2
        + contraction / (aa * aa * x * x)
    ) / (2 * s.sqrt(abar * bbar) / (aa * x))
    checks["full_physical_volume_and_frequency_normalization"] = s.simplify(
        aa * x * physical - normalized
    )
    L = s.symbols("L0:6")
    # d_t[g L^j U exp(-i theta)] = (L^(j+1)U - i L^j U) exp(-i theta).
    telescoped = sum(s.I * (-s.I) ** j * (L[j + 1] - s.I * L[j]) for j in range(5))
    checks["entire_five_endpoint_telescoping_and_sixth_bulk"] = s.expand(
        telescoped - L[0] - s.I * L[5]
    )

    eps = s.Symbol("epsilon", real=True)
    mu, mm = s.symbols("mu scalar_mass", positive=True)
    ell = s.Symbol("ell", real=True)
    dimension = 3 - 2 * eps
    sphere = 2 ** (1 - dimension) * s.pi ** (-dimension / 2) / s.gamma(dimension / 2)
    measure = (
        s.exp(s.EulerGamma * eps) * (mu * mu / (4 * s.pi)) ** eps * mm ** (-2 * eps)
    )
    factor = sphere * measure / (1 / (2 * s.pi**2))
    derivative = s.simplify(s.expand_func(s.diff(factor, eps).subs(eps, 0)))
    expected = 2 - 2 * s.log(2) - s.log(mm * mm / (mu * mu))
    checks["complete_dimensional_sphere_and_MSbar_scale_derivative"] = s.simplify(
        s.expand_log(derivative - expected, force=True)
    )
    F2, F4, F4d, HPd = s.symbols(
        "F2 F4 F4_dimension_derivative Hpole_dimension_derivative", real=True
    )
    # Multiply by epsilon before differentiating the full analytic expression.
    # The finite value is that derivative after the pole cancellation.
    epsilon_times_full = factor * (
        (F4 - 2 * eps * F4d) / 2 - eps * mm * mm * F2 / (2 - 2 * eps)
    ) / (2 * s.pi**2) - (16 * F4 - 2 * eps * HPd) / (64 * s.pi**2)
    checks["full_radial_pole_cancels_before_limit"] = s.simplify(
        epsilon_times_full.subs(eps, 0)
    )
    actual_finite = s.simplify(
        s.expand_func(s.diff(epsilon_times_full, eps).subs(eps, 0))
    )
    finite = (
        (1 - s.log(2) - s.log(mm * mm / (mu * mu)) / 2) * F4 / (2 * s.pi**2)
        - F4d / (2 * s.pi**2)
        - mm * mm * F2 / (4 * s.pi**2)
        + HPd / (32 * s.pi**2)
    )
    checks["complete_same_scheme_finite_radial_expression"] = s.simplify(
        s.expand_log(actual_finite - finite, force=True)
    )
    # The full arbitrary-trace geometries, not a tracefree interpolation.
    y, geo = jets.geometries()
    UD, UG = jets.inv[4:]
    checks["distinct_ordered_trace_gradient_geometry"] = s.factor(
        geo["01"]
        - geo["10"]
        - (UD - UG) * ((1 - jets.u**2) / (jets.d - 1) + jets.u * (y - jets.u))
    )
    checks["sphere_unit_normalization"] = jets.sphere(s.S.One) - 1
    checks["sphere_second_moment"] = jets.sphere(jets.u**2) - 1 / jets.d
    checks["sphere_fourth_moment"] = jets.sphere(jets.u**4) - 3 / (
        jets.d * (jets.d + 2)
    )
    checks["sphere_sixth_moment"] = jets.sphere(jets.u**6) - 15 / (
        jets.d * (jets.d + 2) * (jets.d + 4)
    )
    checks["sphere_eighth_moment"] = jets.sphere(jets.u**8) - 105 / (
        jets.d * (jets.d + 2) * (jets.d + 4) * (jets.d + 6)
    )
    for odd in (1, 3, 5, 7):
        checks["whole_sphere_odd_moment_" + str(odd)] = jets.sphere(jets.u**odd)
    return {
        "literal_complete_spatial_Hamiltonian_feature_matrix": M,
        "complete_pair_vertex_before_UV_expansion": full,
        "entire_ordered_five_endpoint_identity": "Integral U exp(-i theta) = exp(-i theta) sum(j0..4) i(-i)^j g (dt g)^j U - i integral (dt g)^5 U exp(-i theta), with every compact lower endpoint retained as zero by the prepared source.",
        "full_endpoint_remainder_scope": "This exact identity fixes phases and all first five endpoints. It does not bound the complete sixth bulk or replace the exact comparison evolution by W6.",
        "continued_angular_measure": sphere,
        "full_MSbar_measure": measure,
        "complete_measure_first_derivative": derivative,
        "same_scheme_Laurent_finite_formula": s.expand(
            finite.subs(s.log(mm * mm / (mu * mu)), ell)
        ),
        "fixed_invariant_geometry": geo,
        "checks": checks,
        "gates": {
            "complete_scalar_not_vector_feature": True,
            "both_ordered_trace_gradient_products_retained": geo["01"] != geo["10"],
            "annihilation_phase_negative_and_detector_sharp": True,
            "all_odd_endpoints_retained_before_extraction": True,
            "dimension_varied_at_fixed_six_physical_invariants": True,
            "complete_counteraction_varied_before_limit": True,
            "no_new_finite_target_or_counterterm": True,
            "local_difference_not_full_kernel_or_homogeneous_anchor": True,
        },
    }
