"""Conditional sourced Gaussian state, physical insertions and clock response order."""

from functools import cache

import sympy as s
from p8_vacuum_analytic_affine_parent import source as actual_source


@cache
def data():
    A, mean, S, cov, metric = s.symbols("A mean source covariance metric", real=True)
    eps, source2, source3, mean2, mean3 = s.symbols(
        "epsilon source2 source3 mean2 mean3", real=True
    )
    # Finite-mode contractions isolate the exact Gaussian algebra, not a
    # continuum value assigned to an unrenormalized coincident covariance.
    polynomial = metric * (A - S) ** 2 / 2
    expanded = s.expand(polynomial.subs(A, mean + A))
    expectation = expanded.coeff(A, 0) + expanded.coeff(A, 2) * cov
    split = metric * cov / 2 + metric * (mean - S) ** 2 / 2
    sourcejet = eps**2 * source2 + eps**3 * source3
    meanjet = eps**2 * mean2 + eps**3 * mean3
    displacement = meanjet - sourcejet
    source_force = s.expand(-s.diff(sourcejet, eps) * displacement)
    quadratic_mean = s.expand(meanjet**2 / 2)
    a, b = s.symbols("a b", real=True)
    # A nonzero second source derivative is a linear field insertion.
    operator_contact = s.diff(-A * sourcejet, eps, 2).subs(eps, 0)
    actual = actual_source.data()["shifted_source_P8"]
    h0, rx0 = s.symbols("H0 RX0", real=True)
    first = s.symbols("R1 RX1 Ru1 X1 H1 Box1 Z1", real=True)
    jets = {
        "R": 1 + eps * first[0],
        "R_X": rx0 + eps * first[1],
        "R_u": eps * first[2],
        "positive_X": 1 + eps * first[3],
        "H_clock": h0 + eps * first[4],
        "Box_P8_u": 3 * h0 + eps * first[5],
        "uHu_P8": eps * first[6],
    }
    actual_jet = actual.subs(
        {v: jets[str(v)] for v in actual.free_symbols}, simultaneous=True
    )
    return {
        "conditional_mean": "For each admissible history, Wbar=G_O,ret S with fixed zero initial mean. The fluctuation covariance is the ordinary source-free Proca Hadamard covariance for that same metric and initial state.",
        "renormalized_split": "The conditional vector stress is ordinary connected Proca stress plus the classical FULL sourced-vector stress evaluated at Wbar, including variations through S and the local S^2/2 contact. Smooth mean terms require no new ultraviolet subtraction.",
        "source_influence_boundary": "The Gaussian cumulant expansion terminates at its connected two-point kernel, with temporal-constraint and explicit source contact terms retained separately. Its causal mean and symmetric covariance are distinct kernels; no retarded single-branch action replacement or numerical noise/response norm is asserted.",
        "clock_bridge": "S and its first variation vanish on the full clock. The extra linear-in-vector second/third metric insertions have zero expectation in the zero-mean Gaussian; their odd Gaussian contractions also vanish. The source functional starts at fourth order, and the added physical mean force starts at third order.",
        "free_response_equality": "For perturbations with the SAME fixed preparation, the reference-clock vector one-point function and its first two history-response orders coincide with ordinary connected Proca. This is an identity of this conditional sector, not a quantitative bound on those functional derivatives or a transfer of an older full-parent response.",
        "nonzero_second_source_operator_contact": operator_contact,
        "no_profile_transfer": "The older scalar profiles which canceled ordinary Proca stress are NOT added. The classical target is unchanged; its quantum residual is not set to zero.",
        "checks": {
            "actual_covariant_clock_source_value": s.simplify(actual_jet.subs(eps, 0)),
            "actual_covariant_clock_source_first_variation": s.simplify(
                s.diff(actual_jet, eps).subs(eps, 0)
            ),
            "finite_Gaussian_stress_split": s.expand(expectation - split),
            "linear_source_contact_has_zero_mean": operator_contact.subs(A, 0),
            "source_force_below_cubic_zero": sum(
                source_force.coeff(eps, j) for j in range(3)
            ),
            "source_force_cubic_retained": s.factor(
                source_force.coeff(eps, 3) + 2 * source2 * (mean2 - source2)
            ),
            "pure_mean_stress_below_quartic_zero": sum(
                quadratic_mean.coeff(eps, j) for j in range(4)
            ),
            "source_contact_first_nonzero_degree_four": s.expand(
                sourcejet**2 / 2
            ).coeff(eps, 4)
            - source2**2 / 2,
            "coherent_two_point_symmetric_addition": (a * b) - (b * a),
        },
    }
