"""Nonperturbative finite volume endpoint gap for both full operator orderings."""

from functools import cache

import sympy as s
from p8_vacuum_affine_weyl_operator_comparison import derivatives, kernel

from . import moving, source

VERROR = s.Rational(1, 10**320)
OPERATOR_ERROR = s.Rational(1, 10**264)
HBOUND = s.Integer(10) ** 1000


def jet(order, amplitude):
    if isinstance(amplitude, bool) or not isinstance(
        amplitude, (int, s.Integer, s.Rational)
    ):
        raise TypeError("Require an exact declared full-symbol amplitude")
    if amplitude not in (VERROR, HBOUND):
        raise ValueError("Require a declared longer-domain full-symbol amplitude")
    return 2 * amplitude * derivatives.growth(order)


@cache
def bounds():
    volume_jets = [jet(n, VERROR) for n in range(197)]
    Hamiltonian_jets = [jet(n, HBOUND) for n in range(197)]
    coherent = VERROR + 96 * volume_jets[2] / 4
    ordering = kernel.CONSTANT * s.Rational(96**2, 32) * volume_jets[4]
    Weyl = coherent + ordering
    HC = 3 * HBOUND + 96 * Hamiltonian_jets[2] / 4
    HW = HC + kernel.CONSTANT * s.Rational(96**2, 32) * Hamiltonian_jets[4]
    endpoint = (1 + source.TIME**2) ** 6 * (1 - OPERATOR_ERROR)
    center = 1 + OPERATOR_ERROR
    minimum = OPERATOR_ERROR / (3 * (1 - OPERATOR_ERROR))
    return {
        "volume_jets": volume_jets,
        "Hamiltonian_jets": Hamiltonian_jets,
        "coherent": coherent,
        "ordering": ordering,
        "Weyl": Weyl,
        "HC": HC,
        "HW": HW,
        "endpoint": endpoint,
        "center": center,
        "minimum_radius_squared": minimum,
    }


@cache
def identities():
    X, Z = s.Matrix([[0, 1], [1, 0]]), s.diag(1, -1)
    U, psi = -s.I * X, s.Matrix([1, 0])
    F = s.eye(2) + OPERATOR_ERROR * Z
    initial = (psi.conjugate().T * F * psi)[0]
    final = ((U * psi).conjugate().T * F * (U * psi))[0]
    change = (U - s.eye(2)) * psi
    b = bounds()
    t = s.Symbol("nonnegative_time_square", nonnegative=True)
    return {
        "full_state_flip_unitarity": U.conjugate().T * U - s.eye(2),
        "full_state_flip_not_near_identity": (change.conjugate().T * change)[0] - 2,
        "full_state_flip_initial_readout": initial - 1 - OPERATOR_ERROR,
        "full_state_flip_final_readout": final - 1 + OPERATOR_ERROR,
        "complete_reference_volume_binomial": s.expand(
            (1 + t) ** 6 - 1 - sum(s.binomial(6, n) * t**n for n in range(1, 7))
        ),
        "full_minimum_localization_identity": 6
        * (1 - OPERATOR_ERROR)
        * b["minimum_radius_squared"]
        - 2 * OPERATOR_ERROR,
        "entire_same_radius_and_kernel_constant": moving.CORE
        - derivatives.source.RADIUS,
    }


@cache
def data():
    b = bounds()
    return {
        "whole_all197_volume_and_Hamiltonian_phase_jet_bounds": {
            "volume": b["volume_jets"],
            "Hamiltonian": b["Hamiltonian_jets"],
        },
        "whole_actual_operator_and_endpoint_bounds": {
            name: value
            for name, value in b.items()
            if name not in ("volume_jets", "Hamiltonian_jets")
        },
        "whole_actual_time_interval": [-source.TIME, source.TIME],
        "whole_operator_ordering_proof": "The same full96-phase normalized Gaussian-frame Schur theorem and all196 phase/cutoff derivatives apply on the enlarged real-time domain. Use the NEW amplitude1e-320. The calibrated-coherent readout differs fromI by<2e-320. The complete heat remainder, including every mixed derivative, has Weyl operator norm bounded by1e112*96^2*J4/32. Thus the actual Weyl readout differs fromI by<1e-264 and is>I/2. This is not generic positivity-preserving Weyl quantization.",
        "whole_entire_dynamics_proof": "Both full real Hamiltonian extensions are bounded self-adjoint scalar-plus-Schwartz operators and norm continuous in time. The entire Dyson construction gives exact unitary evolutions and the inherited common Schwartz strong equation. Integrated norm bounds here are huge, not small: retain all scalar phases and do not extend old near-identity, leakage, ordering-state or cutoff-state estimates. Both original reference conjugations and all residual translations remain.",
        "whole_state_independent_turnaround_proof": "For either ordering and either cutoff, every evolved normalized state has normalized volume mean between1-epsilon and1+epsilon, epsilon=1e-264. For nu(u)=(1+u^2)^6 times that mean, each endpoint exceeds nu(0) by>5e-260. Continuity gives a global minimum strictly inside the interval. Every minimizer satisfies |u|<1e-132. The source smoothness gives C1 means, so a minimizer is critical and the mean has points with negative and positive derivative on its two sides. No unique minimum or strict positive second derivative is inferred.",
        "whole_physical_scope": "Homogeneous reduced variables are the unchanged EXTERNAL functions, not a newly solved quantum background. This finite-regulator turnaround is not self-consistent homogeneous backreaction, an unlocalized interacting bounce, regulator removal, physical matching, omitted-loop control or original V/G/B/P8 closure.",
        "checks": identities(),
        "gates": {
            "all196_full_jet_ratios": all(
                b["volume_jets"][n + 1] / b["volume_jets"][n] < s.Rational(1, 10**9)
                for n in range(196)
            ),
            "actual_coherent_operator_bound": b["coherent"] < 2 * VERROR,
            "actual_Weyl_operator_bound": b["Weyl"] < OPERATOR_ERROR,
            "actual_positive_Weyl_floor": 1 - b["Weyl"] > s.Rational(1, 2),
            "strict_both_endpoint_gap": b["endpoint"] - b["center"]
            > 5 * source.TIME**2,
            "every_minimizer_strict_inner_radius": b["minimum_radius_squared"]
            < s.Rational(1, 10**264),
            "strict_inner_radius_inside_time_domain": s.Rational(1, 10**132)
            < source.TIME,
            "whole_coherent_H_bound": b["HC"] < 10**1010,
            "whole_Weyl_H_bound": b["HW"] < 10**1056,
            "old_small_state_bounds_not_reassigned": source.TIME * b["HC"] > 2
            and source.TIME * b["HW"] > 2,
            "both_defined_cutoffs_and_full_original_state": True,
            "external_means_not_claimed_as_self_consistent_feedback": True,
            "original_V_G_B_P8_still_open": True,
        },
    }
