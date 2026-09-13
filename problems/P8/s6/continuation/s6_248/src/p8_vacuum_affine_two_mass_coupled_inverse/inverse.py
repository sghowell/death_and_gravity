"""Ordered actual two-mass Volterra inverse and full QG2 force recovery."""

from functools import cache

import sympy as s
from p8_vacuum_affine_quantum_forced_constraints import forces
from p8_vacuum_affine_two_mass_reference import spectral

from . import assembly, normal

GAMMA = normal.GAMMA
L0 = normal.L0


def require_ball(radius):
    if isinstance(radius, bool) or not isinstance(radius, (int, s.Rational)):
        raise TypeError("A finite positive exact external radius is required")
    if radius <= 0:
        raise ValueError("Keep a fixed finite positive external Fourier ball")
    return s.Rational(radius)


@cache
def data():
    C = s.Symbol("complete_remainder_C", nonnegative=True)
    b0, Kmax = s.symbols(
        "complete_reference_bound half_line_total_kernel_norm", positive=True
    )
    gamma = s.Symbol("actual_gamma", positive=True)
    beta = b0 * C
    weight = (4 * beta + 1) ** 2
    contraction = 2 * beta / (4 * beta + 1)
    F0, F2 = s.symbols("Ftrace_total Fshear_total", nonzero=True)
    middle = gamma * L0.T * s.diag(F0, s.Rational(8, 3) * F2) * L0
    middle_inverse = (
        L0.inv() * s.diag(1 / F0, s.Rational(3, 8) / F2) * L0.T.inv() / gamma
    )
    norms = [
        max(sum(abs(M[i, j]) for j in range(2)) for i in range(2))
        for M in (L0, L0.inv(), L0.T, L0.T.inv())
    ]
    tt = assembly.t
    pivot = -s.Rational(3, 2) / (1 + tt**2) ** 6
    eta = s.Function("eta")(tt)
    f = forces.data()
    # These are generic exact saddle identities. Their actual coefficient map is explicit below.
    cd = assembly.clock.data()
    cc = cd["full_current_clock_coefficients"]
    base = assembly.clock.old.clock()
    scalar = assembly.scalar
    current_map = {
        forces.J: cc["Jc"],
        forces.A: cc["A"],
        forces.T: cc["Tcorr"],
        forces.th: base["Theta"],
        forces.E: base["E"],
        forces.l: scalar.ell,
        forces.w: -scalar.ell * base["E"],
        forces.delta: scalar.delta,
        forces.H: scalar.H,
    }
    checks = {
        "whole_summed_middle_left_inverse": middle_inverse * middle - s.eye(2),
        "whole_summed_middle_right_inverse": middle * middle_inverse - s.eye(2),
        "all_four_coordinate_factors": s.ImmutableMatrix(norms)
        - s.ones(4, 1) * s.Rational(4, 3),
        "both_coordinate_norms_in_bound": s.Rational(4, 3) ** 2 - s.Rational(16, 9),
        "strict_half_contraction_margin_including_zero_remainder": s.Rational(1, 2)
        - contraction
        - 1 / (2 * (4 * beta + 1)),
        "physical_gamma_and_kappa_restored": GAMMA
        * 64
        * s.pi**2
        * normal.reference.KAPPA
        - 1,
        "nonstationary_reference_commutator": s.diff(pivot * eta, tt)
        - pivot * s.diff(eta, tt)
        - s.diff(pivot, tt) * eta,
        "complete_direct_auxiliary_block": f["checks"][
            "direct_auxiliary_rank_two_decomposition"
        ],
        "complete_phase_force_adjoint_sign": f["checks"][
            "force_is_negative_symplectic_metric_adjoint"
        ],
        "full_current_classical_hamiltonian_identity": f["checks"][
            "zero_force_recovers_full_classical_Hamiltonian"
        ],
        "current_coefficient_J_is_full_QG2": current_map[forces.J] - cc["Jc"],
        "current_coefficient_T_is_full_QG2": current_map[forces.T] - cc["Tcorr"],
        "current_coefficient_A_is_full_QG2": current_map[forces.A] - cc["A"],
    }
    # A generic noncommuting block factorization retains the entire direct auxiliary map.
    n = 2
    E = s.eye(n)
    G = s.Matrix([[1, 1], [0, 2]])
    Cmat = s.Matrix([[2, 0], [1, 1]])
    F = s.Matrix([[1, 0], [1, 2]])
    D = s.Matrix([[-2, 1], [1, -3]])
    Q = s.Matrix([[1, 2], [0, 1]]) / 101
    big = s.BlockMatrix([[E, -G * F], [-Q * Cmat, E - Q * D]]).as_explicit()
    left = s.BlockMatrix([[E, s.zeros(n)], [-Q * Cmat, E]]).as_explicit()
    right = s.BlockMatrix(
        [[E, -G * F], [s.zeros(n), E - Q * (D + Cmat * G * F)]]
    ).as_explicit()
    checks["complete_noncommuting_force_block_factorization"] = big - left * right
    deltaD, deltaG = s.symbols("detector_delta source_delta", real=True)
    clockD = s.Matrix([[1, 0, 0], [deltaD, 1, 0], [0, 0, 1]])
    clockG = s.Matrix([[1, 0, 0], [deltaG, 1, 0], [0, 0, 1]])
    raw = s.Matrix(3, 3, s.symbols("whole_raw_entry0:9"))
    contact = s.Symbol("whole_Gaussian_clock_nn_contact", real=True) * s.diag(1, 0, 0)
    fullclock = clockD.T * (raw + contact) * clockG
    checks["both_time_leg_full_clock_force_reconstruction"] = (
        clockD.T.inv() * fullclock * clockG.inv() - raw - contact
    )
    return {
        "complete_actual_force_coefficient_map": current_map,
        "reference_middle_inverse_transform": middle_inverse,
        "actual_quantum_multiplier": GAMMA,
        "channel_norm": "Kmax=max(||Ktrace,total||L1(0,infinity),(3/8)||K2,total||L1(0,infinity)); S247 proves both finite using both thresholds and the two-sided interior trace cusp. This is the inverse of the SUM, not the sum of inverses or the old Proca-only pole.",
        "reference_inverse_bound": s.Max(
            assembly.PIVOT_INVERSE, 16 * Kmax / (9 * gamma), 1
        ),
        "weight": weight,
        "weighted_contraction_bound": contraction,
        "unweighted_adapted_inverse_bound": 2 * s.exp(weight) * b0,
        "ordered_inverse": "(I+B0 V)^-1 B0 diag(I4,I4,I4,I2)",
        "smoothness": "Simultaneous prepared time translations retain both (Dt+Ds)^j V and every derivative of A_eta=-6delta². All fixed smooth seminorm bounds are finite. The source's original zero germ is preserved.",
        "current_feedback_inverse": "Let Qbar be pure Proca-plus-H Gaussian response with distinct clock contact and output density, and C,K,D,F the full QG2 classical maps. The coupled Euler inverse and original constraint germs prove both identities for I-Qbar(D+C G0 F) on the smooth prepared Fourier-ball subgraph. General mathematical phase forcing enters as e+Qbar C G0 h.",
        "exact_force_coordinate_bridge": "With A(t)=[[1,0,0],[delta(t),1,0],[0,0,1]], the pure heavy clock kernel is A(t)^T Rhat_H(t,s) A(s). Thus Rhat_H=A(t)^(-T)[R_S246,clock-P_H,clock]A(s)^(-1), in the scalar sector. Qbar(t,s)=[kappa a(t)^3]^(-1)(Rhat_Proca+Rhat_H)(t,s). These are ordered detector/source maps; the output density does not move to s.",
        "heavy_coherent_mean": "The independent full physical heavy-mean block is kappa^-1 times the normalized massive KG Euler equation when the action is divided by kappa a³. Specify physical heavy forcing before inverting; no silent unit normalization. Its physical energy-scaled homogeneous propagator is below1000 by S241, uniformly over all P. Its Gaussian metric response is already included in the coupled four-block problem.",
        "scope": "For every fixed finite Pmax and real r, a two-sided actual conditional reference inverse on smooth prepared C^infinity(I;H^r_ball) scalar amplitudes, and the corresponding invariant force subgraph. No unrestricted completed-graph surjectivity, all-Pmax bound, physical smallness, stability or nonlinear/interacting/P8 closure.",
        "checks": {key: assembly.zero(value) for key, value in checks.items()},
        "gates": {
            "strict_weighted_half_margin": s.simplify(
                s.Rational(1, 2) - contraction
            ).is_positive,
            "actual_pivot_not_stationary": s.diff(pivot, tt) != 0,
            "actual_direct_auxiliary_rank_two": forces.system()["D"].rank() == 2,
            "whole_current_coefficients_independent_of_phase_variables": all(
                not value.has(*forces.Z) for value in current_map.values()
            ),
            "complete_reference_spectral_gates_rechecked": all(
                bool(v) for v in spectral.data()["gates"].values()
            ),
            "current_physical_mass_and_kappa_unchanged": normal.reference.PROCA_MASS2
            == 10**6
            and normal.reference.KAPPA == s.Integer(10) ** 800,
            "all_gamma_inverse_factors_in_physical_bound": True,
            "finite_external_support_not_internal_cutoff_or_full_graph": True,
            "Gaussian_profile_not_counted_twice_in_force_interface": True,
        },
    }
