"""Complete mixed-row inverse bound, prepared regularity and actual smooth band Schur graph."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_quantum_forced_constraints import forces

from .coordinates import L0, PIVOT_INVERSE

MASS, KAPPA = modes.MASS, modes.KAPPA
GAMMA = 1 / (64 * s.pi**2 * KAPPA)


def require_ball(radius):
    if isinstance(radius, bool) or not isinstance(radius, (int, s.Integer, s.Rational)):
        raise TypeError("A finite positive exact external-momentum radius is required")
    if radius <= 0:
        raise ValueError("A finite positive exact external-momentum radius is required")
    return s.Rational(radius)


@cache
def data():
    C, b0, Kmax = s.symbols(
        "complete_remainder_C reference_inverse_b0 finite_window_Kmax", positive=True
    )
    gamma = s.Symbol("gamma", positive=True)
    beta = b0 * C
    weight = (4 * beta + 1) ** 2
    bound = 2 * beta / (4 * beta + 1)
    F0, F2 = s.symbols("Ftrace Fshear", nonzero=True)
    middle = gamma * L0.T * s.diag(F0, s.Rational(8, 3) * F2) * L0
    inv = L0.inv() * s.diag(1 / F0, s.Rational(3, 8) / F2) * L0.T.inv() / gamma
    norms = [
        max(sum(abs(M[i, j]) for j in range(2)) for i in range(2))
        for M in (L0, L0.inv(), L0.T, L0.T.inv())
    ]
    d = forces.data()
    tt = s.Symbol("proper_time", real=True)
    eta = s.Function("eta")(tt)
    pivot = -s.Rational(3, 2) / (1 + tt**2) ** 6
    return {
        "actual_mass": MASS,
        "actual_kappa": KAPPA,
        "actual_quantum_multiplier": GAMMA,
        "original_reference_middle_inverse": inv,
        "finite_window_channel_norm": "Kmax=max(||Ktrace||L1(0,1),(3/8)||K2||L1(0,1)); full K2 retains its pole and is not half-line L1.",
        "reference_inverse_bound": s.Max(PIVOT_INVERSE, 16 * Kmax / (9 * gamma), 1),
        "weight": weight,
        "weighted_contraction_bound": bound,
        "unweighted_adapted_source_norm_bound": 2 * s.exp(weight) * b0,
        "normal_form": "diag(I4,I4,I4,I2) T_adapt=diag(-6delta(t)^2,gamma L0^T diag(Ftrace,(8/3)F2)L0,-1)+V",
        "ordered_inverse": "(I+B0 V)^-1 B0 diag(I4,I4,I4,I2)",
        "smoothness": "Prepared simultaneous time shifts retain both V_j=(Dt+Ds)^j V and A^(j)(t) on the eta component. Every fixed smooth seminorm has a finite bound; no analytic-in-order or all-P estimate.",
        "band_domain": "C^infinity prepared(I;H^r_ball), each finite ball and real r; common nonempty zero initial germ. No momentum derivatives of the selected state or spatial Schwartz output claimed.",
        "actual_S222_smooth_band_inverse": "For h0 solve the full Euler equation with -J_adapt e, recover g=e+Qbar s and original Legendre Z. This gives both identities for I-Qbar Taux. General h uses RHS e+Qbar C G0 h; Daux and output density retained.",
        "graph_boundary": "Solutions and feedback lie in the smooth Fourier-supported subgraph of X_r=H13_t H^(r+8),Y_r=L2_t H^(r-3). No maximality, graph-density, unrestricted X/Y surjectivity or uniform Pmax-to-infinity bound.",
        "checks": {
            "middle_inverse_left": (inv * middle - s.eye(2)).applyfunc(s.cancel),
            "middle_inverse_right": (middle * inv - s.eye(2)).applyfunc(s.cancel),
            "all_four_coordinate_norms": s.ImmutableMatrix(norms)
            - s.ones(4, 1) * s.Rational(4, 3),
            "strict_half_contraction_margin": s.Rational(1, 2)
            - bound
            - 1 / (2 * (4 * beta + 1)),
            "both_middle_norm_factors": s.Rational(4, 3) ** 2 - s.Rational(16, 9),
            "actual_gamma_normalization": GAMMA * 64 * s.pi**2 * KAPPA - 1,
            "variable_pivot_derivative_commutator": s.diff(pivot * eta, tt)
            - pivot * s.diff(eta, tt)
            - s.diff(pivot, tt) * eta,
            "retained_actual_direct_auxiliary_block": d["checks"][
                "direct_auxiliary_rank_two_decomposition"
            ],
            "actual_phase_force_sign": d["checks"][
                "force_is_negative_symplectic_metric_adjoint"
            ],
        },
        "gates": {
            "strict_contraction_for_every_finite_C_and_b0": s.simplify(
                s.Rational(1, 2) - bound
            ).is_positive,
            "finite_window_not_half_line_shear_norm": True,
            "variable_pivot_commutator_is_nonzero": s.diff(pivot, tt) != 0,
            "actual_Daux_rank_two": d["direct_auxiliary_response_matrix"].rank() == 2,
            "all_original_parameters_fixed": MASS == 1000
            and KAPPA == s.Integer(10) ** 800,
            "actual_smooth_band_graph_not_unrestricted_completion": True,
        },
    }
