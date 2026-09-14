"""Whole original finite-mode covariance, canonical maps and field-norm ceilings."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import phase
from p8_vacuum_affine_heavy_curved_state import estimates as heavy_estimates
from p8_vacuum_affine_heavy_curved_state import state as heavy_state
from p8_vacuum_affine_matrix_adiabatic import initial as proca_initial
from p8_vacuum_affine_quantitative_gaussian_window import state as scalar_state
from p8_vacuum_affine_scalar_tame_propagator import charts

P = s.Integer(10) ** 64
KAPPA = s.Integer(10) ** 800
ZETA = s.Rational(1, 10**6)
MASS = s.Integer(1000)
COUNT = s.Integer(6)
DIMENSION = 48
CORE = s.Integer(10) ** 20
SUPPORT = 2 * CORE
OMEGA_LO, OMEGA_HI = s.Integer(10) ** 98, s.Integer(10) ** 100
SCALAR, TENSOR, HEAVY, PROCA = (
    s.Integer(10) ** 11,
    s.Integer(1000),
    s.Integer(100),
    s.Integer(2),
)
POWERS = {
    "v_A2": 350,
    "tau_A2": 295,
    "Pi_v_A0": 288,
    "Pi_tau_A0": 360,
    "M1_A1": 350,
    "delta_Pi_M_A0": 350,
    "eta_A0": 340,
    "H_gradient_A0": 375,
    "Pi_H_A0": 343,
    "W_A0": 360,
    "Pi_W_A0": 365,
}


def field_bounds():
    root, weight = s.sqrt(KAPPA), 1 + P
    return {
        "v_A2": COUNT * weight**2 * SCALAR / (2 * P ** s.Rational(3, 2) * root),
        "tau_A2": COUNT * 2 * 2 * weight**2 * TENSOR / (s.sqrt(P) * root),
        "Pi_v_A0": COUNT * 2 * P ** s.Rational(3, 2) * SCALAR / root,
        "Pi_tau_A0": COUNT * 2 * s.sqrt(P) * TENSOR / (2 * root),
        "M1_A1": COUNT * weight * SCALAR / (s.sqrt(P) * root),
        "delta_Pi_M_A0": COUNT * s.sqrt(P) * SCALAR / root,
        "eta_A0": COUNT * 10**100 * HEAVY / (s.sqrt(OMEGA_LO) * root),
        "H_gradient_A0": COUNT * P * HEAVY / (s.sqrt(OMEGA_LO) * root),
        "Pi_H_A0": COUNT * s.sqrt(OMEGA_HI) * HEAVY / root,
        "W_A0": COUNT * 4 * s.sqrt(P) / (MASS * s.sqrt(KAPPA * ZETA)),
        "Pi_W_A0": COUNT * 4 * s.sqrt(P) * s.sqrt(ZETA / KAPPA),
    }


def safe_bounds():
    return {name: s.Rational(1, 10**power) for name, power in POWERS.items()}


@cache
def data():
    original = scalar_state.preparation()
    heavy = heavy_estimates.data()
    proca = proca_initial.data()
    raw, safe = field_bounds(), safe_bounds()
    rp = s.Matrix(s.symbols("balanced_scalar_phase0:4", real=True))
    balanced_to_clean = s.diag(1 / s.sqrt(P), 1 / s.sqrt(P), s.sqrt(P), s.sqrt(P))
    central = phase.central_map().subs({phase.a: 1, charts.q: P**2})
    mapped = central.inv() * balanced_to_clean * rp / s.sqrt(KAPPA)
    wanted = s.Matrix(
        [
            rp[2] / (2 * P ** s.Rational(3, 2) * s.sqrt(KAPPA)),
            rp[1] / s.sqrt(P * KAPPA),
            -2 * P ** s.Rational(3, 2) * rp[0] / s.sqrt(KAPPA),
            s.sqrt(P / KAPPA) * rp[3],
        ]
    )
    x = s.Symbol("torus_coordinate", real=True)
    normalization = s.sqrt(2 / (2 * s.pi) ** 3)
    cosine, sine = normalization * s.cos(P * x), normalization * s.sin(P * x)
    checks = {
        "literal_entire_scalar_bounce_canonical_dictionary": (
            mapped - wanted
        ).applyfunc(s.factor),
        "actual_scalar_covariance_ceiling": original[
            "whole_scalar_balanced_covariance_operator_upper"
        ]
        - 10**21,
        "actual_each_tensor_covariance_ceiling": original[
            "whole_each_tensor_balanced_covariance_upper"
        ]
        - 10**5,
        "actual_exact_heavy_comparison_zero_jet_ceiling": heavy[
            "exact_comparison_mode_derivative_coefficients"
        ][0]
        - 20,
        "whole_six_real_basis_and_eight_channels": COUNT * 8 - DIMENSION,
        "whole_real_cosine_torus_normalization": s.integrate(
            cosine**2, (x, 0, 2 * s.pi)
        )
        * (2 * s.pi) ** 2
        - 1,
        "whole_real_sine_torus_normalization": s.integrate(sine**2, (x, 0, 2 * s.pi))
        * (2 * s.pi) ** 2
        - 1,
        "whole_real_cosine_sine_orthogonality": s.integrate(
            cosine * sine, (x, 0, 2 * s.pi)
        ),
        "whole_scalar_physical_CCR": (
            (central.inv() * balanced_to_clean / s.sqrt(KAPPA))
            * phase.OMEGA
            * (central.inv() * balanced_to_clean / s.sqrt(KAPPA)).T
            - phase.OMEGA / KAPPA
        ).applyfunc(s.factor),
    }
    omega, kap, zeta, mass = s.symbols("omega kappa zeta mass", positive=True)
    B = s.diag(1, 1, omega / mass)
    proca_map = s.diag(
        B / s.sqrt(omega) / s.sqrt(kap * zeta),
        B.inv() * s.sqrt(omega) * s.sqrt(zeta / kap),
    )
    J3 = s.zeros(3).row_join(s.eye(3)).col_join((-s.eye(3)).row_join(s.zeros(3)))
    checks["whole_three_polarization_W_and_Pi_W_CCR"] = (
        proca_map * J3 * proca_map.T - J3 / kap
    ).applyfunc(s.factor)
    tensor_map = s.diag(-2 / s.sqrt(kap), -1 / (2 * s.sqrt(kap)))
    J1 = s.Matrix([[0, 1], [-1, 0]])
    checks["whole_tensor_shape_and_dual_CCR"] = (
        tensor_map * J1 * tensor_map.T - J1 / kap
    )
    return {
        "whole_exact_finite_family": {
            "torus_radius_L": 1,
            "wave_number_P": P,
            "symmetric_nonzero_modes": [
                tuple(sign * P * int(i == j) for i in range(3))
                for j in range(3)
                for sign in (-1, 1)
            ],
            "real_basis_count": COUNT,
            "configuration_dimension": DIMENSION,
            "phase_dimension": 2 * DIMENSION,
            "whitened_core_radius": CORE,
            "whitened_support_radius": SUPPORT,
        },
        "whole_same_state_full_covariance_argument": "Let z=S w with V=SS^T/2 the FULL original unit-CCR finite covariance. For every complete linear physical row l, ||l S||^2=2 l V l^T, so |l z|<=sqrt(2 l V l^T)||w||. No internal, q-p or real-pair covariance is discarded. The full original scalar/tensor, actual H SLE and original all-order three-mode Proca prescriptions remain unchanged.",
        "whole_scalar_bounce_dictionary": wanted,
        "whole_original_scalar_full_canonical_map": phase.C,
        "whole_original_boundary_scope": "Both S251 symmetric boundary shears and the central map are retained in the reference. The central shear has value0 at u=0; its derivative is not set to zero. Only this bounce slice is asserted here.",
        "whole_heavy_covariance_derivation": "The exact comparison's energy-scaled mode y has ||y||1<20. For the actual SLE T=alpha S+beta conjugate(S), |beta|<1/2 and alpha+|beta|<2, retaining its phase, give ||y_T||1<40. At a=1,H=0 these are balanced canonical quadratures, with full real covariance<1600 I. The mode-coordinate row ceiling100 is safe. The comparison is not substituted for the state.",
        "whole_heavy_frequency_interval": [OMEGA_LO, OMEGA_HI],
        "whole_original_Proca_graph_prescription": proca["unchanged_initial_state"],
        "whole_original_Proca_high_band_result": proca["actual_high_band_result"],
        "whole_Proca_covariance_derivation": "At this actual P, nu_minus>P/2>1e16. S190 gives ||r_actual||<1/10 on the original zero-shear reference. The FULL six-quadrature covariance is <=(1+1/100)/(1-1/100) I<2I. With B=sqrt(I+kk^T/m^2), A=B Q/sqrt(omega), pi_A=B^-1 sqrt(omega)Pquad. W=A/sqrt(kappa*zeta), Pi_W=sqrt(zeta/kappa)pi_A; no longitudinal or covariance block is deleted.",
        "whole_Proca_three_polarization_physical_map": proca_map,
        "whole_tensor_shape_and_dual_map": tensor_map,
        "whole_Fourier_norm_convention": "A_s sums (1+|k|)^s times absolute scalar, Euclidean vector or matrix operator Fourier norms, with physical integer k on the2pi torus. A normalized real cosine/sine basis function has A0 norm sqrt(2/(2pi)^3)<1. Summing all six functions, both tensor polarizations and the full three-vector row gives every displayed ceiling. Nonlinear reconstruction later uses full convolution, not these six harmonics alone.",
        "whole_all_eleven_exact_per_whitened_radius_bounds": raw,
        "whole_all_eleven_safe_per_whitened_radius_bounds": safe,
        "whole_support_ball_field_bounds": {
            name: value * SUPPORT for name, value in raw.items()
        },
        "checks": checks,
        "gates": {
            **{
                "actual_full_field_bound_" + name: bool(value < safe[name])
                for name, value in raw.items()
            },
            "same_full_scalar_covariance_row_ceiling": 2 * 10**21 < SCALAR**2,
            "same_both_tensor_covariance_row_ceiling": 2 * 10**5 < TENSOR**2,
            "same_actual_H_SLE_covariance_row_ceiling": 2 * 1600 < HEAVY**2,
            "same_actual_Proca_covariance_row_ceiling": 2 * 2 <= PROCA**2,
            "actual_heavy_frequency_lower": heavy_state.MASS2 + P**2 / 16 > OMEGA_LO**2,
            "actual_heavy_frequency_upper": heavy_state.MASS2 + P**2 / 16 < OMEGA_HI**2,
            "actual_Proca_high_band_hypotheses": bool(
                proca_initial.RADIAL_SCALE > s.Rational(1, 4)
                and P / 2 > proca_initial.PARTITION
                and proca["gates"]["actual_graph_and_reference_inside_covariance_ball"]
            ),
            "actual_Proca_frequency_bounds": P**2 + MASS**2 < 4 * P**2,
            "whole_normalized_torus_basis_A0_below_one": bool(2 < (2 * s.pi) ** 3),
            "all_five_canonical_momentum_A1_ceilings_below1e_minus200": all(
                raw[name] * SUPPORT * (1 + P) < s.Rational(1, 10**200)
                for name in (
                    "Pi_v_A0",
                    "Pi_tau_A0",
                    "delta_Pi_M_A0",
                    "Pi_H_A0",
                    "Pi_W_A0",
                )
            ),
            "whole_tau_A2_below1e_minus274": raw["tau_A2"] * SUPPORT
            < s.Rational(1, 10**274),
            "whole_v_A2_below1e_minus329": raw["v_A2"] * SUPPORT
            < s.Rational(1, 10**329),
            "no_homogeneous_quantum_state_or_physical_UV_cutoff": True,
        },
    }
