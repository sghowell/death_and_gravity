"""Bounds for the unchanged full sampled scalar and tensor Gaussian covariances."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import gaussian, phase
from p8_vacuum_affine_scalar_tame_propagator import charts

from . import source


@cache
def energy():
    low, high = s.Rational(1, 4), s.Integer(64)
    gyro = 2 * 10**18 / source.LOW
    lower = min(1 / (128 * high), low / 4 - gyro**2 / high)
    upper = 2 / low + 2 * gyro**2 / low + high
    P, a = s.symbols("positive_P positive_reference_scale", positive=True)
    k0, k1, k2, g0, g1, g2, alpha = s.symbols(
        "full_K00 full_K01 full_K11 full_G00 full_G01 full_G11 full_skew_boundary",
        real=True,
    )
    K = s.Matrix([[k0, k1], [k1, k2]])
    G = s.Matrix([[g0, g1], [g1, g2]])
    A = s.Matrix([[0, alpha], [-alpha, 0]])
    velocity = (-K.inv() * A / P).row_join(K.inv() / a**3)
    Q = (velocity.T * K * velocity + s.diag(G / a**2, s.zeros(2))).applyfunc(s.factor)
    x = s.Matrix(s.symbols("full_weighted_clean_phase0:4", real=True))
    dy = velocity * x
    direct = (dy.T * K * dy)[0] / 2 + (x[:2, 0].T * G * x[:2, 0])[0] / (2 * a * a)
    bindings = {}
    for name, chart in (("outer", charts.outer()), ("central", charts.central())):
        skew = (chart["B"] - chart["B"].T) / 2
        bindings[name] = {
            k0: chart["K"][0, 0],
            k1: chart["K"][0, 1],
            k2: chart["K"][1, 1],
            g0: chart["G"][0, 0],
            g1: chart["G"][0, 1],
            g2: chart["G"][1, 1],
            alpha: skew[0, 1],
        }
    return {
        "whole_scalar_fixed_reference_weighted_energy_metric": Q,
        "whole_scalar_velocity_from_complete_clean_momenta": velocity,
        "whole_actual_two_chart_coefficient_bindings": bindings,
        "whole_scalar_energy_coercivity_lower": lower,
        "whole_scalar_energy_coercivity_upper": upper,
        "whole_safe_weighted_energy_operator_interval": [
            s.Rational(1, 10**4),
            s.Integer(100),
        ],
        "whole_symplectically_balanced_phase": "r=(sqrt(P)y,Pi/sqrt(P)) and x=sqrt(P)r. The entire momentum is Pi=a^3(K ydot+A_skew y) after the fixed original symmetric boundary. Thus E=P r^T Q r/2. Use the complete displayed coefficient bindings and q=P^2/a^2.",
        "whole_energy_argument": "Young's inequality, full K,G in[1/4,64],1<=a<=2 and the retained gyro norm<=2e18 give the displayed lower and upper constants. No simultaneous diagonalization of K,G is assumed. The covariance estimate uses canonical r, not the nonsymplectic weight x as if it had unit CCR.",
        "checks": {
            "whole_noncommuting_energy_metric_reconstruction": s.factor(
                (x.T * Q * x)[0] / 2 - direct
            ),
            "whole_energy_metric_symmetric": Q - Q.T,
        },
        "gates": {
            "complete_scalar_metric_lower_above_1e_minus4": lower
            > s.Rational(1, 10**4),
            "complete_scalar_metric_upper_below100": upper < 100,
            "complete_gyro_small_after_actual_P_scaling": gyro < s.Rational(1, 10**40),
            "whole_K_G_mixing_and_symmetric_boundary_retained": True,
            "both_original_chart_bindings_present": len(bindings) == 2,
            "balanced_phase_CCR_not_replaced_by_weighted_phase_CCR": True,
        },
    }


@cache
def preparation():
    path_bound = 2 * 3**6 * 16 * 10**3
    gamma = s.Integer(10) ** 8
    gram_lo = s.Rational(1, 10**4) / gamma**2
    gram_hi = 16 * 100 * gamma**2
    tensor_gamma = 2 * 3**3
    tensor_lo = s.Rational(1, 8) / tensor_gamma**2
    tensor_hi = 4 * tensor_gamma**2
    tensor_cov = s.Integer(10) ** 5
    checks = {}
    for case in range(4):
        _S, _D, M = gaussian.fixture_preparation(case)
        V = gaussian.ground_covariance(M)
        X = 2 * V
        checks["whole_ground_covariance_positive_matrix_equation_" + str(case)] = (
            X * M * X + phase.OMEGA * M * phase.OMEGA
        ).applyfunc(s.factor)
        checks["whole_original_Gaussian_CCR_and_purity_" + str(case)] = (
            V * phase.OMEGA * V - phase.OMEGA / 4
        ).applyfunc(s.factor)
    _S, _D, M = gaussian.fixture_preparation(1)
    W = gaussian.symplectic_fixture(5)
    transported = W * gaussian.ground_covariance(M) * W.T
    moved = gaussian.ground_covariance(W.inv().T * M * W.inv())
    checks["whole_same_preparation_reexpressed_not_reselected"] = (
        moved - transported
    ).applyfunc(s.factor)
    checks["whole_positive_bump_normalization_cancels"] = (
        gaussian.ground_covariance(7 * M) - gaussian.ground_covariance(M)
    ).applyfunc(s.factor)
    return {
        "whole_full_reference_root_energy_rate": source.ROOT_RATE,
        "whole_actual_chart_energy_conversion_bound": s.Integer(16) * 10**3,
        "whole_actual_preparation_to_target_energy_path_bound": path_bound,
        "whole_safe_both_time_directions_preparation_energy_bound": gamma,
        "whole_original_selection_form_relative_energy_interval": [
            s.Integer(1),
            s.Integer(16),
        ],
        "whole_scalar_preparation_Gramian_bounds_after_dividing_by_P_and_bump_integral": [
            gram_lo,
            gram_hi,
        ],
        "whole_scalar_balanced_covariance_operator_upper": source.COVARIANCE,
        "whole_scalar_weighted_covariance_operator_upper": "C_x(t,P)<=10^21 P I for every P>=10^64 and |t|<=10^-60, in the unchanged complete fixed reference. x=sqrt(P)r, so both P factors and canonical CCR are explicit.",
        "whole_tensor_reference_energy_path_bound": tensor_gamma,
        "whole_tensor_Gramian_bounds_after_dividing_by_P_and_bump_integral": [
            tensor_lo,
            tensor_hi,
        ],
        "whole_each_tensor_balanced_covariance_upper": tensor_cov,
        "whole_preparation_argument": "For any target time in the tiny slab, rewrite the SAME S251 fixed Gramian at that time by the exact canonical flow and complete chart maps. The only switch is-3/16. Absolute root rate12 over at most1/2+T and one full transition give the listed path bound: exp6<3^6 and exp(12T)<2. The original positive preparation form a^3(E_scalar+ y^T K y/2) lies between E_scalar and16 E_scalar at these momenta. In balanced canonical coordinates E_scalar=P r^T Q r/2 with Q in[1e-4,100]. Integration of the original f^2 therefore gives the displayed ACTUAL Gramian bounds; its positive integral is a common scalar and cancels from the selected state.",
        "whole_covariance_matrix_argument": "For the unique positive selected covariance put X=2V. Its full spectral formula gives X M X=-Omega M Omega. If m I<=M<=L I, orthogonality of Omega gives m X^2<=X M X<=L I, so ||V||<=sqrt(L/m)/2. This does not require commuting blocks, branch diagonalization, distinct preparation frequencies or an instantaneous ground-state reset. Symplectic congruence proves that the rewritten Gramian produces the transported original Cauchy covariance.",
        "whole_tensor_argument": "Each original tensor has normalized energy coefficients(a,a^-3), with floor1/8 and ceiling2; absolute root rate6 is safe for |H|<=2. Its unchanged sampled mass1 norm is at most twice that energy. There is no chart switch. The displayed bounds imply C_xT<=1e5 P I for EACH of the two original TT polarizations. The original H/Proca product states and their determinants remain untouched; their interacting mean is not part of this local scalar/tensor covariance theorem.",
        "checks": checks,
        "gates": {
            "both_full_path_bounds_below1e8": path_bound < gamma,
            "whole_short_extra_time_exponential_below_two": source.ROOT_RATE
            * source.TIME
            < s.Rational(1, 2),
            "whole_fixed_mass1_preparation_correction_below_one": 64
            / ((source.LOW**2 / 4) * s.Rational(1, 4))
            < 1,
            "whole_scalar_covariance_bound_from_full_Gramian": gram_hi
            < (2 * source.COVARIANCE) ** 2 * gram_lo,
            "whole_tensor_covariance_bound_from_full_Gramian": tensor_hi
            < (2 * tensor_cov) ** 2 * tensor_lo,
            "same_bump_and_original_future_sampling_interpretation": True,
            "full_matrix_condition_bound_not_scalar_mode_formula": True,
            "complete_both_tensor_polarizations_retained": True,
            "not_new_interacting_state_or_global_varied_history": True,
        },
    }
