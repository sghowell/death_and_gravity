"""Residual translations and the unchanged pure Gaussian mode prescription."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import phase
from p8_vacuum_affine_heavy_state_response import response as heavy
from p8_vacuum_affine_matrix_adiabatic import frame as proca
from p8_vacuum_affine_spatial_gauge import gauge


def rotation(n):
    return s.zeros(n).row_join(s.eye(n)).col_join((-s.eye(n)).row_join(s.zeros(n)))


def gaussian_graph(n):
    av = s.symbols(
        "symmetric_A_" + str(n) + "_0:" + str(n * (n + 1) // 2), complex=True
    )
    bv = s.symbols(
        "antisymmetric_B_" + str(n) + "_0:" + str(n * (n - 1) // 2), complex=True
    )
    A = s.zeros(n)
    B = s.zeros(n)
    count = 0
    for i in range(n):
        for j in range(i, n):
            A[i, j] = A[j, i] = av[count]
            count += 1
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            B[i, j] = bv[count]
            B[j, i] = -bv[count]
            count += 1
    return A.row_join(B).col_join((-B).row_join(A))


@cache
def data():
    checks = {}
    graphs = []
    for n in (1, 2, 3, 8):
        Z = gaussian_graph(n)
        S = rotation(n)
        q = s.Matrix(s.symbols("real_mode_" + str(n) + "_0:" + str(2 * n), real=True))
        checks["whole_gaussian_graph_symmetry_" + str(n)] = Z - Z.T
        checks["whole_gaussian_translation_covariance_" + str(n)] = Z * S - S * Z
        checks["whole_gaussian_translation_generator_annihilation_" + str(n)] = (
            s.expand((S * q).dot(Z * q))
        )
        checks["whole_configuration_rotation_zero_divergence_" + str(n)] = s.trace(S)
        graphs.append(
            {
                "physical_channel_count": n,
                "whole_allowed_graph": Z,
                "configuration_rotation": S,
            }
        )
    a, b, pc, ps = s.symbols(
        "real_cos real_sin canonical_cos_p canonical_sin_p", real=True
    )
    x, k = s.symbols("periodic_angle nonzero_wave_number", real=True)
    q = a * s.cos(x) + b * s.sin(x)
    p = 2 * pc * s.cos(x) + 2 * ps * s.sin(x)
    charge = k * (b * pc - a * ps)
    checks["whole_real_Fourier_pair_normalization"] = s.expand(
        s.integrate(p * k * s.diff(q, x), (x, 0, 2 * s.pi)) / (2 * s.pi) - charge
    )
    # Three generators with non-collinear modes, not a one-axis angular example.
    qs = s.Matrix(s.symbols("canonical_pair_q0:6", real=True))
    pp = s.Matrix(s.symbols("canonical_pair_p0:6", real=True))
    vectors = [(1, 0, 0), (0, 1, 0), (1, -1, 2)]
    Ji = s.Matrix(
        [
            sum(
                vectors[m][i] * (qs[2 * m + 1] * pp[2 * m] - qs[2 * m] * pp[2 * m + 1])
                for m in range(3)
            )
            for i in range(3)
        ]
    )
    bracket = s.Matrix(
        3,
        3,
        lambda i, j: sum(
            s.diff(Ji[i], qs[r]) * s.diff(Ji[j], pp[r])
            - s.diff(Ji[i], pp[r]) * s.diff(Ji[j], qs[r])
            for r in range(6)
        ),
    )
    checks["whole_three_residual_translation_Poisson_brackets"] = bracket
    phase_data = phase.data()
    maps = {
        "physical_kappa_density_map": phase.C,
        "central_canonical_swap": phase.central_map(),
    }
    for name in ("outer", "central"):
        maps[name + "_complete_symmetric_boundary"] = phase.chart_phase(name)["shear"]
    for name, T in maps.items():
        doubled = s.diag(T, T)
        S = rotation(T.rows)
        checks["whole_original_map_translation_equivariance_" + name] = (
            doubled * S - S * doubled
        )
    old_heavy = heavy.data()
    old_proca = proca.data()
    checks["whole_fixed_heavy_Bogoliubov_purity"] = old_heavy["checks"][
        "full_fixed_Bogoliubov_CCR"
    ]
    checks["whole_original_Proca_complete_covariance_purity"] = old_proca["checks"][
        "exact_noncommuting_graph_covariance_purity"
    ].applyfunc(s.cancel)
    mixed = s.diag(s.Rational(1, 2), s.Rational(1, 2))
    Jmixed = s.diag(1, -1)
    checks["whole_mixed_state_invariance_is_insufficient"] = (
        mixed * Jmixed - Jmixed * mixed
    )
    checks["whole_mixed_state_zero_charge_mean"] = s.trace(mixed * Jmixed)
    checks["whole_mixed_state_nonzero_charge_variance"] = s.trace(mixed * Jmixed**2) - 1
    oldzero = gauge.zero_modes()
    checks["mean_zero_gauge_bracket_retains_translation"] = oldzero["checks"][
        "whole_nonzero_translation_bracket"
    ]
    return {
        "whole_residual_translation_constraint": "J_i=integral[Pi_v partial_i v+Pi_tau:partial_i tau+Pi_W dot partial_i W+Pi_M partial_i M1+Pi_H partial_i H]. All three J_i=0 remain after solving the nonzero spatial gauge/constraint pairs.",
        "whole_real_pair_charge": charge,
        "whole_independent_three_direction_diagnostic": Ji,
        "whole_general_translation_invariant_pure_Gaussian_graphs": graphs,
        "whole_original_scalar_canonical_maps": maps,
        "whole_original_fixed_sampling": {
            "initial_time": phase_data["initial_time"],
            "support": phase_data["fixed_sampling_support"],
            "amplitude": phase_data["fixed_sampling_amplitude"],
        },
        "whole_original_H_preparation": old_heavy["unchanged_preparation"],
        "whole_original_Proca_preparation": old_proca["actual_pure_graph"],
        "whole_pure_Gaussian_kernel_argument": "For every finite nonzero real Fourier pair, S=[0,I;-I,0]. A homogeneous pure Gaussian has graph Z=Z^T, Re Z>0 and ZS=SZ, equivalently Z=[A,B;-B,A], A^T=A, B^T=-B. Its Schwartz wavefunction exp(-q^T Z q/2) is annihilated by J=-i hbar(Sq).partial_q because (Sq)^T Zq=0. Odd internal cross blocks are not deleted. The same pure reference prescriptions and full canonical maps preserve this symmetry.",
        "whole_group_and_operator_domain": "The coordinate rotation acts by ordinary pullback on finite-dimensional L2. It is a periodic strongly continuous unitary representation of the spatial torus; its three strongly commuting self-adjoint generators are the closures of the displayed complete first-order operators. The Gaussian is a Schwartz common zero-charge vector. Group averaging therefore leaves this nonzero-mode factor unchanged, without conditioning or renormalizing it.",
        "whole_regulator_and_zero_mode_boundary": "A symmetric finite torus-mode evaluation of the original continuum mode prescriptions is an explicit IR/UV regulator family, not literally the original R3 state. Only nonzero modes are assigned this Gaussian factor. Homogeneous tracefree shape has five modes and is not assigned a new state; translations act trivially on any separate homogeneous factor. The nonlinear shape chart exists only on its stated classical neighborhood; an unbounded Gaussian is not claimed supported there or nonlinearly transported there.",
        "whole_quantum_dynamic_boundary": "For a specified translation-invariant finite Weyl Hamiltonian on a suitable common domain, quadratic J has no higher Moyal terms, so the commutator is i hbar times the quantized Poisson bracket. This algebraic identity alone does not construct that full Hamiltonian, ordering, domain, evolution, continuum regulator or interacting mean.",
        "checks": checks,
        "gates": {
            "all_eight_physical_nonzero_channels_included": graphs[-1][
                "physical_channel_count"
            ]
            == 8,
            "odd_internal_Gaussian_cross_blocks_allowed": gaussian_graph(2)[:2, 2:]
            != s.zeros(2),
            "same_fixed_reference_initial_time": phase_data["initial_time"]
            == -s.Rational(1, 2),
            "same_fixed_reference_sampling_support": phase_data[
                "fixed_sampling_support"
            ]
            == (-s.Rational(15, 32), -s.Rational(13, 32)),
            "same_full_heavy_and_Proca_reference": all(old_heavy["gates"].values())
            and all(old_proca["gates"].values()),
            "mixed_state_counterexample_has_positive_charge_variance": s.trace(
                mixed * Jmixed**2
            )
            > 0,
            "mean_zero_descriptors_not_a_Lie_subalgebra": oldzero["gates"][
                "mean_zero_complement_not_a_Lie_subalgebra"
            ],
            "translation_group_not_phase_space_oscillator_rotation": True,
            "no_homogeneous_state_selection_or_spurious_two_TT_zero_modes": True,
            "no_nonlinear_Gaussian_support_or_interacting_evolution_claim": True,
        },
    }
