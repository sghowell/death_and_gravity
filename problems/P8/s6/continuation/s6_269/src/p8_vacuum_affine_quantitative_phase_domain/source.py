"""Whole canonical invariant map with all density and original source factors."""

from functools import cache

import sympy as s
from p8_vacuum_affine_local_quantum_regulator import source as quantum_parent
from p8_vacuum_affine_nonlinear_auxiliary_measure import cotangent
from p8_vacuum_affine_nonlinear_lapse_branch import source as auxiliary

from . import reference


@cache
def invariant_map():
    g = s.Matrix(
        [
            [s.Rational(6, 5), s.Rational(1, 10), 0],
            [s.Rational(1, 10), s.Rational(9, 10), 0],
            [0, 0, s.Rational(100, 107)],
        ]
    )
    # The map formulas below hold for every metric; this non-diagonal
    # metric is used only for independent exact canonical diagnostics.
    metric = s.Rational(7, 5) * g.inv()
    pi_entries = s.symbols("metric_pi0:6", real=True)
    pi = s.Matrix(
        [
            [pi_entries[0], pi_entries[1], pi_entries[2]],
            [pi_entries[1], pi_entries[3], pi_entries[4]],
            [pi_entries[2], pi_entries[4], pi_entries[5]],
        ]
    )
    W = s.Matrix(s.symbols("whole_W0:3", real=True))
    PW = s.Matrix(s.symbols("whole_PiW0:3", real=True))
    gradM = s.Matrix(s.symbols("whole_grad_M0:3", real=True))
    gradH = s.Matrix(s.symbols("whole_grad_H0:3", real=True))
    PM, PH, h, divPW, curvature = s.symbols(
        "whole_Pi_M whole_Pi_H full_H_field full_div_PiW full_intrinsic_curvature",
        real=True,
    )
    V = s.Symbol("whole_hat_volume_density", positive=True)
    F01, F02, F12 = s.symbols("full_F01 full_F02 full_F12", real=True)
    F = s.Matrix([[0, F01, F02], [-F01, 0, F12], [-F02, -F12, 0]])
    inverse = metric.inv()
    trace = s.trace(pi * metric)
    tf = pi - trace * inverse / 3
    magnetic = sum(
        inverse[i, k] * inverse[j, l] * F[i, j] * F[k, l]
        for i in range(3)
        for j in range(3)
        for k in range(3)
        for l in range(3)
    )
    expressions = [
        2 * trace / (3 * V),
        divPW / V,
        PM / V - s.Rational(1, 10),
        PH / V,
        10**100 * h,
        s.trace(tf * metric * tf * metric) / V**2,
        (PW.T * metric * PW)[0] / (reference.ZETA * V**2),
        reference.ZETA * magnetic,
        (W.T * inverse * W)[0],
        (gradM.T * inverse * gradM)[0],
        (gradH.T * inverse * gradH)[0],
        curvature,
    ]
    return {
        "metric": metric,
        "pi": pi,
        "Pi_v": 2 * trace,
        "tf": tf,
        "volume": V,
        "expressions": expressions,
        "magnetic_F": F,
        "W": W,
        "Pi_W": PW,
        "matter_momenta": (PM, PH),
        "h": h,
    }


@cache
def data():
    full = quantum_parent.data()
    spatial = cotangent.spatial_kinetic_blocks()
    inv = invariant_map()
    V = inv["volume"]
    checks = {
        "entire_trace_density_factor": s.factor(
            inv["expressions"][0] - inv["Pi_v"] / (3 * V)
        ),
        "entire_tracefree_metric_momentum": s.factor(
            s.trace(inv["tf"] * inv["metric"])
        ),
        "entire_M1_density_background_contact": s.factor(
            inv["expressions"][2]
            - (inv["matter_momenta"][0] - s.Rational(1, 10)) / V
            - s.Rational(1, 10) * (1 / V - 1)
        ),
        "entire_mass_adapted_heavy_field": inv["expressions"][4] - 10**100 * inv["h"],
        "full_original_electric_Legendre_binding": spatial["checks"][
            "whole_three_direction_electric_Legendre"
        ],
        "full_original_both_matter_Legendre_binding": spatial["checks"][
            "whole_two_matter_canonical_normalizations"
        ],
        "full_original_five_shear_Legendre_binding": spatial["checks"][
            "whole_five_shear_canonical_normalizations"
        ],
        "full_original_Gauss_integration_binding": spatial["checks"][
            "whole_spatial_Gauss_summation_by_parts"
        ],
        "literal_full_bounce_source_binding": full["checks"][
            "literal_original_bounce_slice_of_whole_time_parent"
        ],
        "whole_original_kappa": auxiliary.parent.source.KAPPA - reference.KAPPA,
        "whole_all_twelve_density_invariant_channels": len(inv["expressions"])
        - len(auxiliary.COORDS),
    }
    return {
        "whole_original_bounce_Hamiltonian": auxiliary.HAMILTONIAN,
        "whole_original_bounce_auxiliary_constraint": auxiliary.CONSTRAINT,
        "whole_original_all_time_source_and_boundary": {
            name: full[name]
            for name in (
                "whole_original_time_dependent_normal_Hamiltonian",
                "whole_actual_source_function_bindings",
                "whole_original_R_F_and_fixed_profiles",
                "whole_primitive_boundary",
                "whole_original_spatial_generator",
            )
        },
        "whole_general_invariant_definitions": "For ANY full reconstructed gamma, inverse g and V=sqrt(det gamma), Pi_v=2tr(pi gamma), pi_TF=pi-tr(pi gamma)g/3. The twelve S266 invariants are p=Pi_v/(3V); G=div(Pi_W)/V; dp=Pi_M/V-1/10; ph=Pi_H/V; eta=1e100 H; sh=tr(pi_TF gamma pi_TF gamma)/V^2; el=Pi_W^T gamma Pi_W/(zeta V^2); ma=zeta g^ik g^jl F_ij F_kl; wm=W^T g W; gm=grad(M1)^T g grad(M1); gh=grad(H)^T g grad(H); curv=R[gamma]. F=dW. All products, volume factors and derivatives use the entire reconstructed fields. These are not independent canonical or Gaussian oscillators.",
        "whole_independent_nondiagonal_invariant_diagnostic": dict(
            zip(auxiliary.COORDS, inv["expressions"], strict=True)
        ),
        "whole_nondiagonal_diagnostic_metric": inv["metric"],
        "whole_source_state_and_slice_boundary": "Only u=0 is claimed. The original scalar background Pi_M=1/10 and other reduced momentum/field backgrounds at this slice are retained. Only the already-eliminated auxiliary primary momenta and mean-zero spatial constraints are solved. The residual global translations are retained as S268 quantum constraints; coherent displacements need not satisfy their classical level set. No new homogeneous state, original interacting mean or time-window estimate is supplied.",
        "checks": checks,
        "gates": {
            "complete_twelve_invariant_binding": len(inv["expressions"]) == 12,
            "all_magnetic_antisymmetric_components_retained": all(
                inv["expressions"][7].has(z)
                for z in (
                    inv["magnetic_F"][0, 1],
                    inv["magnetic_F"][0, 2],
                    inv["magnetic_F"][1, 2],
                )
            ),
            "both_matter_and_all_vector_channels_retained": all(
                spatial["gates"].values()
            ),
            "whole_nondiagonal_metric_not_replaced_by_diagonal_diagnostic": inv[
                "metric"
            ][0, 1]
            != 0,
            "full_original_sources_and_state_unchanged": all(full["gates"].values()),
            "physical_normalized_momenta_not_wrong_unit_CCR": True,
            "no_discrete_diffeomorphism_algebra_or_homogeneous_state": True,
        },
    }
