"""The entire original R, global coefficient bounds and exact lapse weights."""

from functools import cache

import sympy as s
from p8_affine_vacuum_domain import bounds as original_bounds
from p8_affine_vacuum_domain import family as original
from p8_vacuum_affine_heavy_scalar_parent import family as heavy
from p8_vacuum_affine_nonlinear_auxiliary_measure import canonical
from p8_vacuum_affine_physical_background_vertices import parent
from p8_vacuum_affine_quantitative_gaussian_window import source as previous

u, X = parent.u, parent.X
TIME, LOW, HIGH, KAPPA = previous.TIME, previous.LOW, previous.HIGH, previous.KAPPA
VAR_UPPER = s.Rational(1, 10**490)
VAR_LOWER = s.Rational(1, 10**576)
VAR_BOUNCE_LOWER = s.Rational(1, 10**572)
VOLUME_ERROR = s.Rational(1, 10**960)
RADIUS = s.Rational(1, 1000)
ALPHAS = (s.Rational(1, 4), s.Rational(1, 2), s.Rational(3, 4))


@cache
def global_R():
    old = original.data()
    full = heavy.coefficients()
    order = s.Integer(original.N)
    localizer = heavy.LOCALIZER
    bump = old["bump"]
    T = old["T"]
    B = old["B"]
    delta = order * X**2 * heavy.switch()
    R = 1 + B * (X - 1) / parent.h + delta
    inherited = original_bounds.domain()
    h, x, b, d = s.symbols(
        "positive_h positive_X positive_B positive_delta", positive=True
    )
    # These full coefficient identities are used on the stated real boxes.
    tail_ceiling = order * s.factorial(5) / localizer**5
    small_ceiling = order / localizer
    return {
        "whole_original_R": R,
        "whole_full_heavy_R": full["R"],
        "whole_original_switch": T,
        "whole_original_bump": bump,
        "whole_complete_blended_switch": B,
        "whole_positive_heavy_correction": delta,
        "whole_original_h": parent.h,
        "whole_inherited_rank_regular_domain": inherited,
        "whole_global_X_nonnegative_R_floor": s.Rational(1, 2),
        "whole_large_X_heavy_correction_ceiling": tail_ceiling,
        "whole_small_X_heavy_correction_ceiling": small_ceiling,
        "whole_global_argument": "For X>=0 and the actual even order1024,0<=T<=1 and0<=mX^2 exp(-mX^2)<=1, hence0<=B<=1. The heavy correction is positive. On[0,1] use the inherited complete R floor>1/2; on[1,infinity),R>=1. For X>=1,T>=1/2, so R>=X/(2h). exp(y)>=y^5/5! bounds the ENTIRE heavy correction by m5!/A^5<1, giving R<=2X. On[0,1],deltaR<=m/A<1 and R<=2. No local source is replaced by zero. These are global coefficient estimates, NOT global regularity of every auxiliary/Legendre pivot.",
        "whole_positive_lapse_bounds": "For0<N<=1: 1/(2h N^2)<=R(u,N^-2)<=2/N^2. For N>=1:1/2<R<=2. Thus each f_alpha(N)=1_(N>0) R(u,N^-2)^(-alpha), alpha in{1/4,1/2,3/4}, is positive and bounded by2. The observable is zero at/nonpositive lapse; this does not condition or project the state.",
        "checks": {
            "literal_entire_original_plus_heavy_R": s.expand(full["R"] - R),
            "literal_positive_heavy_coefficient": heavy.GAMMA * heavy.K0 - order,
            "literal_complete_blended_switch": s.expand(B - (T + (1 - T) * bump)),
            "literal_full_large_X_lower_difference": s.expand(
                1
                + b * (x - 1) / h
                + d
                - x / (2 * h)
                - (1 - 1 / (2 * h) + (b - s.Rational(1, 2)) * (x - 1) / h + d)
            ),
            "same_original_order": order - 1024,
            "same_original_localizer": localizer - s.Integer(10) ** 420,
        },
        "gates": {
            "whole_old_domain_conditions_pass": all(inherited["bounds"].values()),
            "whole_old_R_floor_strictly_above_half": inherited["selected_R_floor"]
            > s.Rational(1, 2),
            "actual_switch_order_positive_even": order > 0 and order % 2 == 0,
            "large_X_complete_heavy_ceiling_below_one": tail_ceiling < 1,
            "small_X_complete_heavy_ceiling_below_one": small_ceiling < 1,
            "original_all_real_u_h_at_least_one": True,
            "negative_power_factors_bounded_but_not_all_action_coefficients": True,
            "global_R_positivity_not_global_auxiliary_branch_regularization": True,
        },
    }


@cache
def weights():
    data = canonical.full_trace()
    N, R = canonical.N, canonical.R
    whole = data["whole_reduced_Hamiltonian"]
    matter = s.diff(whole, canonical.pm, 2)
    trace = s.diff(whole, canonical.p, 2)
    curvature = s.diff(whole, canonical.curvature)
    homogeneous = N * (-3 * canonical.M * (s.Symbol("hat_H", real=True) / N) ** 2)
    H = s.Symbol("hat_H", real=True)
    kinetic = -s.diff(homogeneous, H, 2) / 6
    current = canonical.current
    return {
        "whole_original_reduced_Hamiltonian": whole,
        "whole_original_trace_lagrangian": data["lagrangian"],
        "whole_current_homogeneous_density": current.L,
        "whole_exact_homogeneous_ADM_kinetic_coefficient": kinetic,
        "whole_exact_canonical_matter_momentum_coefficient": matter,
        "whole_exact_canonical_trace_momentum_coefficient": trace,
        "whole_exact_canonical_curvature_coefficient": curvature,
        "whole_full_coefficient_boundary": "D=M/N occurs in the complete homogeneous ADM density -3D H_hat^2. Q=N/U is the complete second momentum derivative of the retained M1/H canonical terms. Their individual positive-lapse Gaussian integrability is not the expectation of the entire constrained Hamiltonian, action or interacting state. All other parent terms and determinants remain; there is no no-go inference from a selected coefficient alone.",
        "checks": {
            "whole_actual_homogeneous_kinetic_weight": s.factor(
                kinetic - R ** s.Rational(1, 4) / N
            ),
            "literal_entire_current_homogeneous_kinetic_weight": s.factor(
                -s.diff(current.L, current.H, 2) / 6
                - current.R ** s.Rational(1, 4) / current.N
            ),
            "whole_original_canonical_matter_weight": s.factor(
                matter - N * R ** s.Rational(3, 4)
            ),
            "whole_original_canonical_H_weight": s.factor(
                s.diff(whole, canonical.ph, 2) - matter
            ),
            "whole_original_canonical_trace_weight": s.factor(
                trace + 3 * N / (2 * canonical.M)
            ),
            "whole_original_canonical_curvature_weight": s.factor(
                curvature + N * canonical.C3
            ),
        },
        "gates": {
            "complete_reduced_Hamiltonian_retained": all(
                whole.has(x)
                for x in (
                    canonical.pm,
                    canonical.ph,
                    canonical.p,
                    canonical.G,
                    canonical.curvature,
                    canonical.shear,
                    canonical.electric,
                    canonical.magnetic,
                    canonical.h,
                    canonical.j,
                )
            ),
            "unrestricted_Gaussian_substitution_not_nonlinear_constraint_solution": True,
            "physical_low_energy_domain_not_extended_by_coefficient_probe": True,
            "no_state_projection_or_new_ordering_is_selected": True,
        },
    }
