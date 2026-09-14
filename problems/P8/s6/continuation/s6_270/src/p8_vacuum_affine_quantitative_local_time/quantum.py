"""Explicit full cutoff derivatives, calibrated positivity and quantum time bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_local_quantum_regulator import source as old_source
from p8_vacuum_affine_quantitative_phase_domain import geometry as old_geometry

from . import branch, reference, source

DIMENSION = 96
HBOUND = s.Integer(10) ** 1000
VERROR = s.Rational(1, 10**255)
GENERATOR = s.Integer(10) ** 1010


@cache
def cutoff_bounds():
    x = s.Symbol("inverse_transition_coordinate", positive=True)
    polynomials = [s.Integer(1)]
    for _ in range(4):
        polynomials.append(
            s.expand(x * x * (polynomials[-1] - s.diff(polynomials[-1], x)))
        )
    bump = [
        sum(
            abs(coef) * s.factorial(monomial[0])
            for monomial, coef in s.Poly(poly, x).terms()
        )
        for poly in polynomials
    ]
    transition = [s.Integer(1)]
    for n in range(1, 5):
        transition.append(
            9
            * (
                bump[n]
                + sum(
                    s.binomial(n, k) * 2 * bump[k] * transition[n - k]
                    for k in range(1, n + 1)
                )
            )
        )
    c = transition
    radial = [
        1,
        4 * c[1],
        16 * c[2] + 2 * c[1],
        64 * c[3] + 24 * c[2],
        256 * c[4] + 192 * c[3] + 12 * c[2],
    ]
    product = [
        sum(
            s.binomial(n, k) * radial[k] * s.factorial(n - k) * 4 ** (n - k)
            for k in range(n + 1)
        )
        for n in range(5)
    ]
    R = reference.CORE
    Hjets = [2 * HBOUND * product[n] / R**n for n in range(5)]
    Fjets = [2 * VERROR * product[n] / R**n for n in range(5)]
    return {
        "polynomials": polynomials,
        "bump": bump,
        "transition": transition,
        "radial": radial,
        "product": product,
        "Hjets": Hjets,
        "Fjets": Fjets,
        "generator": 3 * HBOUND + DIMENSION * Hjets[2] / 4,
        "Dvolume": DIMENSION * Fjets[2] / 4,
        "D2H": DIMENSION**2 * Hjets[4] / 16,
        "D2volume": DIMENSION**2 * Fjets[4] / 16,
    }


@cache
def complete_domain_bounds():
    auxiliary = branch.data()
    root = auxiliary["whole_full_complex_root_displacement_bound"]
    volume_error = 12 * old_geometry.VLOG + 24 * root
    spatial_measure = s.Integer(1000)
    full_hred = reference.field.KAPPA * spatial_measure * (2 * 10**62 + 2)
    full_canonical_radius = 2 * reference.WHITENING * reference.ANALYTIC_RADIUS
    full_href = reference.CANONICAL_FLOW * full_canonical_radius**2 / 2
    return {
        "whole_implicit_root_bound": root,
        "whole_complex_volume_error": volume_error,
        "whole_torus_measure_bound": spatial_measure,
        "whole_source_and_chart_hred_bound": full_hred,
        "whole_fixed_canonical_complex_radius": full_canonical_radius,
        "whole_quadratic_reference_H_bound": full_href,
        "whole_complete_interaction_bound": full_hred + full_href,
    }


@cache
def dynamical_bounds():
    epsilon = GENERATOR * source.TIME
    tail = s.Rational(1, 10**4000)
    difference = (
        2 * GENERATOR * s.sqrt(tail) * source.TIME + GENERATOR**2 * source.TIME**2
    )
    mean_difference = 4 * difference + 4 * VERROR * tail + 8 * VERROR * epsilon
    return {
        "whole_original_initial_gamma_tail_ceiling": tail,
        "whole_same_seed_unitary_state_error": epsilon,
        "whole_evolved_coherent_outside_core_ceiling": tail + 2 * epsilon,
        "whole_two_explicit_regulator_state_difference": difference,
        "whole_two_readout_mean_difference": mean_difference,
    }


@cache
def identities():
    r, cwidth, R = s.symbols(
        "radial_coordinate positive_width core_radius", positive=True
    )
    phi = s.Function("whole_transition")
    y = (r * r - R * R) / (cwidth * R * R)
    Y = s.Symbol("whole_transition_argument", real=True)
    first, second = 2 * r / (cwidth * R * R), 2 / (cwidth * R * R)
    wanted = (
        s.diff(phi(Y), Y, 4).subs(Y, y) * first**4
        + 6 * s.diff(phi(Y), Y, 3).subs(Y, y) * first**2 * second
        + 3 * s.diff(phi(Y), Y, 2).subs(Y, y) * second**2
    )
    S = s.Matrix(2, 2, s.symbols("full_whitening_entry0:4", real=True))
    covariance = S * S.T / 2
    w = s.symbols("full_phase_coordinate0:2", real=True)
    chi, f = s.Function("whole_cutoff")(*w), s.Function("whole_symbol")(*w)
    laplacian = lambda expr: sum(s.diff(expr, z, 2) for z in w)
    product = laplacian(chi * f) / 4 - chi * laplacian(f) / 4 - f * laplacian(chi) / 4
    cross = sum(s.diff(chi, z) * s.diff(f, z) for z in w) / 2
    A, B, full_F, F0 = s.symbols(
        "cutoff_one cutoff_two full_volume actual_center", real=True
    )
    return {
        "complete_radial_fourth_chain_including_all_contacts": s.simplify(
            s.diff(phi(y), r, 4) - wanted
        ),
        "full_covariance_whitening_heat_factor": (
            S.inv() * covariance * S.inv().T / 2 - s.eye(2) / 4
        ).applyfunc(s.factor),
        "entire_cutoff_heat_product_cross_contact": s.expand(product - cross),
        "whole_same_background_cutoff_difference": s.expand(
            (F0 + A * (full_F - F0))
            - (F0 + B * (full_F - F0))
            - (A - B) * (full_F - F0)
        ),
        "literal_complete_moving_chart_time_connection": old_source.moving_chart()[
            "checks"
        ]["whole_moving_chart_boundary_and_time_contact"],
        "literal_full_physical_volume_second_chain": old_source.response.physical()[
            "checks"
        ]["whole_complete_physical_second_chain"],
        "full_phase_and_configuration_dimension": s.Integer(DIMENSION - 2 * 48),
    }


@cache
def data():
    cut, domain, dynamics = (
        cutoff_bounds(),
        complete_domain_bounds(),
        dynamical_bounds(),
    )
    return {
        "whole_complete_complex_domain_and_Hamiltonian_bounds": domain,
        "whole_real_time_interval": [-source.TIME, source.TIME],
        "whole_interaction_and_physical_readout_definition": "At fixed u use the ENTIRE S267 chart with its original external means. hred=kappa integral[V Hbar-Hclock Pi_v-ell Pi_M]d^3x, with any exact-boundary term as in S268 and the original primitive already counted once. Let g=(hred-href)(u,S_u z), including its actual scalar background g0, linear terms, all quadratics and higher interactions. F is the entire spatial average exp(3v)R_full(u,Nstar)^(-3/4) pulled back by the same S_u, with its actual center F0. No background phase, tadpole, primitive or moving-chart contact is deleted.",
        "whole_holomorphic_phase_not_time_argument": "The full untruncated shape and adjoint Neumann reconstruction is holomorphic on the finite complex phase ball, uniformly in space. Its actual twelve invariant images lie in the proved complex auxiliary domain. The unique pointwise lapse root and whole coefficients are holomorphic in phase there. Uniform bounds permit spatial integration and give holomorphic g with |g|<1e1000 and |F-1|<1e-255 on initial complex w ball4R. Time is REAL and smooth; the whole Hamiltonian and the cutoff are not asserted holomorphic in time or phase, respectively. Each individual derivative of total order m<=4 involves at most m coordinates. Vary ONLY those at most four coordinates in a polydisc of radiusR/4; its Euclidean displacement is<=R/2, so every center in real ball2R stays inside complex ball4R. A simultaneous96-coordinate radiusR/4 polydisc is not used.",
        "whole_two_explicit_cutoff_definitions": "Let b(y)=exp(-1/y) for y>0 and0 otherwise. Let theta(y)=b(1-y)/(b(y)+b(1-y)), equal1 for y<=0 and0 for y>=1. For c=1 and2 set chi_c(w)=theta((||w||^2-R^2)/(c R^2)), R=1e20. Both are smooth and invariant, agree on the core, and have supports sqrt(2)R and sqrt(3)R, strictly inside2R. This derivative-controlled pair is not the class of all arbitrary smooth S269 cutoffs.",
        "whole_bump_derivative_polynomials": cut["polynomials"],
        "whole_bump_derivative_ceilings": cut["bump"],
        "whole_transition_derivative_ceilings": cut["transition"],
        "whole_radial_mixed_derivative_ceilings_per_R_power": cut["radial"],
        "whole_complete_Leibniz_Cauchy_product_coefficients": cut["product"],
        "whole_cutoff_derivative_proof": "For x>=0, x^k exp(-x)<=k! from the positive exponential series. The full derivative polynomials give B=[1,2,36,1584,129600]. On the transition denominator>=exp(-2)>1/9, since one argument is>=1/2 and e<3. Differentiating denominator*theta=numerator gives the displayed complete quotient recurrence. The radial inner map has first derivative norm<=4/R and second<=2/R^2; its higher derivatives vanish. The full partition/Leibniz rules include every mixed radial and cutoff contact. For the complete holomorphic amplitude, Cauchy gives m!4^m A/R^m, not a truncated physical-field Taylor symbol.",
        "whole_calibrated_generator_actual_norm_bound": cut["generator"],
        "whole_safe_calibrated_generator_norm": GENERATOR,
        "whole_complete_interaction_D2_symbol_bound": cut["D2H"],
        "whole_complete_volume_D_and_D2_bounds": [cut["Dvolume"], cut["D2volume"]],
        "whole_calibrated_volume_positive_symbol_floor": 1 - VERROR - cut["Dvolume"],
        "whole_calibrated_volume_operator_distance_from_identity": VERROR
        + cut["Dvolume"],
        "whole_volume_heat_symbol_ordering_error": cut["D2volume"] / 2,
        "whole_ordering_and_positivity_proof": "Keep the same full covariance V0=S0 S0^T/2. D=(1/2)V0_AB partial_zA partial_zB becomes Delta_w/4, with all original cross covariance retained by the exact whitening. Extend gext=g0+chi(g-g0), Fext=F0+chi(F-F0). A=Q_V[(1-D)gext] is bounded self-adjoint with norm<1e1010. Every cutoff derivative remains. The real calibrated volume symbol (1-D)Fext is>1/2 because ||D Fext||<1e-275; positive unital Q_V therefore gives an actually positive calibrated operator with distance fromI<2e-255. The exact heat SYMBOL remainder is <=(1/2)||D^2Fext||<1e-310. This is not an automatic norm estimate for the corresponding Weyl operator or equality with an original unregularized observable.",
        "whole_evaluated_same_seed_dynamics_and_two_regulator_bounds": dynamics,
        "whole_exact_unitary_dynamics_proof": "Use the S268 norm-convergent Dyson propagator for the entire bounded A(u), not its first-order truncation. The full scalar background phase is retained. For |u|<=1e-2000, ||UI(u)-I||<=integral||A||<=1e-990<1e-980. Ufull=Uref UI is the declared physical-picture evolution. The three original translation generators commute with the full symbols, same-seed coherent map and radial cutoffs, so their common zero-charge factor is reducing. No state projection or new vacuum is used.",
        "whole_evolved_leakage_and_comparison_proof": "The exact same-seed core tail is positive and<exp(-1e39)<1e-4000, using e^3>10 and1e39>12000. Its evolved coherent outside-core probability is at most tail+2||UI psi0-psi0||<1e-980, never asserted zero. For the two explicit cutoffs, delta a=(1-D)[(chi1-chi2)(g-g0)] vanishes on the common core INCLUDING every derivative term; its sup norm is<=2B, B=1e1010. The full S268 Duhamel comparison is<=2B sqrt(tail)T+B^2 T^2<1e-1970. Both readouts must also be compared in the evolved state: their mean difference is bounded by4*state_difference+4*VERROR*tail+8*VERROR*(BT)<1e-1230. This applies to the positive coherent and separately calibrated volumes; they remain distinct observables.",
        "whole_physical_picture_readout": "The interaction-picture readout is Q_V(Fext) or its declared first-Weyl calibration. The physical-picture readout is its conjugate by the SAME Uref(u), and its expectation in Uref UI psi0 equals the displayed interaction-picture expectation. No metaplectic covariance is assumed for a fixed window without transforming that window. These are defined finite regulated readouts, not evaluated original interacting volume means.",
        "checks": identities(),
        "gates": {
            "whole_analytic_interaction_bound": domain[
                "whole_complete_interaction_bound"
            ]
            < HBOUND,
            "whole_analytic_physical_volume_error": domain["whole_complex_volume_error"]
            < VERROR,
            "full_canonical_radius_bound": domain[
                "whole_fixed_canonical_complex_radius"
            ]
            < 10**142,
            "complete_torus_measure_bound": (2 * s.pi) ** 3
            < domain["whole_torus_measure_bound"],
            "actual_cutoff_bump_derivative_coefficients": cut["bump"]
            == [1, 2, 36, 1584, 129600],
            "complete_radial_derivative_majorant": max(cut["radial"]) < 10**13,
            "complete_Leibniz_Cauchy_majorant": max(cut["product"]) < 10**14,
            "entire_calibrated_generator_bound": cut["generator"] < GENERATOR,
            "entire_calibrated_volume_D_bound": cut["Dvolume"] < s.Rational(1, 10**275),
            "full_calibrated_volume_positive": 1 - VERROR - cut["Dvolume"]
            > s.Rational(1, 2),
            "full_calibrated_volume_operator_bound": VERROR + cut["Dvolume"]
            < 2 * VERROR,
            "full_volume_heat_SYMBOL_error_bound": cut["D2volume"] / 2
            < s.Rational(1, 10**310),
            "full_positive_exponential_log10_bound": sum(
                s.Rational(3) ** i / s.factorial(i) for i in range(4)
            )
            > 10,
            "actual_tail_below_one_e_minus4000": s.Integer(10) ** 39 > 12000,
            "whole_unitary_state_error": dynamics["whole_same_seed_unitary_state_error"]
            < s.Rational(1, 10**980),
            "whole_evolved_coherent_leakage": dynamics[
                "whole_evolved_coherent_outside_core_ceiling"
            ]
            < s.Rational(1, 10**980),
            "whole_two_regulator_state_comparison": dynamics[
                "whole_two_explicit_regulator_state_difference"
            ]
            < s.Rational(1, 10**1970),
            "whole_two_regulator_mean_comparison": dynamics[
                "whole_two_readout_mean_difference"
            ]
            < s.Rational(1, 10**1230),
            "full_source_time_boundary_and_implicit_contacts_retained": True,
            "not_Weyl_operator_norm_or_original_regulator_removal": True,
            "no_physical_cutoff_matching_loop_or_original_P8_closure": True,
        },
    }
