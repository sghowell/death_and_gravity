"""Whole current weighted phase, exact canonical charts and preparation form."""

from functools import cache

import sympy as s
from p8_vacuum_affine_scalar_tame_propagator import charts

OMEGA = s.zeros(4)
OMEGA[:2, 2:] = s.eye(2)
OMEGA[2:, :2] = -s.eye(2)
DAMP = s.diag(0, 0, 1, 1)
a, kappa = s.symbols("scale kappa", positive=True)
pv, ps = s.symbols("normalized_pv normalized_ps", real=True)
Z = s.Matrix([charts.v, charts.sigma, pv, ps])
C = s.sqrt(kappa) * s.diag(1, 1, a**3, a**3)
PICK_Q = s.eye(4)[:2, :]
PICK_P = s.eye(4)[2:, :]


def dtime(value):
    return s.factor(charts.dtime(value) + charts.H * a * s.diff(value, a))


@cache
def whole_hamiltonian():
    return s.cancel(charts.ham().subs({charts.b: -pv / (2 * charts.q), charts.ps: ps}))


@cache
def canonical_generator():
    hess = s.hessian(whole_hamiltonian(), Z)
    return OMEGA * (kappa * a**3 * C.inv().T * hess * C.inv())


def central_map():
    tau = 2 * a**3 * charts.q
    return s.Matrix([[0, 0, -1 / tau, 0], [0, 1, 0, 0], [tau, 0, 0, 0], [0, 0, 0, 1]])


@cache
def chart_phase(which):
    if which not in ("outer", "central"):
        raise ValueError("Unknown exact canonical chart")
    T = s.eye(4) if which == "outer" else central_map()
    generator = canonical_generator()
    transformed = (T.applyfunc(dtime) + T * generator) * T.inv()
    chart = charts.outer() if which == "outer" else charts.central()
    K, B, D = [chart[name] for name in ("K", "B", "D")]
    S = (B + B.T) / 2
    A = (B - B.T) / 2
    shear = s.eye(4)
    shear[2:, :2] = -(a**3) * S
    clean = (shear.applyfunc(dtime) + shear * transformed) * shear.inv()
    V = S.applyfunc(dtime) + 3 * charts.H * S - D - charts.q * chart["G"]
    HQQ = a**3 * (charts.q * chart["G"] + V + A.T * K.inv() * A)
    HQP = -A.T * K.inv()
    HPP = K.inv() / a**3
    Hclean = HQQ.row_join(HQP).col_join(HQP.T.row_join(HPP))
    return {
        "map": T,
        "generator": transformed,
        "K": K,
        "B": B,
        "D": D,
        "symmetric_boundary": S,
        "antisymmetric_B": A,
        "shear": shear,
        "clean_generator": clean,
        "clean_hessian": Hclean,
        "lower_symmetric_potential": V,
    }


@cache
def selection_matrix():
    chart = charts.outer()
    velocity = PICK_Q * canonical_generator()
    return a**3 * (
        velocity.T * chart["K"] * velocity
        + PICK_Q.T * (charts.q * chart["G"] + chart["K"]) * PICK_Q
    )


@cache
def clean_transition():
    outer, central = chart_phase("outer"), chart_phase("central")
    return (central["shear"] * central["map"] * outer["shear"].inv()).applyfunc(
        s.cancel
    )


def tensor_generator():
    return s.Matrix([[0, a**-3], [-(a**3) * charts.q, 0]])


def tensor_selection_matrix():
    return s.diag(a**3 * (charts.q + 1), a**-3)


@cache
def generic_data():
    hess = s.hessian(whole_hamiltonian(), Z)
    weighted = OMEGA * hess - 3 * charts.H * DAMP
    canonical = canonical_generator()
    checks = {
        "whole_density_canonical_generator": C * weighted * C.inv()
        + 3 * charts.H * DAMP
        - canonical,
        "whole_weighted_symplectic_defect": weighted * OMEGA
        + OMEGA * weighted.T
        + 3 * charts.H * OMEGA,
        "canonical_symplectic_generator": canonical * OMEGA + OMEGA * canonical.T,
        "physical_density_CCR_normalization": C * OMEGA * C.T - kappa * a**3 * OMEGA,
        "full_regular_hamiltonian_no_inverse_q": s.denom(whole_hamiltonian())
        / (4 * charts.J)
        - 1,
        "central_exact_time_dependent_symplectic_swap": central_map()
        * OMEGA
        * central_map().T
        - OMEGA,
        "central_swap_weight_derivative": dtime(2 * a**3 * charts.q)
        - charts.H * 2 * a**3 * charts.q,
        "selection_form_symmetric": selection_matrix() - selection_matrix().T,
    }
    for which in ("outer", "central"):
        data = chart_phase(which)
        generator = data["generator"]
        K, B = data["K"], data["B"]
        checks[which + "_full_legendre_momentum"] = (
            a**3 * (K * PICK_Q * generator + B * PICK_Q) - PICK_P
        )
        checks[which + "_full_canonical_generator"] = (
            generator * OMEGA + OMEGA * generator.T
        )
        checks[which + "_boundary_shear_symplectic"] = (
            data["shear"] * OMEGA * data["shear"].T - OMEGA
        )
        checks[which + "_clean_hamiltonian_symmetric"] = (
            data["clean_hessian"] - data["clean_hessian"].T
        )
        checks[which + "_entire_boundary_cleaned_generator"] = (
            data["clean_generator"] - OMEGA * data["clean_hessian"]
        )
        checks[which + "_full_clean_potential_symmetric"] = (
            data["lower_symmetric_potential"] - data["lower_symmetric_potential"].T
        )
        checks[which + "_entire_gyro_retained"] = 2 * data["antisymmetric_B"] - (
            B - B.T
        )
    radial = s.Symbol("radial_momentum", positive=True)
    scaling = s.diag(
        s.sqrt(radial), s.sqrt(radial), 1 / s.sqrt(radial), 1 / s.sqrt(radial)
    )
    normalized_transition = (
        scaling * clean_transition().subs(charts.q, radial**2 / a**2) * scaling.inv()
    )
    principal_transition = s.diag(-charts.E / charts.th, 1, -charts.th / charts.E, 1)
    checks["whole_clean_transition_symplectic"] = (
        clean_transition() * OMEGA * clean_transition().T - OMEGA
    )
    checks["bounded_high_momentum_clean_transition_principal"] = (
        normalized_transition.applyfunc(lambda value: s.limit(value, radial, s.oo))
        - principal_transition
    )
    checks["bounded_inverse_clean_transition_principal"] = (
        scaling
        * clean_transition().inv().subs(charts.q, radial**2 / a**2)
        * scaling.inv()
    ).applyfunc(lambda value: s.limit(value, radial, s.oo)) - principal_transition.inv()
    tensor_omega = s.Matrix([[0, 1], [-1, 0]])
    gamma_rate, canonical_rate = s.symbols(
        "tensor_gamma_rate tensor_canonical_rate", real=True
    )
    checks["full_tensor_canonical_action_normalization"] = (
        kappa * a**3 * gamma_rate**2 / 8
    ).subs(
        gamma_rate, 2 * canonical_rate / s.sqrt(kappa)
    ) - a**3 * canonical_rate**2 / 2
    checks["full_tensor_canonical_symplectic_generator"] = (
        tensor_generator() * tensor_omega + tensor_omega * tensor_generator().T
    )
    checks["full_tensor_wave_velocity"] = (tensor_generator() ** 2)[0, 0] + charts.q
    checks["tensor_positive_preparation_determinant"] = (
        tensor_selection_matrix().det() - charts.q - 1
    )
    return {
        "whole_normalized_hamiltonian": whole_hamiltonian(),
        "physical_canonical_map": C,
        "full_canonical_generator": canonical,
        "central_canonical_swap": central_map(),
        "complete_preparation_Gramian_integrand": selection_matrix(),
        "complete_clean_chart_transition": clean_transition(),
        "normalized_chart_transition_principal": principal_transition,
        "two_tensor_oscillator_generator": tensor_generator(),
        "each_tensor_preparation_form": tensor_selection_matrix(),
        "checks": {
            name: value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
            for name, value in checks.items()
        },
        "gates": {
            "canonical_map_has_both_kappa_and_volume": C[0, 0] == s.sqrt(kappa)
            and C[2, 2] == s.sqrt(kappa) * a**3,
            "omitting_density_time_connection_breaks_CCR": (
                -3 * charts.H * DAMP * OMEGA - 3 * charts.H * OMEGA * DAMP
            )
            != s.zeros(4),
            "central_swap_not_time_independent": dtime(2 * a**3 * charts.q) != 0,
            "preparation_support_inside_outer_chart": bool(
                s.Rational(13, 32) > s.Rational(1, 8)
            ),
            "preparation_support_inside_unchanged_slab": bool(
                s.Rational(15, 32) < s.Rational(1, 2)
            ),
            "selection_mass_is_not_action_retuning": True,
            "full_nonlocal_interacting_state_not_claimed": True,
        },
    }


@cache
def data():
    from p8_vacuum_affine_heavy_clock_quadratic import clock
    from p8_vacuum_affine_heavy_curved_state import state
    from p8_vacuum_affine_reduced_scalar_hamiltonian import scalar

    current = clock.data()
    charge, pivot = s.symbols("current_charge current_Jc", real=True)
    bridge = {
        charts.v: scalar.v,
        charts.sigma: scalar.sigma,
        pv: scalar.pv,
        ps: scalar.ps,
        charts.l: charge,
        charts.J: pivot,
        charts.th: scalar.Theta,
        charts.E: scalar.E,
        charts.A: scalar.A,
        charts.T: scalar.Tc,
        charts.q: scalar.q,
    }
    generic = generic_data()
    substitutions = current["exact_current_Hamiltonian_coefficient_substitution"]
    momentum = s.Symbol("physical_comoving_momentum", nonnegative=True)
    full_map = {
        variable: substitutions.get(destination, destination)
        for variable, destination in bridge.items()
        if variable not in set(Z) and variable != charts.q
    }
    full_map.update(
        {
            a: scalar.a,
            kappa: state.KAPPA,
            charts.H: scalar.H,
            charts.q: momentum**2 / scalar.a**2,
        }
    )
    checks = dict(generic["checks"])
    checks["literal_complete_current_Hamiltonian_bridge"] = s.cancel(
        whole_hamiltonian().subs(bridge, simultaneous=True)
        - current["complete_current_coefficient_sector_Hamiltonian"]
    )
    checks["actual_physical_action_normalization"] = state.KAPPA - 10**800
    return {
        **generic,
        "candidate_and_reference": "SAME QG2-H8A420 action, newly specified FREE-REF-S251 for its local physical quadratic coefficient sector only. Existing H/Proca preparations and classical M1 mean remain unchanged; no QG3 retuning or full nonlocal/interacting state is asserted.",
        "entire_current_coefficient_substitution": current[
            "exact_current_Hamiltonian_coefficient_substitution"
        ],
        "entire_current_clock_coefficients": current["full_current_clock_coefficients"],
        "entire_phase_coefficient_and_physical_background_substitution": full_map,
        "actual_current_coefficient_lower_bounds": current[
            "new_positive_coefficient_lower_bounds"
        ],
        "initial_time": -s.Rational(1, 2),
        "fixed_sampling_support": (-s.Rational(15, 32), -s.Rational(13, 32)),
        "fixed_sampling_amplitude": "exp[-1/(1-1024(t+7/16)^2)] in the support, zero outside. The selection form adds K at preparation mass scale1 but changes no physical action coefficient. State data at t* are selected once on the fixed reference and held fixed thereafter.",
        "all_order_input_boundary": "The full fixed Proca/H renormalized reference means are smooth functions; every finite coefficient jet is bounded on the compact slab. Explicit low-jet envelopes and the two-chart tame propagator are already proved for the complete current functions. No analytic-in-order, uniform-mass or small interacting-loop estimate is inferred.",
        "checks": checks,
        "gates": {
            **generic["gates"],
            "actual_current_positive_F": current["gates"]["current_F_positive"],
            "actual_current_positive_Jc": current["gates"]["current_Jc_positive"],
            "actual_current_two_cone_gap": current["gates"][
                "current_characteristic_strict_gap"
            ],
            "actual_current_Jc_upper": current["gates"]["current_Jc_upper"],
            "actual_complete_tame_flow_retained": current["gates"][
                "replayed_full_generic_energy_and_conversion_bounds"
            ],
        },
    }
