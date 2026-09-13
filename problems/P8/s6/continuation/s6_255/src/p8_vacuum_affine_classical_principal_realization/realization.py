"""Actual current constraint data, strict margins and on-shell map contacts."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_principal_obstruction import background as old
from p8_vacuum_affine_coupled_principal_obstruction import coupled as c
from p8_vacuum_affine_coupled_principal_obstruction import principal
from p8_vacuum_affine_physical_background_vertices import parent

from . import homogeneous as bg
from . import rolling

PROFILE_BOUND = s.Integer(2) / 10**400


@cache
def data():
    u, N = bg.u, bg.N
    rows = parent.lapse_jets()["rows"]
    U = rows["U"]
    Z0, ZN, ZNN = U[0], U[1] - U[0], U[2] - 2 * U[1] + 2 * U[0]
    F0, F1, F2 = rows["Fhat"][:3]
    m = s.Symbol("initial_M1_rate", real=True)
    Qn = (F0 + F1 + ZN * m**2 / 2).subs(u, 0)
    Cnn = (2 * F1 + F2 + ZNN * m**2 / 2).subs(u, 0) / 2
    J = Cnn - (ZN.subs(u, 0) * m) ** 2 / (2 * Z0.subs(u, 0))
    rho, pressure = parent.RHO.subs(u, 0), parent.PRESSURE.subs(u, 0)
    T = rho - 3 * pressure / 2
    root_square = s.Rational(1, 100) - 4 * T
    root_J = s.factor(J.subs(m**2, root_square))
    root_Cnn = s.factor(Cnn.subs(m**2, root_square))
    # Full off-clock initial root, not the root of a truncated clock action.
    formal_U = bg.R ** -s.Rational(3, 4)
    primitive_N = old.coefficients()["primitive_lapse_derivative"]
    initial_geometry = s.diff(
        N * formal_U * bg.F + 9 * formal_U * s.diff(bg.R, u) ** 2 / (16 * bg.R * N), N
    ) - s.diff(primitive_N, u)
    initial_ZN = s.diff(formal_U / N, N)
    exact_root = s.sqrt(-2 * initial_geometry / initial_ZN)
    direct_geometry = bg.action()["complete_lapse_constraint"].subs(
        {c.H: 0, bg.hbar: 0, bg.mc: 0, bg.mh: 0}, simultaneous=True
    )
    primitive_rule = {s.diff(bg.I, u, N): s.diff(primitive_N, u)}
    refJ = s.Rational(243, 160) + rho - 7 * pressure / 8
    refC = s.Rational(38, 25) + rho - 7 * pressure / 8
    Tbound = 5 * PROFILE_BOUND / 2
    lower_square = s.Rational(1, 100) - 4 * Tbound
    lower_J = s.Rational(243, 160) - 15 * PROFILE_BOUND / 8 - 3 * Tbound / 4
    lower_Cnn = s.Rational(38, 25) - 15 * PROFILE_BOUND / 8 - Tbound / 4
    checks = {
        "entire_current_central_Z": Z0.subs(u, 0) - 1,
        "entire_current_central_ZN": ZN.subs(u, 0) - s.Rational(1, 2),
        "entire_current_central_ZNN": ZNN.subs(u, 0) + s.Rational(1, 4),
        "whole_constraint_at_old_fixed_M1_rate": Qn.subs(m, s.Rational(1, 10)) - T,
        "whole_constraint_at_exact_classical_root": Qn.subs(m**2, root_square),
        "exact_current_root_pivot_J": root_J - refJ - 3 * T / 4,
        "exact_current_root_pivot_Cnn": root_Cnn - refC - T / 4,
        "full_offclock_initial_constraint_primitive_elimination": (
            direct_geometry - initial_geometry
        ).subs(primitive_rule, simultaneous=True),
        "full_offclock_initial_positive_root_solves_constraint": initial_geometry
        + initial_ZN * exact_root**2 / 2,
        "actual_current_R_lapse_transversality": parent.lapse_jets()[
            "R_lapse_jets_zero_through_four"
        ][1].subs(u, 0)
        + 2,
        "actual_current_L2_at_classical_constraint_center": (
            4 * (rows["C3"][0] + rows["C3"][1])
        ).subs(u, 0)
        + 1,
    }
    fast = principal.fast_data()["positive_fast_rate_fourth_power"]
    onshell = rolling.coefficient_map()
    clean_contacts = {c.L0: 0, c.Vvv: 0, c.vs[0]: 0, c.vs[1]: 0}
    checks["entire_fast_rate_survives_all_classical_Euler_contacts"] = (
        fast.subs(clean_contacts, simultaneous=True) - fast
    )
    return {
        "whole_current_lapse_constraint_at_N1_u0": Qn,
        "actual_central_classical_M1_rate_squared": root_square,
        "actual_central_classical_pivot_J": root_J,
        "actual_central_classical_pivot_Cnn": root_Cnn,
        "entire_offclock_initial_geometry_constraint": initial_geometry,
        "entire_offclock_initial_ZN": initial_ZN,
        "entire_offclock_positive_M1_constraint_root": exact_root,
        "initial_root_evaluation": "Apply the full current R,F bindings and their derivatives first, and then setu0,N=1+epsilon,Hhat0,hbar0,mh0. The primitive derivative is eliminated exactly by I_N=3U R_u R_N/(4RN). Use a_hat1 and the old M1 field value, with its RATE set to this positive root. This is a comparison initial datum, not a replacement for the fixed reference preparation.",
        "entire_current_source_bindings": old.data()[
            "entire_current_function_bindings"
        ],
        "whole_fixed_profile_bindings": old.data()["whole_fixed_profile_bindings"],
        "actual_profile_C5_bound": PROFILE_BOUND,
        "central_positive_root_square_lower_bound": lower_square,
        "central_positive_J_lower_bound": lower_J,
        "central_positive_Cnn_lower_bound": lower_Cnn,
        "central_M1_rate_difference_upper_bound": 40 * Tbound,
        "whole_rolling_background_parameter_map": onshell[
            "entire_rolling_four_mode_parameter_map"
        ],
        "unforced_realization": "The actual fixed profiles obey the prior normalized C5 bound2*10^-400. AtN1 the constraint root has m^2>9/1000,J>3/2,Cnn>3/2 and |m-1/10|<10^-397. Continuity of the FULL root gives admitted initial data at every sufficiently small fixed nonzero epsilon. R_N=-2 andL2=-1 at the center imply r=R-1 nonzero and L2 nonzero at those data. The regular homogeneous Euler system supplies a local unforced classical solution, and all strict margins persist on a compact subinterval of positive length. Apply the full S254 fast theorem with the complete rolling map. This proves an unbounded-momentum Gaussian/Jacobi instability of those local CLASSICAL comparison solutions, not of the original fixed quantum mean or below a controlled cutoff.",
        "proximity_boundary": "At epsilon0 the comparison rate contains the explicitly retained profile correction, bounded as listed; it has not been proved either zero or nonzero. Do not claim an unforced family arbitrarily close to the old local off-shell reference without resolving that mismatch. The family is arbitrarily close to its constraint-satisfying classical center and quantitatively within the stated profile-sized correction of the old central data.",
        "checks": checks,
        "gates": {
            "actual_positive_root_square_strict_margin": bool(
                lower_square > s.Rational(9, 1000)
            ),
            "actual_J_strict_margin_at_constraint_center": bool(
                lower_J > s.Rational(3, 2)
            ),
            "actual_Cnn_strict_margin_at_constraint_center": bool(
                lower_Cnn > s.Rational(3, 2)
            ),
            "actual_M1_rate_correction_bound": bool(
                40 * Tbound < s.Integer(1) / 10**397
            ),
            "whole_root_not_a_clock_polynomial_substitution": exact_root.has(
                bg.R, bg.F
            ),
            "all_lower_Euler_contacts_leave_positive_fast_coefficient": not any(
                fast.has(key) for key in clean_contacts
            ),
            "local_classical_comparison_not_quantum_mean": True,
            "no_uniform_heavy_mass_lifetime_or_cutoff_onset": True,
        },
    }


@cache
def map_contact_data():
    # Independent finite-dimensional polynomial re-entry of the full chain
    # rule. Its functional version uses adjoints and compactly supported
    # variations; written proof specifies regular differential maps.
    y = s.Matrix(s.symbols("field_chart0:4", real=True))
    x = s.Matrix(s.symbols("original_field0:4", real=True))
    epsilon = s.Symbol("chart_curvature", real=True)
    field_map = s.Matrix(
        [
            y[0] + epsilon * y[1] * y[2],
            y[1] + epsilon * y[0] ** 2,
            y[2] + epsilon * y[3] ** 2,
            y[3] + epsilon * y[0] * y[1],
        ]
    )
    tadpole = s.Matrix(s.symbols("Euler_onepoint0:4", real=True))
    Q = s.Matrix(4, 4, s.symbols("quadratic0:16", real=True))
    Q = (Q + Q.T) / 2
    action = tadpole.dot(x) + (x.T * Q * x)[0] / 2
    zero = dict.fromkeys(y, 0)
    actual = s.hessian(action.subs(dict(zip(x, field_map)), simultaneous=True), y).subs(
        zero
    )
    Jacobian = field_map.jacobian(y).subs(zero)
    contact = sum(
        (tadpole[i] * s.hessian(field_map[i], y).subs(zero) for i in range(4)),
        s.zeros(4),
    )
    return {
        "whole_example_nonlinear_field_map": field_map,
        "whole_example_background_Jacobian": Jacobian,
        "whole_example_onepoint_contact": contact,
        "on_shell_bridge": "For the full local classical action S and an admissible nonlinear field map F, delta^2(S o F)=DF* delta^2 S DF+delta S D^2F, with the functional adjoint and the retained boundaries. On the constructed unforced background every full Euler expression vanishes, so the second term vanishes. Within regular local charts whose physical phase maps and inverses have polynomial momentum bounds, the exp(cP^(3/2)) propagator obstruction is preserved. This does not cover an arbitrary nonlocal map, changed UV action, or the nonstationary full quantum effective action.",
        "checks": {
            "whole_nonlinear_field_map_second_variation_chain": actual
            - Jacobian.T * Q * Jacobian
            - contact,
            "whole_on_shell_map_contact_zero": contact.subs(
                dict.fromkeys(tadpole, 0), simultaneous=True
            ),
            "whole_example_linear_map_invertible": Jacobian.det() - 1,
        },
        "gates": {
            "off_shell_contact_is_genuinely_nonzero": contact != s.zeros(4),
            "on_shell_contact_removed_by_Euler_not_assumption": True,
            "regular_physical_phase_map_class_explicit": True,
            "full_quantum_effective_Euler_not_assumed_zero": True,
        },
    }
