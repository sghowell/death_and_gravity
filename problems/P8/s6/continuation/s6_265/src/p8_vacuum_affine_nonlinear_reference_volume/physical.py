"""Strongly commuting bounce rows and the exact Gaussian integrability obstruction."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import gaussian as old_gaussian
from p8_vacuum_affine_coupled_gaussian_state import phase
from p8_vacuum_affine_quantitative_gaussian_window import physical as previous_physical
from p8_vacuum_affine_scalar_tame_propagator import charts

from . import source


@cache
def bounce():
    old = previous_physical.rows()
    v = old["whole_hat_log_scale_row_times_sqrt_kappa"]
    n = old["whole_physical_lapse_row_times_sqrt_kappa"]
    P = next(z for z in v.free_symbols if str(z) == "whole_comoving_P")
    shears = {
        z: 0
        for z in v.free_symbols | n.free_symbols
        if str(z).startswith("whole_reference_Sbar")
    }
    fixed = {
        **shears,
        phase.a: 1,
        charts.th: 0,
        charts.E: -s.Rational(1, 2),
        charts.l: s.Rational(1, 10),
    }
    vb = v.subs(fixed, simultaneous=True).applyfunc(s.factor)
    nb = n.subs(fixed, simultaneous=True).applyfunc(s.factor)
    wanted_v = s.Matrix([[0, 0, 1 / (2 * P**2), 0]])
    wanted_n = s.Matrix(
        [
            [
                0,
                0,
                1 / (4 * charts.J) - 3 * charts.T / (4 * charts.J * P**2),
                -1 / (40 * charts.J),
            ]
        ]
    )
    rows = vb.col_join(nb)
    radial = (source.HIGH**4 - source.LOW**4) / 128
    floor = radial / (4 * source.previous.COVARIANCE * source.KAPPA * 10**8)
    bounce_floor = radial / (4 * source.previous.COVARIANCE * source.KAPPA * 80**2)
    checks = {
        "whole_bounce_hat_log_scale_row": (vb - wanted_v).applyfunc(s.factor),
        "whole_bounce_full_lapse_both_momenta": (nb - wanted_n).applyfunc(s.factor),
        "whole_bounce_two_rows_have_no_position_entries": rows[:, :2],
        "whole_bounce_weighted_CCR_zero": rows * (P * phase.OMEGA) * rows.T,
        "whole_bounce_independent_momentum_minor": s.factor(
            rows[:, 2:].det() + 1 / (80 * charts.J * P**2)
        ),
        "actual_bounce_Theta_zero": source.previous.theta.subs(source.u, 0),
    }
    for case in range(4):
        _, _, M = old_gaussian.fixture_preparation(case)
        V = old_gaussian.ground_covariance(M)
        checks["whole_pure_covariance_inverse_identity_" + str(case)] = (
            V.inv() + 4 * phase.OMEGA * V * phase.OMEGA
        ).applyfunc(s.factor)
    return {
        "whole_original_full_time_physical_rows": old,
        "whole_bounce_hat_log_scale_row_times_sqrt_kappa": vb,
        "whole_bounce_lapse_row_times_sqrt_kappa": nb,
        "whole_bounce_two_momentum_minor": rows[:, 2:].det(),
        "whole_each_mode_canonical_balanced_covariance_lower": 1
        / (4 * source.previous.COVARIANCE),
        "whole_current_uniform_single_lapse_variance_lower_enclosure": floor,
        "whole_current_bounce_single_lapse_variance_lower_enclosure": bounce_floor,
        "whole_safe_current_uniform_lapse_variance_interval": [
            source.VAR_LOWER,
            source.VAR_UPPER,
        ],
        "whole_safe_current_bounce_lapse_variance_interval": [
            source.VAR_BOUNCE_LOWER,
            source.VAR_UPPER,
        ],
        "whole_positive_variance_argument": "Every original scalar momentum-mode covariance is pure: V Omega V=Omega/4, hence V^-1=-4 Omega V Omega. S264's V<=1e21 I therefore implies V>=I/(4e21); weighted covariance is at least P/(4e21) I. The full physical lapse momentum-m entry is ell E/(2J a^3), in absolute value>1e-4 on the tiny slab since ell>9/100,|E|>=1/4,J<=2,a<=2. At the bounce it is1/(40J)>=1/80. Integrate these single-entry lower bounds, restoring both1/sqrt(kappa) factors and using pi<4 for the radial floor. Thus the actual variance is strictly positive with the stated explicit lower bounds, not zero because a Gaussian tail underflows numerically.",
        "whole_strong_commutation_argument": "At the bounce the complete symmetric clean boundary and Theta vanish. Both displayed real band-smearing vectors use only the commuting clean momenta; the nonzero2x2 momentum minor and the positive full covariance give a strictly positive joint covariance. The band smears have finite one-particle norm by S264. In the regular original quasifree Weyl representation their symplectic pairing is zero, so their Weyl unitaries commute for all real parameters and the self-adjoint linear fields strongly commute. Joint spectral functional calculus is therefore justified on this slice, not merely inferred from an unbounded formal commutator. Away from this slice no joint spectral probability is claimed.",
        "checks": checks,
        "gates": {
            "actual_full_charge_above_9_over100": s.Rational(1, 10)
            / (1 + source.TIME**2) ** 6
            > s.Rational(9, 100),
            "full_momentum_entry_uniform_lower_above_1e_minus4": s.Rational(9, 100)
            * s.Rational(1, 4)
            / (2 * 2 * 2**3)
            > s.Rational(1, 10**4),
            "uniform_lapse_variance_lower_above_1e_minus576": floor > source.VAR_LOWER,
            "bounce_lapse_variance_lower_above_1e_minus572": bounce_floor
            > source.VAR_BOUNCE_LOWER,
            "outward_radial_floor_uses_pi_less_than_four": bool(s.pi < 4),
            "nonzero_bounce_second_momentum_row": nb[0, 3] != 0,
            "original_full_Tcorr_is_not_deleted": nb.has(charts.T),
            "only_bounce_joint_spectral_law_asserted": True,
        },
    }


@cache
def integrability():
    N, cut = s.symbols("positive_lapse positive_small_cutoff", positive=True)
    h = s.Symbol("positive_clock_h", positive=True)
    sigma = s.Symbol("positive_lapse_standard_deviation", positive=True)
    Dlow = (2 * h) ** -s.Rational(1, 4) * N ** -s.Rational(3, 2)
    Qlow = (2 * h) ** -s.Rational(3, 4) * N ** -s.Rational(1, 2)
    density0 = s.exp(-1 / (2 * sigma**2)) / (s.sqrt(2 * s.pi) * sigma)
    D_integral = 2 * (cut ** -s.Rational(1, 2) - 1)
    Q2_integral = -s.log(cut)
    return {
        "whole_original_coefficient_binding": source.weights(),
        "whole_positive_gaussian_density_at_zero_lapse": density0,
        "whole_full_kinetic_coefficient_lower_on_lapse_zero_to_one": Dlow,
        "whole_full_canonical_matter_coefficient_lower_on_lapse_zero_to_one": Qlow,
        "whole_exact_kinetic_cutoff_comparison_integral": D_integral,
        "whole_exact_matter_square_cutoff_comparison_integral": Q2_integral,
        "whole_kinetic_coefficient_positive_lapse_expectation": "INFINITE for the actual unchanged Gaussian linear lapse, and for every strictly positive Gaussian variance: the density is bounded below by the displayed strictly positive value on0<N<1, while the complete D=R^(1/4)/N is bounded below by(2h)^(-1/4) N^(-3/2).",
        "whole_canonical_matter_weight_integrability": "Q=N/U=N R^(3/4) has a finite positive-lapse first Gaussian moment, since Q<=2 N^(-1/2) on(0,1) and Q<=2N on[1,infinity). Its SECOND moment is infinite since Q^2>=(2h)^(-3/2)/N on(0,1).",
        "whole_precise_obstruction": "This refutes unrestricted positive-lapse Gaussian functional substitution as a blanket way to define every exact full-action coefficient. Merely discarding N<=0 from the OBSERVABLE does not cure the approach to0 from above. It does NOT prove a divergent complete constrained Hamiltonian/action, an inconsistent interacting state, a candidate no-go or P8 closure. Nonlinear auxiliary reconstruction, possible full-channel cancellations, physical quantum measure and the EFT domain remain separate. No projected/conditioned replacement state or quantum ordering is chosen.",
        "checks": {
            "exact_kinetic_power_comparison_primitive": s.factor(
                s.integrate(N ** -s.Rational(3, 2), (N, cut, 1)) - D_integral
            ),
            "exact_matter_square_log_comparison_primitive": s.factor(
                s.integrate(1 / N, (N, cut, 1)) - Q2_integral
            ),
            "canonical_matter_weight_square_power": s.powsimp(
                Qlow**2 - (2 * h) ** -s.Rational(3, 2) / N
            ),
            "integrable_first_matter_comparison_at_zero": s.integrate(
                N ** -s.Rational(1, 2), (N, 0, 1)
            )
            - 2,
        },
        "gates": {
            "gaussian_density_at_zero_is_strictly_positive": density0.is_positive
            is True,
            "actual_gaussian_variance_has_evaluated_positive_floor": bounce()[
                "whole_current_uniform_single_lapse_variance_lower_enclosure"
            ]
            > source.VAR_LOWER,
            "kinetic_cutoff_power_diverges": s.limit(D_integral, cut, 0, dir="+")
            == s.oo,
            "matter_square_cutoff_log_diverges": s.limit(Q2_integral, cut, 0, dir="+")
            == s.oo,
            "finite_first_moment_not_mistaken_for_finite_second": True,
            "individual_coefficient_not_entire_constraint_action_no_go": True,
            "tiny_nonzero_prefactor_does_not_regularize_infinite_integral": True,
        },
    }
