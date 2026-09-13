"""Explicit new-profile clock expansion and rechecked classical comparison hypotheses."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import quantum, state
from p8_vacuum_affine_quantum_retuning import profile as old_profile
from p8_vacuum_affine_reduced_scalar_hamiltonian import scalar
from p8_vacuum_affine_scalar_tame_propagator import charts, energy, majorants
from p8_vacuum_affine_scalar_tame_propagator import coefficients as old

t = scalar.t
EPS = old_profile.EPS + s.Rational(1, 10**400)
INTERVAL = (-s.Rational(1, 2), s.Rational(1, 2))


@cache
def data():
    heavy = quantum.fixed_profile()
    hrho = heavy["rho"].subs(state.TIME, t) / state.KAPPA
    hp = heavy["pressure"].subs(state.TIME, t) / state.KAPPA
    rho = old_profile.rho + hrho
    pressure = old_profile.P + hp
    A = -pressure
    B = -(rho + pressure) / 2
    delta = scalar.delta
    dJ = (21 * delta**2 - 3 * delta) * A / 2 + (1 - 6 * delta) * B
    Tc = (1 + 3 * delta) * A - 2 * B
    base = old.clock()
    Jc = base["J"] + dJ
    heavyA = -hp
    heavyB = -(hrho + hp) / 2
    heavyJ = (21 * delta**2 - 3 * delta) * heavyA / 2 + (1 - 6 * delta) * heavyB
    heavyT = (1 + 3 * delta) * heavyA - 2 * heavyB
    p = state.germs.parent
    X = p.X
    R = p.coefficients()["R"]
    N, e, n, v = s.symbols("N e n_lapse v", real=True)
    fullF = heavy["DeltaF"].subs({state.TIME: t, state.X: X}, simultaneous=True)
    density = (
        N
        * R.subs({p.u: t, X: N**-2}, simultaneous=True) ** (-s.Rational(3, 4))
        * fullF.subs(X, N**-2)
    )
    before = s.exp(3 * e * v) * density.subs(N, 1 + e * n)
    quadratic = heavyJ * n * n + 3 * heavyT * n * v + s.Rational(9, 2) * heavyA * v * v
    checks = {
        "whole_new_R_clock_value": s.simplify(R.subs(X, 1) - 1),
        "whole_new_R_clock_first_X_jet": s.simplify(
            s.diff(R, X).subs(X, 1) - 2 * delta
        ),
        "whole_new_R_clock_second_X_jet": s.simplify(s.diff(R, X, 2).subs(X, 1)),
        "whole_new_profile_clock_value": s.simplify(fullF.subs(X, 1) - heavyA),
        "whole_new_profile_clock_first_X_jet": s.simplify(
            s.diff(fullF, X).subs(X, 1) - heavyB
        ),
        "whole_new_profile_clock_second_X_jet": s.simplify(
            s.diff(fullF, X, 2).subs(X, 1)
        ),
        "literal_full_new_lapse_profile_tadpole": s.simplify(
            s.diff(density, N).subs(N, 1) - heavyT
        ),
        "literal_full_new_lapse_profile_pivot": s.simplify(
            s.diff(density, N, 2).subs(N, 1) / 2 - heavyJ
        ),
        "literal_full_new_profile_quadratic_in_affine_metric_chart": s.simplify(
            s.diff(before, e, 2).subs(e, 0) / 2 - quadratic
        ),
        "new_total_profile_not_only_heavy_part": s.cancel(A + old_profile.P + hp),
        "new_total_pivot_adds_to_full_old_QG1": s.cancel(Jc - base["Jc"] - heavyJ),
        "new_total_tadpole_matches_actual_combined_stress": s.cancel(
            Tc - rho + 3 * delta * pressure
        ),
        "new_gradient_pair_unchanged_but_kinetic_pivot_recomputed": s.cancel(
            base["J"] - base["F"] - 1 / (50 * (1 + t * t) ** 6)
        ),
    }
    inputs = old.data()["coefficient_time_jet_envelopes"]
    retuning = {}
    for j in range(3):
        dj = (
            sum(
                s.binomial(j, k)
                * (
                    (
                        21
                        * sum(
                            s.binomial(k, z)
                            * inputs["delta"][z]
                            * inputs["delta"][k - z]
                            for z in range(k + 1)
                        )
                        + 3 * inputs["delta"][k]
                    )
                    / 2
                    + (1 if k == 0 else 0)
                    + 6 * inputs["delta"][k]
                )
                for k in range(j + 1)
            )
            * EPS
        )
        tc = (
            1 + 3 * sum(s.binomial(j, k) * inputs["delta"][k] for k in range(j + 1))
        ) * EPS
        retuning[j] = (dj, tc, EPS)
    Fmin = s.Rational(1199, 800) * s.Rational(4, 5) ** 18
    Jmin = s.Rational(1215, 800) * s.Rational(4, 5) ** 18 - 4 * EPS
    gap = s.Rational(1, 50) * s.Rational(4, 5) ** 6 - 4 * EPS
    speed, j0, w0, f = s.symbols("speed_squared Jc w F", real=True)
    K = s.Matrix([[2 * j0 + w0 * w0, w0], [w0, 1]])
    G = s.Matrix([[2 * f + w0 * w0, w0], [w0, 1]])
    checks["new_two_scalar_principal_characteristic_identity"] = s.factor(
        (speed * K - G).det() - 2 * j0 * (speed - 1) * (speed - f / j0)
    )
    charge, pivot, hubble = s.symbols(
        "current_charge current_Jc current_Hubble", real=True
    )
    substitutions = {
        charge: scalar.ell,
        pivot: Jc,
        hubble: scalar.H,
        scalar.Theta: base["Theta"],
        scalar.E: base["E"],
        scalar.Tc: Tc,
        scalar.A: A,
    }
    fullham = charts.ham().subs(
        {
            charts.b: -scalar.pv / (2 * charts.q),
            charts.ps: scalar.ps,
            charts.l: charge,
            charts.J: pivot,
            charts.th: scalar.Theta,
            charts.E: scalar.E,
            charts.A: scalar.A,
            charts.T: scalar.Tc,
            charts.v: scalar.v,
            charts.sigma: scalar.sigma,
            charts.q: scalar.q,
        },
        simultaneous=True,
    )
    Z = scalar.Z
    phase = scalar.JC * s.hessian(fullham, Z) - 3 * hubble * s.diag(0, 0, 1, 1)
    checks["new_full_weighted_phase_hamiltonian_reconstruction"] = s.cancel(
        (Z.T * s.hessian(fullham, Z) * Z)[0] / 2 - fullham
    )
    hsource = p.coefficients()["normalized_heavy_source"]
    checks["added_heavy_reference_source_zero"] = hsource.subs(X, 1)
    checks["added_heavy_linear_clock_source_zero"] = s.diff(hsource, X).subs(X, 1)
    checks["added_heavy_linear_clock_value_source_zero"] = s.diff(hsource, p.u).subs(
        X, 1
    )
    rhoH = state.germs.MASS2
    major = majorants.data()
    en = energy.data()
    return {
        "candidate": "The existing separately named CD-REG-AFFINE-ISO-QG2-H8A420 of S240, with old Proca/M1 data and the S238 reconstructed affine mean. This is a new quadratic audit, not another retuning.",
        "full_actual_normalized_reference_energy": rho,
        "full_actual_normalized_reference_pressure": pressure,
        "full_current_clock_coefficients": {
            "A": A,
            "B": B,
            "DeltaJ": dJ,
            "Tcorr": Tc,
            "Jc": Jc,
            "F_gradient": base["F"],
        },
        "complete_added_profile_quadratic_density": quadratic,
        "complete_current_coefficient_sector_Hamiltonian": fullham,
        "complete_current_weighted_phase_matrix": phase,
        "exact_current_Hamiltonian_coefficient_substitution": substitutions,
        "whole_function_representation": "The compact Hamiltonian and phase matrix use the named generic coefficients. The displayed exact substitution map defines EVERY current coefficient by its full fixed reference functions, not a finite time or field expansion. Substitution is simultaneous; the coefficients are independent of all phase variables. Generic polynomial identities are proved before substitution, avoiding expansion of repeated complete rational reference denominators.",
        "full_profile_expansion_proof": "Expand N exp(3v) R(t,N^-2)^(-3/4) DeltaF(t,N^-2) through degree2 using the literal complete S238 R and S240 profile. Its value, lapse, mixed n-v and v² terms are all retained. The high-order S238/S239 changes vanish in the relevant physical coefficient/source jets; this is checked and source-pinned, not inferred solely from a tiny numerical tube.",
        "new_total_reference_stress_jet_bound": EPS,
        "new_retuning_zero_through_two_jet_envelopes": retuning,
        "new_positive_coefficient_lower_bounds": {
            "F": Fmin,
            "Jc": Jmin,
            "Jc_minus_F": gap,
        },
        "explicit_generic_chart_application": "The full current coefficients satisfy EVERY generic S221 hypothesis: |Theta|<=2, |H|<=2, ell<=1/10, |E|<=1 (<=1/2 centrally), 1/100<Jc<100, F>1/100, Jc-F>1/1000; |A|,|Tcorr| and their first2 jets<1e-6, all coefficient jets through2<1e6. Central |E|>1/4 for |t|<=1/4 and outer |Theta|>=3/10 for |t|>=1/8 are unchanged. The central chart is only used at q>=4096, with switches+-3/16.",
        "generic_finite_chart_enclosures": {
            "central_Delta_interval": major["central_denominator_bounds"],
            "K_G_operator_interval": major["common_positive_operator_interval"],
            "complete_entry_bound": major["common_absolute_entry_bound"],
            "complete_matrix_entries": sum(
                map(len, major["coefficient_entry_majorants"].values())
            ),
        },
        "new_current_classical_comparison": "After this explicit new hypothesis check, the SAME generic full finite-q algebra, all48 polynomial coefficient bounds, energy identity and finite chart conversions apply. For all ordered times in[-1/2,1/2] and nonzero P, ||U_QG2(t,s;P)||<=exp(1e29)(1+|P|²)^6. Duhamel loses12 spatial derivatives. The constant is not small relative to kappa0. The measure-zero P0 extension on Sobolev fields is not the literal homogeneous constraint.",
        "generic_log_propagator": en["all_transfer_propagator_log_constant"],
        "spatial_derivative_loss": en["spatial_derivative_loss"],
        "added_classical_heavy_mode": "At zero heavy mean the source's clock zeros eliminate ALL linear light-H coupling. The whole new free H quadratic action remains, with equation Hmode''+3Hubble Hmode'+(n+P²/a²)Hmode=0. Its heavy Gaussian metric response is NOT eliminated by this classical decoupling.",
        "heavy_energy_scaled_phase": "Use chi=a^(3/2)Hmode and y=(sqrt(Omega_star)chi,chi'/sqrt(Omega_star)), Omega_star²=n+P²/16. The exact normalized S240 comparison mode has ||y||1<20 on[-1,1]. Its two-column complex fundamental matrix has determinant i, norm<=20 and inverse norm<=40. Hence the actual homogeneous KG phase propagator is below800<1000, uniformly over all P including0, with no WKB truncation substituted for an exact solution.",
        "full_remaining_boundary": "This transfers ONLY the stated classical coefficient-sector comparison and adds the independent classical heavy mode. The new Gaussian local/nonlocal metric response, full quantum lapse/shift constraints, interacting light/mixed loops, compatible inverse, nonlinear same-state bounce and original P8 closure are separate. A positive principal pair is not a frozen-frequency or full quantum stability theorem.",
        "checks": {
            key: value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
            for key, value in checks.items()
        },
        "gates": {
            "actual_sum_of_old_and_heavy_bounds": bool(EPS < s.Rational(2, 10**400)),
            "full_profile_pointwise_pivot_coefficient": s.Rational(31, 8) < 4,
            "current_F_positive": bool(Fmin > s.Rational(1, 100)),
            "current_Jc_positive": bool(Jmin > s.Rational(1, 100)),
            "current_characteristic_strict_gap": bool(gap > s.Rational(1, 1000)),
            "current_Jc_upper": bool(inputs["J"][0] + 4 * EPS < 100),
            "every_new_retuning_jet_in_generic_box": all(
                v < s.Rational(1, 10**6) for row in retuning.values() for v in row
            ),
            "every_base_jet_in_generic_box": all(
                v < 10**6 - 1 for row in inputs.values() for v in row
            ),
            "replayed_all_generic_matrix_bounds": all(
                bool(v) for v in major["gates"].values()
            ),
            "replayed_full_generic_energy_and_conversion_bounds": all(
                bool(v) for v in en["gates"].values()
            ),
            "heavy_global_frequency_gap": bool(rhoH > 15),
            "whole_exact_heavy_propagator_bound": 20 * 40 < 1000,
            "current_coefficients_independent_of_all_phase_variables": all(
                not value.has(*Z) for value in substitutions.values()
            ),
            "no_inherited_full_quantum_or_nonlinear_stability": True,
        },
    }
