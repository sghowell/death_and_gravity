"""Independent nonstationary, derivative-loss and compact-time controls."""

from functools import cache

import sympy as s

from . import chart


@cache
def data():
    N, vhat, c = s.symbols("N vhat c", real=True)
    h = s.symbols("h", positive=True)
    omega = -s.log(1 + (N**-2 - 1) / h) / 4
    density = c * N * s.exp(3 * (vhat + omega))
    Hess = s.Matrix(
        [
            [s.diff(density, x, y).subs({N: 1, vhat: 0}) for y in (N, vhat)]
            for x in (N, vhat)
        ]
    )
    q = chart.data()
    subst = {
        sym: (-c if str(sym) == "rho" else c)
        for sym in q["complete_nonstationary_clock_contact"].free_symbols
        if str(sym) in ("rho", "pressure")
    }
    expected = q["complete_nonstationary_clock_contact"].subs(subst).subs(chart.u, 0)
    t = s.symbols("t", real=True)
    H = 4 * t / (1 + t * t)
    kappa = s.symbols("kappa", positive=True)
    freq = s.symbols("frequency", positive=True)
    # A smooth cutoff equal to one near t=0 has identical local jets.
    test = s.cos(freq * t) / freq**4
    volterra = s.Matrix([[0, 0], [1, 0]])
    action_gradient = (volterra + volterra.T) / 2
    return {
        "uncancelled_cosmological_density_clock_Hessian_at_bounce": s.ImmutableMatrix(
            Hess.subs(h, 1)
        ),
        "nonzero_second_map_cosmological_contact_at_bounce": -3 * c / 2,
        "prepared_norm_diagnostic": "For a fixed smooth compact cutoff equal to one near the bounce, f_frequency=frequency^-4 cutoff*cos(frequency*t) has bounded C4 norm but its tenth local derivative is -frequency^6. This diagnoses the stated C10 comparison proof; it does NOT prove that the actual retarded stress operator cannot satisfy a sharper estimate.",
        "flat_local_density": "The fixed mu=m local Proca prescription contains c=5m^4/(128pi^2), rho=-c and P=c. Its physical stress variation vanishes but its uncancelled clock-current response does not.",
        "finite_time_boundary": "The displayed ratio of this one local term to kappa H^2 grows quadratically in either tail. This is not a total-state asymptotic estimate or a no-go, and no vacuum/Newton counterterm is changed here.",
        "causal_boundary": "Mean-current pullback preserves the retarded kernel plus local contacts. It is not the gradient of half a single-branch source times its retarded inverse; no symmetric-noise bound follows.",
        "checks": {
            "literal_cosmological_density_contact_matrix": s.ImmutableMatrix(
                (Hess.subs(h, 1) - expected).applyfunc(s.simplify)
            ),
            "nonzero_cosmological_NN_bounce": s.simplify(
                Hess[0, 0].subs(h, 1) - 15 * c / 4
            ),
            "nonzero_cosmological_cross_bounce": s.simplify(
                Hess[0, 1].subs(h, 1) - 15 * c / 2
            ),
            "physical_constant_stress_response_zero": s.diff(c, N),
            "cutoff_plateau_fourth_jet": s.diff(test, t, 4).subs(t, 0) - 1,
            "cutoff_plateau_tenth_jet": s.diff(test, t, 10).subs(t, 0) + freq**6,
            "specified_local_tail_ratio_leading": s.limit(
                (c / (kappa * H * H)) / t**2, t, s.oo
            )
            - c / (16 * kappa),
            "causal_vs_single_branch_control": action_gradient[0, 1] - s.Rational(1, 2),
        },
        "controls": {
            "second_map_contact_is_nonzero": True,
            "C10_proof_not_relabelled_C4": True,
            "local_tail_not_total_quantum_asymptotic": True,
            "retarded_response_not_single_branch_Hessian": True,
        },
    }
