"""Unchanged coupled reference rows and their actual regular bounce symbols."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import gaussian, phase, uv
from p8_vacuum_affine_physical_background_vertices import coupled, parent
from p8_vacuum_affine_scalar_tame_propagator import charts


def transported_covariance(gramian, flow):
    if not isinstance(flow, s.MatrixBase) or flow.shape != (4, 4):
        raise TypeError("Require the whole exact four-phase transport")
    if any(not value.is_Rational for value in flow):
        raise TypeError("The finite diagnostic accepts exact rational transport only")
    if flow * phase.OMEGA * flow.T != phase.OMEGA:
        raise ValueError("The whole fixed-state transport must be symplectic")
    return (flow * gaussian.ground_covariance(gramian) * flow.T).applyfunc(s.simplify)


@cache
def reference_rows():
    a, kappa = phase.a, phase.kappa
    J, th, ell, E, T, q = charts.J, charts.th, charts.l, charts.E, charts.T, charts.q
    physical_map = phase.C
    normalized_v = s.Matrix([[1, 0, 0, 0]])
    normalized_n = s.Matrix([[-(2 * E * q + 3 * T), 3 * th * ell, th, ell * E]]) / (
        2 * J
    )
    vr = normalized_v * physical_map.inv()
    nr = normalized_n * physical_map.inv()
    substitutions = {
        coupled.D: 1,
        coupled.Z: 1,
        coupled.r: 0,
        coupled.dH: 0,
        coupled.J: J,
        coupled.th: th,
        coupled.w: -ell * E,
        coupled.c: ell,
        coupled.Lnv: 2 * E * q + 3 * T,
        coupled.v: charts.v,
        coupled.sigma: charts.sigma,
        coupled.pv: phase.pv,
        coupled.ps: phase.ps,
    }
    entire_n = (coupled.lapse_numerator() / (2 * coupled.J)).subs(
        substitutions, simultaneous=True
    )
    covariance = s.Matrix(
        4,
        4,
        lambda i, j: s.Symbol(
            "fixed_symmetric_covariance_" + str(min(i, j)) + str(max(i, j)), real=True
        ),
    )
    Cvv = (vr * covariance * vr.T)[0]
    Cvn = (vr * covariance * nr.T)[0]
    bracket = (vr * phase.OMEGA * nr.T)[0]
    return {
        "whole_original_physical_canonical_map": physical_map,
        "whole_normalized_lapse_row": normalized_n,
        "whole_physical_log_volume_row": vr,
        "whole_physical_lapse_row": nr,
        "whole_symmetric_canonical_covariance": covariance,
        "whole_Cvv_from_fixed_coupled_covariance": Cvv,
        "whole_Cvn_from_fixed_coupled_covariance": Cvn,
        "whole_lapse_volume_equal_time_bracket": bracket,
        "whole_actual_covariance_definition": "Use S(t,t*;P) V(M(P)) S(t,t*;P)^T from S251, with its unchanged full preparation Gramian and every actual coefficient, kappa and endpoint map. Both first-momentum and cross entries remain. The other tensor/Proca/H product blocks are unchanged, not integrated out or re-prepared.",
        "finite_state_definition": "For each selected nonzero real cosine/sine Fourier pair, use two copies of this SAME four-phase covariance. Their Wick/Weyl second moments give 2/Vol times each mode row contraction. This explicitly named finite periodic diagnostic samples the fixed radial Cauchy covariance; it does not replace the original R3 state or specify its homogeneous constraint sector. On R3 a bounded radial momentum band restricts the original two-point function without new preparation.",
        "checks": {
            "literal_entire_reference_lapse_constraint_row": s.factor(
                entire_n - (normalized_n * phase.Z)[0]
            ),
            "whole_physical_volume_row_normalization": vr
            - s.Matrix([[1, 0, 0, 0]]) / s.sqrt(kappa),
            "whole_phase_density_lapse_bracket": s.factor(
                bracket - th / (2 * J * kappa * a**3)
            ),
            "whole_bounce_linear_lapse_volume_commutator_zero": bracket.subs(th, 0),
            "whole_symmetric_cross_covariance": s.factor(
                Cvn - (nr * covariance * vr.T)[0]
            ),
        },
        "gates": {
            "lapse_cross_covariance_has_canonical_momenta": Cvn.has(covariance[0, 2])
            and Cvn.has(covariance[0, 3]),
            "general_reduced_lapse_and_volume_do_not_commute": bracket != 0,
            "bounce_zero_commutator_not_nonlinear_ordering_proof": True,
            "old_complete_reference_not_new_instantaneous_vacuum": True,
        },
    }


@cache
def bounce_symbols():
    J = s.Symbol("whole_positive_current_bounce_J", positive=True)
    cs = s.Symbol("whole_positive_current_scalar_speed", positive=True)
    E = s.Symbol("nonzero_central_E", real=True)
    ell = s.Symbol("whole_current_matter_charge", real=True)
    a, k, kap = s.symbols(
        "positive_scale positive_comoving_k actual_kappa", positive=True
    )
    q = k * k / (a * a)
    F = s.Symbol("whole_positive_current_gradient", positive=True)
    R = s.Matrix([[E / s.sqrt(2 * J), 0], [-ell * E / s.sqrt(2 * J), 1]])
    K0 = s.Matrix([[2 * J / E**2 + ell**2, ell], [ell, 1]])
    G0 = s.Matrix([[2 * F / E**2 + ell**2, ell], [ell, 1]])
    vrow = R.inv().T[0, :] / (2 * s.sqrt(kap) * a ** s.Rational(3, 2) * q)
    nrow = (ell * E * R.inv().T[1, :] - E * R.inv().T[0, :]) / (
        2 * J * s.sqrt(kap) * a ** s.Rational(3, 2)
    )
    pi_cov = s.diag(k * cs / a, k / a) / 2
    vv = s.factor((vrow * pi_cov * vrow.T)[0])
    vn = s.factor((vrow * pi_cov * nrow.T)[0])
    nn = s.factor((nrow * pi_cov * nrow.T)[0])
    central = charts.central()
    Bbounce = central["B"].subs({charts.th: 0, charts.H: 0}).applyfunc(s.factor)
    u = parent.u
    delta = 1 / (2 * parent.h)
    Ec = 1 - 3 * delta
    thc = parent.H - u / (1 + u * u) ** 4
    Fc = (
        thc * s.diff(Ec, u)
        - Ec * s.diff(thc, u)
        + parent.H * Ec * thc
        - thc**2
        - parent.ell**2 * Ec**2 / 2
    )
    # These are whole current reference coefficients, including the original
    # fixed profile functions; no evaluation of an interacting mean enters.
    Jc = (
        Fc
        + 1 / (50 * (1 + u * u) ** 6)
        - (21 * delta**2 - 3 * delta) * parent.PRESSURE / 2
        - (1 - 6 * delta) * (parent.RHO + parent.PRESSURE) / 2
    )
    Jbounce = s.factor(Jc.subs(u, 0))
    whole_current = parent.reference_data()["whole_reference_coefficients"]
    old_f, old_j, old_pairs = uv.principal_pairs()
    old_central = old_pairs["central"]
    old_symbols = set().union(*(value.free_symbols for value in old_central.values()))
    bridge = {old_f: F, old_j: J, phase.a: a}
    bridge.update(
        {
            symbol: E if symbol.name == "nonzero_E" else ell
            for symbol in old_symbols
            if symbol.name in ("nonzero_E", "ell")
        }
    )
    T = s.Symbol("whole_current_Tcorr", real=True)
    p_b, p_s = s.symbols(
        "clean_central_momentum_b clean_central_momentum_sigma", real=True
    )
    exact_v = p_b / (2 * s.sqrt(kap) * a**3 * q)
    exact_n = (ell * E * p_s - (E + 3 * T / (2 * q)) * p_b) / (
        2 * J * s.sqrt(kap) * a**3
    )
    Lambda, k0 = s.symbols("UV_endpoint positive_fixed_lower_endpoint", positive=True)
    vv_log = s.factor(vv * k**3 / (2 * s.pi**2))
    vn_quad = s.factor(vn * k / (4 * s.pi**2))
    nn_quartic = s.factor(nn / k / (8 * s.pi**2))
    return {
        "whole_exact_bounce_central_boundary_B": Bbounce,
        "whole_current_central_principal_K": K0,
        "whole_current_central_principal_G": G0,
        "whole_current_principal_diagonalizing_map": R,
        "whole_exact_bounce_metric_coordinate": exact_v,
        "whole_exact_bounce_lapse_coordinate_including_Tcorr": exact_n,
        "whole_leading_bounce_metric_momentum_row": vrow,
        "whole_leading_bounce_lapse_momentum_row": nrow,
        "whole_leading_fixed_reference_momentum_covariance": pi_cov,
        "whole_fixed_reference_bounce_Cvv_principal": vv,
        "whole_fixed_reference_bounce_Cvn_principal": vn,
        "whole_fixed_reference_bounce_Cnn_principal": nn,
        "whole_actual_bounce_J_with_unchanged_profiles": Jbounce,
        "whole_actual_bounce_F": s.factor(Fc.subs(u, 0)),
        "whole_actual_bounce_E": Ec.subs(u, 0),
        "whole_actual_bounce_charge": parent.ell.subs(u, 0),
        "whole_R3_Cvv_log_coefficient": vv_log,
        "whole_R3_Cvn_Lambda_squared_coefficient": vn_quad,
        "whole_R3_Cnn_Lambda_fourth_coefficient": nn_quartic,
        "whole_gauge_mean_v_log_coefficient": -s.Rational(9, 4) * vv_log,
        "whole_gauge_mean_N_Lambda_squared_coefficient": -s.Rational(9, 4) * vn_quad,
        "ultraviolet_argument": "The unchanged S251 two-cone classical-symbol/WKB proof gives pi-pi covariance diag(k cs/a,k/a)/2+O(1) in the normalized cleaned central frame. At H=Theta=0 its ENTIRE symmetric boundary B is zero; exact Qv=Pb/(2a^3 q), including the kappa field factor, and the entire lapse constraint give the stated symbols. Thus Cvv=A/k^3+O(k^-4), Cvn=B/k+O(k^-2), Cnn=D k+O(1). The R3 radial measure is k^2 dk/(2pi^2). For a fixed positive lower endpoint the respective integrals have log Lambda+O(1), Lambda^2+O(Lambda), and Lambda^4+O(Lambda^3) behavior. These are asymptotic unrenormalized coordinate covariances, not an actual EFT cutoff or quantitative interacting loop estimate.",
        "checks": {
            "literal_entire_fixed_current_J_coefficient": s.factor(
                whole_current["Jc"] - Jc
            ),
            "literal_existing_central_principal_R": (
                old_central["R"].subs(bridge, simultaneous=True) - R
            ).applyfunc(s.factor),
            "literal_existing_central_principal_K": (
                old_central["K"].subs(bridge, simultaneous=True) - K0
            ).applyfunc(s.factor),
            "literal_existing_central_principal_G": (
                old_central["G"].subs(bridge, simultaneous=True) - G0
            ).applyfunc(s.factor),
            "whole_principal_canonical_K_diagonalization": (
                R.T * K0 * R - s.eye(2)
            ).applyfunc(s.factor),
            "whole_principal_two_speeds": (R.T * G0 * R - s.diag(F / J, 1)).applyfunc(
                s.factor
            ),
            "whole_bounce_boundary_vanishes_without_dropping_Tcorr": Bbounce,
            "whole_reference_metric_cubic_decay": s.factor(
                vv - (2 * J * cs / E**2 + ell**2) / (8 * kap * k**3)
            ),
            "whole_reference_lapse_cross_linear_decay": s.factor(
                vn + cs / (4 * E * kap * a * a * k)
            ),
            "whole_reference_lapse_variance_linear_growth": s.factor(
                nn - cs * k / (4 * J * kap * a**4)
            ),
            "whole_entire_lapse_Tcorr_contact": s.factor(
                exact_n
                - (ell * E * p_s / (s.sqrt(kap) * a**3) - (2 * E * q + 3 * T) * exact_v)
                / (2 * J)
            ),
            "actual_current_bounce_F": s.factor(Fc.subs(u, 0) - s.Rational(1199, 800)),
            "actual_current_bounce_J": s.factor(
                Jbounce
                - s.Rational(243, 160)
                - parent.RHO.subs(u, 0)
                + s.Rational(7, 8) * parent.PRESSURE.subs(u, 0)
            ),
            "whole_radial_Cvv_coefficient": s.factor(
                s.integrate(k * k * vv / (2 * s.pi**2), (k, k0, Lambda))
                - vv_log * s.log(Lambda / k0)
            ),
            "whole_radial_Cvn_coefficient": s.factor(
                s.integrate(k * k * vn / (2 * s.pi**2), (k, k0, Lambda))
                - vn_quad * (Lambda * Lambda - k0 * k0)
            ),
            "whole_radial_Cnn_coefficient": s.factor(
                s.integrate(k * k * nn / (2 * s.pi**2), (k, k0, Lambda))
                - nn_quartic * (Lambda**4 - k0**4)
            ),
        },
        "gates": {
            "full_bounce_Tcorr_source_contact_kept": exact_n.has(T),
            "fixed_pressure_and_energy_profiles_both_kept": Jbounce.has(
                parent.PRESSURE.subs(u, 0)
            )
            and Jbounce.has(parent.RHO.subs(u, 0)),
            "whole_matter_momentum_mix_kept_in_metric_symbol": vv.has(ell),
            "continuum_Cvv_has_positive_log_coefficient": True,
            "actual_bounce_Cvn_positive_E_negative": True,
            "full_current_symbol_not_a_cutoff_bound_or_quantum_mean": True,
        },
    }
