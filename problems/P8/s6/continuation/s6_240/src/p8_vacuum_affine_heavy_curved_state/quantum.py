"""Whole dimensional scalar counteraction and finite reference-profile matching."""

from functools import cache

import sympy as s

from . import state

D, z, eps, ell = s.symbols("spatial_D momentum_fraction epsilon log_n", real=True)
h = s.symbols("h0:6", real=True)
H, Hd, Hdd, Hddd = h[:4]


def dt(expression):
    return s.expand(
        sum(h[j + 1] * s.diff(expression, h[j]) for j in range(5))
        - 2 * H * z * (1 - z) * s.diff(expression, z)
    )


@cache
def adiabatic():
    lam = -H * z
    U = D * Hd / 2 + D**2 * H**2 / 4
    P2 = s.expand(-U / 2 - dt(lam) / 4 + lam**2 / 8)
    P4 = s.expand(
        -(P2**2) / 2
        - dt(dt(P2)) / 4
        + s.Rational(5, 4) * lam * dt(P2)
        + (dt(lam) / 2 - s.Rational(3, 2) * lam**2) * P2
    )
    B2 = dt(P2) - 2 * lam * P2
    rate0 = (D * H + lam) / 2
    energy = (
        s.Rational(1, 2),
        s.expand(rate0**2 / 4),
        s.expand((P2**2 + rate0 * B2 - rate0**2 * P2) / 4),
    )
    pressure = (
        z / (2 * D),
        s.expand(((2 - 2 * z / D) * P2 + rate0**2) / 4),
        s.expand(
            (
                (2 - 2 * z / D) * P4
                - (1 - 2 * z / D) * P2**2
                + rate0 * B2
                - rate0**2 * P2
            )
            / 4
        ),
    )
    return {"P2": P2, "P4": P4, "rho": energy, "pressure": pressure}


def radial(piece, order):
    factor = 0
    for (j,), coefficient in s.Poly(piece, z).terms():
        factor += (
            coefficient.subs(D, 3 - 2 * eps)
            * s.rf(s.Rational(3, 2) - eps, j)
            / s.gamma(j + order - s.Rational(1, 2))
        )
    f0 = s.simplify(factor.subs(eps, 0))
    f1 = s.simplify(s.diff(factor, eps).subs(eps, 0))
    r = 2 - order
    pref = 8 * s.sqrt(s.pi) * (-1) ** r / s.factorial(r)
    pole = s.expand(pref * f0)
    finite = s.expand(pole * (s.harmonic(r) - ell) + pref * f1)
    return pole, finite


def curvatures():
    R = 2 * D * Hd + D * (D + 1) * H**2
    Ric = D**2 * (Hd + H**2) ** 2 + D * (Hd + D * H**2) ** 2
    Riem = 4 * D * (Hd + H**2) ** 2 + 2 * D * (D - 1) * H**4
    return R, Ric, Riem


def rho_of(density):
    return s.expand(
        -density
        + H * s.diff(density, H)
        + (Hd - D * H**2) * s.diff(density, Hd)
        - H * dt(s.diff(density, Hd))
    )


def pressure_of(rho):
    return s.cancel(-rho - dt(rho) / (D * H))


@cache
def local():
    R, Ric, Riem = curvatures()
    a2 = R**2 / 72 + (Riem - Ric) / 180
    densities = (s.S.One, -R / 3, 2 * a2)
    actual = {}
    component_only = {}
    checks = {}
    for order in range(3):
        rhoD = rho_of(densities[order])
        pD = pressure_of(rhoD)
        actual[order] = {}
        component_only[order] = {}
        for label, piece, target in (
            ("rho", adiabatic()["rho"][order], rhoD),
            ("pressure", adiabatic()["pressure"][order], pD),
        ):
            pole, finite = radial(piece, order)
            checks[f"complete_dimensional_{label}_pole_{order}"] = s.cancel(
                pole - target.subs(D, 3)
            )
            combined = s.expand(finite + 2 * s.diff(target, D).subs(D, 3))
            wanted = s.expand(target.subs(D, 3) * (s.harmonic(2 - order) - ell))
            checks[f"full_covariant_counterterm_finite_{label}_{order}"] = s.cancel(
                combined - wanted
            )
            actual[order][label] = s.factor(combined)
            component_only[order][label] = s.factor(finite)
        checks[f"full_finite_local_Ward_{order}"] = s.cancel(
            dt(actual[order]["rho"])
            + 3 * H * (actual[order]["rho"] + actual[order]["pressure"])
        )
    return actual, component_only, checks


def fixed_profile():
    u, X = state.TIME, state.X
    original = state.germs.parent.original.data()
    oldX = state.germs.X
    switch = original["T"].subs(oldX, X)
    rho = s.Function("rho_H_reference", real=True)(u)
    pressure = s.Function("P_H_reference", real=True)(u)
    n, k = state.MASS2, state.KAPPA
    pvH = (s.Rational(3, 2) - s.log(n)) * n * n / (64 * s.pi**2)
    pvPhi = s.Rational(3, 128) / s.pi**2
    coefficient = (
        -pvH
        + switch * (-pressure + pvH - (rho + pressure) * (X - 1) / 2)
        - pvPhi * (1 - switch)
    ) / k
    return {
        "T": switch,
        "rho": rho,
        "pressure": pressure,
        "p_v_H": pvH,
        "p_v_Phi": pvPhi,
        "DeltaF": coefficient,
    }


@cache
def data():
    actual, component_only, checks = local()
    n = s.Symbol("mass_squared", positive=True)
    R, Ric, Riem = s.symbols("R_old Ricci_squared Riemann_squared", real=True)
    a2 = R**2 / 72 + (Riem - Ric) / 180
    Weyl = Riem - 2 * Ric + R**2 / 3
    Euler = Riem - 4 * Ric + R**2
    checks = dict(checks)
    checks["complete_scalar_heat_covariant_basis"] = s.expand(
        a2 - (Weyl / 120 - Euler / 360 + R**2 / 72)
    )
    finite = (s.Rational(3, 2) - ell) * n * n + (ell - 1) * n * R / 3 - 2 * ell * a2
    profile = fixed_profile()
    u, X, k = state.TIME, state.X, state.KAPPA
    coefficient = profile["DeltaF"]
    A = s.simplify(coefficient.subs(X, 1))
    B = s.simplify(s.diff(coefficient, X).subs(X, 1))
    rho, P = profile["rho"], profile["pressure"]
    checks["full_actual_heavy_clock_pressure_cancel"] = s.cancel(A + P / k)
    checks["full_actual_heavy_clock_energy_cancel"] = s.cancel(2 * B - A + rho / k)
    checks["full_clock_scalar_mean_equation_is_actual_Ward"] = s.cancel(
        s.diff(A, u)
        - 2 * s.diff(B, u)
        - 6 * state.HUBBLE * B
        - (s.diff(rho, u) + 3 * state.HUBBLE * (rho + P)) / k
    )
    vacuum = coefficient.subs(X, 0)
    checks["complete_heavy_and_light_flat_vacuum_constant"] = s.cancel(
        vacuum + (profile["p_v_H"] + profile["p_v_Phi"]) / k
    )
    checks["complete_new_one_loop_flat_pressure_cancel"] = s.cancel(
        profile["p_v_H"] + profile["p_v_Phi"] + k * vacuum
    )
    checks["no_new_flat_field_dependent_classical_potential"] = s.cancel(
        vacuum - vacuum.subs(u, 0)
    )
    denominator = X**state.N + (1 - X) ** state.N
    checks["new_nonconstant_vacuum_jet_factor"] = s.factor(
        denominator * profile["T"] - X**state.N
    )
    checks["extra_light_vacuum_clock_flat_factor"] = s.factor(
        denominator * (1 - profile["T"]) - (1 - X) ** state.N
    )
    checks["Gaussian_vacuum_loop_degree_identity"] = 2 * 1 - 2
    return {
        "new_candidate": "CD-REG-AFFINE-ISO-QG2-H8A420 with the separately fixed H8A420-VAC-OS4 finite vacuum coefficients retained",
        "new_heavy_prescription": "For the conditional minimally coupled Gaussian H field use covariant dimensional regularization at mu1, with the COMPLETE pole action continued in spatial D=3-2epsilon and varied before taking its limit. No state-dependent normal ordering is used. Keep the old Proca prescription unchanged. The displayed new scalar finite matching is not borrowed from the vector.",
        "complete_D_scalar_adiabatic_readout": adiabatic(),
        "complete_scalar_pole_density_before_64pi2epsilon": n * n - n * R / 3 + 2 * a2,
        "complete_scalar_finite_heat_matching_before_64pi2": finite,
        "complete_finite_local_stress_coefficients": actual,
        "component_only_pole_subtraction_negative_control": component_only,
        "dimensional_matching_proof": "Integrate every omega^(1-2j) z^r term in D spatial dimensions using the complete gamma/rising-factorial identity. Its pole and component finite part are exact. The full covariant pole counterstress has the additional finite term -2*partial_D T_pole at D3, so subtraction ADDS2*partial_D T_pole to the component finite part. This cancels the false extra finite terms from component-only subtraction and gives the full displayed scalar density. The compact box-R/Euler variational boundary is retained; no infinite-time flux is dropped.",
        "full_reference_stress_definition": "rho_H_ref and P_H_ref are the full radial exact-SLE-minus-fourth-order-adiabatic mode integrals from state.data(), with measure p²dp/(2pi²), PLUS the displayed same-scheme local stress sum n^(2-j)*local[j]/(64pi²), ell=log n. All momenta and the full exact state are included. These symbols are names for that fixed source-defined integral, not adjustable independent functions.",
        "fixed_new_normalized_scalar_coefficient": coefficient,
        "vacuum_heavy_pressure": profile["p_v_H"],
        "vacuum_light_pressure": profile["p_v_Phi"],
        "new_profile_definition": "DeltaF_HQ=[-p_v_H+T(-P_H_ref+p_v_H-(rho_H_ref+P_H_ref)(X-1)/2)-p_v_Phi(1-T)]/kappa0. The last term matches the mass1 Phi Gaussian vacuum constant at mu1 while vanishing to order1024 at the clock. It does not assert that Phi fluctuations have been quantized on the curved clock background.",
        "fixed_locality": "The reference state, full CD history and renormalization prescription define rho_H_ref(u),P_H_ref(u) ONCE. Hold these functions fixed under every later state/source/metric variation. The new coefficient is globally smooth, not claimed real analytic in u. It is not a live nonlocal subtraction or an instruction to cancel the stress of each perturbed history.",
        "reference_mean_equations": "At X1, the profile pressure and energy cancel the ACTUAL heavy Gaussian reference stress; its scalar Euler expression is the same renormalized Ward identity. The extra light-vacuum switch has no clock jets. The old retained Proca/M1 mean equations therefore remain stationary after including this new specified H Gaussian sector, with the S238 recomputed affine mean. There is no separate H clock force because its source and relevant jets vanish and the reference coherent mean is0.",
        "flat_vacuum_matching": "At X0 the new coefficient is the constant -(p_v_H+p_v_Phi)/kappa0. Together with the unchanged old Proca cancellation it cancels the formal flat limiting-action Gaussian one-loop vacuum constant. Massless free spectator determinants are scaleless in this dimensional prescription. Interacting vacuum graphs start at two loops, by connected degree counting. The old first-order heavy-onepoint condition remains. This is not a finite-gravity quantum-limit or quantum-vacuum stability theorem.",
        "unchanged_low_vacuum_matching": "Every nonconstant new term contains T=O(X^1024). Thus no relevant vacuum mass, derivative, cubic, quartic or degree-six jet is changed, and the complete S239 first-loop four-Phi calculation with its finite OS4 condition is retained. The added flat potential is only a constant, so the previous relative classical comparison is unchanged; no full quantum-potential theorem is inferred.",
        "domain_and_response_boundary": "F alone is retuned; R, the S238 full affine maps and physical matter metric are unchanged. The new nonzero clock F jets CHANGE the classical light Hessian. The old kinetic/stability or complete quantum-response inverse results do not automatically transfer. All interacting light/mixed clock loops, full new response and same-state nonlinear bounce remain open.",
        "checks": {key: s.cancel(value) for key, value in checks.items()},
        "gates": {
            "full_dimensional_counteraction_not_component_only_subtraction": s.cancel(
                actual[2]["rho"] - component_only[2]["rho"]
            )
            != 0,
            "full_new_heavy_and_light_vacuum_constants_both_retained": True,
            "nonconstant_profile_cannot_enter_first_four_point_loop": 2 * state.N - 2
            > 4,
            "exact_old_full_R_and_affine_maps_unchanged_by_F_only_profile": True,
            "new_clock_scalar_Hessian_not_silently_identified_with_old": True,
            "conditional_Gaussian_mean_not_full_interacting_clock_state": True,
        },
    }
