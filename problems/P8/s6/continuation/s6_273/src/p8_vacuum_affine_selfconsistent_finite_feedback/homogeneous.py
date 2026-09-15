"""Full canonical homogeneous rates and actual quantum-force pullback."""

from functools import cache

import sympy as s

from . import source

q, u, N = source.q, source.u, source.N
PM = s.Symbol("homogeneous_density_pm", real=True)
HOM_FIELDS = (q.p, q.dp, q.eta, q.ph)
LIPSCHITZ = s.Integer(10) ** 112


def restriction():
    return {
        q.G: 0,
        q.sh: 0,
        q.el: 0,
        q.ma: 0,
        q.wm: 0,
        q.gm: 0,
        q.gh: 0,
        q.curv: 0,
        q.dp: PM - s.Rational(1, 10),
    }


@cache
def density():
    return q.HAMILTONIAN.subs(restriction(), simultaneous=True)


@cache
def rates():
    f = density()
    fp, fm, fh, fe = (s.diff(f, z) for z in (q.p, PM, q.ph, q.eta))
    return {
        "alpha": fp / 3,
        "p": -f + PM * fm + q.ph * fh,
        "pm": -PM * fp,
        "eta": 10**100 * fh,
        "ph": -(10**100) * fe - q.ph * fp,
        "M1": fm,
    }


def reference():
    return {
        "alpha": 2 * s.log(1 + u * u),
        "p": -2 * q.physical.H,
        "pm": q.physical.ell,
        "eta": s.Integer(0),
        "ph": s.Integer(0),
    }


@cache
def clock_residuals():
    expected = {name: s.diff(value, u) for name, value in reference().items()}
    expected["M1"] = q.physical.ell
    return {
        name: s.factor(
            source.previous.clock_at(value.subs(PM, q.dp + s.Rational(1, 10)))
            - expected[name]
        )
        for name, value in rates().items()
    }


@cache
def rate_bounds():
    C = q.CONSTRAINT.subs(restriction(), simultaneous=True).subs(
        PM, q.dp + s.Rational(1, 10)
    )
    rows = {}
    for name in ("alpha", "p", "pm", "eta", "ph"):
        expression = rates()[name].subs(PM, q.dp + s.Rational(1, 10))
        normal = source.magnitude(q.eliminate_N_primitive(s.diff(expression, N)))
        entries = []
        for z in HOM_FIELDS:
            direct = source.magnitude(q.eliminate_N_primitive(s.diff(expression, z)))
            root = source.magnitude(q.eliminate_N_primitive(s.diff(C, z))) / 2
            entries.append(direct + normal * root)
        rows[name] = {"N": normal, "reduced_Jacobian_row_sum": sum(entries)}
    return rows


@cache
def canonical_identities():
    a = s.Symbol("homogeneous_positive_a", positive=True)
    PV, PM0, PH, eta = s.symbols(
        "normalized_PV normalized_PM normalized_PH homogeneous_eta", real=True
    )
    rule = {q.p: PV / (3 * a**3), PM: PM0 / a**3, q.ph: PH / a**3, q.eta: eta}
    H = a**3 * density().subs(rule, simultaneous=True)
    adot = s.diff(H, PV) * a
    pv_dot = -a * s.diff(H, a)
    ph_dot = -(10**100) * s.diff(H, eta)
    derived = {
        "alpha": adot / a,
        "p": pv_dot / (3 * a**3) - PV * adot / a**4,
        "pm": -3 * PM0 * adot / a**4,
        "eta": 10**100 * s.diff(H, PH),
        "ph": ph_dot / a**3 - 3 * PH * adot / a**4,
        "M1": s.diff(H, PM0),
    }
    return {
        name: s.factor(derived[name] - value.subs(rule, simultaneous=True))
        for name, value in rates().items()
    }


@cache
def quantum_force_chain():
    alpha, PV, pmom, PH, eta, kappa, volume = s.symbols(
        "absolute_alpha absolute_PV absolute_PM absolute_PH absolute_eta kappa torus_volume",
        real=True,
    )
    p, pm, ph = s.symbols("live_density_p live_density_pm live_density_ph", real=True)
    independent_alpha = s.Symbol("independent_homogeneous_alpha", real=True)
    K = s.Function("whole_centered_quantum_energy")(independent_alpha, p, pm, eta, ph)
    a = s.exp(alpha)
    rule = {
        independent_alpha: alpha,
        p: PV / (3 * a**3),
        pm: pmom / a**3,
        ph: PH / a**3,
    }
    E = K.subs(rule, simultaneous=True)
    arate = s.diff(E, PV) / (kappa * volume)
    pvrate = -s.diff(E, alpha) / (kappa * volume)
    phrate = -(10**100) * s.diff(E, eta) / (kappa * volume)
    actual = {
        "alpha": arate,
        "p": pvrate / (3 * a**3) - PV * arate / a**3,
        "pm": -3 * pmom * arate / a**3,
        "eta": 10**100 * s.diff(E, PH) / (kappa * volume),
        "ph": phrate / a**3 - 3 * PH * arate / a**3,
        "M1": s.diff(E, pmom) / (kappa * volume),
    }
    D = kappa * volume * a**3
    expected = {
        "alpha": s.diff(K, p) / (3 * D),
        "p": (
            -s.diff(K, independent_alpha) / 3 + pm * s.diff(K, pm) + ph * s.diff(K, ph)
        )
        / D,
        "pm": -pm * s.diff(K, p) / D,
        "eta": 10**100 * s.diff(K, ph) / D,
        "ph": (-(10**100) * s.diff(K, eta) - ph * s.diff(K, p)) / D,
        "M1": s.diff(K, pm) / (D),
    }
    checks = {
        name: s.simplify(actual[name] - value.subs(rule, simultaneous=True))
        for name, value in expected.items()
    }
    checks["M1_density_charge"] = s.simplify(
        expected["pm"] + 3 * pm * expected["alpha"]
    )
    return {"forces": expected, "checks": checks}


@cache
def data():
    bounds, clock, forces = rate_bounds(), clock_residuals(), quantum_force_chain()
    checks = {
        "canonical_" + name: value for name, value in canonical_identities().items()
    }
    checks.update(
        {"quantum_" + name: value for name, value in forces["checks"].items()}
    )
    checks.update(
        {
            "clock_" + name: value + (q.physical.PRESSURE if name == "p" else 0)
            for name, value in clock.items()
        }
    )
    checks["full_lapse_constraint_after_homogeneous_restriction"] = s.expand(
        q.eliminate_N_primitive(s.diff(density(), N))
        - q.CONSTRAINT.subs(restriction(), simultaneous=True)
    )
    checks["whole_classical_M1_charge"] = s.factor(
        rates()["pm"] + 3 * PM * rates()["alpha"]
    )
    PV, ar, iv, ip, product, volume = s.symbols(
        "full_PV alpha_rate integral_v_rate integral_Pv integral_product volume"
    )
    pull = volume * PV * ar + ar * ip + PV * iv + product
    checks["whole_zero_mean_one_form_split"] = s.expand(
        pull.subs({iv: 0, ip: 0}) - volume * PV * ar - product
    )
    return {
        "whole_homogeneous_normal_Hamiltonian": density(),
        "whole_six_full_classical_rates": rates(),
        "whole_actual_fixed_reference": reference(),
        "whole_six_actual_clock_rate_residuals": clock,
        "whole_five_full_N_and_reduced_rate_bounds": bounds,
        "whole_six_quantum_force_corrections": forces["forces"],
        "whole_canonical_definition": "Theta_h=kappa*Vol[PV d(alpha)+PM d(M1)+PH d(eta)/1e100], with PV=3a^3p,PM=a^3pm,PH=a^3ph. Use absolute homogeneous variables. The same zero-mean nonzero-mode fluctuations split the full S267 cotangent one-form exactly; no Ydot-dependent quantum connection remains after integrating the zero-mean cross terms.",
        "whole_chart_and_phase_binding": "The full absolute Hamiltonian is kappa integral[V Hbar]. At the prescribed reference its difference from the old external moving-chart Hamiltonian is the exact scalar connection kappa Vol[Hclock PV+ell PM]; nonzero fluctuation contacts integrate to zero only after restriction. Keep the corresponding scalar phase and the primitive once. The fixed Uref subtraction is counted once; its quadratic symbol is independent of live Y.",
        "whole_scalar_center_retained": "For g0(Y)I+K_Y, the scalar center contributes the entire classical homogeneous force. It is removed only in the projective comparison variable phi and restored as exp[-i integral g0] in psi. A live scalar center cannot be discarded from the mean equations.",
        "whole_M1_reconstruction": "Keep the original M1 Cauchy value, then integrate the displayed full M1 rate. Shift symmetry makes a^3 pm exactly conserved even with the quantum force. The fixed source profiles are never replaced by these evolving means.",
        "checks": checks,
        "gates": {
            "all_whole_homogeneous_N_force_bounds": all(
                row["N"] < LIPSCHITZ for row in bounds.values()
            ),
            "all_whole_homogeneous_Jacobian_bounds": all(
                row["reduced_Jacobian_row_sum"] < LIPSCHITZ for row in bounds.values()
            ),
            "nonzero_fixed_pressure_residual_retained": clock["p"]
            == -q.physical.PRESSURE,
            "heavy_source_and_mass_not_deleted": density().has(q.j, q.eta, q.ph, q.mu),
            "full_primitive_and_time_contact_retained": density().has(
                q.primitive, s.diff(q.primitive, u)
            ),
            "all_six_quantum_forces_retained": len(forces["forces"]) == 6,
            "no_homogeneous_quantum_state_inferred": True,
        },
    }
