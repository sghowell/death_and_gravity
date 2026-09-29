"""Adopted S0: canonical domain and source-free Gaussian reference bridge.

This checks the new candidate, not a rewrite of the old sourced parent.
No full light/metric covariant measure, causal response norm, physical
matching result, gravity verdict or P8 closure is inferred.
"""

import copy
import hashlib
import importlib.util
import json

import sympy as s
from p8_match1_rate_input import REPO
from p8_source1 import inputs as previous_inputs

SPEC = "docs/candidates/p8-rate4-covzero-s0.json"
OLD_SPEC = "docs/candidates/p8-rate4-covzero-1.json"
COUPLED = "problems/P8/s6/continuation/s6_253/src/p8_vacuum_affine_physical_background_vertices/coupled.py"
CLOCK = "problems/P8/s6/continuation/s6_174/certificates/polynomial-vacuum-analytic-affine-parent.json"
VECTOR = "problems/P8/s6/continuation/s6_176/certificates/polynomial-vacuum-affine-proca-gaussian.json"
EXTRA = {
    "scripts/p8_source1.py": "dbe3e7d8a0f7892289c832d60ffe19ae782df6e40ff434925e5c1251e287abfb",
    "docs/assessment-2026-09-28-p8-source1-feasibility.md": "4d651b907b14a67464990033cf1e48a9a4f951f935e792b92adf56a6f3e3d255",
    "docs/candidates/p8-source-1-study.json": "8095e223de1a734daadd73bad39f34d38f2cd5a52aa229fbe81be71022ce69e2",
    "docs/validation/p8-source1-2026-09-28.json": "f6869c26a6fcb0d1425512d27a2bfd6da660081f00275db265a518e5a4da080a",
}


def inputs():
    manifest = previous_inputs()
    for name, expected in EXTRA.items():
        assert hashlib.sha256((REPO / name).read_bytes()).hexdigest() == expected
        assert name not in manifest or manifest[name] == expected
        manifest[name] = expected
    assert all(name in manifest for name in (COUPLED, CLOCK, VECTOR, OLD_SPEC))
    return manifest


def validate_spec(spec):
    expected = {
        "candidate": "QG2-H8A420-RATE4-COVZERO-S0 v1",
        "source_change_adopted": True,
        "physical_boundary_rule_adopted": True,
        "original_parent_files_modified": False,
        "common_quantum_reference_complete": False,
        "physical_matching_established": False,
        "physical_gravity_verdict_established": False,
        "original_P8_open": True,
        "order": "first quantum order only",
    }
    if any(
        spec.get(k) != v or type(spec.get(k)) is not type(v)
        for k, v in expected.items()
    ):
        raise ValueError("Keep adoption separate from completion")
    old = json.loads((REPO / OLD_SPEC).read_text())
    if spec["additional_finite_functional"] != old["additional_finite_functional"]:
        raise ValueError("Do not reset the adopted finite boundary")
    for key in (
        "physical_metric_signature",
        "clock_norm",
        "states",
        "matching_targets",
        "detector",
        "physical_cutoff",
    ):
        if spec["retained"][key] != old["retained"][key]:
            raise ValueError("Retained physical data changed")
    new = spec["classical_change"]
    if (
        new["source"] != "0"
        or new["source_in_regulator_dimension"]
        != "0 at every d in the selected continuation"
    ):
        raise ValueError("The entire source is zero, including its continuation")
    if new["zeta"] != "1/1000000" or new["canonical_vector_mass"] != "1000":
        raise ValueError("The vector remains massive with its old normalization")
    if new["new_vector_action"] != "kappa*sqrt(-g_phys)*[-zeta*F(W)^2/4+W^2/2]":
        raise ValueError("Do not retain a partial source coupling")
    if spec["reference_seed"]["source_in_hat_variables"] != "0":
        raise ValueError("Zero source is not R=1")
    for key in (
        "dimension",
        "mu",
        "A1_A2",
        "A3",
        "A4",
        "A5",
        "metric_map_C",
        "proca_legacy_transport",
    ):
        if spec["reference_seed"][key] != old["reference_seed"][key]:
            raise ValueError("Retain the approved continuation and legacy transport")


def full_trace_and_clock_domain():
    N = s.Symbol("N", positive=True)
    R = s.Function("R", positive=True)(N)
    B, F, E = (s.Function(name)(N) for name in ("B", "Fhat", "E_other"))
    U, M = R ** -s.Rational(3, 4), R ** s.Rational(1, 4)
    aa = -M / 3
    K, T, p, G = s.symbols("K T p G")
    lag = aa * K**2 + B * K + F + U * T**2 / 2
    Kstar = (p - B) / (2 * aa)
    raw = N * ((p * K - lag).subs(K, Kstar) - T * G + E)
    Tstar = -G / U
    reduced = N * ((p - B) ** 2 / (4 * aa) - F + G**2 / (2 * U) + E)
    assert s.simplify(s.diff(lag, K).subs(K, Kstar) - p) == 0
    assert s.simplify(raw.subs(T, Tstar) - reduced) == 0
    assert s.simplify(s.diff(raw, T).subs(T, Tstar)) == 0
    D = s.hessian(raw, (N, T))
    assert D[1, 1] == -N * U
    assert s.simplify(D[0, 1].subs(T, Tstar) + D[1, 1] * s.diff(Tstar, N)) == 0
    schur = (D[0, 0] - D[0, 1] ** 2 / D[1, 1]).subs(T, Tstar)
    assert s.simplify(schur - s.diff(reduced, N, 2)) == 0
    assert s.diff(s.diff(reduced, N, 2), s.diff(E, N, 2)) == N
    clock = json.loads((REPO / CLOCK).read_text())[
        "actual_clock_vacuum_and_constraint_blocks"
    ]["clock"]
    u = s.Symbol("u", real=True)
    J = s.sympify(clock["J"], locals={"u": u})
    polynomial = s.Poly(s.cancel(800 * (1 + u**2) ** 18 * J), u)
    assert polynomial.degree() == 34
    assert all(c > 0 for c in polynomial.all_coeffs() if c != 0)
    assert polynomial.nth(0) == 1215
    assert J.subs(u, 0) == s.Rational(243, 160)
    # For every |u|<=Tmax, positivity of the numerator supplies this bound.
    Tmax = s.Symbol("Tmax", nonnegative=True)
    bound = s.Rational(243, 160) / (1 + Tmax**2) ** 18
    return {
        "T_star": str(Tstar),
        "temporal_pivot": "-N*R^(-3/4)",
        "lapse_pivot": "d_N^2 H_red at fixed canonical data and spatial jets",
        "compact_clock_J_lower_bound": str(bound),
        "all_time_uniform_radius_proved": False,
        "vacuum_unitary_extension_used": False,
    }


def gaussian_source_free_reduction():
    # This pinned module imports only functools and sympy and defines algebra.
    # Do not invoke its ancestor acceptance, data(), or cached report builders.
    module_spec = importlib.util.spec_from_file_location(
        "_p8_s0_pinned_coupled", REPO / COUPLED
    )
    c = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(c)
    lag = (
        c.base_lagrangian()
        + c.K * (c.Ad - c.P * c.A0) ** 2 / 2
        + c.Z * c.A0**2 / 2
        - c.Yv * c.longitudinal**2 / 2
    )
    rates, momenta = s.Matrix([c.vd, c.sd, c.Ad]), s.Matrix([c.pv, c.ps, c.pA])
    kinetic = s.hessian(lag, rates)
    assert kinetic == s.diag(-6 * c.D, c.Z, c.K)
    zero_rates = dict.fromkeys(rates, 0)
    linear = s.Matrix([s.diff(lag, x).subs(zero_rates) for x in rates])
    solution = kinetic.inv() * (momenta - linear)
    raw = s.expand(
        (momenta.dot(rates) - lag).subs(
            dict(zip(rates, solution, strict=True)), simultaneous=True
        )
    )
    auxiliary = (c.n, c.b, c.A0)
    Daux = s.hessian(raw, auxiliary).applyfunc(s.cancel)
    assert Daux == s.diag(-2 * c.J, -2 * c.D / 3, -c.Z)
    aux_solution = s.solve([s.diff(raw, x) for x in auxiliary], auxiliary)
    reduced = s.cancel(raw.subs(aux_solution, simultaneous=True))
    inherited = c.reduced_hamiltonian().subs({c.r: 0, c.rN: 0})
    assert s.cancel(reduced - inherited) == 0
    H = s.hessian(reduced, c.PHASE)
    light_indices, vector_indices = [0, 1, 3, 4], [2, 5]
    assert H.extract(light_indices, vector_indices) == s.zeros(4, 2)
    expected_vector = s.diag(c.Yv, 1 / c.K + c.P**2 / c.Z)
    assert (H.extract(vector_indices, vector_indices) - expected_vector).applyfunc(
        s.cancel
    ) == s.zeros(2)
    assert s.factor(Daux.det()) == -4 * c.D * c.Z * c.J / 3
    assert s.factor(kinetic.det() * Daux.det()) == 8 * c.D**2 * c.Z**2 * c.J * c.K
    # Full lower constraint bracket, not a fictitious zero lower block.
    constraints = s.Matrix([s.diff(raw, x) for x in auxiliary])
    gradients = constraints.jacobian(c.PHASE)
    lower = (gradients * c.OMEGA * gradients.T).applyfunc(s.cancel)
    assert lower != s.zeros(3)
    zero = s.zeros(3)
    dirac = zero.row_join(-Daux).col_join(Daux.T.row_join(lower))
    invD = Daux.inv()
    inverse = (invD.T * lower * invD).row_join(invD.T).col_join((-invD).row_join(zero))
    assert (dirac * inverse - s.eye(6)).applyfunc(s.cancel) == s.zeros(6)
    assert s.factor(dirac.det() - Daux.det() ** 2) == 0
    physical_constraint_brackets = s.zeros(6, 3).row_join(c.OMEGA * gradients.T)
    correction = physical_constraint_brackets * inverse * physical_constraint_brackets.T
    assert correction.applyfunc(s.cancel) == s.zeros(6)
    # At the common finite regulator the same canonical source-free measure
    # follows after integrating delta(p_aux) delta(chi) |det(Daux)|.
    # Track the entire old source difference, including its local contact.
    numerator = c.th * (c.pv + 3 * c.c * c.sigma) / c.D - c.w * c.ps / c.Z - c.Lnv * c.v
    source_mix = c.r * c.th / c.D + c.rN * c.dH
    old_minus_new = c.P * c.pA * (
        c.r * (c.pv + 3 * c.c * c.sigma) / (2 * c.D)
        - 3 * numerator * source_mix / (2 * c.J)
    ) + c.P**2 * c.pA**2 * (9 * source_mix**2 / (4 * c.J) - 3 * c.r**2 / (4 * c.D))
    assert s.cancel(c.reduced_hamiltonian() - reduced - old_minus_new) == 0
    epsilon, alpha, beta = s.symbols("epsilon alpha beta")
    source_hessian = s.hessian(old_minus_new, c.PHASE).subs(
        {c.r: epsilon * alpha, c.dH: epsilon * beta}, simultaneous=True
    )
    first = source_hessian.diff(epsilon).subs(epsilon, 0)
    second = source_hessian.diff(epsilon, 2).subs(epsilon, 0)
    assert first.extract(light_indices, light_indices) == s.zeros(4)
    assert first.extract(vector_indices, vector_indices) == s.zeros(2)
    assert first.extract(light_indices, vector_indices) != s.zeros(4, 2)
    assert second.extract(light_indices, light_indices) == s.zeros(4)
    assert second.extract(vector_indices, vector_indices) != s.zeros(2)
    return {
        "coupled_phase_split": "four-dimensional coupled light/M1 plus two-dimensional longitudinal Proca for all admitted homogeneous background jets",
        "auxiliary_determinant": "-4*D*Z*J/3",
        "configuration_density_squared": "8*D^2*Z^2*J*K",
        "nonzero_lower_constraint_bracket_retained": True,
        "reduced_physical_Dirac_bracket_is_canonical": True,
        "old_source_first_probe_vertex_only_off_diagonal": True,
        "old_source_second_contact_not_zero": True,
        "common_finite_regulator_metric_light_mean_preserved": True,
        "source_off_does_not_set_R_or_its_derivatives_to_constants": True,
    }


def physical_vector_vertices():
    R, N, a, k, zeta = s.symbols("R N a k zeta", positive=True)
    ahat = R ** s.Rational(1, 4) * a
    U, Cchi = R ** -s.Rational(3, 4), R ** -s.Rational(1, 4)
    K, Z, Yv = zeta * Cchi / (N * ahat**2), U / N, N * Cchi / ahat**2
    HL = s.diag(ahat**3 * Yv / zeta, zeta * (1 / K + k**2 / Z) / ahat**3).applyfunc(
        s.simplify
    )
    expected_L = N * s.diag(a / zeta, 1 / a + zeta * k**2 / a**3)
    assert (HL - expected_L).applyfunc(s.cancel) == s.zeros(2)
    assert HL.diff(R) == s.zeros(2)
    HT = N * s.diag(a / zeta + k**2 / a, 1 / a)
    eta, v, e = s.symbols("eta v epsilon", real=True)
    rows = {}
    for name, matrix in (("L", HL), ("T", HT)):
        H0 = matrix.subs(N, 1)
        Hv = (a * H0.diff(a)).applyfunc(s.cancel)
        Hvv = (a * Hv.diff(a)).applyfunc(s.cancel)
        perturbed = matrix.subs(
            {N: 1 + e * eta, a: a * s.exp(e * v)}, simultaneous=True
        )
        assert (perturbed.diff(e).subs(e, 0) - eta * H0 - v * Hv).applyfunc(
            s.simplify
        ) == s.zeros(2)
        assert (
            perturbed.diff(e, 2).subs(e, 0) - 2 * eta * v * Hv - v**2 * Hvv
        ).applyfunc(s.simplify) == s.zeros(2)
        if name == "T":
            assert (Hvv - H0).applyfunc(s.cancel) == s.zeros(2)
        else:
            assert (Hvv - s.diag(a / zeta, 1 / a + 9 * zeta * k**2 / a**3)).applyfunc(
                s.cancel
            ) == s.zeros(2)
        rows[name] = {
            "H": str(H0),
            "H_v": str(Hv),
            "H_vv": str(Hvv),
            "H_eta_eta": "0",
            "H_eta_v": "H_v",
            "multiplicity": 1 if name == "L" else 2,
        }
    # The same endpoint shear and scaling act on the old state, not a reset.
    g, d = s.symbols("g d", positive=True)
    transform = s.Matrix([[1 / s.sqrt(g), 0], [-s.sqrt(g) * d, s.sqrt(g)]])
    omega = s.Matrix([[0, 1], [-1, 0]])
    assert transform.T * omega * transform == omega
    # gL=a/(1+zeta*q), gT=a; logarithmic derivative at fixed comoving k.
    Hubble = s.Symbol("Hubble")
    gL, gT = a / (1 + zeta * k**2 / a**2), a
    z = zeta * k**2 / (a**2 + zeta * k**2)
    assert (
        s.cancel(
            Hubble * a * s.diff(gL, a) / (2 * gL) - Hubble * (s.Rational(1, 2) + z)
        )
        == 0
    )
    assert Hubble * a * s.diff(gT, a) / (2 * gT) == Hubble / 2
    # Homogeneous metric Ward identity before the mode sum: the Hamiltonian
    # evolution term in E' vanishes; E'=Hubble*O_v=-3*Hubble*a^3*pressure.
    Q, Pi = s.symbols("Q Pi")
    phase = s.Matrix([Q, Pi])
    for matrix in (HL.subs(N, 1), HT.subs(N, 1)):
        energy = (phase.T * matrix * phase)[0] / 2
        evolution = omega * matrix * phase
        canonical_term = (
            s.diff(energy, Q) * evolution[0] + s.diff(energy, Pi) * evolution[1]
        )
        pressure = -s.diff(energy, a) / (3 * a**2)
        assert s.cancel(canonical_term) == 0
        assert (
            s.cancel(Hubble * a * s.diff(energy, a) + 3 * Hubble * a**3 * pressure) == 0
        )
    return rows


def proca_measure_and_legacy():
    controls = []
    # Finite derivative complexes; d1*d0=0 is the only determinant hypothesis.
    d0 = s.Matrix([1, 2, -1, 0])
    d1 = s.Matrix([[-2, 1, 0, 0], [1, 0, 1, 0], [0, 0, 0, 1]])
    assert d1 * d0 == s.zeros(3, 1)
    for mass2 in (s.Integer(1), s.Integer(7), s.Integer(10) ** 6):
        P = d1.T * d1 + mass2 * s.eye(4)
        K1 = P + d0 * d0.T
        K0 = d0.T * d0 + mass2 * s.eye(1)
        assert P * d0 == mass2 * d0
        assert P.det() * K0.det() == K1.det() * mass2**d0.cols
        controls.append(
            {"mass_squared": str(mass2), "identity": "det(P)=det(K1)*m^(2*n0)/det(K0)"}
        )
    # Recheck the already derived COVZERO/S176 conversion, not a new loop result.
    d, Ro, Ric2, Riem2, mass, ell = s.symbols(
        "d R_old Ricci_squared Riemann_squared positive_mass log_mass_ratio"
    )
    scalar_heat = Ro**2 / 72 + (Riem2 - Ric2) / 180
    h0, h1 = d - 1, (d - 7) * Ro / 6
    h2 = (d - 1) * scalar_heat + Ric2 / 2 - Ro**2 / 6 - Riem2 / 12
    pole_d = mass**4 * h0 - 2 * mass**2 * h1 + 2 * h2
    pole = pole_d.subs(d, 4)
    evanescent = -2 * s.diff(pole_d, d)
    finite = (
        mass**4 * (s.Rational(3, 2) - ell) * h0
        - 2 * mass**2 * (1 - ell) * h1
        - 2 * ell * h2
    ).subs(d, 4) + evanescent
    legacy = json.loads((REPO / VECTOR).read_text())[
        "specified_vector_stress_and_reference_clock_response"
    ]["prescription"]
    assert (
        s.expand(pole - s.sympify(legacy["ordinary_pole_density_before_common_factor"]))
        == 0
    )
    assert (
        s.expand(
            finite
            - s.sympify(
                legacy["ordinary_finite_heat_matching_density_before_64_pi_squared"]
            )
        )
        == 0
    )
    assert (
        s.expand(
            finite.subs(ell, s.log(mass**2))
            + s.log(mass**2) * pole
            - finite.subs(ell, 0)
        )
        == 0
    )
    assert (
        s.expand(evanescent + 2 * (mass**4 - mass**2 * Ro / 3 + 2 * scalar_heat)) == 0
    )
    return {
        "finite_complex_controls": controls,
        "ordinary_vector_minus_scalar_determinant": True,
        "known_evanescent_finite_term_already_in_S176_not_added_twice": True,
        "known_mu1000_to_mu1_transport_preserved": True,
        "full_light_metric_measure_equivalence_proved": False,
    }


def parity_and_clock_frame():
    # Canonical W parity, including temporal pair; even rest sector untouched.
    omega = s.zeros(4).row_join(s.eye(4)).col_join((-s.eye(4)).row_join(s.zeros(4)))
    parity = -s.eye(8)
    assert parity.T * omega * parity == omega and parity.det() == 1
    # Coupled light/M1 covariance need not factor into two scalar states.
    Vlight = s.Matrix(
        [
            [2, 1, 0, 0],
            [1, 2, 0, 0],
            [0, 0, 1, s.Rational(1, 3)],
            [0, 0, s.Rational(1, 3), 1],
        ]
    )
    Vvector = s.diag(s.Rational(1, 2), 2)
    covariance = s.diag(Vlight, Vvector)
    mode_parity = s.diag(1, 1, 1, 1, -1, -1)
    assert mode_parity * covariance * mode_parity.T == covariance
    even = s.diag(2, 3, 4, 5, 6, 7)
    odd = s.zeros(6)
    odd[0, 4] = odd[4, 0] = 1
    assert mode_parity * odd * mode_parity == -odd
    assert s.trace(odd * covariance) == 0
    assert s.trace(even * covariance * odd * covariance.T) == 0
    # A tempting fixed-spatial-metric Horndeski frame fails on the actual clock.
    N = s.Symbol("N", positive=True)
    R = s.Function("R", positive=True)(N)
    F = s.Function("F", positive=True)(N)  # transformed lapse, not scalar density
    A4, B4 = -(R ** s.Rational(1, 4)) / 2, R ** s.Rational(3, 4) / 2
    Anew, Bnew = F * A4 / N, N * B4 / F
    condition = s.factor(Anew + Bnew + F * s.diff(Bnew, N) / s.diff(F, N))
    rhs = -N * (B4 + N * s.diff(B4, N)) / A4
    assert (
        s.simplify(condition - A4 * (F * s.diff(F, N) - rhs) / (N * s.diff(F, N))) == 0
    )
    RX, X = s.symbols("RX X")
    assert (
        s.simplify(
            rhs.subs(s.diff(R, N), -2 * X * RX / N)
            - N * s.sqrt(R) * (1 - 3 * X * RX / (2 * R))
        )
        == 0
    )
    t = s.Symbol("t", real=True)
    Eclock = 1 - s.Rational(3, 2) / (1 + t**2) ** 3
    t2 = (s.Rational(3, 2)) ** s.Rational(1, 3) - 1
    assert s.simplify(Eclock.subs(t**2, t2)) == 0
    assert 0 < t2 < 1
    return {
        "zero_odd_vector_expectations_in_parity_preserving_reference": True,
        "all_vector_quantum_response_zero": False,
        "restricted_Horndeski_lapse_map_equation": "F*F_N=N*sqrt(R)*(1-3X*R_X/(2R))",
        "singular_on_clock_at": "u^2=(3/2)^(1/3)-1",
        "scope": "fixed hat spatial metric, lapse-only map; not every field redefinition or a failure of the original regular canonical chart",
    }


def audit():
    before = inputs()
    spec = json.loads((REPO / SPEC).read_text())
    validate_spec(spec)
    rejected = []
    for path, value in [
        (("source_change_adopted",), False),
        (("common_quantum_reference_complete",), True),
        (("physical_matching_established",), True),
        (("physical_gravity_verdict_established",), True),
        (("original_P8_open",), False),
        (("classical_change", "source"), "S_old"),
        (("classical_change", "zeta"), "1"),
        (("retained", "physical_cutoff"), "1000"),
        (("retained", "states"), "new vacuum"),
        (("additional_finite_functional", "additional_complement"), "reset zeros"),
        (("reference_seed", "source_in_hat_variables"), "(R-1)*K"),
    ]:
        bad = copy.deepcopy(spec)
        target = bad
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        try:
            validate_spec(bad)
        except ValueError:
            rejected.append(".".join(path))
        else:
            raise AssertionError("Unapproved change accepted")
    result = {
        "milestone": "S0-REFERENCE.CANONICAL-AND-VECTOR-BRIDGE",
        "outcome": "S0_ADOPTED_CANONICAL_DOMAIN_AND_GAUSSIAN_VECTOR_BRIDGE_ESTABLISHED_FULL_REFERENCE_IN_PROGRESS",
        "protected_input_files": len(before),
        "source_change_adopted": True,
        "full_trace_and_compact_clock": full_trace_and_clock_domain(),
        "full_gaussian_reduction": gaussian_source_free_reduction(),
        "physical_Proca_probe_vertices": physical_vector_vertices(),
        "Proca_measure_and_inherited_finite_conversion": proca_measure_and_legacy(),
        "parity_and_restricted_frame_admission": parity_and_clock_frame(),
        "rejected_inputs": rejected,
        "common_quantum_reference_complete": False,
        "renormalized_response_norm_bound_established": False,
        "physical_matching_established": False,
        "physical_gravity_verdict_established": False,
        "original_P8_open": True,
    }
    assert inputs() == before
    return result


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
