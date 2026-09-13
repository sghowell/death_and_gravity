"""Independent full clock, ADM, tensor, contact, graph and scope checks."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_heavy_clock_quadratic import audit, clock, local, tensor


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_exact_residual_including_all_matrix_entries(name):
    value = audit.residuals()[name]
    assert all(e == 0 for e in value) if isinstance(value, s.MatrixBase) else value == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_every_complete_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_exact_complete_counts_and_frontier():
    assert (
        len(audit.residuals()),
        audit.scalar_entry_count(),
        len(audit.gates()),
        len(audit.controls()),
        audit.rejected_inputs(),
        len(audit.matching()),
    ) == (49, 64, 41, 9, 284, 97)
    assert len(audit.frontier()) == 9
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.validate_scope(audit.frontier(), audit.matching())


@pytest.mark.parametrize(
    "time,momentum,block",
    [
        (-s.Rational(1, 2), 1, "clock"),
        (s.Rational(1, 2), 10**300, "clock"),
        (-1, 0, "added_classical_heavy"),
        (1, 10**300, "added_classical_heavy"),
        (-s.Rational(1, 2), 0, "finite_local_scalar"),
        (s.Rational(1, 2), 0, "finite_local_tensor"),
    ],
)
def test_full_closed_stated_domains(time, momentum, block):
    assert audit.require_domain(time, momentum, block) == (time, momentum, block)


@pytest.mark.parametrize(
    "label,time,momentum",
    [
        ("central", 0, 100),
        ("central", s.Rational(1, 4), 100),
        ("outer", s.Rational(1, 8), 100),
        ("outer", -s.Rational(1, 2), 100),
        ("central", s.Rational(3, 16), 100),
        ("outer", s.Rational(3, 16), 100),
        ("regular_low", 0, s.Rational(1, 2)),
    ],
)
def test_actual_finite_chart_boundaries(label, time, momentum):
    assert audit.require_chart(label, time, momentum) == (label, time, momentum)


@pytest.mark.parametrize(
    "label,order",
    [
        ("physical_scalar_two_jets", 2),
        ("same_clock_scalar_two_jets", 2),
        ("local_tensor_total_order", 4),
        ("coefficient_phase_spatial_loss", 12),
        ("heavy_energy_scaled_phase", 0),
    ],
)
def test_actual_input_graph_orders(label, order):
    assert audit.require_graph(label, order) == (label, order)


@pytest.mark.parametrize(
    "stage",
    [
        "literal_full_QG2_clock_quadratic",
        "explicit_new_generic_chart_hypotheses",
        "current_coefficient_sector_twelve_derivative_comparison",
        "independent_classical_heavy_energy_phase",
        "full_physical_ADM_finite_heat_Hessian",
        "full_covariant_vacuum_constant_cancellation",
        "full_fixed_heavy_profile_Hessian_bound",
        "complete_weighted_local_tensor_Euler_graph",
        "complete_same_clock_local_Hessian",
    ],
)
def test_precise_allowed_stage(stage):
    assert audit.require_stage(stage) == stage


def test_exact_actual_parameters_and_nonzero_new_stress_budget():
    assert audit.require_parameters(local.state.MASS2, local.state.KAPPA) == (
        local.state.MASS2,
        local.state.KAPPA,
    )
    assert s.Rational(1, 10**400) < clock.EPS < s.Rational(2, 10**400)
    assert clock.EPS - clock.old_profile.EPS != 0
    p = clock.data()
    assert p["new_positive_coefficient_lower_bounds"]["Jc"] > s.Rational(1, 100)
    assert p["new_positive_coefficient_lower_bounds"]["Jc_minus_F"] > s.Rational(
        1, 1000
    )
    assert set(p["new_retuning_zero_through_two_jet_envelopes"]) == {0, 1, 2}
    assert all(
        0 < b < s.Rational(1, 10**6)
        for row in p["new_retuning_zero_through_two_jet_envelopes"].values()
        for b in row
    )


def test_independent_complete_current_Legendre_reduction():
    p = clock.data()
    Hreport = p["complete_current_coefficient_sector_Hamiltonian"]
    byname = {str(v): v for v in Hreport.free_symbols}
    J, l = byname["current_Jc"], byname["current_charge"]
    theta, E, A, T = clock.scalar.Theta, clock.scalar.E, clock.scalar.A, clock.scalar.Tc
    n, v, vd, sigma, sd, b, q = (
        clock.scalar.n,
        clock.scalar.v,
        clock.scalar.vd,
        clock.scalar.sigma,
        clock.scalar.sd,
        clock.scalar.b,
        clock.scalar.q,
    )
    pv, ps = clock.scalar.pv, clock.scalar.ps
    w = -l * E
    L = (
        -3 * vd * vd
        + (J + w * w / 2 - 3 * theta * theta) * n * n
        + 6 * theta * n * vd
        + sd * sd / 2
        + w * n * sd
        - 3 * l * vd * sigma
        + 2 * b * (vd - theta * n)
        + l * b * sigma
        + q * v * v
        + 2 * E * q * n * v
        - q * sigma * sigma / 2
        + 3 * T * n * v
        + s.Rational(9, 2) * A * v * v
    )
    velocities = {vd: (-pv + 6 * theta * n - 3 * l * sigma + 2 * b) / 6, sd: ps - w * n}
    fullH = s.expand((pv * vd + ps * sd - L).subs(velocities, simultaneous=True))
    shifted = s.factor(fullH.subs(b, pv / 2))
    numerator = theta * pv - w * ps + 3 * theta * l * sigma - (2 * E * q + 3 * T) * v
    assert s.cancel(s.diff(shifted, n) + 2 * J * n - numerator) == 0
    assert s.cancel(shifted.subs(n, numerator / (2 * J)) - Hreport) == 0
    assert all(
        not value.has(*clock.scalar.Z)
        for value in p["exact_current_Hamiltonian_coefficient_substitution"].values()
    )


def test_complete_weighted_phase_damping_and_substitution():
    p = clock.data()
    names = {
        str(v): v for v in p["complete_current_weighted_phase_matrix"].free_symbols
    }
    h = names["current_Hubble"]
    expected = clock.scalar.JC * s.hessian(
        p["complete_current_coefficient_sector_Hamiltonian"], clock.scalar.Z
    ) - 3 * h * s.diag(0, 0, 1, 1)
    assert (p["complete_current_weighted_phase_matrix"] - expected).applyfunc(
        s.cancel
    ) == s.zeros(4)
    assert p["exact_current_Hamiltonian_coefficient_substitution"][h] == clock.scalar.H
    assert p["generic_log_propagator"] == 10**29 and p["spatial_derivative_loss"] == 12
    assert p["generic_finite_chart_enclosures"]["complete_matrix_entries"] == 48


@pytest.mark.parametrize("ratio", (1, 2, 4))
def test_independent_normalized_heavy_fundamental_matrix_bound(ratio):
    with mp.workdps(80):
        Omega = mp.mpf(100)
        W = Omega * ratio
        chi = mp.exp(mp.j * mp.mpf(".73")) / mp.sqrt(2 * W)
        chip = (-mp.mpf(".1") - mp.j * W) * chi
        y = (mp.sqrt(Omega) * chi, chip / mp.sqrt(Omega))
        F = mp.matrix([[y[0], mp.conj(y[0])], [y[1], mp.conj(y[1])]])
        assert abs(mp.det(F) - mp.j) < mp.mpf("1e-70")
        norm = lambda matrix: max(
            sum(abs(matrix[i, j]) for i in range(2)) for j in range(2)
        )
        assert norm(F) < 20 and norm(F**-1) < 40
        assert norm(F) * norm(F**-1) < 800
    assert (
        "Gaussian metric response is NOT eliminated"
        in clock.data()["added_classical_heavy_mode"]
    )


def test_independent_nonlinear_ADM_identity_for_all_four_metric_jets():
    raw = local.full()
    n, zz = raw["n"], raw["zeta"]
    N = 1 + n
    invN = local.jet(1 - n + n * n)
    invN2 = local.jet(1 - 2 * n + 3 * n * n)
    shift = (
        s.diff(local.ed * local.Bd * local.pd + local.eg * local.Bg * local.pg, local.x)
        / local.a**2
    )
    H = s.diff(local.a, local.t) / local.a
    rate = local.jet(H + s.diff(zz, local.t) - shift * s.diff(zz, local.x))
    divergence = s.diff(shift, local.x)
    K = local.jet((3 * rate - divergence) * invN)
    KK = local.jet(((rate - divergence) ** 2 + 2 * rate**2) * invN2)
    spatial = (
        -local.expjet(-2 * zz)
        * (4 * s.diff(zz, local.x, 2) + 2 * s.diff(zz, local.x) ** 2)
        / local.a**2
    )
    lapse_laplacian = (
        local.expjet(-2 * zz)
        * (s.diff(N, local.x, 2) + s.diff(zz, local.x) * s.diff(N, local.x))
        / local.a**2
    )
    independent = local.jet(
        spatial
        + KK
        + K * K
        + 2 * invN * (s.diff(K, local.t) - shift * s.diff(K, local.x))
        - 2 * invN * lapse_laplacian
    )
    assert s.simplify(local.jet(raw["Rold"] - independent)) == 0


def test_complete_first_geometry_under_time_coordinate_change():
    p = local.data()
    shift = s.Function("time_shift")(local.t)
    H = s.diff(local.a, local.t) / local.a
    change = {local.ng: -s.diff(shift, local.t), local.zg: -H * shift, local.Bg: shift}
    transformed = (
        p["complete_first_curvature_source"].subs(change, simultaneous=True).doit()
    )
    assert (
        s.simplify(transformed + shift * s.diff(p["actual_reference_R"], local.t)) == 0
    )
    assert (
        s.simplify(
            p["complete_scalar_Weyl_source"].subs(change, simultaneous=True).doit()
        )
        == 0
    )


def test_complete_local_Hessian_has_symmetry_and_retained_lower_orders():
    p = local.data()
    full = p["complete_matched_local_heavy_Hessian_before_64pi2"]
    swap = {
        local.nd: local.ng,
        local.ng: local.nd,
        local.zd: local.zg,
        local.zg: local.zd,
        local.Bd: local.Bg,
        local.Bg: local.Bd,
    }
    assert s.expand(full - full.subs(swap, simultaneous=True).doit()) == 0
    coefficients = p["complete_finite_scalar_coefficients_before_64pi2"]
    N = next(
        v for v in coefficients["vacuum"].free_symbols if str(v) == "heavy_mass_squared"
    )
    assert s.degree(s.expand(full), N) == 1
    assert p["complete_polynomial_term_count"] == 73
    assert set(p["complete_bilinear_and_spatial_orders"]) == {
        (1, 1, 0),
        (1, 1, 2),
        (1, 1, 4),
    }
    raw = local.full()
    Rd = raw["Rold"].coeff(local.ed, 1).coeff(local.eg, 0) / local.pd
    Rg = raw["Rold"].coeff(local.eg, 1).coeff(local.ed, 0) / local.pg
    principal = (
        2 * coefficients["R_old_squared"] * local.a**3 * Rd * Rg
        + coefficients["Weyl_squared"] * p["complete_Weyl_mixed_density"]
    )
    remainder = s.expand(full - principal)
    assert remainder != 0
    actual = (
        remainder.subs(local.k, 0).coeff(s.diff(local.zd, local.t, 2)).coeff(local.zg)
    )
    target = (
        18
        * local.a**3
        * (
            coefficients["R_old"]
            + 2 * coefficients["R_old_squared"] * p["actual_reference_R"]
        )
    )
    assert s.simplify(actual - target) == 0


def test_independent_full_state_profile_second_chart_contact():
    p = local.data()
    rs, ps = s.symbols(
        "fixed_normalized_state_integral_energy fixed_normalized_state_integral_pressure",
        real=True,
    )
    d = s.Symbol("delta", real=True)
    u = local.t
    vd, vg = s.Function("vD")(u), s.Function("vG")(u)
    eps, eta = s.symbols("eps eta")
    n = eps * local.nd + eta * local.ng
    v = eps * vd + eta * vg
    N = 1 + n
    X = N**-2
    density = (
        N
        * s.exp(3 * v)
        * (1 + 2 * d * (X - 1)) ** (-s.Rational(3, 4))
        * (-ps - (rs + ps) * (X - 1) / 2)
    )
    mixed = s.diff(density, eps, eta).subs({eps: 0, eta: 0})
    actual = s.expand(local.a**3 * mixed.subs(d, 1 / (2 * (1 + u * u) ** 3)))
    assert s.expand(actual - p["retained_exact_fixed_state_profile_clock_Hessian"]) == 0
    contact = p["retained_nonzero_state_profile_second_clock_contact"]
    assert contact != 0
    assert s.simplify(
        contact.subs(
            {u: 0, local.a: 1, ps: 1, local.nd: 1, local.ng: 1}, simultaneous=True
        ).doit()
    ) == s.Rational(3, 2)


def test_complete_same_clock_graph_has_both_sources_of_terms():
    p = local.data()
    n, kappa = local.state.MASS2, local.state.KAPPA
    state = s.Rational(10**310, 1) / (n * kappa)
    assert (
        p["complete_same_clock_local_response_bound"]
        == 144 * p["complete_normalized_matched_local_Hessian_bound"] + 128 * state
    )
    assert 0 < p["complete_same_clock_local_response_bound"] < s.Rational(1, 10**580)
    assert (
        0
        < p["full_fixed_profile_plus_finite_heat_Hessian_bound"]
        < s.Rational(1, 10**580)
    )
    assert s.Rational(5, 4) ** 6 * 32 < 128
    assert 1 + s.Rational(1, 2) + s.Rational(3, 2) + s.Rational(33, 4) < 12
    assert "NOT applied to the state-profile piece alone" in p["same_clock_pullback"]
    assert "No quantum lapse/shift constraint is eliminated" in p["same_clock_pullback"]


def test_tensor_full_graph_and_wrong_adjoint_negative_control():
    p = tensor.data()
    assert 0 < p["full_normalized_local_tensor_graph_bound"] < s.Rational(1, 10**597)
    coefficients = p["complete_differential_coefficients"]
    assert (4, 0) in coefficients and (0, 2) in coefficients
    assert all(j + power <= 4 for j, power in coefficients)
    assert all(key in coefficients for key in ((3, 0), (2, 2), (1, 2), (0, 4)))
    assert s.simplify(coefficients[4, 0] - tensor.ell / 120) == 0
    assert s.simplify(coefficients[0, 2] - tensor.A / tensor.a**2) == 0
    f = s.Function("adjoint_test")(tensor.t)
    D = (
        s.diff(f, tensor.t, 2)
        + tensor.H * s.diff(f, tensor.t)
        + tensor.k**2 * f / tensor.a**2
    )
    dadj = (
        s.diff(tensor.a**3 * f, tensor.t, 2)
        - s.diff(tensor.a**3 * tensor.H * f, tensor.t)
        + tensor.a * tensor.k**2 * f
    ) / tensor.a**3
    assert s.simplify(dadj - D) != 0


def test_full_proper_conformal_tensor_conversion():
    f, g = s.Function("f")(tensor.t), s.Function("g")(tensor.t)

    def conformal(value):
        return (
            tensor.a * s.diff(tensor.a * s.diff(value, tensor.t), tensor.t)
            + tensor.k**2 * value
        )

    def proper(value):
        return (
            s.diff(value, tensor.t, 2)
            + tensor.H * s.diff(value, tensor.t)
            + tensor.k**2 * value / tensor.a**2
        )

    assert (
        s.simplify(
            conformal(f) * conformal(g) / tensor.a - tensor.a**3 * proper(f) * proper(g)
        )
        == 0
    )


def test_full_scope_preserves_remaining_response_and_no_inverse_claim():
    obs = audit.observable()
    assert "remaining exact SLE/subtraction-dependent" in obs["remaining"]
    assert "full quantum constraints/inverse" in obs["remaining"]
    assert "12 spatial derivatives lost" in obs["clock"]
    assert "Gaussian connected metric response is not zero" in obs["heavy"]
    assert "second-chart contact" in local.data()["same_clock_pullback"]
    assert (
        "physical pole of the full quantum theory" in tensor.data()["response_boundary"]
    )


def test_independent_full_256_component_Weyl_contraction_and_compact_boundary():
    t, x, y, z, k = s.symbols("t x y z k", real=True)
    coords = (t, x, y, z)
    f, g = s.Function("f")(t), s.Function("g")(t)
    metric = s.diag(1, -1, -1, -1)

    @cache
    def linear_weyl(field, sign):
        amp = field * s.exp(sign * s.I * k * x)
        perturb = s.diag(0, 0, -amp / s.sqrt(2), amp / s.sqrt(2))
        riem = {}
        for a in range(4):
            for b in range(4):
                for c in range(4):
                    for d in range(4):
                        riem[a, b, c, d] = (
                            s.diff(perturb[a, d], coords[b], coords[c])
                            + s.diff(perturb[b, c], coords[a], coords[d])
                            - s.diff(perturb[a, c], coords[b], coords[d])
                            - s.diff(perturb[b, d], coords[a], coords[c])
                        ) / 2
        ric = s.Matrix(
            4, 4, lambda b, d: sum(metric[a, a] * riem[a, b, a, d] for a in range(4))
        )
        scalar = sum(metric[a, a] * ric[a, a] for a in range(4))
        return {
            (a, b, c, d): s.expand(
                riem[a, b, c, d]
                - (
                    metric[a, c] * ric[b, d]
                    - metric[a, d] * ric[b, c]
                    - metric[b, c] * ric[a, d]
                    + metric[b, d] * ric[a, c]
                )
                / 2
                + scalar
                * (metric[a, c] * metric[b, d] - metric[a, d] * metric[b, c])
                / 6
            )
            for a, b, c, d in riem
        }

    left, right = linear_weyl(f, -1), linear_weyl(g, 1)
    mixed = s.expand(
        2
        * sum(
            metric[a, a]
            * metric[b, b]
            * metric[c, c]
            * metric[d, d]
            * value
            * right[a, b, c, d]
            for (a, b, c, d), value in left.items()
        )
    )
    wanted = (s.diff(f, t, 2) + k * k * f) * (s.diff(g, t, 2) + k * k * g)
    boundary = -2 * k * k * s.diff(f * g, t, 2)
    assert s.expand(mixed - wanted - boundary) == 0
    assert s.expand(mixed - wanted) != 0


def test_independent_literal_unimodular_tensor_curvature_and_heat_kinetic():
    t, k, a = local.t, local.k, local.a
    amp = local.ed * local.zd * local.pd + local.eg * local.zg * local.pg
    metric = s.diag(
        1,
        -a * a,
        -a * a * local.expjet(amp / s.sqrt(2)),
        -a * a * local.expjet(-amp / s.sqrt(2)),
    )
    inverse = s.diag(
        1,
        -1 / (a * a),
        -local.expjet(-amp / s.sqrt(2)) / (a * a),
        -local.expjet(amp / s.sqrt(2)) / (a * a),
    )
    Gamma = {}
    for aa in range(4):
        for bb in range(4):
            for cc in range(4):
                Gamma[aa, bb, cc] = local.jet(
                    inverse[aa, aa]
                    * (
                        s.diff(metric[aa, cc], local.coords[bb])
                        + s.diff(metric[aa, bb], local.coords[cc])
                        - s.diff(metric[bb, cc], local.coords[aa])
                    )
                    / 2
                )
    Ric = s.Matrix(
        4,
        4,
        lambda bb, dd: local.jet(
            sum(
                s.diff(Gamma[aa, dd, bb], local.coords[aa])
                - s.diff(Gamma[aa, aa, bb], local.coords[dd])
                + sum(
                    Gamma[aa, aa, rr] * Gamma[rr, dd, bb]
                    - Gamma[aa, dd, rr] * Gamma[rr, aa, bb]
                    for rr in range(4)
                )
                for aa in range(4)
            )
        ),
    )
    R = -local.jet(sum(inverse[aa, aa] * Ric[aa, aa] for aa in range(4)))
    H = s.diff(a, t) / a
    R0 = 6 * (s.diff(H, t) + 2 * H * H)
    assert s.simplify(R.subs({local.ed: 0, local.eg: 0}) - R0) == 0
    assert s.expand(R.coeff(local.ed, 1).coeff(local.eg, 0)) == 0
    assert s.expand(R.coeff(local.eg, 1).coeff(local.ed, 0)) == 0
    mixed_R = local.mixed(R)
    kinetic = s.diff(local.zd, t) * s.diff(
        local.zg, t
    ) - k * k * local.zd * local.zg / (a * a)
    assert s.simplify(mixed_R - kinetic / 2) == 0
    assert s.expand(local.mixed(R * R) - 2 * R0 * mixed_R) == 0
    coefs = local.data()["complete_finite_scalar_coefficients_before_64pi2"]
    A = coefs["R_old"] / 4 + coefs["R_old_squared"] * R0 / 2
    full = a**3 * (
        coefs["R_old"] * mixed_R + coefs["R_old_squared"] * local.mixed(R * R)
    )
    assert s.simplify(full - 2 * a**3 * A * kinetic) == 0
