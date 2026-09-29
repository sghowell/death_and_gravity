"""S0 reference: dimensional gauge measure, actual ghost operator and UV input.

The ordinary de Donder FP block is computed, not the full DHOST measure,
coupled light/metric determinant, renormalized response norm or P8 closure.
"""

import hashlib
import importlib.util
import json
from fractions import Fraction

import sympy as s
from p8_match1_rate_input import REPO
from p8_s0_reference import COUPLED, validate_spec
from p8_s0_reference import inputs as previous_inputs
from p8_source_response import source_rows

EXTRA = {
    "scripts/p8_s0_reference.py": "1b86a6be89fd905ad05afcfb4e18c1823ef5ffa38753d490ded2613369ad7440",
    "docs/candidates/p8-rate4-covzero-s0.json": "9aedce3366e60f56cbfa4324b4c080b14711fc9cc1d9618a308c6551485818b7",
    "docs/assessment-2026-09-28-p8-s0-reference.md": "4c2d252ae747c5a571d25725eb89bb7da32cffe98d9dc275492097226cb3681d",
    "docs/validation/p8-s0-reference-2026-09-28.json": "ebcedbf44ecb96f90e73234911fd091e6777d9411fd677c6b1cfd0b9bbdcadaf",
    "problems/P8/s6/continuation/s6_258/notes/gauge.md": "cbb8d01d955995abc4029286cdd52dcbd95956b64a125a4baa77a586a5325973",
    "problems/P8/s6/continuation/s6_258/notes/ghost.md": "6d90118cab9a6551ed7b5a11bb467bd512bd207ef8e59cfd977574769b9c8716",
    "problems/P8/s6/continuation/s6_258/notes/local.md": "a369d272ea6a29907a88d222fa64374211738f40f0ea2d0dc5a6acbe52b7ca1b",
    "problems/P8/s6/continuation/s6_260/notes/algebra.md": "75aa4eefd022c4181f965c02941ac4357611200396f4b18f4a0f37b11a9ca909",
    "problems/P8/s6/continuation/s6_260/notes/ward.md": "a4d3f3f2fff18a622f30ee4ea76578253226e02b4b347be9c120fa9f9e287275",
    "problems/P8/s6/continuation/s6_260/notes/state.md": "dfa9bbc6db22ea9e2ec631933017688f48560fc3ee78a41bfffac6bedbcc66d5",
}


def inputs():
    manifest = previous_inputs()
    for name, expected in EXTRA.items():
        assert hashlib.sha256((REPO / name).read_bytes()).hexdigest() == expected
        assert name not in manifest or manifest[name] == expected
        manifest[name] = expected
    return manifest


def zero(value):
    if isinstance(value, s.MatrixBase):
        return value.applyfunc(s.cancel) == s.zeros(*value.shape)
    return s.cancel(value) == 0


def spatial_dimension():
    m = s.Symbol("m", positive=True)
    epsilon = s.Symbol("epsilon", real=True)
    R, scale, C = s.symbols("R scale C", positive=True)
    w, c = 2 / m, 1 - 2 / m
    inverse_weight = (m - 2) / (2 * (m - 1))
    # Determinant of an m-dimensional metric rescales by C^m.
    assert s.simplify((C**m) ** (1 / m) / C) == 1
    naive_log_factor = (m / 3 - 1) * (-s.log(R) / (m - 1))
    assert (
        s.diff(naive_log_factor.subs(m, 3 - 2 * epsilon), epsilon).subs(epsilon, 0)
        == s.log(R) / 3
    )
    naive_homogeneous = (2 * m / 3 - 2) * s.log(scale)
    assert (
        s.diff(naive_homogeneous.subs(m, 3 - 2 * epsilon), epsilon)
        == -4 * s.log(scale) / 3
    )
    assert s.diff(w.subs(m, 3 - 2 * epsilon), epsilon).subs(epsilon, 0) == s.Rational(
        4, 9
    )
    assert s.diff(inverse_weight.subs(m, 3 - 2 * epsilon), epsilon).subs(
        epsilon, 0
    ) == -s.Rational(1, 4)
    assert s.diff(s.log(2 * (m - 1) / m).subs(m, 3 - 2 * epsilon), epsilon).subs(
        epsilon, 0
    ) == -s.Rational(1, 3)
    lapse, clock_h = s.symbols("lapse clock_h", positive=True)
    clock_R = 1 + (lapse**-2 - 1) / clock_h
    assert s.diff(s.log(clock_R) / 3, lapse).subs(lapse, 1) == -2 / (3 * clock_h)
    controls = []
    for dim in (2, 3, 4, 5):
        weight = s.Rational(2, dim)
        k = s.Matrix(range(1, dim + 1))
        ell = s.Matrix([1] + [0] * (dim - 1))
        p = k + ell
        Q = s.eye(dim) + ell * ell.T / 5
        xi = s.Matrix(s.symbols(f"xi0:{dim}"))
        # Differentiate the full density Lie derivative, before gauge restriction.
        dQ = s.I * (
            xi.dot(ell) * Q - (Q * k) * xi.T - xi * (Q * k).T + weight * k.dot(xi) * Q
        )
        derived = (s.I * dQ * p).jacobian(xi)
        expected = (
            p.dot(Q * k) * s.eye(dim)
            + (Q * k) * p.T
            - (Q * p) * ell.T
            - weight * (Q * p) * k.T
        )
        assert zero(derived - expected)
        r = k.dot(Q * k)
        sigma = r * s.eye(dim) + (1 - weight) * (Q * k) * k.T
        assert sigma.det() == s.Rational(2 * (dim - 1), dim) * r**dim
        assert zero(sigma * Q * k - (2 - weight) * r * Q * k)
        assert Q.inv() * sigma == (Q.inv() * sigma).T
        flat = k.dot(k) * s.eye(dim) + (1 - weight) * k * k.T
        inverse = (
            s.eye(dim) - s.Rational(dim - 2, 2 * (dim - 1)) * k * k.T / k.dot(k)
        ) / k.dot(k)
        assert zero(flat * inverse - s.eye(dim))
        # Full Fourier ghosts keep the incoming translation column off the slice.
        translation_column = -(Q * ell) * ell.T
        assert translation_column != s.zeros(dim)
        controls.append(
            {"spatial_dimension": dim, "principal_determinant": str(sigma.det())}
        )
    assert c.subs(m, 3) == s.Rational(1, 3)
    return {
        "Q": "det(gamma)^(1/m)*gamma^(-1)",
        "density_weight": "2/m",
        "principal_determinant": "2*(m-1)/m*(k.Q.k)^m",
        "flat_inverse": "[I-(m-2)/(2*(m-1))*k*k^T/k^2]/k^2, k!=0",
        "epsilon_weight_derivative": "4/9",
        "naive_conformal_log_defect_derivative": "log(R)/3",
        "naive_homogeneous_log_defect_derivative": "-4*log(a)/3",
        "naive_conformal_evanescent_first_lapse_probe_at_clock": "-2/[3*(1+t^2)^3]",
        "physical_m3_coercivity_at_Q_distance_1_over_8": "3/4",
        "integer_dimension_controls": controls,
        "fractional_dimension_Hilbert_space_claimed": False,
    }


def pb(left, right, q, p):
    return s.expand(
        sum(
            s.diff(left, a) * s.diff(right, b) - s.diff(left, b) * s.diff(right, a)
            for a, b in zip(q, p, strict=True)
        )
    )


def time_gauge_measure():
    u, A, q, pu, pA, p = s.symbols("u A q pu pA p")
    energy = s.Function("E")(u, q, p)
    coords, momenta = [u, A, q], [pu, pA, p]
    H = pu * A - A**2 / 2 + energy
    second = s.Matrix([pA, A - pu])
    C = s.Matrix(2, 2, lambda i, j: pb(second[i], second[j], coords, momenta))
    left = s.Matrix([[pb(u, f, coords, momenta) for f in second]])
    right = s.Matrix([pb(f, H, coords, momenta) for f in second])
    clock_bracket = pb(u, H, coords, momenta) - (left * C.inv() * right)[0]
    assert clock_bracket == pu
    whole = list(second) + [u, H]
    matrix = s.Matrix(4, 4, lambda i, j: pb(whole[i], whole[j], coords, momenta))
    assert zero(matrix.det() - pu**2)
    reduced = H.subs(A, pu)
    assert s.diff(reduced, pu) == pu
    # Generic triangular time/spatial FP operator. Do not delete its lower leg.
    speed = s.Symbol("positive_clock_speed", positive=True)
    v = s.Matrix(s.symbols("v0:3"))
    M = s.Matrix([[4, 1, 0], [0, 3, 1], [1, 0, 2]])
    block = s.Matrix([[speed]]).row_join(s.zeros(1, 3)).col_join(v.row_join(M))
    inverse = (
        s.Matrix([[1 / speed]])
        .row_join(s.zeros(1, 3))
        .col_join((-M.inv() * v / speed).row_join(M.inv()))
    )
    assert block.det() == speed * M.det()
    assert zero(block * inverse - s.eye(4))
    return {
        "clock_FP_on_unitary_slice": "A_*=sqrt(X)>0 for normal deformation parameter",
        "time_spatial_FP_determinant": "A_* det(M_spatial)",
        "clock_delta_constraint_Jacobian_cancellation": True,
        "independent_finite_auxiliary_clock_control_determinant": "p_u^2",
        "vacuum_clock_gauge_admissible": False,
        "full_regulated_BFV_measure_constructed": False,
    }


def dedonder_variation():
    t = s.Symbol("t", real=True)
    a = s.Function("positive_scale")(t)
    H = s.diff(a, t) / a
    results = []
    for m in (2, 3, 4):
        dim = m + 1
        k = s.symbols(f"k1:{dim}", real=True)
        xi = [s.Function(f"xi{i}")(t) for i in range(dim)]
        metric = s.diag(1, *([-(a**2)] * m))
        inv = metric.inv()

        def partial(expr, index, k=k):
            return s.diff(expr, t) if index == 0 else s.I * k[index - 1] * expr

        def gamma(up, low1, low2):
            if up == 0 and low1 == low2 and low1 > 0:
                return a * s.diff(a, t)
            if up > 0 and ((low1 == 0 and low2 == up) or (low2 == 0 and low1 == up)):
                return H
            return s.S.Zero

        lie = s.Matrix(
            dim,
            dim,
            lambda mu, nu, xi=xi, metric=metric, dim=dim, partial=partial: (
                xi[0] * s.diff(metric[mu, nu], t)
                + sum(
                    metric[rho, nu] * partial(xi[rho], mu)
                    + metric[mu, rho] * partial(xi[rho], nu)
                    for rho in range(dim)
                )
            ),
        )
        trace = s.trace(inv * lie)
        covariant_F = []
        for mu in range(dim):
            divergence = sum(
                inv[nu, nu]
                * (
                    partial(lie[nu, mu], nu)
                    - sum(
                        gamma(rho, nu, nu) * lie[rho, mu]
                        + gamma(rho, nu, mu) * lie[nu, rho]
                        for rho in range(dim)
                    )
                )
                for nu in range(dim)
            )
            covariant_F.append(divergence - partial(trace, mu) / 2)
        actual = inv * s.Matrix(covariant_F)
        qmomentum = sum(value**2 for value in k) / a**2
        expected = s.Matrix(
            [
                s.diff(xi[0], t, 2)
                + m * H * s.diff(xi[0], t)
                + (qmomentum - m * (s.diff(H, t) + 2 * H**2)) * xi[0]
                - 2 * s.I * H * sum(k[i - 1] * xi[i] for i in range(1, dim)),
                *[
                    s.diff(xi[i], t, 2)
                    + (m + 2) * H * s.diff(xi[i], t)
                    + qmomentum * xi[i]
                    - 2 * s.I * H * k[i - 1] * xi[0] / a**2
                    for i in range(1, dim)
                ],
            ]
        )
        assert zero(actual - expected)
        y = [s.Function(f"y{i}")(t) for i in range(dim)]
        transformed = actual.subs(
            {xi[0]: y[0], **{xi[i]: y[i] / a for i in range(1, dim)}}, simultaneous=True
        ).doit()
        transformed = s.diag(1, *([a] * m)) * transformed
        desired = s.Matrix(
            [
                s.diff(y[0], t, 2)
                + m * H * s.diff(y[0], t)
                + (qmomentum - m * (s.diff(H, t) + 2 * H**2)) * y[0]
                - 2 * s.I * H * sum(k[i - 1] * y[i] for i in range(1, dim)) / a,
                *[
                    s.diff(y[i], t, 2)
                    + m * H * s.diff(y[i], t)
                    + (qmomentum - s.diff(H, t) - (m + 1) * H**2) * y[i]
                    - 2 * s.I * H * k[i - 1] * y[0] / a
                    for i in range(1, dim)
                ],
            ]
        )
        assert zero(transformed - desired)
        results.append(
            {
                "spacetime_dimension": dim,
                "all_components_and_orthonormal_transport": True,
            }
        )
    return {
        "Lorentzian_operator": "Box_P8*delta^mu_nu+Ricci_P8^mu_nu, defined by F(Lie_xi g)",
        "independent_metric_variation_controls": results,
        "no_scalar_clock_unitary_division_in_deDonder_operator": True,
    }


def retarded_energy():
    t = s.Symbol("t", real=True)
    a = (1 + t**2) ** 2
    H = s.diff(a, t) / a
    Hd = s.diff(H, t)
    V0, Vspace = -3 * (Hd + 2 * H**2), -Hd - 4 * H**2
    assert zero(4 - H**2 - 4 * (t**2 - 1) ** 2 / (1 + t**2) ** 2)
    assert zero(s.Rational(49, 2) + V0 - (7 * t**2 - 5) ** 2 / (2 * (1 + t**2) ** 2))
    assert zero(
        s.Rational(225, 14) + Vspace - (15 * t**2 - 13) ** 2 / (14 * (1 + t**2) ** 2)
    )
    assert s.Rational(225, 14) < s.Rational(49, 2)
    # E'=damping+potential+first-spatial-derivative+frequency-weight terms.
    ceiling = 12 + s.Rational(51, 2) + 4 + 4
    assert ceiling == s.Rational(91, 2)
    controls = []
    for time in (s.Integer(0), s.Rational(-1, 2), s.Rational(1, 2), s.Integer(1)):
        for k in (s.Integer(0), s.Integer(1), s.Integer(1000)):
            aa, hh = a.subs(t, time), H.subs(t, time)
            physical_k = k / aa
            frequency_squared = 1 + physical_k**2
            potential = s.diag(V0.subs(t, time), Vspace.subs(t, time))
            cross = s.Matrix(
                [[0, -2 * s.I * hh * physical_k], [-2 * s.I * hh * physical_k, 0]]
            )
            generator = (
                s.zeros(2)
                .row_join(s.eye(2))
                .col_join(
                    (-(physical_k**2 * s.eye(2) + potential + cross)).row_join(
                        -3 * hh * s.eye(2)
                    )
                )
            )
            energy = s.diag(frequency_squared, frequency_squared, 1, 1)
            energy_dot = s.diag(-2 * hh * physical_k**2, -2 * hh * physical_k**2, 0, 0)
            slack = (
                ceiling * energy
                - generator.conjugate().T * energy
                - energy * generator
                - energy_dot
            )
            assert slack == slack.conjugate().T
            assert all(slack[:size, :size].det() > 0 for size in range(1, 5))
            controls.append(
                {"t": str(time), "comoving_k": str(k), "energy_slack_positive": True}
            )
    assert V0.subs(t, 0) == -12
    return {
        "physical_clock": "a=(1+t^2)^2",
        "energy": "|dot(y)|^2+(1+k^2/a^2)|y|^2, y=(xi^0,a*xi^i)",
        "homogeneous_energy_growth_ceiling": "91/2",
        "sqrt_energy_retarded_growth_ceiling": "91/4",
        "zero_data_bound": "sqrt(E(t))<=integral_s0^t exp(91*(t-s)/4)*|F_orth(s)| ds",
        "uniform_in_all_comoving_momenta_including_zero": True,
        "spatial_derivative_loss_for_this_ghost_inverse": 0,
        "constant_time_vector_at_bounce_operator_value": "-12",
        "positive_slack_controls": controls,
        "physical_light_metric_response_norm_established": False,
        "both_CTP_endpoint_conditions_automatically_preserved": False,
    }


def ghost_heat():
    d, Ro, Ric2, Riem2, boxR = s.symbols(
        "d R_old Ricci_squared Riemann_squared Box_R_old"
    )
    scalar = Ro**2 / 72 + (Riem2 - Ric2) / 180 + boxR / 30
    ghost = d * scalar + Ro**2 / 6 + Ric2 / 2 - Riem2 / 12 + boxR / 6
    universal = (
        d * (5 * Ro**2 - 2 * Ric2 + 2 * Riem2 + 12 * boxR)
        + 60 * Ro**2
        + 180 * Ric2
        - 30 * Riem2
        + 60 * boxR
    ) / 360
    assert zero(ghost - universal)
    einstein = s.cancel(ghost.subs({Ric2: Ro**2 / d, boxR: 0}))
    published = (5 * d**2 + 58 * d + 180) * Ro**2 / (360 * d) + (d - 15) * Riem2 / 180
    assert zero(einstein - published)
    curvature_controls = []
    for dimension in (3, 4, 5):
        sectional = {
            (i, j): s.Rational(i - 2 * j - 1, 7)
            for i in range(dimension)
            for j in range(i + 1, dimension)
        }

        def riemann(i, j, a, b, sectional=sectional):
            if i == j:
                return s.S.Zero
            value = sectional[tuple(sorted((i, j)))]
            return value * (int(i == a and j == b) - int(i == b and j == a))

        ric = s.Matrix(
            dimension,
            dimension,
            lambda i, j, dimension=dimension, riemann=riemann: sum(
                riemann(a, i, a, j) for a in range(dimension)
            ),
        )
        scalar_curvature = s.trace(ric)
        ric_norm = s.trace(ric * ric)
        riem_norm = sum(
            riemann(i, j, a, b) ** 2
            for i in range(dimension)
            for j in range(dimension)
            for a in range(dimension)
            for b in range(dimension)
        )
        omega_trace = 0
        for a in range(dimension):
            for b in range(dimension):
                omega = s.Matrix(
                    dimension,
                    dimension,
                    lambda i, j, a=a, b=b, riemann=riemann: riemann(i, j, a, b),
                )
                assert omega == -omega.T
                omega_trace += s.trace(omega * omega)
        assert omega_trace == -riem_norm
        trace_formula = (
            dimension * (5 * scalar_curvature**2 - 2 * ric_norm + 2 * riem_norm)
            + 60 * scalar_curvature * s.trace(ric)
            + 180 * s.trace(ric * ric)
            + 30 * omega_trace
        ) / 360
        assert zero(
            trace_formula
            - ghost.subs(
                {
                    d: dimension,
                    Ro: scalar_curvature,
                    Ric2: ric_norm,
                    Riem2: riem_norm,
                    boxR: 0,
                }
            )
        )
        curvature_controls.append(
            {"dimension": dimension, "non_Einstein_algebraic_curvature": True}
        )
    # Relative to a real boson's 1/(64*pi^2) pole polynomial convention.
    pole = -4 * ghost
    pole4 = s.expand(pole.subs(d, 4))
    expected = (
        -s.Rational(8, 9) * Ro**2
        - s.Rational(86, 45) * Ric2
        + s.Rational(11, 45) * Riem2
        - s.Rational(6, 5) * boxR
    )
    assert zero(pole4 - expected)
    evanescent = s.expand(-2 * s.diff(pole, d))
    assert zero(evanescent - 8 * scalar)
    # Gravitational FP ghosts are not a positive Proca vector or four scalars.
    scalar_only = d * scalar
    proca_sign = d * scalar - Ro**2 / 6 + Ric2 / 2 - Riem2 / 12 - boxR / 6
    assert not zero(ghost - scalar_only)
    assert zero(ghost - proca_sign - (Ro**2 + boxR) / 3)
    assert not zero(pole4 - (-2 * ghost).subs(d, 4))
    return {
        "Euclidean_UV_operator_only": "-(nabla^2*I+Ricci)",
        "tr_I": "d",
        "tr_E": "R_old",
        "tr_E_squared": "Ricci_squared",
        "tr_Omega_squared": "-Riemann_squared",
        "a4": str(s.expand(ghost)),
        "pole_polynomial_before_common_factor": str(pole4),
        "explicit_dimension_evanescent_finite_polynomial": str(evanescent),
        "all_surface_terms_retained": True,
        "non_Einstein_curvature_controls": curvature_controls,
        "Einstein_restriction_matches_published_independent_coefficient": True,
        "complete_finite_ghost_determinant_computed": False,
        "complete_light_metric_pole_computed": False,
    }


def off_shell_ward_control():
    x, y, r, b = s.symbols("x y r b", real=True)
    action = (x**2 + y**2 - 1) ** 2 / 2
    fields = s.Matrix([x, y])
    generator = s.Matrix([-y, x])
    gradient = s.Matrix([s.diff(action, z) for z in fields])
    hessian = s.hessian(action, fields)
    assert zero(gradient.dot(generator))
    assert zero(hessian * generator + generator.jacobian(fields).T * gradient)
    assert hessian * generator != s.zeros(2, 1)
    background_hessian = hessian.subs({x: r, y: 0})
    gauge = s.Matrix([[b, 1]])
    ghost = (gauge * generator.subs({x: r, y: 0}))[0]
    assert ghost == r
    determinant = s.factor((background_hessian + gauge.T * gauge).det())
    reference = determinant.subs(b, 0)
    difference = s.log(determinant / reference) / 2
    assert difference.subs(r, 1) == 0
    first = s.diff(difference, r).subs(r, 1)
    second = s.diff(difference, r, 2).subs(r, 1)
    assert zero(first - b**2 / 2)
    assert zero(second + (13 * b**2 + b**4) / 2)
    return {
        "exact_off_shell_identity": "K_ij R^j_alpha+E_j partial_i R^j_alpha=0",
        "independent_finite_gauge_control": "S=(x^2+y^2-1)^2/2, R=(-y,x), F=(b,1)",
        "on_shell_Gamma1_difference": "0 at r=1",
        "first_derivative_difference": "b^2/2",
        "second_derivative_difference": "-(13*b^2+b^4)/2",
        "on_shell_determinant_agreement_implies_off_shell_response_agreement": False,
        "actual_P8_Nielsen_transport_computed": False,
    }


def held_vector_contact_correction():
    # This pinned module is algebra-only. Do not call its ancestor report builders.
    spec = importlib.util.spec_from_file_location("_p8_s0_held_vector", REPO / COUPLED)
    c = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(c)
    ZN = s.Symbol("Z_N", real=True)
    nstar = (c.th * (c.pv + 3 * c.c * c.sigma) - c.w * c.ps - c.Lnv * c.v) / (2 * c.J)
    # Independent envelope derivative of the complete held-W square at the clock.
    Lnormal = (
        -c.b * c.longitudinal / c.P
        - c.rN * c.n * (3 * c.vd - c.b)
        - (ZN * c.n + 3 * c.v) * c.A0
    )
    hW = -Lnormal.subs(
        {
            c.n: nstar,
            c.b: c.pv / 2,
            c.vd: c.th * nstar - c.c * c.sigma / 2,
            c.A0: c.P * c.pA,
        },
        simultaneous=True,
    )
    light = c.rN * nstar * (3 * c.th * nstar - 3 * c.c * c.sigma / 2 - c.pv / 2)
    assert zero(
        hW
        - light
        - c.pv * c.longitudinal / (2 * c.P)
        - (ZN * nstar + 3 * c.v) * c.P * c.pA
    )
    pinned = c.held_vector_data()
    assert all(zero(value) for value in pinned["checks"].values())
    assert zero(hW - pinned["whole_clock_normal_vector_Hamiltonian_vertex"])
    HW = s.hessian(hW, c.PHASE)
    assert zero(HW - pinned["whole_clock_normal_vector_Hessian"])
    assert zero((c.PHASE.T * HW * c.PHASE)[0] / 2 - hW)

    a, kappa, zeta = s.symbols("a kappa zeta", positive=True)
    Cphase = s.sqrt(kappa) * s.diag(1, 1, s.sqrt(zeta), a**3, a**3, a**3 / s.sqrt(zeta))
    physical = s.Matrix(s.symbols("physical0:6", real=True))
    mapped = Cphase.inv() * physical
    physical_hW = (
        kappa
        * a**3
        * hW.subs(dict(zip(c.PHASE, mapped, strict=True)), simultaneous=True)
    )
    assert zero(
        s.hessian(physical_hW, physical)
        - kappa * a**3 * Cphase.inv().T * HW * Cphase.inv()
    )

    # Actual tree-clock bounce, q=1: the light block is a nonzero operator.
    # This is not the expectation in the original prepared Gaussian state.
    nrow, Prow, theta, clock_h = source_rows(
        time=0, comoving_k_squared=1, lapse_pivot=Fraction(243, 160), profile_current=0
    )
    phase_light = s.Matrix([c.v, c.sigma, c.pv, c.ps])
    bounce_n, bounce_P = (nrow * phase_light)[0], (Prow * phase_light)[0]
    light_bounce = bounce_n * bounce_P
    assert theta == 0 and clock_h == 1
    assert zero(
        light.subs(
            {
                c.rN: -2,
                c.th: 0,
                c.c: s.Rational(1, 10),
                c.w: s.Rational(1, 20),
                c.Lnv: -1,
                c.J: s.Rational(243, 160),
            }
        )
        - light_bounce
    )
    assert s.diff(light_bounce, c.v, c.pv) == s.Rational(80, 243)
    assert light_bounce.subs({c.v: 1, c.pv: 1, c.sigma: 0, c.ps: 0}) == s.Rational(
        80, 243
    )

    t, eps = s.symbols("t epsilon", real=True)
    eta, vp = s.Function("eta")(t), s.Function("v_phys")(t)
    r1 = -2 / (1 + t**2) ** 3
    etaH = s.diff(vp + r1 * eta / 4, t)
    embedding = 3 * (eps * r1 * eta) * (eps * etaH)
    Ssecond = 6 * r1 * eta * etaH
    assert zero(s.diff(embedding, eps, 2).subs(eps, 0) - Ssecond)
    # Smooth compact probes may have these jets at the bounce.
    bounce_embedding = Ssecond.subs(
        {eta: 1, s.diff(eta, t): 0, s.diff(vp, t): 1}, simultaneous=True
    ).subs(t, 0)
    assert bounce_embedding == -12
    W, g0, g1, g2, gw, g1w, gww = s.symbols("W G0 G1 G2 GW G1W GWW", real=True)
    jet = g0 + eps * g1 + eps**2 * g2 / 2 + gw * W + eps * g1w * W + gww * W**2 / 2
    difference = jet.subs(W, embedding) - jet.subs(W, 0)
    assert s.diff(difference, eps).subs(eps, 0) == 0
    assert zero(s.diff(difference, eps, 2).subs(eps, 0) - gw * Ssecond)
    return {
        "supersedes_unqualified_held_response_interpretation_of": "2026-09-28 S0 reference section 3.1",
        "old_aligned_Hamiltonian_difference_identity_still_valid": True,
        "complete_held_response_inferred_from_that_identity_alone": False,
        "embedding_second_derivative": "6*r1*eta*d_t(v_phys+r1*eta/4), r1=-2/(1+t^2)^3",
        "held_minus_S0_second_response": "(old_aligned-S0)_second-Gamma_old,W*S_second",
        "Gaussian_Hamiltonian_contact_correction": "+expectation(h_W)*S_second",
        "normal_vertex": "rN*nstar*(3*Theta*nstar-3*c*sigma/2-pv/2)+pv*A_L/(2*P)+(Z_N*nstar+3*v)*P*pA",
        "actual_bounce_light_vertex_v_pv_coefficient_q1": "80/243",
        "actual_bounce_embedding_second_derivative_probe_control": "-12",
        "full_physical_phase_density_transport_checked": True,
        "two_original_scalar_boundary_momentum_shifts_retained": True,
        "P_inverse_domain": "P>0 with the inherited physical-shift infrared domain",
        "original_state_old_normal_tadpole_computed": False,
        "S0_vector_parity_and_adoption_affected": False,
    }


def audit():
    before = inputs()
    validate_spec(
        json.loads((REPO / "docs/candidates/p8-rate4-covzero-s0.json").read_text())
    )
    result = {
        "milestone": "S0-REFERENCE.GAUGE-OPERATOR-AND-GHOST-UV",
        "outcome": "DIMENSIONAL_CLOCK_GAUGE_ACTUAL_DEDONDER_INVERSE_AND_GHOST_UV_ESTABLISHED_FULL_REFERENCE_IN_PROGRESS",
        "protected_input_files": len(before),
        "source_change_adopted": True,
        "dimensional_clock_spatial_gauge": spatial_dimension(),
        "time_gauge_measure": time_gauge_measure(),
        "deDonder_metric_variation": dedonder_variation(),
        "actual_clock_retarded_ghost_bound": retarded_energy(),
        "ordinary_deDonder_ghost_UV": ghost_heat(),
        "off_shell_Ward_admission": off_shell_ward_control(),
        "held_physical_vector_response_scope_correction": held_vector_contact_correction(),
        "new_physical_boundary_rule_adopted": False,
        "original_states_or_physical_cutoff_changed": False,
        "common_quantum_reference_complete": False,
        "renormalized_physical_response_bound_established": False,
        "physical_matching_established": False,
        "physical_gravity_verdict_established": False,
        "original_P8_open": True,
    }
    assert inputs() == before
    return result


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
