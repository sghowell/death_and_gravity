"""Authorized COVZERO boundary, dimensional seed and legacy scale transport.

Exact controls accompany the written derivations. A regular dimensional
kinetic seed and a known Gaussian scale conversion are not a completed
constrained quantum measure, one-loop reference, gravity test or P8 solution.
"""

import copy
import hashlib
import json
from fractions import Fraction as F

import sympy as s
from p8_match1_rate_input import REPO
from p8_rate4_d8 import matrix, rate4_d8_coefficients
from p8_response_closure import inputs as previous_inputs

SPEC = "docs/candidates/p8-rate4-covzero-1.json"
EXTRA = {
    "scripts/p8_response_closure.py": "06b7c8689063c922cd0dfad95f89faae85a34a8bd621aa1bda302042e071b1be",
    "docs/assessment-2026-09-24-p8-response-closure-blockers.md": "832ffed57b81ad49d7466bd3e5791a9bbaee0d08a16d8e03fe26b60a84621426",
}
VECTOR_REPORT = "problems/P8/s6/continuation/s6_176/certificates/polynomial-vacuum-affine-proca-gaussian.json"


def inputs():
    manifest = previous_inputs()
    for name, expected in EXTRA.items():
        assert hashlib.sha256((REPO / name).read_bytes()).hexdigest() == expected
        assert name not in manifest or manifest[name] == expected
        manifest[name] = expected
    return manifest


def validate_spec(spec):
    """Reject silent changes to the authorized scope, not certify a theory."""
    required = {
        "candidate": "QG2-H8A420-RATE4-COVZERO-1 v1",
        "physical_boundary_rule_adopted": True,
        "common_quantum_reference_complete": False,
        "physical_matching_established": False,
        "physical_gravity_verdict_established": False,
        "original_P8_open": True,
        "order": "first quantum order only",
    }
    if any(
        spec.get(k) != v or type(spec.get(k)) is not type(v)
        for k, v in required.items()
    ):
        raise ValueError("Candidate identity, authority or completion scope changed")
    if spec["retained"]["physical_metric_signature"] != "+---":
        raise ValueError("The inherited physical metric has signature +---")
    if spec["retained"]["clock_norm"] != "X=g_phys^{mu nu}*partial_mu(u)*partial_nu(u)":
        raise ValueError("X=N^-2 on the clock; no extra minus sign")
    if spec["retained"]["physical_cutoff"] is not None:
        raise ValueError("No physical cutoff was authorized")
    if spec["additional_finite_functional"]["basis"] != ["1", "S", "U", "S^2"]:
        raise ValueError("Keep the four authorized contacts")
    if spec["reference_seed"]["mu"] != "1":
        raise ValueError("Transport legacy schemes to the declared common scale")
    if len(spec["open_reference_obligations"]) != 5:
        raise ValueError("Do not erase the outstanding reference obligations")


def kinetic_hessian(m, *, a_star, b_squared, R, RX, continued=True):
    """Actual covariant scalar-metric kinetic form on a tilted ADM slice.

    Finite integer-dimensional algebra fixtures only. The written invariant
    Schur proof, not these ranks, supplies the analytic-in-d statement.
    """
    if type(m) is not int or not 2 <= m <= 5:
        raise ValueError("Finite dimension controls use integer m=2 through 5")
    a_star, b_squared, R, RX = map(s.Rational, (a_star, b_squared, R, RX))
    X = a_star**2 - b_squared
    if X == 0 or R <= 0 or b_squared < 0:
        raise ValueError("Use a regular non-null fixture with positive R")
    velocity = s.Symbol("V_star")
    entries = s.symbols(f"k0:{m * (m + 1) // 2}")
    K = s.zeros(m)
    index = 0
    for i in range(m):
        for j in range(i, m):
            K[i, j] = K[j, i] = entries[index]
            index += 1
    A3 = RX / X
    factor = s.Rational(3 * m - 2, 2 * (m - 1)) if continued else s.Rational(7, 4)
    A4, A5 = -RX / X - factor * RX**2 / R, RX**2 / (R * X)
    Z = a_star * velocity + b_squared * K[0, 0]
    literal = (
        R * (s.trace(K * K) - s.trace(K) ** 2) / 2
        - 2 * RX * s.trace(K) * Z
        + A3 * a_star * (velocity + a_star * s.trace(K)) * Z
        + (A4 + a_star**2 * A5) * Z**2
    )
    return s.hessian(literal, (*entries, velocity))


def extended_chart_control(covector):
    """Regular disformal point Jacobian with independent gradient A_mu.

    R=(1+b*X^2)^2 is a generic-R control, not a replacement parent. Each
    contracted point-map Hessian is symmetric, so the phase check uses an
    arbitrary symmetric lower block. It is not a quantum ordering theorem.
    """
    metric = s.diag(1, -1, -1, -1)
    A = s.Matrix(covector)
    raised = metric.inv() * A
    X = (A.T * raised)[0]
    u = s.Rational(1, 3)
    b, bu = (1 + u**2) / 7, 2 * u / 7
    denominator = 1 + b * X**2
    C, D = 1 / denominator, b * X / denominator
    CX, DX = -2 * b * X / denominator**2, b * (1 - b * X**2) / denominator**2
    Cu, Du = -bu * X**2 / denominator**2, bu * X / denominator**2
    pairs = [(i, j) for i in range(4) for j in range(i, 4)]
    dX = s.Matrix([-(1 if i == j else 2) * raised[i] * raised[j] for i, j in pairs])
    direction = s.Matrix([CX * metric[i, j] + DX * A[i] * A[j] for i, j in pairs])
    Jg = C * s.eye(10) + direction * dX.T
    JA = s.Matrix(
        10,
        4,
        lambda r, k: (
            direction[r] * 2 * raised[k]
            + D
            * (
                int(pairs[r][0] == k) * A[pairs[r][1]]
                + int(pairs[r][1] == k) * A[pairs[r][0]]
            )
        ),
    )
    Ju = s.Matrix([Cu * metric[i, j] + Du * A[i] * A[j] for i, j in pairs])
    J = Jg.row_join(JA).row_join(Ju).col_join(s.zeros(5, 10).row_join(s.eye(5)))
    assert s.cancel(J.det() - C**9) == 0
    physical = C * metric + D * A * A.T
    assert s.cancel((A.T * physical.inv() * A)[0] - X) == 0
    assert s.cancel(physical.det() / metric.det() - C**3) == 0
    H = s.diag(*range(1, 16))
    H[0, 14] = H[14, 0] = s.Rational(2, 7)
    inverse_transpose = J.inv().T
    phase = J.row_join(s.zeros(15)).col_join(
        (-inverse_transpose * H).row_join(inverse_transpose)
    )
    omega = s.zeros(15).row_join(s.eye(15)).col_join((-s.eye(15)).row_join(s.zeros(15)))
    assert phase.T * omega * phase == omega
    assert phase.det() == 1
    return str(X)


def audit():
    before = inputs()
    spec = json.loads((REPO / SPEC).read_text())
    validate_spec(spec)
    m, R, RX, Ru, X = s.symbols("m R R_X R_u X", positive=True)
    lapse_speed, K, KK, V = s.symbols("s K KK V", real=True)
    factor = (3 * m - 2) / (2 * (m - 1))
    A3, A4, A5 = RX / X, -RX / X - factor * RX**2 / R, RX**2 / (R * X)
    D = s.factor(X * (A3 + A4) + X**2 * A5)
    assert D == -m * RX**2 * X / (2 * R * (m - 1))
    omega_X, omega_u = -RX / (2 * (m - 1) * R), -Ru / (2 * (m - 1) * R)
    shift = lapse_speed * omega_u + 2 * lapse_speed * omega_X * V
    physical_K, physical_KK = K + m * shift, KK + 2 * shift * K + m * shift**2
    raw = (
        -R * (physical_K**2 - physical_KK) / 2
        - lapse_speed * RX * physical_K * V
        + D.subs(X, lapse_speed**2) * V**2
        - lapse_speed * Ru * physical_K
    )
    expected = (
        -R * (K**2 - KK) / 2
        - lapse_speed * Ru * K / 2
        + 3 * m * lapse_speed**2 * Ru**2 / (8 * (m - 1) * R)
        + m * lapse_speed**2 * Ru * RX * V / ((m - 1) * R)
    )
    assert s.factor(raw - expected) == 0
    assert s.factor(s.diff(raw, K, V)) == s.factor(s.diff(raw, V, 2)) == 0
    delta = X * RX / ((m - 1) * R)
    E = -2 * X * RX - X**2 * A4
    gradient = (
        E
        + (m - 1) * (m - 2) * R * delta**2 / 2
        + 2 * (m - 1) * (R / 2 - X * RX) * delta
    )
    assert s.factor(gradient) == 0

    # Invariant Schur complement at ANY tilt b^2/X, not just unitary gauge.
    tilt = s.Symbol("tilt", real=True)
    norm = 1 + (m - 1) * (tilt - 1) ** 2
    trace = (m - 1) * tilt - m
    beta = 1 + tilt - factor
    assert s.factor(beta - (norm - trace**2 / (m - 1)) / 2) == 0
    naive = s.factor(1 + tilt - s.Rational(7, 4) - (norm - trace**2 / (m - 1)) / 2)
    assert s.cancel(naive + (m - 3) / (4 * (m - 1))) == 0
    ranks = []
    for dimension in (2, 3, 4, 5):
        for b2 in (F(0), F(1, 2)):
            args = {"a_star": 1, "b_squared": b2, "R": F(7, 6), "RX": F(1, 9)}
            new = kinetic_hessian(dimension, **args)
            old = kinetic_hessian(dimension, **args, continued=False)
            assert new.rank() == new.rows - 1
            assert old.rank() == old.rows - int(dimension == 3)
            ranks.append([dimension, str(b2), new.rank(), old.rank()])
    # The no-division-by-A_star branch: scalar gradient purely spatial.
    zero_time = kinetic_hessian(3, a_star=0, b_squared=F(1, 8192), R=1, RX=F(-1, 4096))
    assert zero_time.rank() == zero_time.rows - 1

    # The full source square keeps a regular trace pivot near m=3.
    assert (
        s.cancel(s.Rational(1, 2) - (R - 1) ** 2 / R + (R - 2) * (2 * R - 1) / (2 * R))
        == 0
    )
    assert s.Rational(5, 2) / (s.Rational(5, 2) - 1) == s.Rational(5, 3)
    assert 1 - s.Rational(5, 3) / 2 == s.Rational(1, 6)

    eps = s.Symbol("epsilon")
    C = R ** (-1 / (m - 1))
    U = R ** (-m / (2 * (m - 1)))
    M = R * U
    curvature = R * U / (2 * C)
    assert s.diff(A4.subs(m, 3 - 2 * eps), eps).subs(eps, 0) == -(RX**2) / (4 * R)
    for expression, logarithmic_derivative in (
        (C, -s.log(R) / 2),
        (U, -s.log(R) / 4),
        (M, -s.log(R) / 4),
        (curvature, s.log(R) / 4),
    ):
        assert (
            s.simplify(
                s.diff(expression.subs(m, 3 - 2 * eps), eps).subs(eps, 0)
                / expression.subs(m, 3)
                - logarithmic_derivative
            )
            == 0
        )
    assert A4.subs(m, 3) == -RX / X - 7 * RX**2 / (4 * R)

    # Null-gradient regularity is inherited from R-1=X^2*A, not 1/X itself.
    a0, a1, box_u, quad_u, H = s.symbols("a0 a1 box_u quad_u H")
    germ = 1 + a0 * X**2 + a1 * X**3
    Cg = s.series(germ ** (-1 / (m - 1)), X, 0, 3).removeO()
    assert s.simplify(s.diff((1 - Cg) / X, X).subs(X, 0) - a0 / (m - 1)) == 0
    source_scalar = (a0 + a1 * X) * (X * box_u - quad_u - m * H * X**2)
    assert source_scalar.subs(X, 0) == -a0 * quad_u
    assert (
        s.cancel(source_scalar - (germ - 1) * (box_u / X - quad_u / X**2 - m * H)) == 0
    )
    source_clock = s.simplify(
        ((germ - 1) * (box_u / X - quad_u / X**2 - m * H) * s.sqrt(X)).subs(
            {box_u: V + s.sqrt(X) * K, quad_u: X * V}
        )
    )
    assert s.simplify(source_clock - (germ - 1) * (K - m * H * s.sqrt(X))) == 0

    extended_norms = [
        extended_chart_control(values)
        for values in (
            [s.Rational(1, 4), 0, 0, 0],
            [0, s.Rational(1, 128), 0, 0],
            [s.Rational(1, 4), s.Rational(1, 4), 0, 0],
        )
    ]
    assert extended_norms == ["1/16", "-1/16384", "0"]

    # Independent minimal-vector-minus-scalar heat coefficient computation.
    d, Ro, Ric2, Riem2, mass, ell = s.symbols(
        "d R_old Ricci_squared Riemann_squared positive_mass log_mass_ratio"
    )
    scalar_heat = Ro**2 / 72 + (Riem2 - Ric2) / 180
    h0, h1 = d - 1, (d - 7) * Ro / 6
    h2 = (d - 1) * scalar_heat + Ric2 / 2 - Ro**2 / 6 - Riem2 / 12
    vector_heat = -(Ro**2) / 8 + s.Rational(29, 60) * Ric2 - Riem2 / 15
    assert s.expand(h2.subs(d, 4) - vector_heat) == 0
    pole_d = mass**4 * h0 - 2 * mass**2 * h1 + 2 * h2
    pole = s.expand(pole_d.subs(d, 4))
    evanescent = -2 * s.diff(pole_d, d)
    finite = s.expand(
        (
            mass**4 * (s.Rational(3, 2) - ell) * h0
            - 2 * mass**2 * (1 - ell) * h1
            - 2 * ell * h2
        ).subs(d, 4)
        + evanescent
    )
    legacy = json.loads((REPO / VECTOR_REPORT).read_text())[
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
    assert s.diff(h2, d) == scalar_heat
    transport = s.log(mass**2) * pole
    assert (
        s.expand(finite.subs(ell, s.log(mass**2)) + transport - finite.subs(ell, 0))
        == 0
    )

    # Bulk compact-variation stress check. R_old=+6(H'+2H^2), NOT R_P8.
    ht, hp, hpp, hppp = s.symbols("H Hprime Hsecond Hthird")

    def dt(expr):
        return s.diff(expr, ht) * hp + s.diff(expr, hp) * hpp + s.diff(expr, hpp) * hppp

    r_flat = 6 * hp + 12 * ht**2
    ric_flat = 9 * (hp + ht**2) ** 2 + 3 * (hp + 3 * ht**2) ** 2
    riem_flat = 12 * (hp + ht**2) ** 2 + 12 * ht**4
    f = s.expand(pole.subs({Ro: r_flat, Ric2: ric_flat, Riem2: riem_flat}))
    rho = s.expand(
        -f
        + ht * s.diff(f, ht)
        + (hp - 3 * ht**2) * s.diff(f, hp)
        - ht * dt(s.diff(f, hp))
    )

    # Independently vary q=log(a): H=q', H'=q''. The Euler derivative
    # of a^3*f is not defined using the continuity identity being tested.
    def volume_dt(expr):
        return dt(expr) + 3 * ht * expr

    pressure = s.expand(
        f - volume_dt(s.diff(f, ht)) / 3 + volume_dt(volume_dt(s.diff(f, hp))) / 3
    )
    pressure_from_continuity = s.cancel(-rho - dt(rho) / (3 * ht))
    assert s.expand(pressure - pressure_from_continuity) == 0
    assert (
        s.expand(
            rho
            + 3 * mass**4
            + 6 * mass**2 * ht**2
            - hp**2
            + 6 * ht**2 * hp
            + 2 * ht * hpp
        )
        == 0
    )
    assert s.expand(dt(rho) + 3 * ht * (rho + pressure)) == 0
    assert rho.subs({ht: 0, hp: 4}) == 16 - 3 * mass**4
    assert pressure.is_polynomial(ht, hp, hpp, hppp, mass)

    # Common-reference contacts stay symbolic until COMPLETE inputs arrive.
    born = (F(2), F(3), F(5), F(7))
    hard, conversion = (F(1), F(-2), F(4), F(-3)), (F(1, 3), F(1, 5), F(-1, 7), F(2, 9))
    coefficients = rate4_d8_coefficients(
        born=born, other_hard_real=hard, rate_conversion=conversion
    )
    for i, row in enumerate(matrix()):
        assert (
            sum(x * c for x, c in zip(row, coefficients, strict=True))
            + hard[i]
            + born[i] * conversion[i] / 2
            == 0
        )

    rejected = 0
    for field, value in (
        ("common_quantum_reference_complete", True),
        ("physical_gravity_verdict_established", True),
        ("physical_boundary_rule_adopted", False),
    ):
        changed = copy.deepcopy(spec)
        changed[field] = value
        try:
            validate_spec(changed)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("Unsupported status accepted")
    for key, value in (
        ("physical_metric_signature", "-+++"),
        ("clock_norm", "-g*du*du"),
        ("physical_cutoff", "1e100"),
    ):
        changed = copy.deepcopy(spec)
        changed["retained"][key] = value
        try:
            validate_spec(changed)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("Changed physical input accepted")
    for bad in (True, 1, 6, 3.0):
        try:
            kinetic_hessian(bad, a_star=1, b_squared=0, R=1, RX=1)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("Unsupported dimension fixture accepted")
    assert inputs() == before
    return {
        "milestone": "RATE4-COVZERO-1.REFERENCE-DIMENSIONAL-SEED-AND-LEGACY-TRANSPORT",
        "outcome": "AUTHORIZED_BOUNDARY_RULE_WITH_REGULAR_DIMENSIONAL_SEED_NOT_COMPLETE_QUANTIZATION",
        "physical_boundary_rule_adopted": True,
        "common_quantum_reference_complete": False,
        "protected_input_files": len(before),
        "candidate_spec_sha256": hashlib.sha256((REPO / SPEC).read_bytes()).hexdigest(),
        "dimension_dependent_A4": "-R_X/X-(3m-2)R_X^2/[2(m-1)R], m=d-1",
        "pure_scalar_metric_kinetic_degeneracy": "invariant arbitrary-tilt Schur identity and nine finite rank controls; not a full quantum measure theorem",
        "tilted_rank_controls": ranks,
        "zero_time_gradient_rank_control": True,
        "naive_four_dimensional_coefficient_schur_defect": "-(m-3)R_X^2/[4(m-1)R] in Z coordinates; vanishes only at m=3 or R_X=0",
        "full_unitary_ADM_cancellations": "K*V, V^2 and spatial lapse-gradient terms cancel; time primitive retained",
        "source_trace_pivot": ">1/6 for 5/2<=m<=7/2 and 1/2<R<6/5; local trace block only",
        "null_gradient_regular_chart_and_source": True,
        "extended_gradient_chart_controls": extended_norms,
        "extended_measure_scope": "independent A_mu=du auxiliary presentation makes the disformal map a regular point map; full cotangent phase Jacobian one, including nonzero null gradients. Not the full quantum gauge/constraint or ordering proof.",
        "A4_evanescent_derivative": "-R_X^2/(4R) at epsilon=0; not discarded before pole subtraction",
        "proca_heat_and_finite_coefficients": "independently reproduce pinned S176 covariant bulk coefficients, including -2*mV^4+2*mV^2*R_old/3-4*a4s finite evanescence",
        "proca_mu1000_to_common_mu1_transport": "+log(1e6)*P_V/(64*pi^2) preserves the old reference; not an adjustable new complementary coupling",
        "transport_negative_control": "omitting it changes bounce energy by log(1e6)*(3e12-16)/(64*pi^2); not canceled by the clock-null RATE4 lift",
        "sign_erratum": "prior closure note's prose X=-g*du*du is a typo: +--- requires X=+g*du*du. Its scripts already used X=N^-2 and results are unchanged.",
        "rate_contact_fixture_equalities": 4,
        "invalid_inputs_rejected": rejected,
        "remaining": spec["open_reference_obligations"] + spec["independent_remaining"],
        "physical_gravity_verdict_established": False,
        "original_P8_open": True,
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
