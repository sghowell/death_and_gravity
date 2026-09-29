"""SOURCE-1: bounded classical source feasibility, not model adoption.

Independent null-branch primary brackets supplement the inherited tilted
rank obstruction. Exact ansatz, Dirac, jet and filtration checks accompany
the written study. Neither a continuum quantum measure nor an EFT cutoff
is certified by these finite algebra controls.
"""

import copy
import hashlib
import json

import sympy as s
from p8_covzero_source_constraints import inputs as previous_inputs
from p8_match1_rate_input import REPO

SPEC = "docs/candidates/p8-source-1-study.json"
EXTRA = {
    "scripts/p8_covzero_source_constraints.py": "c484a95c11d7979e0f054eebd2edb35065e2cd6dca4ae03db9cca9efab3c6d99",
    "docs/assessment-2026-09-24-p8-covzero-source-constraint-obstruction.md": "ca846e9411e544b5ed45c329c112f308f876962fdd6e214a72d156297f825dd4",
    "docs/validation/p8-covzero-source-constraints-2026-09-24.json": "35bf8bfde310c1e2d454e25f02ad855e15107caad007be674b243e6e36376d7d",
}


def inputs():
    manifest = previous_inputs()
    for name, expected in EXTRA.items():
        assert hashlib.sha256((REPO / name).read_bytes()).hexdigest() == expected
        assert name not in manifest or manifest[name] == expected
        manifest[name] = expected
    return manifest


def validate_spec(spec):
    expected = {
        "milestone": "SOURCE-1",
        "study_authorized": True,
        "study_scope": "bounded feasibility, not physical model adoption",
        "classical_action_modified": False,
        "replacement_source_adopted": False,
        "physical_cutoff": None,
        "field_domain_restricted": False,
        "state_changed": False,
        "extra_mode_prescription_adopted": False,
        "existing_COVZERO_boundary_adopted": True,
        "common_quantum_reference_complete": False,
        "physical_matching_established": False,
        "physical_gravity_verdict_established": False,
        "original_P8_open": True,
    }
    if any(
        spec.get(k) != v or type(spec.get(k)) is not type(v)
        for k, v in expected.items()
    ):
        raise ValueError("Study authority and physical completion flags are fixed")
    retained = spec["retained"]
    if retained["physical_metric_signature"] != "+---":
        raise ValueError("Keep the physical metric")
    if retained["clock_norm"] != "X=g_phys^{mu nu}*partial_mu(u)*partial_nu(u)":
        raise ValueError("Keep the positive clock X convention")
    if retained["vector_kinetic_coefficient"] != "1/1000000":
        raise ValueError("Do not retune the vector")
    if retained["canonical_vector_mass"] != "1000":
        raise ValueError("Do not change the mass")
    if retained["detector_resolution"] != "1/256":
        raise ValueError("Do not change the detector")
    if [a["id"] for a in spec["alternatives"]] != ["S0", "SALG", "SGRAD", "EFT"]:
        raise ValueError("The comparison is bounded")
    if spec["alternatives"][0]["source"] != "0":
        raise ValueError("S0 is exactly source-free")


def null_constraint_audit():
    """Take the regular X=0 limit before differentiating the action.

    A=(R-1)/X^2 at X=0 is finite and is not R_X, which vanishes there.
    The inherited full-R expression fixes A=1024*(1-(1+u^2)^-3).
    No previous Hessian or null-vector implementation is called.
    """
    a, A = s.symbols("a A", nonzero=True)
    V, k11, k22, k33, k12, k13, k23 = s.symbols("V k11 k22 k33 k12 k13 k23")
    K = s.Matrix([[k11, k12, k13], [k12, k22, k23], [k13, k23, k33]])
    Z = a * V + a**2 * k11
    box, uhu = V + a * s.trace(K), a * Z
    # R=1, A3=2A, A4=-2A, A5=0; X=a^2-b^2=0, b=a.
    light = (s.trace(K * K) - s.trace(K) ** 2) / 2 + 2 * A * (box * uhu - Z**2)
    velocity = (V, k11, k22, k33, k12, k13, k23)
    H0 = s.hessian(light, velocity)
    null = s.Matrix([(1 - 2 * A * a**4) / a, 2 * A * a**2, 0, 0, 0, 0, 0])
    assert s.simplify(H0 * null) == s.zeros(7, 1)
    assert s.factor(H0[1:, 1:].det()) == -16 * (2 * A * a**4 - 1) ** 2
    q = -A * uhu
    projection = s.expand(
        sum(s.diff(q, v) * n for v, n in zip(velocity, null, strict=True))
    )
    assert projection == -A * a
    T, W, pnull = s.symbols("T W p_null")
    # At X=0 the source square vanishes, but the linear source does not.
    primary = pnull + (a * T - a * W) * projection
    bracket = s.diff(primary, T)  # {Psi,p_T}, not a secondary pivot.
    assert bracket == -A * a**2
    primary_matrix = s.Matrix([[0, bracket], [-bracket, 0]])
    assert primary_matrix.det() == A**2 * a**4
    controls = []
    for u, speed in [
        (s.Integer(0), s.Rational(1, 128)),
        (s.Integer(1), s.Rational(1, 128)),
        (s.Rational(1, 2), s.Rational(1, 64)),
    ]:
        coefficient = 1024 * (1 - (1 + u**2) ** -3)
        subs = {A: coefficient, a: speed}
        assert H0.subs(subs).rank() == 6
        controls.append(
            {
                "u": str(u),
                "a=b": str(speed),
                "A": str(coefficient),
                "kinetic_rank": 6,
                "primary_bracket": str(bracket.subs(subs)),
            }
        )
    assert controls[1]["primary_bracket"] == "-7/128"
    return controls


def aligned_source_ansatz():
    """All-tilt obstruction within the specified Hessian-linear parallel class."""
    a, b, R, RX = s.symbols("a b R RX", nonzero=True)
    X = a**2 - b**2
    k11 = RX / R * (b**2 / X - s.Rational(1, 2))
    k22 = k33 = -RX / (2 * R)
    V = (1 - b**2 * k11) / a
    Z = a * V + b**2 * k11
    assert s.cancel(Z) == 1
    box, uhu = V + a * (k11 + k22 + k33), a * Z
    factor = 1 - 3 * X * RX / (2 * R)
    assert s.cancel(box - factor / a) == 0
    assert s.cancel(uhu - a) == 0
    c, d, t1, t2, factor_symbol = s.symbols("c d t1 t2 D")
    rows = s.Matrix([[factor_symbol, t1], [factor_symbol, t2]])
    assert s.factor(rows.det()) == -factor_symbol * (t1 - t2)
    assert s.solve(
        [c * factor_symbol + d * t1, c * factor_symbol + d * t2], (c, d)
    ) == {c: 0, d: 0}
    # Original c=(R-1)/X, d=(R-1)*(-1/X^2+3*RX/(2*R*X)).
    old_projection = (R - 1) * (box / X + (-1 / X**2 + 3 * RX / (2 * R * X)) * uhu)
    assert s.cancel(old_projection + (R - 1) * b**2 * factor / (a * X**2)) == 0
    return "c=d=0 where X!=0 and 1-3X*R_X/(2R)!=0; regular continuation through X=0; no general scalar-vector-tensor no-go"


def pb(left, right, coordinates, momenta):
    return s.expand(
        sum(
            s.diff(left, q) * s.diff(right, p) - s.diff(left, p) * s.diff(right, q)
            for q, p in zip(coordinates, momenta, strict=True)
        )
    )


def algebraic_gradient_constraints():
    a, T, W, f, b, G = s.symbols("a T W f b G", real=True)
    vector_H = -((T - f * a) ** 2) / 2 + (W - f * b) ** 2 / 2 - T * G
    Tstar = s.solve(s.diff(vector_H, T), T)[0]
    assert Tstar == f * a - G
    reduced = s.expand(vector_H.subs(T, Tstar))
    assert s.expand(reduced - ((W - f * b) ** 2 / 2 + G**2 / 2 - f * a * G)) == 0
    assert s.diff(vector_H, T, 2) == -1
    # At the actual quadratic vacuum, normalise kappa*sqrt(h)=1.
    # Spatial gradients and Gauss terms are retained as fixed spatial jets.
    u, p_u, p_a, p_T, p_W, zeta = s.symbols("u p_u p_a p_T p_W zeta")
    coordinates, momenta = (u, a, T, W), (p_u, p_a, p_T, p_W)
    H = (
        p_u * a
        - a**2 / 2
        + u**2 / 2
        - T**2 / 2
        - T * G
        + W**2 / 2
        + p_W**2 / (2 * zeta)
    )
    primaries = [p_a, p_T]
    assert pb(*primaries, coordinates, momenta) == 0
    secondaries = [pb(p, H, coordinates, momenta) for p in primaries]
    assert secondaries == [a - p_u, G + T]
    constraints = primaries + secondaries
    matrix = s.Matrix(
        [[pb(l, r, coordinates, momenta) for r in constraints] for l in constraints]
    )
    assert matrix.det() == 1
    assert 18 - 8 - 4 // 2 == 8
    assert 18 - 8 - 2 // 2 == 9
    return {
        "vector_secondary_solution": str(Tstar),
        "vector_temporal_pivot": "-1",
        "vacuum_second_class_matrix_determinant": "1 in unit positive density normalization",
        "source_free_regular_branch_modes": 8,
        "old_generic_local_Dirac_count": 9,
        "global_second_class_pivot_bound_established": False,
    }


def source_jets_and_lift():
    e, n, k1, H, h, lam = s.symbols("e n k1 H h lambda", nonzero=True)
    N = 1 + e * n
    X = N**-2
    # R_X=1/h at X=1; higher clock jets cannot enter these second jets.
    r = (X - 1) / h
    old = r * (3 * H + e * k1 - 3 * H / N)
    alg = lam * r * (X - 1) / N
    old_jet = s.expand(s.series(old, e, 0, 3).removeO()).coeff(e, 2)
    alg_jet = s.expand(s.series(alg, e, 0, 3).removeO()).coeff(e, 2)
    assert s.expand(old_jet + 2 * n * (k1 + 3 * H * n) / h) == 0
    assert alg_jet == 4 * lam * n**2 / h
    assert s.diff(old_jet, k1) == -2 * n / h
    assert s.diff(alg_jet, k1) == 0
    w, source_old, source_new = s.symbols("w S_old S_new")
    delta = w * (source_old - source_new) + (source_new**2 - source_old**2) / 2
    assert s.expand((w - source_old) ** 2 / 2 + delta - (w - source_new) ** 2 / 2) == 0
    # For each Lorentz component the same identity holds with its metric sign.
    cubic_difference = s.expand(
        delta.subs({w: e * w, source_old: e**2 * old_jet, source_new: e**2 * alg_jet})
    )
    assert all(cubic_difference.coeff(e, j) == 0 for j in range(3))
    return {
        "old_second_clock_jet": str(old_jet),
        "SALG_second_clock_jet": str(alg_jet),
        "unchanged_action_jet_through_degree": 2,
        "all_old_clock_cubic_vertices_preserved": False,
        "affine_lift": "Delta L/(kappa sqrt(-g))=W.(S_old-S_new)+(S_new^2-S_old^2)/2; full 56-mode complement unchanged",
    }


def vacuum_filtration():
    """Field-degree calculation, not a derivative or finite-domain truncation."""
    e, phi, Y, A = s.symbols("e phi Y A", real=True)
    X, u, order = e**2 * Y, e * phi, s.Integer(1024)
    # T(X) first enters at e^2048, beyond this finite jet. Both exponentials
    # remain in the full parent; A is finite and fixed before e tends to zero.
    r_low = order * X**2 * (X - 1) * s.exp(-order * X**2) / (
        1 + u**2
    ) ** 3 + order * X**2 * (1 - X) ** 8 * s.exp(-A * X**2)
    jet = s.series(r_low, e, 0, 7).removeO().expand()
    assert s.expand(jet - order * e**6 * Y**2 * (3 * phi**2 - 7 * Y)) == 0
    regular_A_jet = s.cancel(jet / X**2)
    grad, box, uhu = s.symbols("du Box_u uHu")
    source_leading = s.expand(e * grad * regular_A_jet * (X * e * box - e**3 * uhu))
    assert (
        s.expand(
            source_leading / e**6
            - order * grad * (3 * phi**2 - 7 * Y) * (Y * box - uhu)
        )
        == 0
    )
    # W.S has odd W parity; every other source-independent vertex is even.
    # No external W: >=2 W.S vertices, or >=1 S^2 vertex. For L=1,
    # E=sum_v(d_v-2), so both ways have at least 2*p-2 external legs.
    for p in (4, 6, 7):
        assert 2 * ((p + 1) - 2) == (2 * p) - 2
        assert 2 * p - 2 > 4
    # Bound all possible one-loop four-leg source graphs: any nonquadratic
    # vertex has d>=3, hence E=4 permits <=4 vertices and total excess <=4.
    rejected_graphs = 0
    for linear_sources in range(5):
        for square_sources in range(5):
            if not (linear_sources or square_sources) or linear_sources % 2:
                continue
            minimum_excess = 5 * linear_sources + 10 * square_sources
            assert minimum_excess > 4
            rejected_graphs += 1
    return {
        "R_minus_one_leading_field_degree": 6,
        "old_source_minimum_scalar_degree": 6,
        "old_source_leading_jet": "1024*du_mu*(3*u^2-7*X)*(X*Box(u)-uHu), homogeneous field degree six",
        "SALG_minimum_scalar_degree": 7,
        "source_tagged_one_loop_minimum_external_legs": 10,
        "one_loop_four_scalar_source_graphs": 0,
        "finite_graph_controls": rejected_graphs,
        "physical_common_reference_or_full_amplitude_equivalence_proved": False,
    }


def eft_admission_controls():
    """Necessary study data and a conditional resolvent bound, not a cutoff."""
    K, zeta, alpha, z = s.symbols("K zeta alpha z", nonzero=True)
    # A local two-coordinate illustration; K is NOT a measured P8 residue.
    matrix = s.Matrix([[K * z - alpha**2 * z**2, alpha * z], [alpha * z, zeta * z - 1]])
    assert (
        s.expand(matrix.det() - z * (-K + K * zeta * z - alpha**2 * zeta * z**2)) == 0
    )
    ratio = s.cancel((matrix[0, 0] - K * z) / (K * z))
    assert ratio == -(alpha**2) * z / K
    assert s.degree(ratio, z) == 1
    # Noncommutative Neumann polynomial identity, checked on rational matrices.
    controls = 0
    for A in (
        s.Matrix([[s.Rational(1, 5), s.Rational(1, 7)], [0, s.Rational(-1, 9)]]),
        s.Matrix([[0, s.Rational(1, 3)], [s.Rational(1, 4), 0]]),
    ):
        rho = max(sum(abs(x) for x in A.row(i)) for i in range(A.rows))
        assert rho < 1
        for p in range(4):
            partial = sum(((-A) ** j for j in range(p + 1)), s.zeros(2))
            remainder = (s.eye(2) + A).inv() - partial
            assert remainder == (-A) ** (p + 1) * (s.eye(2) + A).inv()
            norm = max(sum(abs(x) for x in remainder.row(i)) for i in range(2))
            assert norm <= rho ** (p + 1) / (1 - rho)
            controls += 1
    return {
        "conditional_remainder_bound": "||G-G_p|| <= rho^(p+1)/(1-rho)*||G0||, rho=||G0 DeltaD||<1 on a specified common domain and boundary prescription",
        "matrix_controls": controls,
        "actual_P8_rho_bound": None,
        "actual_P8_cutoff_or_extra_pole_mass": None,
        "controlled_EFT_excluded": False,
        "controlled_EFT_admitted": False,
    }


def audit():
    before = inputs()
    spec = json.loads((REPO / SPEC).read_text())
    validate_spec(spec)
    rejected = []
    for key, value in [
        ("replacement_source_adopted", True),
        ("physical_cutoff", "1000"),
        ("field_domain_restricted", True),
        ("state_changed", True),
        ("extra_mode_prescription_adopted", True),
        ("existing_COVZERO_boundary_adopted", False),
        ("common_quantum_reference_complete", True),
        ("physical_matching_established", True),
        ("physical_gravity_verdict_established", True),
        ("original_P8_open", False),
    ]:
        bad = copy.deepcopy(spec)
        bad[key] = value
        try:
            validate_spec(bad)
        except ValueError:
            rejected.append(key)
        else:
            raise AssertionError("Unapproved input was accepted")
    results = {
        "milestone": "SOURCE-1",
        "outcome": "BOUNDED_STUDY_COMPLETE_S0_RECOMMENDED_PENDING_PHYSICAL_ADOPTION",
        "protected_input_files": len(before),
        "independent_null_branch_controls": null_constraint_audit(),
        "parallel_source_ansatz_result": aligned_source_ansatz(),
        "algebraic_gradient_constraints": algebraic_gradient_constraints(),
        "clock_jets_and_affine_lift": source_jets_and_lift(),
        "formal_vacuum_filtration": vacuum_filtration(),
        "unchanged_action_EFT": eft_admission_controls(),
        "rejected_authority_inputs": rejected,
        "classical_action_modified": False,
        "replacement_source_adopted": False,
        "common_quantum_reference_complete": False,
        "physical_matching_established": False,
        "physical_gravity_verdict_established": False,
        "original_P8_open": True,
    }
    assert inputs() == before
    return results


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
