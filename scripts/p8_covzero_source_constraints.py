"""Coupled-source constraint admission test, not a quantum or P8 certificate.

The literal source lifts the light/metric kinetic null direction away from
the scalar-aligned slice. Exact parent-germ bounds place a spacelike example
inside the retained coefficient domain. This rejects a uniform fixed-count
canonical-measure shortcut, not every perturbative EFT interpretation.
"""

import ast
import hashlib
import json
from fractions import Fraction as F

import sympy as s
from p8_covzero_reference import inputs as previous_inputs
from p8_match1_rate_input import REPO

BASE_SOURCE = (
    "problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/"
    "quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/"
    "dispersion/liouville/response/inverse/curved/frozen/auxiliary/momentum/vertices/"
    "control/neighborhood/alignment/proca/stress/cone/response/nonlocal/dispersion/"
    "curved/nearby/propagation/packets/retarded/compact/quantum/sharp/spectators/"
    "retuned/support/global/hadamard/qei/finite_width/physical_matter/seven_modes/"
    "mean_response/tensor_source/scalar_source/global_response/vacuum/affine_domain/"
    "src/p8_affine_vacuum_domain/family.py"
)
EXTRA = {
    BASE_SOURCE: "4fb427765b60d8b1f19bb95b2637e0eced85dbcef5e089b34de06fc248eab1c7",
    "scripts/p8_covzero_reference.py": "ee114ef55e8a2a5cad89ed094543df6a9b3db7195b407892616553cda97bf0c9",
    "docs/assessment-2026-09-24-p8-covzero-reference-seed.md": "df7b964284f05f4dd26df0919826ba93bb65330167b6b9543b222cc67cfbd131",
    "docs/candidates/p8-rate4-covzero-1.json": "015e18f4f04653bc39ef41d8858927280177aade4981472654956c786bbd9b38",
}
AFFINE_REPORT = "problems/P8/s6/continuation/s6_174/certificates/polynomial-vacuum-analytic-affine-parent.json"
PARENT_REPORT = "problems/P8/s6/continuation/s6_238/certificates/polynomial-vacuum-affine-heavy-scalar-parent.json"


def inputs():
    manifest = previous_inputs()
    for name, expected in EXTRA.items():
        assert hashlib.sha256((REPO / name).read_bytes()).hexdigest() == expected
        assert name not in manifest or manifest[name] == expected
        manifest[name] = expected
    return manifest


def fraction_det(rows):
    """Independent rational elimination for the finite Hessian controls."""
    rows = [[F(x) for x in row] for row in rows]
    determinant = F(1)
    for j in range(len(rows)):
        pivot = next((i for i in range(j, len(rows)) if rows[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            rows[pivot], rows[j] = rows[j], rows[pivot]
            determinant *= -1
        value = rows[j][j]
        determinant *= value
        for i in range(j + 1, len(rows)):
            factor = rows[i][j] / value
            for k in range(j + 1, len(rows)):
                rows[i][k] -= factor * rows[j][k]
    return determinant


def principal_forms():
    """Reconstruct the covariant contractions in a local physical ADM frame."""
    a, b, R, RX = s.symbols("a b R R_X", real=True)
    V = s.Symbol("V_star", real=True)
    entries = s.symbols("k0:6", real=True)
    K = s.Matrix(
        [
            [entries[0], entries[1], entries[2]],
            [entries[1], entries[3], entries[4]],
            [entries[2], entries[4], entries[5]],
        ]
    )
    eta = s.diag(1, -1, -1, -1)
    gradient = s.Matrix([a, b, 0, 0])
    raised = eta * gradient
    X = (gradient.T * raised)[0]
    spatial = gradient[1:, :]
    cross = -K * spatial
    hessian = s.Matrix([[V]]).row_join(cross.T).col_join(cross.row_join(-a * K))
    box = s.trace(eta * hessian)
    contracted = hessian * raised
    uhu = (raised.T * contracted)[0]
    L4 = (contracted.T * eta * contracted)[0]
    Z = a * V + b**2 * K[0, 0]
    assert s.expand(box - V - a * s.trace(K)) == 0
    assert s.expand(uhu - a * Z) == 0
    assert s.expand(L4 - Z**2) == 0
    A3, A4, A5 = RX / X, -RX / X - 7 * RX**2 / (4 * R), RX**2 / (R * X)
    light = (
        R * (s.trace(K * K) - s.trace(K) ** 2) / 2
        - 2 * RX * s.trace(K) * Z
        + A3 * box * uhu
        + A4 * L4
        + A5 * uhu**2
    )
    retained = json.loads((REPO / AFFINE_REPORT).read_text())[
        "regular_classical_parent_and_full_algebraic_reduction"
    ]["source"]["shifted_source_P8"]
    q = s.sympify(
        retained,
        locals={
            "R": R,
            "R_X": RX,
            "R_u": 0,
            "H_clock": 0,
            "positive_X": X,
            "Box_P8_u": box,
            "uHu_P8": uhu,
        },
    )
    # R_u/H_clock terms have no principal velocity; this is not a state reset.
    source = gradient * q
    coupled = light + (source.T * eta * source)[0] / 2
    assert s.expand(coupled - light - X * q**2 / 2) == 0
    return (a, b, R, RX), (*entries, V), light, coupled, q


def parent_binding():
    """Read only the six pinned algebraic assignments, not ancestor imports."""
    X, u, A = s.symbols("X u A", real=True)
    tree = ast.parse((REPO / BASE_SOURCE).read_text())
    function = next(
        n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "data"
    )
    expected = ("T", "bump", "B", "S", "h", "R")
    bindings = {"sp": s, "N": s.Integer(1024), "X": X, "u": u}
    for statement, name in zip(function.body, expected, strict=False):
        assert isinstance(statement, ast.Assign)
        assert len(statement.targets) == 1 and statement.targets[0].id == name
        bindings[name] = eval(
            compile(ast.Expression(statement.value), BASE_SOURCE, "eval"),
            {"__builtins__": {}},
            bindings,
        )
    original = bindings["R"]
    report = json.loads((REPO / PARENT_REPORT).read_text())
    family = report["literal_new_parent_vacuum_and_full_clock_coefficient_tube"][
        "family"
    ]
    assert s.sympify(family["fixed_localizer_parameter"]) == s.Integer(10) ** 420
    addition = s.sympify(
        family["normalized_additions"]["delta_R"],
        locals={"physical_positive_clock_X": X, "positive_localizer_A": A},
    )
    assert s.expand(addition - 1024 * X**2 * (1 - X) ** 8 * s.exp(-A * X**2)) == 0
    N = s.Integer(1024)
    expected_R = (
        1
        + (bindings["T"] + (1 - bindings["T"]) * N * X**2 * s.exp(-N * X**2))
        * (X - 1)
        / (1 + u**2) ** 3
        + addition
    )
    assert s.expand(original + addition - expected_R) == 0
    # T=O(X^1024). Removing T here is a finite jet calculation only.
    low = 1 + N * X**2 * (X - 1) * s.exp(-N * X**2) + addition
    jet = s.series(low, X, 0, 5).removeO().expand()
    assert jet.coeff(X, 2) == 0
    assert jet.coeff(X, 3) == -7 * N
    assert jet.coeff(X, 4) == N * (N + 28 - A)
    return N


def actual_spacelike_interval():
    """Rational side conditions for the written all-epsilon exponential bounds.

    No tiny exponential is numerically rounded to zero. The proof uses
    1-z<=exp(-z)<=1, binomial bounds and T<epsilon^1024<epsilon^4.
    Monotonicity extends the endpoint inequalities to 0<epsilon<=1e-430.
    """
    e, A, N = F(1, 10**430), 10**420, 1024
    assert 0 < e < F(1, 4096)
    assert A * e**2 < 1 and N * e**2 < 1
    assert (1 + e) ** 8 < 2
    assert sum(F(s.binomial(8, k)) * e ** (k - 1) for k in range(2, 9)) < 1
    # The new exponential loses less than epsilon from its lower bound.
    assert A * e * (1 + 8 * e) < 1
    # Lower delta: 6N e^3 - 2e^4 > 5N e^3.
    assert 2 * e < N
    # Upper delta: N e^2 [8e+N e^2(1+e)] < 10N e^3.
    assert N * e * (1 + e) < 2
    # |R_old,X|<5N e and |delta_R,X|<5N e.
    assert e + e**3 / N + 2 * e**2 < 1
    assert 16 * e + 4 * A * e**2 < 1
    # Therefore 2R-3X R_X > 2-30N e^2 > 1.
    assert 30 * N * e**2 < 1
    assert 1 + 10 * N * e**3 < F(6, 5)
    return {
        "domain": "u=0, X=-epsilon, 0<epsilon<=1e-430, actual kappa=kappa0 parent",
        "R_minus_one": "5*1024*epsilon^3 < R-1 < 10*1024*epsilon^3",
        "R_X": "abs(R_X)<10*1024*epsilon",
        "metric_pivot": "2R-3X*R_X>1",
        "scalar_acceleration_hessian": "-(R-1)^2/epsilon < -25*1024^2*epsilon^5 < 0",
        "scope": "actual off-shell coefficient-domain neighborhood, not a demonstrated on-shell background, cutoff, pole mass or instability time",
    }


def admit_fixed_count_measure(*, same_primary_rank, proca_secondary_compensates):
    """The demonstrated rank loss cannot be silently called measure transport."""
    if (
        type(same_primary_rank) is not bool
        or type(proca_secondary_compensates) is not bool
    ):
        raise ValueError("Explicit constraint evidence flags are required")
    if not same_primary_rank and not proca_secondary_compensates:
        raise ValueError(
            "The proposed uniform fixed-count canonical extension is not established"
        )


def audit():
    before = inputs()
    (a, b, R, RX), velocity, light, coupled, q = principal_forms()
    X = a**2 - b**2
    H0, H1 = s.hessian(light, velocity), s.hessian(coupled, velocity)
    gradient = s.Matrix([s.diff(q, v) for v in velocity])
    assert all(s.cancel(v) == 0 for v in H1 - H0 - X * gradient * gradient.T)

    # Null velocity is normalized by Z=1. Written invariant derivation
    # explains why these scalar identities hold for arbitrary nonzero a,X.
    tilt = b**2 / X
    K11 = RX * (tilt - s.Rational(1, 2)) / R
    K22 = -RX / (2 * R)
    Vnull = (1 - b**2 * K11) / a
    null = s.Matrix([K11, 0, 0, K22, 0, K22, Vnull])
    assert all(s.cancel(v) == 0 for v in H0 * null)
    qnull = (R - 1) * b**2 * (-1 + 3 * X * RX / (2 * R)) / (a * X**2)
    assert s.cancel((gradient.T * null)[0] - qnull) == 0
    determinant = -4 * R**4 * (R - 1) ** 2 * b**4 * (2 * R - 3 * X * RX) ** 2 / X**3
    controls = []
    for av, bv, rv, rxv, expected_rank in (
        (F(1), F(0), F(7, 6), F(1, 9), 6),
        (F(1), F(1, 2), F(7, 6), F(1, 9), 7),
        (F(1, 4), F(1, 2), F(7, 6), F(1, 9), 7),
        (F(0), F(1, 128), F(7, 6), F(1, 9), 7),
        (F(1), F(1, 2), F(1), F(1), 6),
    ):
        point = dict(
            zip((a, b, R, RX), map(s.Rational, (av, bv, rv, rxv)), strict=True)
        )
        h0, h1 = H0.subs(point), H1.subs(point)
        assert h0.rank() == 6 and h1.rank() == expected_rank
        exact_det = fraction_det(h1.tolist())
        assert s.Rational(exact_det) == h1.det() == determinant.subs(point)
        controls.append(
            {
                "a": str(av),
                "b": str(bv),
                "R": str(rv),
                "R_X": str(rxv),
                "light_rank": 6,
                "coupled_rank": expected_rank,
                "coupled_determinant": str(exact_det),
            }
        )

    # Pure spatial gradient: no division by the vanishing a is used.
    spatial_light = H0.subs(a, 0)
    assert spatial_light[:, -1] == s.zeros(7, 1)
    metric_det = -4 * R**4 * (2 * R + 3 * b**2 * RX) ** 2
    assert s.factor(spatial_light[:6, :6].det() - metric_det) == 0
    assert s.cancel(H1[-1, -1].subs(a, 0) + (R - 1) ** 2 / b**2) == 0
    assert s.cancel(q.subs(a, 0) + (R - 1) * velocity[-1] / b**2) == 0

    # Full temporal vector check after Legendre transformation. Maxwell
    # gives a Gauss multiplier linear in T, not an extra T^2 cancellation.
    T, W, y, py, G = s.symbols("T W y p_y Gauss", real=True)
    J = a * T - b * W
    source_mass = X * y**2 / 2 - J * y + (T**2 - W**2) / 2
    ysol = (py + J) / X
    reduced = s.cancel((py * y - source_mass).subs(y, ysol) + T * G)
    assert s.cancel(s.diff(reduced, T, 2) - b**2 / X) == 0
    assert s.diff(reduced, T, 2).subs(a, 0) == -1

    assert parent_binding() == 1024
    interval = actual_spacelike_interval()
    rejected = 0
    for rank, compensates in ((False, False), (None, False), (False, None), (0, False)):
        try:
            admit_fixed_count_measure(
                same_primary_rank=rank, proca_secondary_compensates=compensates
            )
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("Unsupported fixed-count measure extension accepted")
    assert inputs() == before
    return {
        "milestone": "RATE4-COVZERO-1.COUPLED-SOURCE-CONSTRAINT-ADMISSION",
        "outcome": "UNIFORM_FIXED_COUNT_CANONICAL_EXTENSION_OBSTRUCTED_BY_RETAINED_CLASSICAL_SOURCE",
        "protected_input_files": len(before),
        "physical_boundary_rule_adopted": True,
        "classical_action_modified": False,
        "common_quantum_reference_complete": False,
        "coupled_kinetic_determinant": "-4*R^4*(R-1)^2*b^4*(2R-3X*R_X)^2/X^3; X=a^2-b^2",
        "generic_rank_controls": controls,
        "pure_spacelike_acceleration_entry": "(R-1)^2/X, nonzero and negative for X<0,R!=1",
        "temporal_Proca_Hamiltonian_Hessian": "b^2/X; nonzero on the rank-lifted branch, -1 at a=0",
        "actual_parent_spacelike_bounds": interval,
        "actual_parent_u0_germ": "R-1=-7168*X^3+1024*(1052-10^420)*X^4+O(X^5); only local jet, not a finite-domain replacement",
        "unchanged": "original clock/unitary constraint count and flat vacuum quadratic spectrum; new pure-light dimensional seed and Proca scale conversion",
        "rejected_inputs": rejected,
        "scope": "local principal and temporal-constraint obstruction to one exact fixed-count canonical measure on the whole regular coefficient domain; not a timed on-shell instability, a no-go for a controlled order-reduced EFT, or exclusion of original P8",
        "needed_before_common_reference": "justify a controlled EFT/order-reduction and domain treatment or investigate a separately authorized classical source completion; no repair, cutoff or field-domain restriction adopted",
        "physical_gravity_verdict_established": False,
        "original_P8_open": True,
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
