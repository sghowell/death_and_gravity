"""Independent tensor contractions, compact functional variations and scope controls."""

import json
from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_affine import dictionary
from p8_affine_vacuum_domain import family
from p8_exceptional_vacuum import analytic
from p8_vacuum_canonical_bounce_gap import audit as previous
from p8_vacuum_target_clock_null import audit, homogeneous, structural, target, verify


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_named_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_unsupported_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_direct_metric_connection_and_all_clock_scalar_contractions():
    a, N = s.symbols("a N", positive=True)
    ad, Nd, v = s.symbols("adot Ndot v", real=True)
    g = s.diag(N * N, -a * a, -a * a, -a * a)
    gi = g.inv()
    dt = lambda f: s.diff(f, a) * ad + s.diff(f, N) * Nd
    partial = lambda f, j: dt(f) if j == 0 else s.Integer(0)
    connection0 = s.Matrix(
        4,
        4,
        lambda i, j: s.simplify(
            sum(
                gi[0, k]
                * (partial(g[k, j], i) + partial(g[k, i], j) - partial(g[i, j], k))
                / 2
                for k in range(4)
            )
        ),
    )
    phi2 = -connection0
    upper = gi * s.Matrix([1, 0, 0, 0])
    contracted = phi2 * upper
    box = s.trace(gi * phi2)
    Z = (upper.T * contracted)[0]
    L3 = s.simplify(Z * box)
    L4 = s.simplify((contracted.T * gi * contracted)[0])
    L5 = s.simplify(Z * Z)
    assert s.simplify(box.subs(ad, a * v) + Nd / N**3 - 3 * v / N**2) == 0
    assert s.simplify(Z + Nd / N**5) == 0
    assert s.simplify(L3.subs(ad, a * v) - Nd * Nd / N**8 + 3 * v * Nd / N**7) == 0
    assert s.simplify(L4 - Nd * Nd / N**8) == 0
    assert s.simplify(L5 - Nd * Nd / N**10) == 0


@pytest.mark.parametrize(
    "time", [-10, -2, -s.Rational(1, 2), 0, s.Rational(1, 2), 2, 10]
)
def test_actual_frozen_witness_and_analytic_clock_jets_at_signed_times(time):
    raw = json.loads(dictionary.WITNESS.read_text())["covariant_functions"]
    u, X = family.u, family.X
    F = s.sympify(raw["F"], locals={"t": u, "X": X})
    d = target.data()
    assert (
        s.simplify(
            F.subs({u: time, X: 1}) - d["full_analytic_target_clock_F0"].subs(u, time)
        )
        == 0
    )
    assert (
        s.simplify(
            s.diff(F, X).subs({u: time, X: 1})
            - d["full_analytic_target_clock_FX"].subs(u, time)
        )
        == 0
    )
    table = d["complete_first_metric_variation_table"]
    for key in ("rho_Euler", "P_Euler", "null_Euler"):
        assert sum(row[key].subs(u, time) for row in table.values()) == 0


@cache
def finite_action_functions():
    u = family.u
    d = target.data()
    functions = {
        "F0": s.lambdify(u, d["full_analytic_target_clock_F0"], "mpmath"),
        "FX": s.lambdify(u, d["full_analytic_target_clock_FX"], "mpmath"),
    }
    table = d["complete_first_metric_variation_table"]
    rows = {
        name: dict(row)
        for name, row in table.items()
        if name not in ("canonical_clock", "scalar_F_minus_canonical")
    }
    rows["scalar_F"] = {
        key: table["canonical_clock"][key] + table["scalar_F_minus_canonical"][key]
        for key in ("rho_Euler", "P_Euler", "null_Euler")
    }
    functions["variations"] = {
        name: {key: s.lambdify(u, value, "mpmath") for key, value in row.items()}
        for name, row in rows.items()
    }
    return functions


@pytest.mark.parametrize("kind", ["lapse", "scale"])
@pytest.mark.parametrize(
    "piece",
    ["Einstein", "scalar_F", "curvature_nonEinstein", "A3", "A4", "A5", "original_M1"],
)
def test_compact_functional_metric_variation_of_complete_clock_jet_packet(kind, piece):
    with mp.workdps(80):
        fn = finite_action_functions()

        def bump(t):
            return (1 - 4 * t * t) ** 4

        def one(t):
            return -32 * t * (1 - 4 * t * t) ** 3

        def two(t):
            return -32 * (1 - 4 * t * t) ** 3 + 768 * t * t * (1 - 4 * t * t) ** 2

        def density(t, epsilon):
            D = 1 + t * t
            h = D**3
            a = D * D * mp.exp(epsilon * bump(t) if kind == "scale" else 0)
            N = 1 + epsilon * bump(t) if kind == "lapse" else mp.mpf(1)
            Nd = epsilon * one(t) if kind == "lapse" else mp.mpf(0)
            v = 4 * t / D + (epsilon * one(t) if kind == "scale" else 0)
            vd = 4 * (1 - t * t) / D**2 + (epsilon * two(t) if kind == "scale" else 0)
            X = 1 / (N * N)
            R = -6 * (vd / N**2 - v * Nd / N**3 + 2 * v * v / N**2)
            L3 = (-Nd / N**5) * (-Nd / N**3 + 3 * v / N**2)
            L4 = Nd * Nd / N**8
            L5 = Nd * Nd / N**10
            # Exact first-clock-jet packet; the arbitrary extra scalar
            # quadratic term exercises insensitivity to unused second jets.
            Rcoef = 1 + (X - 1) / h
            kernels = {
                "Einstein": -R / 2,
                "scalar_F": fn["F0"](t)
                + fn["FX"](t) * (X - 1)
                + mp.mpf(3) / (7 * h * h) * (X - 1) ** 2,
                "curvature_nonEinstein": -(X - 1) * R / (2 * h),
                "A3": L3 / (h * X),
                "A4": (-1 / (h * X) - 7 / (4 * h * h * Rcoef)) * L4,
                "A5": L5 / (h * h * Rcoef * X),
                # Field held fixed during variation, not reset to C/a^3.
                "original_M1": 1 / (200 * D**12 * N * N),
            }
            return N * a**3 * kernels[piece]

        interval = [-mp.mpf(".5"), 0, mp.mpf(".5")]
        e = mp.mpf("1e-12")
        difference = mp.quad(
            lambda t: (density(t, e) - density(t, -e)) / (2 * e), interval
        )
        key = "rho_Euler" if kind == "lapse" else "P_Euler"
        factor = -1 if kind == "lapse" else 3
        expected = factor * mp.quad(
            lambda t: (1 + t * t) ** 6 * fn["variations"][piece][key](t) * bump(t),
            interval,
        )
        assert abs(difference - expected) < mp.mpf("1e-18") * max(1, abs(expected))
        if piece in ("curvature_nonEinstein", "A3") and kind == "lapse":
            assert abs(expected) > mp.mpf(".1")


def test_actual_bounce_source_balance_and_omission_controls():
    table = target.data()["actual_bounce_table"]
    expected = {
        "Einstein": 8,
        "canonical_clock": 1,
        "scalar_F_minus_canonical": -s.Rational(2101, 100),
        "curvature_nonEinstein": 24,
        "A3": -12,
        "A4": 0,
        "A5": 0,
        "original_M1": s.Rational(1, 100),
    }
    assert {name: row["null_Euler"] for name, row in table.items()} == expected
    assert sum(expected.values()) == 0
    assert sum(
        value for name, value in expected.items() if name != "original_M1"
    ) == -s.Rational(1, 100)
    assert (
        sum(
            value for name, value in expected.items() if name != "curvature_nonEinstein"
        )
        == -24
    )
    assert sum(value for name, value in expected.items() if name != "A3") == 12
    extra = (
        expected["scalar_F_minus_canonical"]
        + expected["curvature_nonEinstein"]
        + expected["A3"]
    )
    assert extra == -s.Rational(901, 100)
    assert analytic.KAPPA * extra < -9 * analytic.KAPPA


@pytest.mark.parametrize("h", [1, 2, 8, 100])
def test_zero_background_terms_have_nonzero_constraint_role(h):
    B = -s.Rational(1, 2)
    fx = -s.Rational(1, 2 * h)
    A3 = s.Rational(1, h)
    A4 = -A3 - s.Rational(7, 4 * h * h)
    A5 = s.Rational(1, h * h)
    C = 4 * fx + A3
    D = A3 + A4 + A5
    K, V = s.symbols("K V")
    complete = s.hessian(2 * B * K * K / 3 + C * K * V + D * V * V, (K, V))
    incomplete = s.hessian(2 * B * K * K / 3 + C * K * V + A3 * V * V, (K, V))
    assert complete.rank() == 1 and complete.det() == 0
    assert incomplete.rank() == 2 and incomplete.det() < 0
    assert A4 != 0 and A5 != 0
    # Each deleted term has zero first jet in lapse velocity but nonzero second jet.
    for coeff in (A4, A5):
        assert s.diff(coeff * V * V, V).subs(V, 0) == 0
        assert s.diff(coeff * V * V, V, 2) != 0


@pytest.mark.parametrize("scale,momentum", [("1", "0.1"), ("2", "3"), ("5", "-2")])
def test_independent_matter_variation_and_fixed_momentum_Routhian(scale, momentum):
    with mp.workdps(70):
        a, p = mp.mpf(scale), mp.mpf(momentum)
        w = p / a**3
        original = lambda lapse, b: a**3 * mp.exp(3 * b) * w * w / (2 * lapse)
        direct_rho = -mp.diff(lambda lapse: original(lapse, 0), 1) / a**3
        direct_P = mp.diff(lambda b: original(1, b), 0) / (3 * a**3)
        substituted = lambda lapse, b: lapse * p * p / (2 * a**3 * mp.exp(3 * b))
        wrong_rho = -mp.diff(lambda lapse: substituted(lapse, 0), 1) / a**3
        wrong_P = mp.diff(lambda b: substituted(1, b), 0) / (3 * a**3)
        routh = lambda lapse, b: -substituted(lapse, b)
        right_rho = -mp.diff(lambda lapse: routh(lapse, 0), 1) / a**3
        right_P = mp.diff(lambda b: routh(1, b), 0) / (3 * a**3)
        assert direct_rho > 0 and direct_P > 0
        assert abs(direct_rho - w * w / 2) < mp.mpf("1e-60")
        assert abs(direct_P - w * w / 2) < mp.mpf("1e-60")
        assert abs(wrong_rho + direct_rho) < mp.mpf("1e-60")
        assert abs(wrong_P + direct_P) < mp.mpf("1e-60")
        assert abs(right_rho - direct_rho) < mp.mpf("1e-60")
        assert abs(right_P - direct_P) < mp.mpf("1e-60")


@pytest.mark.parametrize("kappa", [2, 4, 100, 10**800])
def test_vacuum_unit_first_jet_class_separation(kappa):
    k = s.Integer(kappa)
    gap = structural.normalized_first_jet_gap(k)
    assert s.simplify(gap - (s.sqrt(k) - 1) / s.sqrt(k)) == 0
    for unit_derivative in (-1, s.Rational(-1, 3), 0, s.Rational(1, 3), 1):
        assert s.simplify(abs(s.sqrt(k) - unit_derivative) / s.sqrt(k) - gap) >= 0
    if k == analytic.KAPPA:
        assert gap > s.Rational(1, 2)
        assert 4 / k < s.Rational(1, 10**798)
        assert s.sqrt(k) - 1 > 10**399


def test_scope_ancestry_and_all_actual_payloads_exact():
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert "NOT_PROPAGATING_PARENT_OR_QUANTUM_BOUNCE" in audit.matching()[-1]["status"]
    assert (
        "not proposed global replacements" in homogeneous.data()["representative_scope"]
    )
    assert "not independently conserved" in target.data()["interpretation"]
    assert (
        "does not construct a healthy propagating UV parent"
        in target.data()["classical_scope"]
    )
    assert "even on a finite interval" in structural.data()["domain_warning"]
    assert target.data()["normalized_by_physical_kappa"] == 10**800
    for mod in audit.MODULES:
        json.dumps(verify.serialize(verify.payload(mod.data())))
    with pytest.raises(ValueError, match="Inexact or nonfinite"):
        verify.serialize(s.Float(1))
