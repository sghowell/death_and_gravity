"""Exact family, independent local algebra, whole-domain bounds and controls."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_exceptional_vacuum import (
    analytic,
    audit,
    calibration,
    family,
    heavy,
    transition,
    uniform,
    vacuum,
)
from p8_uv import vacuum as original


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_native_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_exact_bound_and_written_proof_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[row[0] for row in calibration.bad_cases()],
)
def test_unsupported_exact_input_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_and_scope_control(name, value):
    assert bool(value), name


def test_actual_vacuum_constant_and_canonical_calibration():
    d = analytic.local_jets()
    assert d["original_retuned_F_at_vacuum"] == -sp.Rational(28, 25)
    assert d["actual_lambda_effective"] == sp.Rational(1, 10**600)
    assert (
        d["leading_canonical_L3_minus_L4_coefficient"]
        == -2
        * d["canonical_switch_order_symbol"]
        / d["canonical_action_normalization_symbol"]
    )
    assert d["canonical_potential_quartic_coefficient"] == -d[
        "canonical_switch_order_symbol"
    ] / (3 * d["canonical_action_normalization_symbol"])


def test_exact_smooth_plateau_has_no_DHOST_coefficients():
    zero = family.constant_branch(0)
    assert zero["F2"] == -sp.Rational(1, 2)
    assert all(zero[key] == 0 for key in ("K", "A1", "A2", "A3", "A4", "A5"))


def test_explicit_transition_jets_and_naive_failure():
    d = transition.data()
    assert [
        d[key]
        for key in (
            "transition_center_R",
            "transition_center_RX",
            "transition_center_correct_A3",
            "transition_center_naive_A3",
        )
    ] == [sp.Rational(3, 4), -sp.Rational(7, 2), -7, 1]


def test_analytic_clock_vanishing_order_is_not_open_tube_equality():
    X = family.X
    assert sp.factor((1 - X * X) / (X - 1)) == -X - 1
    assert (-X - 1).subs(X, 1) == -2
    assert (1 - sp.Rational(9, 10) ** 2) ** analytic.ORDER > 0


@pytest.mark.parametrize("j", range(7))
def test_independent_even_order_switch_majorant_ratio(j):
    def C(n):
        return sum(
            sp.factorial(j)
            * sp.prod(n - i for i in range(j - b))
            * sp.Rational(21, 100) ** (n - j + b)
            * sp.Rational(11, 5) ** (j - 2 * b)
            / (sp.factorial(j - 2 * b) * sp.factorial(b))
            for b in range(j // 2 + 1)
        )

    assert C(1024) == analytic.derivative_upper(j)
    assert sp.Rational(1026, 1024) * C(1026) < C(1024) / 20


def test_two_exact_family_points_hold_canonical_interactions_fixed():
    a, b = calibration.order_point(), calibration.order_point(1026)
    assert a["kappa"] == 10**800
    assert b["kappa"] / a["kappa"] == sp.Rational(1026, 1024)
    assert a["fixed_gamma"] == b["fixed_gamma"] == sp.Rational(1024, 10**800)
    assert a["fixed_lambda"] == b["fixed_lambda"] == sp.Rational(1, 10**600)


def test_massive_contact_at_crossing_symmetric_point():
    d = vacuum.data()
    j = analytic.local_jets()
    contact = d["analytic_candidate_scalar_contact"]
    n, kap = (
        j["canonical_switch_order_symbol"],
        j["canonical_action_normalization_symbol"],
    )
    gamma = sp.Symbol("g", positive=True)
    lam = j["canonical_coupling_symbol"]
    value = sp.factor(
        contact.subs(
            {
                kap: n / gamma,
                original.mass: 1,
                original.s: sp.Rational(4, 3),
                original.transfer: sp.Rational(4, 3),
                original.w: sp.Rational(4, 3),
            },
            simultaneous=True,
        )
    )
    assert sp.factor(value - sp.Rational(8, 3) * lam + sp.Rational(8, 9) * gamma) == 0


def test_heavy_remainder_formula_accepts_exact_fraction_only():
    r = Fraction(7, 3)
    assert calibration.heavy_remainder_bound(r) == calibration.heavy_remainder_bound(
        sp.Rational(7, 3)
    )
    assert (
        0
        < calibration.heavy_remainder_bound(r)
        < calibration.heavy_remainder_bound(heavy.CHANNEL_RADIUS)
    )


def test_scalar_envelope_bounds_a_known_constant():
    rows = uniform.localized_rational_norm(sp.Integer(1))
    assert len(rows) == 15
    assert rows[(0, 0)] >= 1


def test_native_audit_counts_and_claim_boundary():
    assert len(audit.residuals()) == 103
    assert len(audit.gates()) == 53
    assert len(audit.controls()) == 15
    assert audit.rejected_inputs() == 61
