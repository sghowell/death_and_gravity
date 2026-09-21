"""Exact restricted curvature basis, selected UV poles and strict matching scope."""

import pytest
import sympy as s
from p8_vacuum_affine_curvature_contact_basis import (
    audit,
    basis,
    calibration,
    matching,
    source,
    uvpoles,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_residual(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_all_rejected_scopes(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "coefficient", (s.S.Zero, s.S.One, -s.Rational(3, 7), s.sqrt(2), s.log(2))
)
def test_entire_matching_polynomial(coefficient):
    assert (
        matching.extract_local_coefficient(coefficient * basis.canonical) == coefficient
    )


@pytest.mark.parametrize(
    "expression",
    (
        basis.q(0, 1, 0, 1),
        sum(a * a for a in basis.a) * basis.canonical,
        basis.canonical / (1 + basis.a0[0]),
        s.sin(basis.a0[0]),
        s.Symbol("unmatched") * basis.canonical,
    ),
)
def test_out_of_span_or_nonlocal_response_rejected(expression):
    with pytest.raises(ValueError):
        matching.extract_local_coefficient(expression)


def test_all256_multiplicities_and_rank():
    packet = basis.data()
    assert packet["whole_Bose_coefficient_multiplicities"] == {
        0: 136,
        4: 12,
        -2: 48,
        -4: 12,
        2: 48,
    }
    assert packet["whole_symbolic_basis_rank"] == 1
    assert sum(packet["whole_Bose_coefficient_multiplicities"].values()) == 256


def test_soft_and_total_derivative_degrees_are_distinct():
    t = s.Symbol("scaling")
    soft = basis.canonical.subs({a: t * a for a in basis.a0}, simultaneous=True)
    hard = basis.canonical.subs(
        {**{a: t * t * a for a in basis.a0}, **{h: t * t * h for h in basis.h0}},
        simultaneous=True,
    )
    assert s.expand(soft - t * t * basis.canonical) == 0
    assert s.expand(hard - t**6 * basis.canonical) == 0


def test_original_literal_curvature_calibrations():
    assert calibration.data()["whole_original_literal_counts"] == {
        "states": 3,
        "polarizations": 6,
        "mixed_pair_contractions": 126,
        "nonzero_canonical_polarizations": 6,
    }


def test_selected_pole_zero_does_not_assign_finite_chi():
    packet = uvpoles.data()
    assert packet["whole_selected_radiative_UV_residual"] == 0
    assert len(packet["whole_six_original_offshell_local_pole_polynomials"]) == 6
    assert "does not determine" in packet["whole_unassigned_finite_boundary"]
    assert "chi" in audit.observable()["not_established"]


def test_whole_source_and_frontiers_unchanged():
    n, g = source.HEAVY_MASS2, source.CUBIC
    assert source.CONTACT == -g * g * (3 / (n - 2) - 2 / (n - 2) ** 2)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 199
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6


def test_actual_exact_module_counts():
    assert len(basis.data()["checks"]) == 931
    assert len(matching.data()["checks"]) == 8
    assert len(calibration.data()["checks"]) == 135
    assert len(uvpoles.data()["checks"]) == 86
    assert len(source.data()["checks"]) == 275


def test_scope_excludes_higher_and_full_curved_claims():
    boundary = basis.data()["whole_scope"]
    assert "higher derivative" in boundary and "parity-odd" in boundary
    assert "curved actions" in boundary
    assert "internal-graviton" in audit.observable()["not_established"]
