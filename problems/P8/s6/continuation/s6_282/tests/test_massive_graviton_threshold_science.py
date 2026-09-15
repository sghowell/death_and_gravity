"""Independent complete normal-cut, sheet, threshold and pole-anchor checks."""

import pytest
import sympy as s
from p8_affine import verify as base
from p8_vacuum_affine_massive_graviton_threshold_sheet import (
    angular,
    audit,
    endpoint,
    invariant,
    source,
)

ROWS = audit.residuals()
GATES = audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not value.atoms(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_every_unsupported_input(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "ss,tt", ((10, -2), (5, -s.Rational(1, 3)), (2, 1), (4, 0), (1, s.Rational(3, 2)))
)
def test_entire_normal_sheet_including_subthreshold(ss, tt):
    ss, tt = s.sympify(ss), s.sympify(tt)
    assert source.require_real_sheet(ss, tt) == (ss, tt, 4 - ss - tt, s.Integer(1))


@pytest.mark.parametrize("mass", (s.Integer(1), s.Rational(3, 2), s.Integer(7)))
def test_exact_kernel_removable_zero_and_center(mass):
    assert invariant.real_kernel(0, mass) == 1 / (4 * mass)
    assert invariant.real_kernel(2 * mass, mass) == s.pi / (8 * mass)
    assert invariant.kernel(0, mass) == 1 / (4 * mass)
    assert invariant.real_kernel(-4 * mass, mass) == s.atanh(s.sqrt(2) / 2) / (
        4 * s.sqrt(2) * mass
    )


@pytest.mark.parametrize(
    "h,z",
    (
        (s.Integer(2), s.Integer(0)),
        (s.Integer(4), s.Rational(1, 3)),
        (s.Rational(7, 3), s.Rational(-2, 5)),
    ),
)
def test_both_physical_seed_symmetries(h, z):
    j0 = (angular.positive_seed(h, z) + angular.positive_seed(h, -z)) / (2 * h)
    j1 = (angular.positive_seed(h, z) - angular.positive_seed(h, -z)) / 2
    l0 = s.atanh(1 / s.sqrt(h)) / s.sqrt(h)
    assert (
        s.expand(
            angular.complete_average(h, z, j0, j1, l0)
            - angular.complete_average(h, -z, j0, -j1, l0)
        )
        == 0
    )


@pytest.mark.parametrize(
    "tau,z",
    (
        (s.Integer(1), s.Integer(0)),
        (s.Integer(2), s.Rational(1, 3)),
        (s.Rational(3, 2), s.Rational(-2, 5)),
    ),
)
def test_complete_unphysical_seed_symmetries(tau, z):
    h = -tau * tau
    j0 = (angular.negative_seed(tau, z) + angular.negative_seed(tau, -z)) / (2 * h)
    j1 = (angular.negative_seed(tau, z) - angular.negative_seed(tau, -z)) / 2
    l0 = -s.atan(1 / tau) / tau
    assert (
        s.expand(
            angular.complete_average(h, z, j0, j1, l0)
            - angular.complete_average(h, -z, j0, -j1, l0)
        )
        == 0
    )


@pytest.mark.parametrize("tau", (s.Integer(1), s.Integer(2), s.Rational(3, 2)))
def test_whole_unphysical_forward_not_a_truncated_series(tau):
    h = -tau * tau
    l0 = -s.atan(1 / tau) / tau
    j0 = 1 / (2 * h * (h - 1)) + l0 / (2 * h)
    j1 = h * j0 - l0
    actual = angular.complete_average(h, s.Integer(1), j0, j1, l0) / h**2
    assert s.simplify(actual - endpoint.imaginary_beta_forward(1 / tau)) == 0


@pytest.mark.parametrize("mass", (s.Integer(1), s.Rational(3, 2), s.Integer(7)))
def test_regular_external_threshold_at_nonzero_transfer(mass):
    assert s.simplify(
        invariant.regular_threshold_combination(4 * mass, mass).doit()
    ) == s.Rational(29, 40)
    k = s.Integer(11)
    assert (
        s.simplify(
            invariant.whole_cut(4 * mass, 0, mass, k)
            - mass * mass / (16 * s.pi * k * k)
        )
        == 0
    )


@pytest.mark.parametrize("mass", (s.Integer(1), s.Rational(3, 2), s.Integer(7)))
def test_kernel_complementary_continuation(mass):
    assert (
        s.simplify(
            invariant.real_kernel(mass, mass)
            + invariant.real_kernel(3 * mass, mass)
            - s.pi / (2 * s.sqrt(3) * mass)
        )
        == 0
    )


def test_exact_integer_arguments_do_not_introduce_float():
    assert not invariant.whole_cut(4, 0, 1, 11).has(s.Float)
    assert not invariant.regular_threshold_combination(4, 1).has(s.Float)
    assert invariant.kernel(0, 1) == s.Rational(1, 4)


def test_full_source_contract_never_mutates_parent():
    before = dict(source.previous.data()["checks"])
    source.data()
    assert before == source.previous.data()["checks"]
    assert source.data()["checks"] is not source.previous.data()["checks"]
    assert len(before) == 36


def test_whole_packet_is_invariant_under_full_parent_first_order():
    before = base.serialize(audit.packets())
    audit.previous.packets()
    source.data.cache_clear()
    audit.packets.cache_clear()
    audit.residuals.cache_clear()
    audit.gates.cache_clear()
    assert base.serialize(audit.packets()) == before
    assert audit.residuals() == ROWS and audit.gates() == GATES


def test_all_old_frontiers_qualifications_and_rejection_preserved():
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 138 and len(audit.qualifications()) == 6
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert audit.require_parameters(audit.parameters()) == audit.parameters()


@pytest.mark.parametrize("mu,cap,external", ((1, 10, -1), (2, 17, -3), (3, 20, -4)))
def test_pole_running_removes_exact_upper_shell_error(mu, cap, external):
    mu, cap, external = map(s.sympify, (mu, cap, external))
    rho = s.Rational(3, 7)
    derivative = (
        -rho / (s.pi * external)
        + rho * endpoint.compensated_kernel(cap, external) / s.pi
    )
    assert s.factor(derivative - rho / (s.pi * (cap - external))) == 0
    assert s.factor(
        rho * endpoint.compensated_kernel(cap, external) / s.pi
        - rho / (s.pi * (cap - external))
    ) == rho / (s.pi * external)
    v = s.Symbol("v", real=True)
    assert s.diff(1 / (2 * mu + v) + 1 / (2 * mu - v), v, 2).subs(v, 0) / 2 == 1 / (
        4 * mu**3
    )


def test_endpoint_compensation_not_a_value_for_original_pole():
    assert "initial c_L0" in endpoint.data()["required_matching"]
    assert "not the full original amplitude" in endpoint.data()["pole_lift"]
    assert "full_total_discontinuity" not in audit.OBSERVABLES


def test_exact_counts_and_controls():
    assert (len(ROWS), audit.scalar_entry_count(), len(GATES)) == (106, 106, 36)
    assert len(audit.controls()) == 8 and audit.rejected_inputs() == 73
    assert len(audit.frontier()) == 9
