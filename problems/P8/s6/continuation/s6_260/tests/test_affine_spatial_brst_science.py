"""Independent exterior, matrix, actual-integral and scope diagnostics."""

import math
from itertools import combinations

import numpy as np
import pytest
import sympy as s
from p8_affine import verify as certificate
from p8_vacuum_affine_spatial_brst import audit, brst, ward
from p8_vacuum_affine_spatial_brst.exterior import G, Jets
from scipy import integrate, linalg


@pytest.mark.parametrize("packet", tuple(audit.packets()))
def test_complete_exact_packet_and_gates(packet):
    data = audit.packets()[packet]
    assert certificate.certify_residuals(data["checks"])
    assert all(value is True for value in data["gates"].values())


@pytest.mark.parametrize(
    "name,call,args",
    audit.bad_cases(),
    ids=lambda value: value if isinstance(value, str) else None,
)
def test_every_invalid_input_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_exact_full_counts():
    assert len(audit.residuals()) == 45 and audit.scalar_entry_count() == 413
    assert len(audit.gates()) == 34 and all(audit.gates().values())
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 259
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 116
    assert audit.matching()[:-1] == audit.previous.matching()


ODD = ((1, ()), (1, (0,)), (2, (1,)), (3, (2, 3)))


def creation_matrices():
    unit = s.eye(2)
    parity = s.diag(1, -1)
    create = s.Matrix([[0, 0], [1, 0]])
    return [
        s.kronecker_product(*([parity] * i + [create] + [unit] * (3 - i)))
        for i in range(4)
    ]


def matrix_rep(value):
    creators = creation_matrices()
    result = s.zeros(16)
    for mon, coef in value.terms.items():
        part = s.eye(16)
        for generator in mon:
            part = part * creators[ODD.index(generator)]
        result += coef * part
    return result


@pytest.mark.parametrize("seed", range(12))
def test_exterior_product_against_independent_full_Fock_matrices(seed):
    rng = np.random.default_rng(seed)
    subsets = [
        tuple(ODD[i] for i in ids)
        for count in range(5)
        for ids in combinations(range(4), count)
    ]
    first = G({key: s.Rational(int(rng.integers(-3, 4)), 5) for key in subsets})
    second = G({key: s.Rational(int(rng.integers(-3, 4)), 7) for key in subsets})
    assert matrix_rep(first * second) == matrix_rep(first) * matrix_rep(second)


@pytest.mark.parametrize("mu", range(4))
@pytest.mark.parametrize("nu", range(4))
def test_whole_jet_derivatives_commute(mu, nu):
    jet = Jets()
    f = jet.add_tensor("scalar", 0, ())[()]
    W = jet.add_tensor("vector", 1, (1,))
    expression = (f * f + W[1,] * W[2,]) * jet.ghost(1, (0,)) + W[3,] * jet.ghost(2)
    assert not jet.d(jet.d(expression, mu), nu) - jet.d(jet.d(expression, nu), mu)
    assert not jet.brst(jet.d(expression, mu)) - jet.d(jet.brst(expression), mu)


@pytest.mark.parametrize("i", range(1, 4))
def test_graded_Leibniz_and_genuine_ghost_sign(i):
    jet = Jets()
    scalar = jet.add_tensor("scalar", 0, ())[()]
    odd = jet.ghost(i) * scalar
    other = jet.ghost(i % 3 + 1, (0,)) * scalar * scalar
    assert not jet.brst(odd * other) - jet.brst(odd) * other + odd * jet.brst(other)
    assert not (odd * other + other * odd)


def test_omitting_upper_time_connection_term_fails_trace_covariance():
    jet = Jets()
    gamma = jet.add_affine()
    trace = sum((gamma[a, a, 1] for a in range(4)), G())
    correct = jet.brst(trace)
    for a in range(1, 4):
        for b in range(4):
            for c in range(4):
                jet.rules["Gamma" + str(a) + str(b) + str(c)] += (
                    jet.ghost(a, (0,)) * gamma[0, b, c]
                )
    assert jet.brst(trace) - correct


@pytest.mark.parametrize("a", range(1, 4))
def test_nonzero_inhomogeneous_affine_term_cannot_be_omitted(a):
    jet = Jets()
    jet.add_affine()
    rule = jet.rules["Gamma" + str(a) + "00"]
    zero_connection = {symbol: 0 for symbol in jet.reverse}
    remaining = G({key: coef.subs(zero_connection) for key, coef in rule.terms.items()})
    assert not remaining - jet.ghost(a, (0, 0))
    assert remaining


@pytest.mark.parametrize(
    "k", [(1, 0, 0), (1, 2, 3), (-2, 1, 1), (s.Rational(1, 2), -1, 2)]
)
@pytest.mark.parametrize("omega", [-2, 0, s.Rational(3, 2)])
def test_full_seven_ghost_matrix_numerically(k, omega):
    packet = brst.projective_ghost_block()
    sub = dict(zip(packet["whole_three_nonzero_spatial_wavevector"], k, strict=True))
    sub[packet["time_frequency"]] = omega
    matrix = np.array(packet["whole_seven_ghost_gauge_matrix"].subs(sub), dtype=float)
    inverse = np.array(packet["whole_seven_ghost_inverse"].subs(sub), dtype=float)
    assert np.isfinite(matrix).all() and np.isfinite(inverse).all()
    np.testing.assert_allclose(
        linalg.solve(matrix, np.eye(7), check_finite=True),
        inverse,
        rtol=1e-12,
        atol=1e-12,
    )
    expected = 1024 / 3 * sum(float(value) ** 2 for value in k) ** 3
    assert linalg.det(matrix, check_finite=True) == pytest.approx(expected, rel=1e-12)
    audit.require_wavevector(k)


@pytest.mark.parametrize("v", [0.25, 1.0, 2.0])
@pytest.mark.parametrize("t", [0.5, 1.0])
@pytest.mark.parametrize("Jq,Jy", [(0.2, -0.1), (-0.3, 0.1)])
def test_source_dependent_gauge_orbit_by_actual_integrals(v, t, Jq, Jy):
    packet = ward.orbit()

    def weight(q):
        return math.exp(-q * q / (2 * v) + Jq * q + Jy * t * q * q) / math.sqrt(
            2 * math.pi * v
        )

    Z = integrate.quad(weight, -np.inf, np.inf, epsabs=1e-11, epsrel=1e-11)[0]
    qmean = (
        integrate.quad(lambda q: q * weight(q), -np.inf, np.inf, epsabs=1e-11)[0] / Z
    )
    qsecond = (
        integrate.quad(lambda q: q * q * weight(q), -np.inf, np.inf, epsabs=1e-11)[0]
        / Z
    )
    sub = {packet["variance"]: v, packet["gauge_parameter"]: t}
    sub.update(dict(zip(packet["sources"], (Jq, Jy), strict=True)))
    assert math.log(Z) == pytest.approx(
        float(packet["whole_connected_generator"].subs(sub)), rel=1e-10, abs=1e-11
    )
    symbolic = packet["whole_retained_first_and_second_q_moments"].subs(sub)
    assert qmean == pytest.approx(float(symbolic[0]), rel=1e-10, abs=1e-11)
    assert qsecond == pytest.approx(float(symbolic[1]), rel=1e-10, abs=1e-11)
    assert t * qsecond == pytest.approx(
        float(packet["whole_gauge_coordinate_mean"].subs(sub)), rel=1e-10, abs=1e-11
    )


@pytest.mark.parametrize("v,eta", [(0.5, 0.2), (1.0, 0.3), (1.2, 0.2)])
def test_noninvariant_boundary_weight_changes_physical_mean(v, eta):
    def mean(t):
        def w(q):
            return math.exp(-q * q / (2 * v) - eta * t * t * q**4 / 2)

        denominator = integrate.quad(w, -np.inf, np.inf, epsabs=1e-12)[0]
        return (
            integrate.quad(lambda q: q * q * w(q), -np.inf, np.inf, epsabs=1e-12)[0]
            / denominator
        )

    eps = 1e-3
    actual = 2 * (mean(eps) - mean(0)) / eps**2
    assert actual == pytest.approx(-12 * eta * v**3, rel=3e-5, abs=1e-7)


def test_all_actual_source_and_state_boundaries_stay_open():
    packet = brst.projective_ghost_block()
    assert packet["gates"]["original_nonzero_conditional_source_contact_retained"]
    assert "OPEN" in audit.observable()["original_problem"]
    assert ward.differentiated_identity()[
        "whole_off_shell_onepoint_contact"
    ] != s.zeros(2)
    assert audit.require_orbit(s.Rational(1, 2), s.Rational(1, 3), s.Rational(1, 4))


def test_ghost_cache_does_not_keep_discarded_contexts_alive():
    import gc
    import weakref

    jet = Jets()
    jet.sg(1)
    reference = weakref.ref(jet)
    del jet
    gc.collect()
    assert reference() is None


def test_instance_ghost_cache_tracks_rule_replacement_in_diagnostics():
    jet = Jets()
    first, second = jet.field("first"), jet.field("second")
    jet.odd_rules[10] = first
    assert not jet.sg(10) - first
    jet.odd_rules[10] = second
    assert not jet.sg(10) - second
