"""Independent full retarded current, phase, subtraction and norm checks."""

import importlib.util
import sys
from pathlib import Path

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_corrected_spatial_current import assembly, response, subtraction
from p8_vacuum_affine_dimensional_spatial_symbol import density as legacy_density
from p8_vacuum_affine_dimensional_spatial_symbol import geometry as legacy_geometry
from p8_vacuum_affine_dimensional_spatial_symbol import jets as legacy_jets
from p8_vacuum_affine_ordered_scalar_symbol import jets, matching


def two_mode():
    name = "p8_s217_literal_two_mode"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(
            name, Path(__file__).with_name("_two_mode_reference.py")
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


def fixture(which):
    Ds = (
        mp.matrix([[1, 0], [0, 0]]),
        mp.matrix([[0, 1], [2, -1]]),
        mp.matrix([[1, -2], [3, 4]]),
    )
    Gs = (
        mp.matrix([[0, 1], [0, 0]]),
        mp.matrix([[2, -1], [3, 1]]),
        mp.matrix([[-1, 2], [1, 3]]),
    )
    return (
        Ds[which],
        Gs[which],
        mp.mpf(2) / 5,
        mp.mpf(3),
        mp.mpf(1) / 5,
        mp.mpf(7),
        -mp.mpf(1) / 7,
    )


@pytest.mark.parametrize("which", range(3))
def test_exact_unequal_modes_solve_positive_time_dependent_Hamiltonians(which):
    h = two_mode()
    with mp.workdps(65):
        _D, _G, T, w1, b1, w2, b2 = fixture(which)
        J = mp.matrix([[0, 1], [-1, 0]])
        for w, b in ((w1, b1), (w2, b2)):
            for t in (mp.mpf(0), T / 3, T):
                u = h.mode(t, w, b)
                K = mp.exp(2 * b * t)
                V = (w * w + b * b) / K
                assert K > 0 and V > 0
                deriv = mp.matrix(
                    [
                        mp.diff(lambda x, w=w, b=b, j=j: h.mode(x, w, b)[j], t)
                        for j in range(2)
                    ]
                )
                assert mp.norm(deriv - mp.matrix([[0, K], [-V, 0]]) * u) < mp.mpf(
                    "1e-58"
                )
                F = h.fundamental(t, w, b)
                assert mp.norm(F * J * F.T - J) < mp.mpf("1e-58")
                assert abs(2 * mp.im(u[0] * mp.conj(u[1])) - 1) < mp.mpf("1e-58")


@pytest.mark.parametrize("which", range(3))
@pytest.mark.parametrize(
    "fraction", ("0", "0.142857142857142857", "0.666666666666666667", "1")
)
def test_literal_Fock_covariance_and_complex_pair_kernel_at_both_times(which, fraction):
    h = two_mode()
    with mp.workdps(75):
        D, G, T, w1, b1, w2, b2 = fixture(which)
        t = T * mp.mpf(fraction)
        literal = h.fock_integrand(D, G, T, t, w1, b1, w2, b2)
        cov = h.covariance_integrand(D, G, T, t, w1, b1, w2, b2)
        pair = (
            2
            * mp.im(
                mp.conj(h.stripped(D, T, w1, b1, w2, b2))
                * h.stripped(G, t, w1, b1, w2, b2)
                * mp.exp(1j * (w1 + w2) * (T - t))
            )
            * t**6
        )
        assert abs(literal - cov) < mp.mpf("1e-65")
        assert abs(literal - pair) < mp.mpf("1e-65")


@pytest.mark.parametrize("which", range(3))
def test_complete_retarded_integral_all_six_endpoints_and_sixth_bulk(which):
    h = two_mode()
    with mp.workdps(75):
        D, G, T, w1, b1, w2, b2 = fixture(which)
        exact = h.pair_integral(D, G, T, w1, b1, w2, b2)
        endpoints, bulk = h.endpoint_and_bulk(D, G, T, w1, b1, w2, b2)
        cov = mp.quad(
            lambda t: h.covariance_integrand(D, G, T, t, w1, b1, w2, b2), [0, T / 2, T]
        )
        assert len(endpoints) == 6
        assert abs(exact - cov) < mp.mpf("1e-65")
        assert abs(exact - sum(endpoints) - bulk) < mp.mpf("1e-65")
        wrong = h.pair_integral(D, G, T, w1, b1, w2, b2, wrong=True)
        assert abs(wrong - exact) > mp.mpf("1e-10")
        assert abs(bulk) > mp.mpf("1e-10")
        assert abs(endpoints[5]) > mp.mpf("1e-10")
        for j in range(6):
            assert abs(h.source_derivative(G, 0, j, w1, b1, w2, b2)) < mp.mpf("1e-65")


@pytest.mark.parametrize("D", range(4))
@pytest.mark.parametrize("G", range(4))
def test_all_two_mode_source_readout_quadrature_blocks(D, G):
    h = two_mode()
    with mp.workdps(60):
        d = mp.zeros(2)
        g = mp.zeros(2)
        d[D // 2, D % 2] = 1
        g[G // 2, G % 2] = 1
        T, t = mp.mpf(".31"), mp.mpf(".17")
        args = (mp.mpf(5), mp.mpf(".2"), mp.mpf(11), -mp.mpf(".125"))
        literal = h.fock_integrand(d, g, T, t, *args)
        cov = h.covariance_integrand(d, g, T, t, *args)
        assert abs(literal - cov) < mp.mpf("1e-50")


@pytest.mark.parametrize(
    "endpoint,source_jet", ((1, 0), (1, 1), (3, 0), (3, 1), (3, 2), (3, 3))
)
def test_unaveraged_odd_UV_cannot_change_nondecreasing_cutoff_terms(
    endpoint, source_jet
):
    for (geo, j, r), row in legacy_jets.endpoint_products().items():
        if (j, r) != (endpoint, source_jet):
            continue
        for degree, v in enumerate(row):
            value = s.factor(v.subs(legacy_geometry.d, 3))
            assert not value.has(jets.p)
            if j + degree < 4:
                assert value == 0
            if value != 0:
                assert j == 1 and degree == 3


@pytest.mark.parametrize("channel", ("tensor", "vector", "scalar"))
def test_complete_physical_odd_coefficient_is_common_transfer_independent_log(channel):
    expected = jets.m**2 * jets.a**2 * (53 * jets.t**2 + 7) / 15
    for key, v in legacy_density.coefficients().items():
        ch, j, r, degree = key
        if ch != channel or j % 2 == 0:
            continue
        target = expected if (j, r, degree) == (1, 0, 3) else 0
        assert s.factor(v.subs(legacy_geometry.d, 3) - target) == 0


@pytest.mark.parametrize(
    "P,K", (("0.1", "2000"), ("100", "2000"), ("500", "2000"), ("1250", "5000"))
)
def test_exact_odd_log_shell_is_retained_as_decaying_original_mask_error(P, K):
    with mp.workdps(65):
        P, K = mp.mpf(P), mp.mpf(K)
        assert 0 < P <= K / 4

        # A nonconstant anisotropic angular factor checks the unaveraged shell.
        # The actual lower radius along n in the original two-ball intersection:
        def lost(u):
            radius = P * u + mp.sqrt(K * K - P * P * (1 - u * u))
            radius = min(K, radius)
            return (1 + u * u) * mp.log(radius / K)

        value = mp.quad(lost, [-1, 0, P / (2 * K), 1])
        assert value < 0
        assert abs(value) < 2 * mp.mpf(8) / 3 * P / K
        assert abs(value) > mp.mpf("1e-20")


@pytest.mark.parametrize("index", range(3))
def test_corrected_finite_coefficient_keeps_the_required_evanescent_change(index):
    from p8_vacuum_affine_dimensional_spatial_symbol import matching as old

    inv, _g, _ell, rows = matching.finite_input()
    fixed = rows[index].subs(dict.fromkeys(inv[3:], 0))
    difference = s.factor(fixed - old.finite_input()[3][index])
    wanted = (
        jets.p**2 * (3 * jets.t**2 + 1) * (-31 * inv[0] + 30 * inv[1]) / (420 * s.pi**2)
        if index == 0
        else 0
    )
    assert s.factor(difference - wanted) == 0


@pytest.mark.parametrize("K", (2000, 3000, 10000))
def test_actual_original_mask_anchored_subtraction_identity(K):
    m = s.Integer(1000)
    H, KP, K0, LP, L0, err = s.symbols("H KP K0 LP L0 err")
    A3, A1, F2, F4, F = s.symbols("A3 A1 F2 F4 F")
    U0 = s.Symbol("U0")
    UP = U0 + assembly.subtractor(K, m, A3, A1, F2, F4, 0) + err
    actual = H + KP - K0 + UP - U0 - assembly.subtractor(K, m, A3, A1, F2, F4, F)
    limit = assembly.anchored_current(H, LP, L0, F)
    assert s.expand(actual - limit - (KP - LP - (K0 - L0) + err)) == 0


def test_full_normalized_bounds_not_a_same_space_inverse():
    current = assembly.data()
    assert current["current_unrounded_bound"] < assembly.CURRENT
    assert current["tail_unrounded_bound"] < assembly.TAIL
    assert assembly.ward_data()["unrounded_known_bound"] < assembly.WARD_KNOWN
    assert "derivative-losing" in current["canonical_scope"]


@pytest.mark.parametrize(
    "fn",
    (
        response.canonical_data,
        response.decomposition_data,
        subtraction.physical_data,
        subtraction.dimension_data,
        assembly.data,
        assembly.ward_data,
    ),
)
def test_every_packet_exact_residual_and_gate(fn):
    p = fn()
    for key, value in p["checks"].items():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.cancel(v) == 0 for v in entries), key
    assert all(bool(v) for v in p["gates"].values())


@pytest.mark.parametrize(
    "imaginary_dimension",
    (s.Rational(-1, 4), s.Rational(-1, 10), s.Rational(1, 10), s.Rational(1, 4)),
)
def test_holomorphic_branch_pair_not_conjugation_of_dimension(imaginary_dimension):
    d = s.Symbol("d")
    plus = 1 + s.I * (d - 3)
    sharp = 1 - s.I * (d - 3)
    correct = (plus - sharp) / (2 * s.I)
    value = correct.subs(d, 3 + s.I * imaginary_dimension)
    naive = s.im(plus.subs(d, 3 + s.I * imaginary_dimension))
    assert s.simplify(value - s.I * imaginary_dimension) == 0
    assert value != naive


@pytest.mark.parametrize("imaginary_shift", ("-0.025", "-0.01", "0", "0.01", "0.025"))
@pytest.mark.parametrize("direction", (-1, 1))
def test_both_retarded_phase_moduli_have_the_explicit_common_bound(
    imaginary_shift, direction
):
    with mp.workdps(65):
        omega = mp.mpf(31000) + 1j * mp.mpf(imaginary_shift)
        value = abs(mp.exp(direction * 1j * omega))
        assert value < mp.mpf(20) / 19 < 2


from p8_vacuum_affine_corrected_spatial_current import audit


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_every_exact_corrected_current_identity(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(v == 0 for v in entries), name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=lambda v: v if isinstance(v, str) else None
)
def test_every_out_of_scope_or_false_closure_input_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_precise_successor_recovery_not_reendorsement_of_old_formulas():
    before = audit.previous.matching()
    after = audit.matching()
    assert {a["id"] for a, b in zip(before, after) if a != b} == set(audit.CORRECTIONS)
    assert audit.frontier() == audit.previous.frontier()
    assert len(after) == len(before) + 1
    assert "not re-endorsed" in audit.observable()["historical_boundary"]
    assert all(v is True for v in audit.gates().values())


def test_original_cutoff_boundary_and_full_scalar_nonclosure():
    assert audit.require_cutoff(2000) == (2000, 1000)
    assert audit.require_cutoff(s.Rational(4001, 2)) == (s.Rational(4001, 2), 1000)
    assert "three missing" in audit.observable()["Ward_input"]
    assert "NOT_FULL_SCALAR" in audit.ITEM["status"]
