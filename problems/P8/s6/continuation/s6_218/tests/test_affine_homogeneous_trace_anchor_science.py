"""Independent trace ADM, scalar mode, full-volume and response checks."""

import importlib.util
import sys
from functools import cache
from pathlib import Path

import mpmath as mp
import pytest
import sympy as s
from p8_affine.verify import serialize
from p8_vacuum_affine_homogeneous_trace_anchor import (
    anchor,
    comparison,
    homogeneous,
    local,
)
from p8_vacuum_affine_matrix_response_tail import mixed


def literal():
    name = "p8_s218_literal_trace_reference"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(
            name, Path(__file__).with_name("_trace_reference.py")
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


@cache
def vertex_fixture(momentum, epsilon, order):
    ref = literal()
    with mp.workdps(90):
        e = mp.mpf(epsilon)
        direct = mp.diff(lambda x: ref.literal_current_vertex(x, momentum), e, order)
        expected = mp.diff(lambda x: ref.projector_vertex(x, momentum), e, order)
        return direct, expected


@pytest.mark.parametrize("momentum", (0, 1000, 10**12))
@pytest.mark.parametrize("epsilon", ("-.01", "0", ".01"))
@pytest.mark.parametrize("order", range(3))
def test_full_literal_ADM_normalized_trace_vertices_and_all_frame_derivatives(
    momentum, epsilon, order
):
    with mp.workdps(70):
        direct, expected = vertex_fixture(momentum, epsilon, order)
        assert mp.norm(direct - expected) < mp.mpf("1e-40")
        assert mp.norm(direct) < (mp.mpf(3), mp.mpf(1), mp.mpf(1))[order]


@pytest.mark.parametrize("momentum", (0, 1000, 30000))
def test_independent_full_covariance_flow_derivative_and_nondeleted_contact(momentum):
    # A fixed-canonical-data finite IVP verifies the tangent/contact algebra.
    # It is not an alternate all-order preparation or an admitted metric history.
    with mp.workdps(75):
        rebuilt, actual, contact = literal().finite_flow_fixture(momentum)
        assert all(abs(a - b) < mp.mpf("1e-48") for a, b in zip(rebuilt, actual))
        assert abs(contact) > mp.mpf("1e-5")
        assert abs(rebuilt[1] - (actual[1] - contact)) > mp.mpf("1e-5")


@pytest.mark.parametrize("momentum", (0, 1000, 10**9))
@pytest.mark.parametrize("time_order", (0, 1, 4, 8, 11))
@pytest.mark.parametrize("amplitude_order", range(3))
def test_actual_trace_frequency_and_squeeze_mixed_jets_fit_all_old_majorants(
    momentum, time_order, amplitude_order
):
    ref = literal()
    c = mixed.constants()
    # This polynomial tests pointwise finite-jet estimates only; no assertion
    # that it is the theorem's common all-order initial-germ history is made.
    with mp.workdps(75):
        e, t = -mp.mpf(".01"), mp.mpf(1) / 7
        value = ref.point_modes(e, t, momentum)
        derivatives = [
            mp.diff(
                lambda ep, j=j: mp.diff(
                    lambda tt: ref.point_modes(ep, tt, momentum)[j], t, time_order
                ),
                e,
                amplitude_order,
            )
            for j in range(3)
        ]
        bounds = [
            c["omega"][time_order, amplitude_order],
            c["S"][time_order, amplitude_order],
            c["S"][time_order, amplitude_order],
        ]
        derivatives[0] /= value[0]
        for actual, bound in zip(derivatives, bounds):
            limit = mp.mpf(str(s.N(bound, 80)))
            assert abs(actual) <= limit * (1 + mp.mpf("1e-65"))


@pytest.mark.parametrize("order", range(12))
@pytest.mark.parametrize(
    "point",
    (s.Rational(-1, 2), s.Rational(-1, 7), 0, s.Rational(2, 7), s.Rational(1, 2)),
)
def test_each_original_H_time_majorant_with_both_slab_endpoints(order, point):
    value = s.diff(homogeneous.H, homogeneous.t, order).subs(homogeneous.t, point)
    assert abs(value) <= homogeneous.h_majorants()[order]


def test_independent_scale_factor_Euler_current_includes_all_volume_derivatives():
    A, mass, _F, current, rates = literal().literal_euler()
    target = A[0] ** 3 * local.current_polynomial().subs(local.m, mass).subs(
        dict(zip(local.h[:4], rates)), simultaneous=True
    )
    assert s.factor(current - target) == 0


@pytest.mark.parametrize("which", range(2))
def test_literal_integrated_full_chart_Hessian_not_only_pointwise_pole(which):
    ref = literal()
    with mp.workdps(75):
        t, a, G, D, mixed_density = ref.compact_hessian_fixture(which)
        H = 4 * t / (1 + t * t)
        substitute = {local.m: 1000, **{local.h[j]: s.diff(H, t, j) for j in range(4)}}
        paired = s.factor(
            a**3
            * D
            * sum(
                v.subs(substitute) * s.diff(G, t, j)
                for j, v in enumerate(local.trace_linear_coefficients())
            )
        )
        direct = mp.quad(
            s.lambdify(t, mixed_density, "mpmath"), [-mp.mpf(".5"), 0, mp.mpf(".5")]
        )
        actual = mp.quad(
            s.lambdify(t, paired, "mpmath"), [-mp.mpf(".5"), 0, mp.mpf(".5")]
        )
        assert abs(direct - actual) < mp.mpf("1e-48")
        assert abs(direct) > 1


@pytest.mark.parametrize(
    "dimension", (s.Rational(11, 4), 3, s.Rational(25, 8), 4, 5, 6)
)
def test_general_dimension_pole_identity_and_nonzero_volume_boundary(dimension):
    p = local.data()
    dd, H, H1, second, fourth = literal().literal_dimensional_action()
    fixed = {dd: dimension, H: local.H, H1: local.H1}
    assert (
        s.factor(
            second.subs(fixed)
            - p["full_general_d_symbols"]["second"].subs(local.d, dimension)
        )
        == 0
    )
    assert (
        s.factor(
            fourth.subs(fixed)
            - p["full_general_d_symbols"]["fourth"].subs(local.d, dimension)
        )
        == 0
    )
    assert (
        p["full_general_d_symbols"]["fourth_weighted_boundary_coefficient"].subs(
            local.d, dimension
        )
        != 0
    )


@pytest.mark.parametrize("jet", range(5))
def test_full_density_local_trace_Green_relation(jet):
    assert local.data()["checks"]["all_five_full_density_Green_relations"][jet] == 0


@pytest.mark.parametrize("mass", (1000, 2000))
def test_deleted_cosmological_volume_is_a_nonzero_negative_control(mass):
    # Other mass is an algebra-only negative control, not an admitted P8 input.
    flat = (
        local.current_polynomial().subs(dict.fromkeys(local.h, 0)).subs(local.m, mass)
    )
    assert flat == s.Rational(5, 2) * mass**4
    assert flat != 0
    coefficient = (
        local.trace_linear_coefficients()[0]
        .subs(dict.fromkeys(local.h, 0))
        .subs(local.m, mass)
    )
    assert coefficient == flat


@pytest.mark.parametrize(
    "fraction", (s.Rational(1, 100), s.Rational(1, 2), s.Rational(99, 100))
)
def test_deleted_longitudinal_trace_constraint_changes_the_physical_vertex(fraction):
    actual = homogeneous.scalar_vertices()["longitudinal"][1, 1].subs(
        homogeneous.z, fraction
    )
    deleted = -(1 - fraction) / 3
    assert s.factor(actual - deleted + fraction) == 0
    assert actual != deleted


@pytest.mark.parametrize("order", range(3))
def test_complete_local_raw_derivative_majorant_includes_varying_volume(order):
    assert 0 < local.local_bounds()[order] < local.LOCAL
    assert local.local_bounds()[order] > 10**10
    assert anchor.data()["unrounded_trace_current_sums"][order] < anchor.CURRENT[order]


@pytest.mark.parametrize(
    "fn",
    (
        homogeneous.data,
        homogeneous.domination_data,
        local.data,
        comparison.data,
        anchor.data,
    ),
)
def test_every_trace_packet_exact_residual_and_proof_gate(fn):
    packet = fn()
    serialize(packet)
    for key, value in packet["checks"].items():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.cancel(v) == 0 for v in entries), key
    assert all(bool(v) for v in packet["gates"].values())


def test_same_parent_reference_Hessian_not_full_nonlinear_source_free_parent():
    assert "finite-amplitude remainder" in comparison.data()["parent_boundary"]
    assert (
        "not a uniform finite-amplitude" in anchor.data()["finite_amplitude_boundary"]
    )


def test_normalized_trace_anchor_and_orthogonal_full_Hilbert_lift_scope():
    p = anchor.data()
    assert s.Rational(3, 4) * anchor.CURRENT[1] < anchor.FULL_HOMOGENEOUS
    assert "or a full inhomogeneous current" in p["full_homogeneous_Hilbert_lift"]
    assert "not a reduced inverse" in p["canonical_scope"]


from p8_vacuum_affine_homogeneous_trace_anchor import audit


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_every_exact_homogeneous_trace_anchor_identity(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(v == 0 for v in entries), name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=lambda v: v if isinstance(v, str) else None
)
def test_every_unsupported_trace_family_or_false_closure_input(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_original_trace_family_edges_and_full_source_jet_domain():
    assert audit.require_family(s.Rational(-1, 2), s.Rational(-1, 100))[:3] == (
        s.Rational(-1, 2),
        s.Rational(-1, 100),
        1,
    )
    assert audit.require_family(s.Rational(1, 2), s.Rational(1, 100), 0)[:3] == (
        s.Rational(1, 2),
        s.Rational(1, 100),
        0,
    )
    assert audit.require_family(0, 0, s.Rational(2, 3))[2] == s.Rational(2, 3)


def test_complete_homogeneous_anchor_not_a_full_nonzero_transfer_scalar_current():
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == len(audit.previous.matching()) + 1
    assert "NONZERO_TRANSFER_SCALAR_REMAINDERS" in audit.ITEM["status"]
    assert "three nonzero-transfer" in audit.observable()["remaining"]
    assert all(v is True for v in audit.gates().values())
